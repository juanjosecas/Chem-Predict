"""Explicit-scale adapter for externally trained reaction-yield predictors."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import isfinite

from chem_predict.applicability import DomainAssessment, SimilarityDomain, split_reaction_smiles


@dataclass(frozen=True, slots=True)
class YieldPrediction:
    reaction_smiles: str
    yield_percent: float
    model_id: str
    domain: DomainAssessment | None
    warnings: tuple[str, ...]


class YieldPredictor:
    """Wrap a batch callable; standardization parameters come from training.

    No model is downloaded or trained implicitly. Native predictors use the
    supplied callable; from_rxn_yields_model adapts SimpleTransformers' tuple
    return. Input reaction strings are validated but never rewritten/tokenized.
    """

    def __init__(self, predict_batch: Callable[[list[str]], Sequence[float]], *,
                 model_id: str, output_scale: str = "percent",
                 training_mean: float | None = None, training_std: float | None = None,
                 domain: SimilarityDomain | None = None):
        if not model_id.strip():
            raise ValueError("model_id is required for provenance")
        if output_scale not in {"percent", "fraction", "standardized"}:
            raise ValueError("output_scale must be percent, fraction or standardized")
        if output_scale == "standardized":
            if (training_mean is None or training_std is None or
                not isfinite(training_mean) or not isfinite(training_std) or training_std <= 0):
                raise ValueError("Provide finite training_mean and positive training_std in percent units")
        elif training_mean is not None or training_std is not None:
            raise ValueError("training_mean/std are only used with standardized outputs")
        if domain is not None and domain.mode != "reaction":
            raise ValueError("Yield coverage requires a reaction-mode domain")
        self.predict_batch = predict_batch
        self.model_id, self.output_scale = model_id, output_scale
        self.training_mean, self.training_std, self.domain = training_mean, training_std, domain

    @classmethod
    def from_rxn_yields_model(cls, model, **options) -> YieldPredictor:
        """Adapt an already loaded SimpleTransformers regression model."""
        def predict(reactions):
            predictions, _raw_outputs = model.predict(reactions)
            return predictions
        return cls(predict, **options)

    def predict(self, reactions: Sequence[str]) -> list[YieldPrediction]:
        if isinstance(reactions, str):
            raise TypeError("Pass a list of reaction SMILES")
        inputs = list(reactions)
        for text in inputs:
            split_reaction_smiles(text)
        if not inputs:
            return []
        outputs = list(self.predict_batch(inputs))
        if len(outputs) != len(inputs):
            raise ValueError("Predictor returned a different number of yields than inputs")
        results = []
        for text, output in zip(inputs, outputs, strict=True):
            value = float(output)
            if not isfinite(value):
                raise ValueError("Predictor returned a non-finite yield")
            if self.output_scale == "fraction":
                value *= 100
            elif self.output_scale == "standardized":
                value = value * self.training_std + self.training_mean
            warnings = []
            if not 0 <= value <= 100:
                warnings.append("yield_outside_physical_range")
            assessment = self.domain.assess(text) if self.domain is not None else None
            if assessment is None:
                warnings.append("applicability_domain_not_assessed")
            elif not assessment.in_domain:
                warnings.append("outside_structural_reference_domain")
            results.append(YieldPrediction(text, value, self.model_id, assessment, tuple(warnings)))
        return results


__all__ = ["YieldPrediction", "YieldPredictor"]
