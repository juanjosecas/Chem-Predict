"""Rule-based degradation engine.

The engine executes rules; it does not decide whether a rule is scientifically
valid. Rule curation and scoring remain separate concerns.
"""

from __future__ import annotations

from chem_predict.chemistry import apply_reaction, canonicalize_smiles
from chem_predict.core import Conditions, Prediction
from chem_predict.rules import RuleRegistry


class DegradationEngine:
    def __init__(self, registry: RuleRegistry) -> None:
        self.registry = registry

    def predict(
        self,
        parent_smiles: str,
        conditions: Conditions,
        *,
        max_products_per_rule: int = 1000,
    ) -> list[Prediction]:
        parent = canonicalize_smiles(parent_smiles)
        predictions: list[Prediction] = []

        matching_rules = [
            rule
            for rule in self.registry.enabled()
            if rule.when.matches(conditions)
        ]
        matching_rules.sort(key=lambda rule: (-rule.priority, rule.id))

        for rule in matching_rules:
            outcomes = apply_reaction(
                rule.reaction_smarts,
                [parent],
                max_products=max_products_per_rule,
            )
            for products in outcomes:
                predictions.append(
                    Prediction(
                        parent_smiles=parent,
                        rule_id=rule.id,
                        rule_name=rule.name,
                        products=products,
                        priority=rule.priority,
                        source=rule.source,
                        metadata={"rule_metadata": dict(rule.metadata)},
                    )
                )

        return predictions
