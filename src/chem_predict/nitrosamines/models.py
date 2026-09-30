from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any
from math import isfinite


class AmineKind(StrEnum):
    SECONDARY = "secondary_amine"
    TERTIARY = "tertiary_amine"
    DIMETHYL_TERTIARY = "dimethyl_tertiary_amine"
    QUATERNARY = "quaternary_ammonium"


@dataclass(frozen=True, slots=True)
class NitrosamineSite:
    nitrosamine_n_index: int
    nitroso_n_index: int
    oxygen_index: int
    alpha_carbon_indices: tuple[int, int]
    alpha_hydrogens: tuple[int, int]
    cpca_applicable: bool
    exclusion_reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class NitrosatableCenter:
    nitrogen_index: int
    kind: AmineKind
    carbon_neighbor_indices: tuple[int, ...]
    methyl_substituent_count: int
    total_hydrogens: int
    aromatic: bool
    amide_like: bool


@dataclass(frozen=True, slots=True)
class CPCAFeature:
    id: str
    label: str
    score: int


@dataclass(frozen=True, slots=True)
class CPCASiteResult:
    site: NitrosamineSite
    potency_score: int | None
    alpha_hydrogen_score: int | None
    deactivating_features: tuple[CPCAFeature, ...]
    activating_features: tuple[CPCAFeature, ...]
    category: int | None
    ai_ng_per_day: float | None
    rationale: tuple[str, ...]
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class CPCAAssessment:
    smiles: str
    sites: tuple[CPCASiteResult, ...]
    overall_category: int | None
    overall_ai_ng_per_day: float | None
    warnings: tuple[str, ...]
    source: str
    source_version: str


@dataclass(slots=True)
class NitrosationContext:
    nitrite_present: bool | None = None
    nitrosating_agent_present: bool | None = None
    ph: float | None = None
    temperature_c: float | None = None
    storage: bool | None = None
    excipient_nitrite_ppm: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("ph", "temperature_c", "excipient_nitrite_ppm"):
            value = getattr(self, name)
            if value is not None and not isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.temperature_c is not None and self.temperature_c <= -273.15:
            raise ValueError("temperature_c must be above absolute zero")
        if self.excipient_nitrite_ppm is not None:
            if self.excipient_nitrite_ppm < 0:
                raise ValueError("excipient_nitrite_ppm must be non-negative")
            if self.excipient_nitrite_ppm > 0 and self.nitrite_present is False:
                raise ValueError("Positive nitrite ppm contradicts nitrite_present=False")


@dataclass(frozen=True, slots=True)
class NitrosationAssessment:
    centers: tuple[NitrosatableCenter, ...]
    structural_precursor_present: bool
    nitrosating_source_supported: bool | None
    acidic_conditions: bool | None
    flags: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ReactivityEvidence:
    id: str
    nitrosamine_scope: str
    reagent_or_condition: str
    outcome: str
    evidence_type: str
    source_doi: str
    source_url: str
    notes: str | None = None
