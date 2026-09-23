"""Shared domain models.

These classes intentionally contain no RDKit objects so they remain easy to
serialize, test, store, and use from other front ends.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any


class StressType(StrEnum):
    ACID = "acid"
    BASE = "base"
    NEUTRAL = "neutral"
    OXIDATION = "oxidation"
    REDUCTION = "reduction"
    PHOTOLYSIS = "photolysis"
    THERMAL = "thermal"
    HUMIDITY = "humidity"
    EXCIPIENT = "excipient"
    NITROSATION = "nitrosation"
    CUSTOM = "custom"


@dataclass(slots=True)
class Conditions:
    """Conditions supplied for a prediction run.

    Missing values mean "unknown", not zero or false.
    """

    stresses: frozenset[StressType] = field(default_factory=frozenset)
    ph: float | None = None
    temperature_c: float | None = None
    oxygen: bool | None = None
    light: bool | None = None
    duration_h: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ConditionWindow:
    """Applicability constraints attached to a rule."""

    stresses: frozenset[StressType] = field(default_factory=frozenset)
    ph_min: float | None = None
    ph_max: float | None = None
    temperature_min_c: float | None = None
    temperature_max_c: float | None = None
    requires_oxygen: bool | None = None
    requires_light: bool | None = None

    def matches(self, conditions: Conditions) -> bool:
        if self.stresses:
            if not conditions.stresses or self.stresses.isdisjoint(conditions.stresses):
                return False

        if self.ph_min is not None:
            if conditions.ph is None or conditions.ph < self.ph_min:
                return False
        if self.ph_max is not None:
            if conditions.ph is None or conditions.ph > self.ph_max:
                return False

        if self.temperature_min_c is not None:
            if conditions.temperature_c is None or conditions.temperature_c < self.temperature_min_c:
                return False
        if self.temperature_max_c is not None:
            if conditions.temperature_c is None or conditions.temperature_c > self.temperature_max_c:
                return False

        if self.requires_oxygen is not None and conditions.oxygen is not self.requires_oxygen:
            return False
        if self.requires_light is not None and conditions.light is not self.requires_light:
            return False

        return True


@dataclass(slots=True)
class Rule:
    """A transformation rule plus its applicability and provenance."""

    id: str
    name: str
    reaction_smarts: str
    when: ConditionWindow = field(default_factory=ConditionWindow)
    priority: int = 0
    enabled: bool = True
    source: str | None = None
    description: str | None = None
    tags: frozenset[str] = field(default_factory=frozenset)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class Prediction:
    """One rule-derived product set.

    Priority is copied from the rule for deterministic ordering. It is not
    a calibrated probability.
    """

    parent_smiles: str
    rule_id: str
    rule_name: str
    products: tuple[str, ...]
    priority: int = 0
    source: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
