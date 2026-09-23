from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any


class MedchemSeverity(StrEnum):
    DEMERIT = "demerit"
    REJECT = "reject"


@dataclass(frozen=True, slots=True)
class MedchemRule:
    id: str
    name: str
    smarts: str
    severity: MedchemSeverity
    demerit: int = 0
    source_rule: str | None = None
    source_version: str | None = None
    port_status: str = "native_rdkit"
    notes: str | None = None


@dataclass(frozen=True, slots=True)
class MedchemHit:
    rule_id: str
    rule_name: str
    severity: MedchemSeverity
    count: int
    demerit_each: int
    demerit_total: int
    atom_matches: tuple[tuple[int, ...], ...] = ()


@dataclass(frozen=True, slots=True)
class MedchemAssessment:
    smiles: str
    passed: bool
    rejected: bool
    total_demerits: int
    demerit_cutoff: int
    hits: tuple[MedchemHit, ...]
    warnings: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)
