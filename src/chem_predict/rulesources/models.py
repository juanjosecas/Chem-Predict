from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any


class SourceKind(StrEnum):
    REACTION_RULES = "reaction_rules"
    REACTIONS = "reactions"
    TRANSFORMATIONS = "transformations"
    STRUCTURAL_ALERTS = "structural_alerts"


class AccessMode(StrEnum):
    API = "api"
    DOWNLOAD = "download"
    EXPORT = "export"
    EXTERNAL = "external"


@dataclass(frozen=True, slots=True)
class RuleSourceDescriptor:
    id: str
    name: str
    kind: SourceKind
    primary_url: str
    access_mode: AccessMode
    license: str
    redistributable_in_repo: bool
    version_checked: str
    notes: str = ""


@dataclass(frozen=True, slots=True)
class ReactionRuleRecord:
    id: str
    reaction_smarts: str
    source: str
    source_record_ids: tuple[str, ...] = ()
    ec_numbers: tuple[str, ...] = ()
    radius: int | None = None
    score: float | None = None
    valid: bool | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ReactionRecord:
    id: str
    source: str
    reaction_smiles: str | None = None
    parent_smiles: str | None = None
    product_smiles: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)


class SourceFetchError(RuntimeError):
    pass
