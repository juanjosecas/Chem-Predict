from __future__ import annotations

from collections.abc import Iterable

from chem_predict.core import ConditionWindow, Rule
from chem_predict.rules import RuleRegistry
from chem_predict.rulesources.models import ReactionRuleRecord


def to_core_rule(
    record: ReactionRuleRecord,
    *,
    priority: int = 0,
    tags: Iterable[str] = (),
) -> Rule:
    """Convert a sourced reaction template into Chem-Predict's generic Rule."""

    metadata = {
        **record.metadata,
        "source_record_ids": record.source_record_ids,
        "ec_numbers": record.ec_numbers,
        "radius": record.radius,
        "source_score": record.score,
        "source_valid": record.valid,
    }
    return Rule(
        id=record.id,
        name=record.id,
        reaction_smarts=record.reaction_smarts,
        when=ConditionWindow(),
        priority=priority,
        source=record.source,
        tags=frozenset(tags),
        metadata=metadata,
    )


def to_rule_registry(
    records: Iterable[ReactionRuleRecord],
    *,
    priority: int = 0,
    tags: Iterable[str] = (),
) -> RuleRegistry:
    return RuleRegistry(
        to_core_rule(record, priority=priority, tags=tags)
        for record in records
    )
