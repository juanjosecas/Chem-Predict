"""Rule registry and a dependency-free JSON interchange format."""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from chem_predict.core import ConditionWindow, Rule, StressType


class DuplicateRuleError(ValueError):
    pass


class RuleRegistry:
    def __init__(self, rules: Iterable[Rule] | None = None) -> None:
        self._rules: dict[str, Rule] = {}
        if rules is not None:
            self.extend(rules)

    def __len__(self) -> int:
        return len(self._rules)

    def __iter__(self) -> Iterator[Rule]:
        return iter(self._rules.values())

    def __contains__(self, rule_id: str) -> bool:
        return rule_id in self._rules

    def register(self, rule: Rule, *, replace: bool = False) -> None:
        if rule.id in self._rules and not replace:
            raise DuplicateRuleError(f"Rule id already registered: {rule.id}")
        self._rules[rule.id] = rule

    def extend(self, rules: Iterable[Rule], *, replace: bool = False) -> None:
        for rule in rules:
            self.register(rule, replace=replace)

    def get(self, rule_id: str) -> Rule:
        return self._rules[rule_id]

    def enabled(self) -> Iterator[Rule]:
        return (rule for rule in self._rules.values() if rule.enabled)

    @classmethod
    def from_json(cls, path: str | Path) -> "RuleRegistry":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, list):
            raise ValueError("Rule JSON root must be a list")
        return cls(_rule_from_mapping(item) for item in payload)

    def to_json(self, path: str | Path, *, indent: int = 2) -> None:
        payload = [_rule_to_mapping(rule) for rule in self._rules.values()]
        Path(path).write_text(
            json.dumps(payload, indent=indent, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )


def _window_from_mapping(data: dict[str, Any] | None) -> ConditionWindow:
    data = dict(data or {})
    stresses = frozenset(StressType(value) for value in data.pop("stresses", []))
    return ConditionWindow(stresses=stresses, **data)


def _rule_from_mapping(data: dict[str, Any]) -> Rule:
    if not isinstance(data, dict):
        raise TypeError("Each rule must be a JSON object")

    values = dict(data)
    values["when"] = _window_from_mapping(values.get("when"))
    values["tags"] = frozenset(values.get("tags", []))
    return Rule(**values)


def _rule_to_mapping(rule: Rule) -> dict[str, Any]:
    return {
        "id": rule.id,
        "name": rule.name,
        "reaction_smarts": rule.reaction_smarts,
        "when": {
            "stresses": sorted(stress.value for stress in rule.when.stresses),
            "ph_min": rule.when.ph_min,
            "ph_max": rule.when.ph_max,
            "temperature_min_c": rule.when.temperature_min_c,
            "temperature_max_c": rule.when.temperature_max_c,
            "requires_oxygen": rule.when.requires_oxygen,
            "requires_light": rule.when.requires_light,
        },
        "priority": rule.priority,
        "enabled": rule.enabled,
        "source": rule.source,
        "description": rule.description,
        "tags": sorted(rule.tags),
        "metadata": rule.metadata,
    }
