import json

import pytest

from chem_predict.core import ConditionWindow, Conditions, Rule, StressType
from chem_predict.rules import DuplicateRuleError, RuleRegistry


def test_condition_window_requires_known_constrained_values() -> None:
    window = ConditionWindow(stresses=frozenset({StressType.ACID}), ph_max=3.0)
    assert window.matches(Conditions(stresses=frozenset({StressType.ACID}), ph=2.0))
    assert not window.matches(Conditions(stresses=frozenset({StressType.ACID})))
    assert not window.matches(Conditions(stresses=frozenset({StressType.BASE}), ph=2.0))


def test_registry_rejects_duplicate_ids() -> None:
    rule = Rule(id="x", name="x", reaction_smarts="[C:1]>>[C:1]")
    registry = RuleRegistry([rule])
    with pytest.raises(DuplicateRuleError):
        registry.register(rule)


def test_registry_loads_json(tmp_path) -> None:
    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            [
                {
                    "id": "demo",
                    "name": "demo",
                    "reaction_smarts": "[C:1]=[O:2]>>[C:1][O:2]",
                    "when": {"stresses": ["reduction"]},
                }
            ]
        ),
        encoding="utf-8",
    )

    registry = RuleRegistry.from_json(path)
    assert registry.get("demo").when.stresses == frozenset({StressType.REDUCTION})
