from chem_predict.core import ConditionWindow, Conditions, Rule, StressType
from chem_predict.degradation import DegradationEngine
from chem_predict.rules import RuleRegistry


def test_engine_filters_rules_by_conditions() -> None:
    reduction = Rule(
        id="reduction",
        name="reduction",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        when=ConditionWindow(stresses=frozenset({StressType.REDUCTION})),
    )
    acid_only = Rule(
        id="acid",
        name="acid",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        when=ConditionWindow(stresses=frozenset({StressType.ACID})),
    )

    engine = DegradationEngine(RuleRegistry([reduction, acid_only]))
    predictions = engine.predict(
        "CC=O",
        Conditions(stresses=frozenset({StressType.REDUCTION})),
    )

    assert [(item.rule_id, item.products) for item in predictions] == [
        ("reduction", ("CCO",))
    ]
