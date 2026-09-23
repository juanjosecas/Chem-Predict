import pytest

from chem_predict.core import ConditionWindow, Conditions, Rule, StressType
from chem_predict.degradation import MixtureComponent, MixtureDegradationEngine
from chem_predict.rules import RuleRegistry


RULE = Rule(
    id="demo_substitution",
    name="Demonstration substitution",
    reaction_smarts="[CH3:1][Br:2].[OH-:3]>>[CH3:1][O:3].[Br-:2]",
    when=ConditionWindow(stresses=frozenset({StressType.BASE})),
    source="Demonstration only",
)


def test_binary_rule_matches_mixture_and_conditions() -> None:
    engine = MixtureDegradationEngine(RuleRegistry([RULE]))
    components = [
        MixtureComponent("CBr", role="substrate"),
        MixtureComponent("[OH-]", role="reagent"),
        MixtureComponent("CCO", role="solvent"),
    ]
    assert engine.predict(components, Conditions()) == []
    results = engine.predict(components, Conditions(stresses=frozenset({StressType.BASE})))
    assert len(results) == 1
    assert results[0].reactant_indices == (0, 1)
    assert results[0].reactant_roles == ("substrate", "reagent")
    assert set(results[0].products) == {"C[O-]", "[Br-]"}
    assert results[0].source == "Demonstration only"


def test_same_entry_cannot_fill_two_slots_by_default() -> None:
    rule = Rule("dimer", "dimer", "[C:1].[C:2]>>[C:1].[C:2]")
    engine = MixtureDegradationEngine(RuleRegistry([rule]))
    assert engine.predict([MixtureComponent("C")], Conditions()) == []
    assert len(engine.predict([MixtureComponent("C")], Conditions(), allow_self_reaction=True)) == 1


def test_combination_cap_fails_explicitly() -> None:
    engine = MixtureDegradationEngine(RuleRegistry([RULE]))
    components = [MixtureComponent("CBr"), MixtureComponent("[OH-]"), MixtureComponent("[OH-]")]
    with pytest.raises(ValueError, match="exceeds max_combinations"):
        engine.predict(components, Conditions(stresses=frozenset({StressType.BASE})), max_combinations_per_rule=1)
