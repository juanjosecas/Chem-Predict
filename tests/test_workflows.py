import json
from math import exp

import pytest
from rdkit.Chem import Descriptors

from chem_predict import ConditionWindow, Conditions, DegradationEngine, Rule, RuleRegistry
from chem_predict.applicability import SimilarityDomain
from chem_predict.chemistry import mol_from_smiles
from chem_predict.cli import main
from chem_predict.degradation import ThermalStage, compare_scenarios, first_order_exposure
from chem_predict.nitrosamines import (
    NitrosationContext, assess_nitrosation_context, enumerate_secondary_amine_nitrosation_products,
    excipient_nitrite_mass, nitrite_limited_bound, nitrosation_speciation,
)
from chem_predict.yields import YieldPredictor


def rules():
    return RuleRegistry([
        Rule("substitution", "Substitution", "[C:1][Br:2].[O;H1;-1:3]>>[C:1][O;+0:3].[Br-:2]", source="demo"),
        Rule("oxidation", "Oxidation", "[CH2:1][OH:2]>>[CH:1]=[O:2]"),
    ])


def test_mixture_retains_coreactants_and_follows_products():
    result = DegradationEngine(rules()).predict_mixture(
        {"substrate": "CCBr", "reagent": "[OH-]"}, Conditions(), max_depth=2)
    assert result.inputs == {"substrate": "CCBr", "reagent": "[OH-]"}
    assert result.species_depth["CCO"] == 1
    assert result.species_depth["CC=O"] == 2
    assert result.steps[0].reactants == ("CCBr", "[OH-]")
    assert set(result.steps[0].products) == {"CCO", "[Br-]"}
    assert result.steps[0].source == "demo"
    graph = result.to_graph_dict()
    assert all(e["source"] in {n["id"] for n in graph["nodes"]} for e in graph["edges"])
    assert len([e for e in graph["edges"] if e["target"] == "r0"]) == 2
    json.dumps(result.to_dict())
    assert result.to_rows()[0]["reactants"] == "CCBr.[OH-]"


def test_network_limits_cycles_and_same_species_reuse():
    registry = RuleRegistry([
        Rule("reduce", "Reduce", "[C:1]=[O:2]>>[C:1][O:2]"),
        Rule("oxidize", "Oxidize", "[CH2:1][OH:2]>>[CH:1]=[O:2]"),
    ])
    engine = DegradationEngine(registry)
    result = engine.predict_mixture(["CC=O", "O=CC"], Conditions(), max_depth=8)
    assert len(result.inputs) == 2 and len(result.species_depth) == 2
    assert len(result.steps) == 2
    limited = engine.predict_mixture(["CC=O"], Conditions(), max_species=1)
    assert limited.truncated and "species_limit" in limited.warnings
    limited = engine.predict_mixture(["CC=O"], Conditions(), max_steps=1)
    assert limited.truncated and len(limited.steps) == 1
    coupling = DegradationEngine(RuleRegistry([
        Rule("couple", "Couple", "[CH3:1].[CH3:2]>>[CH2:1][CH2:2]")]))
    assert coupling.predict_mixture(["C"], Conditions()).steps == []
    # A pool species can occupy both matching slots.
    dimer = DegradationEngine(RuleRegistry([
        Rule("join", "Join", "[CH4:1].[CH4:2]>>[CH3:1][CH3:2]")]))
    assert dimer.predict_mixture(["C"], Conditions(), max_depth=1).steps[0].products == ("CC",)


def test_network_combination_and_raw_product_limits():
    engine = DegradationEngine(RuleRegistry([
        Rule("reduce", "Reduce", "[C:1]=[O:2]>>[C:1][O:2]")]))
    limited = engine.predict_mixture(["CC=O", "CCC=O"], Conditions(), max_combinations_per_rule=1)
    assert limited.truncated and "combination_limit:reduce" in limited.warnings
    # Equivalent matches collapse after generation, but still exhaust raw cap.
    symmetric = engine.predict_mixture(["O=CCC=O"], Conditions(), max_depth=1, max_products_per_rule=1)
    assert symmetric.truncated and "product_limit:reduce" in symmetric.warnings


def test_network_collects_invalid_products_but_strict_by_default():
    engine = DegradationEngine(RuleRegistry([
        Rule("invalid", "Invalid", "[CH4:1]>>[CH4:1]C")]))
    with pytest.raises(ValueError):
        engine.predict_mixture(["C"], Conditions())
    result = engine.predict_mixture(["C"], Conditions(), on_error="collect")
    assert not result.steps and result.warnings[0].startswith("product_error:invalid")


def test_scenarios_and_duration_roundtrip(tmp_path):
    registry = RuleRegistry([Rule("reduce", "Reduce", "[C:1]=[O:2]>>[C:1][O:2]",
                                  when=ConditionWindow(duration_min_h=2, duration_max_h=5))])
    path = tmp_path / "rules.json"
    registry.to_json(path)
    restored = RuleRegistry.from_json(path)
    results = compare_scenarios(restored, ["CC=O"], {
        "short": Conditions(duration_h=1), "long": Conditions(duration_h=3),
        "unknown": Conditions()})
    assert not results["short"].steps and not results["unknown"].steps
    assert results["long"].steps
    assert not DegradationEngine(rules()).predict("CCBr", Conditions())


def test_domain_preserves_reaction_direction_and_ignores_maps():
    domain = SimilarityDomain(["CCO", "c1ccccc1"], threshold=0.9)
    result = domain.assess("OCC")
    assert result.in_domain and result.nearest_similarity == 1
    assert result.neighbors[0].reference_index == 0
    assert not domain.assess("[Na+]").in_domain
    rxns = SimilarityDomain(["CC=O>>CCO"], mode="reaction", threshold=1)
    assert rxns.assess("[CH3:1][CH:2]=[O:3]>>[CH3:1][CH2:2][OH:3]").in_domain
    assert not rxns.assess("CCO>>CC=O").in_domain
    with pytest.raises(ValueError):
        SimilarityDomain([])
    with pytest.raises(ValueError):
        rxns.assess("CCO")


def test_yield_adapter_scaling_validation_and_domain():
    class Model:
        def predict(self, inputs):
            return [0.5] * len(inputs), []
    adapter = YieldPredictor.from_rxn_yields_model(Model(), model_id="test",
        output_scale="standardized", training_mean=50, training_std=20,
        domain=SimilarityDomain(["CC=O>>CCO"], mode="reaction"))
    result = adapter.predict(["CC=O>>CCO"])[0]
    assert result.yield_percent == 60 and not result.warnings
    assert YieldPredictor(lambda _: [1.2], model_id="fraction", output_scale="fraction").predict(
        ["C>>CC"])[0].yield_percent == 120
    with pytest.raises(ValueError):
        YieldPredictor(lambda _: [], model_id="missing-scale", output_scale="standardized")
    with pytest.raises(ValueError):
        YieldPredictor(lambda _: [], model_id="bad").predict(["C>>CC"])
    with pytest.raises(ValueError):
        YieldPredictor(lambda _: [float("nan")], model_id="bad").predict(["C>>CC"])


def test_process_mass_balance_speciation_and_amide_exclusion():
    assert enumerate_secondary_amine_nitrosation_products("CNC(C)=O") == []
    result = assess_nitrosation_context("CNC(C)=O", NitrosationContext(nitrite_present=True))
    assert not result.structural_precursor_present
    assert "amide_like_centers_require_separate_assessment" in result.flags
    assert assess_nitrosation_context("CCNCC", NitrosationContext(excipient_nitrite_ppm=1)).nitrosating_source_supported
    with pytest.raises(ValueError):
        NitrosationContext(nitrite_present=False, excipient_nitrite_ppm=1)
    species = nitrosation_speciation(3.2, amine_pka=10, nitrous_acid_pka=3.2)
    assert species.nitrous_acid_fraction == 0.5
    assert species.free_amine_fraction < 1e-6
    nitrite = excipient_nitrite_mass(excipient_mass_mg=100, nitrite_ppm=10)
    assert nitrite == 0.001
    budget = nitrite_limited_bound("CCNCC", amine_mass_mg=100, nitrite_mass_mg=nitrite)
    assert budget.maximum_product_mmol == budget.nitrite_mmol
    mw = Descriptors.MolWt(mol_from_smiles("CCN(CC)N=O"))
    assert budget.product_bounds_mg["CCN(CC)N=O"] == pytest.approx(budget.nitrite_mmol * mw)
    zero = nitrite_limited_bound("CCNCC", amine_mass_mg=0, nitrite_mass_mg=1)
    assert zero.maximum_product_mmol == 0


def test_kinetics_known_solution_and_temperature():
    cold = first_order_exposure([ThermalStage(25, 2)], k_reference_per_h=0.1,
                               reference_temperature_c=25, activation_energy_kj_mol=50)
    assert cold.remaining_fraction == pytest.approx(exp(-0.2))
    assert cold.converted_fraction + cold.remaining_fraction == pytest.approx(1)
    hot = first_order_exposure([ThermalStage(50, 2)], k_reference_per_h=0.1,
                              reference_temperature_c=25, activation_energy_kj_mol=50)
    assert hot.converted_fraction > cold.converted_fraction
    with pytest.raises(ValueError):
        first_order_exposure([ThermalStage(25, -1)], k_reference_per_h=0.1,
                             reference_temperature_c=25, activation_energy_kj_mol=50)


def test_cli_json_and_friendly_errors(tmp_path, capsys):
    path = tmp_path / "rules.json"
    rules().to_json(path)
    output = tmp_path / "network.json"
    assert main(["mixture", str(path), "CCBr", "[OH-]", "--output", str(output)]) == 0
    assert len(json.loads(output.read_text())["steps"]) == 2
    assert main(["nitrosation", "CCNCC", "--nitrite"]) == 0
    assert json.loads(capsys.readouterr().out)["direct_products"] == ["CCN(CC)N=O"]
    assert main(["normalize", "invalid"]) == 2
    assert "chem-predict:" in capsys.readouterr().err


@pytest.mark.parametrize("options", [{"max_depth": 0}, {"max_steps": -1}, {"max_species": True}])
def test_invalid_network_limits(options):
    with pytest.raises(ValueError):
        DegradationEngine(rules()).predict_mixture(["CCBr"], Conditions(), **options)


@pytest.mark.parametrize("kwargs", [{"ph": float("nan")}, {"temperature_c": -274}, {"duration_h": -1}])
def test_conditions_reject_invalid_process_values(kwargs):
    with pytest.raises(ValueError):
        Conditions(**kwargs)


def test_window_rejects_inverted_limits():
    with pytest.raises(ValueError):
        ConditionWindow(ph_min=8, ph_max=2)
    with pytest.raises(ValueError):
        ConditionWindow(duration_max_h=-1)


def test_synkit_rejects_unmapped_and_duplicate_maps(monkeypatch):
    from chem_predict.integrations import synkit
    # Validation must fail before dispatching to the optional backend.
    monkeypatch.setattr(synkit, "_load_synkit_io", lambda: object())
    for reaction in ("CC=O>>CCO", "[CH3:1][CH:1]=[O:3]>>[CH3:1][CH2:2][OH:3]"):
        with pytest.raises(ValueError, match="mapped atoms"):
            synkit.reaction_to_its(reaction)
