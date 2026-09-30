"""Run with python examples/modular_workflows.py; no network or ML model needed."""

from chem_predict import Conditions, DegradationEngine, Rule, RuleRegistry
from chem_predict.applicability import SimilarityDomain
from chem_predict.degradation import ThermalStage, first_order_exposure
from chem_predict.nitrosamines import excipient_nitrite_mass, nitrite_limited_bound


registry = RuleRegistry([
    Rule("demo_substitution", "Demo substitution",
         "[C:1][Br:2].[O;H1;-1:3]>>[C:1][O;+0:3].[Br-:2]", source="API demonstration"),
    Rule("demo_oxidation", "Demo oxidation",
         "[CH2:1][OH:2]>>[CH:1]=[O:2]", source="API demonstration"),
])
network = DegradationEngine(registry).predict_mixture(
    {"substrate": "CCBr", "reagent": "[OH-]"}, Conditions(), max_depth=2)
for row in network.to_rows():
    print(row["depth"], row["rule_id"], row["reactants"], "->", row["products"])
print("Warnings:", network.warnings)

domain = SimilarityDomain(["CC=O>>CCO", "CCC=O>>CCCO"], mode="reaction")
print("Coverage:", domain.assess("CCCC=O>>CCCCO"))

nitrite = excipient_nitrite_mass(excipient_mass_mg=100, nitrite_ppm=10)
print("Nitrite bound:", nitrite_limited_bound("CCNCC", amine_mass_mg=100, nitrite_mass_mg=nitrite))

# Illustrative kinetic parameters; replace with measured values.
print("First-order exposure:", first_order_exposure(
    [ThermalStage(25, 24), ThermalStage(40, 2)], k_reference_per_h=0.001,
    reference_temperature_c=25, activation_energy_kj_mol=50))
