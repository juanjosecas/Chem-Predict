"""Process calculations with supplied constants; no inferred formation yields."""

from dataclasses import dataclass
from math import isfinite

from rdkit import Chem
from rdkit.Chem import Descriptors

from chem_predict.chemistry import mol_from_smiles
from chem_predict.nitrosamines.formation import enumerate_secondary_amine_nitrosation_products


def _finite(value: float, name: str, *, nonnegative: bool = False) -> None:
    if not isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f"{name} must be finite" + (" and non-negative" if nonnegative else ""))


def _fraction(exponent: float) -> float:
    if exponent >= 0:
        ratio = 10.0 ** (-exponent)
        return ratio / (1 + ratio)
    return 1 / (1 + 10.0 ** exponent)


@dataclass(frozen=True, slots=True)
class NitrosationSpeciation:
    ph: float
    nitrous_acid_fraction: float
    free_amine_fraction: float
    amine_pka: float
    nitrous_acid_pka: float
    assumptions: tuple[str, ...]


def nitrosation_speciation(ph: float, *, amine_pka: float,
                          nitrous_acid_pka: float) -> NitrosationSpeciation:
    """Aqueous single-site equilibrium fractions, not rates or formation yields.

    Supply pKa values appropriate to the medium and temperature. Activity
    corrections and coupled protonation equilibria are not modeled.
    """
    for name, value in (("ph", ph), ("amine_pka", amine_pka),
                        ("nitrous_acid_pka", nitrous_acid_pka)):
        _finite(value, name)
    return NitrosationSpeciation(
        ph, _fraction(ph - nitrous_acid_pka), _fraction(amine_pka - ph),
        amine_pka, nitrous_acid_pka,
        ("aqueous_equilibrium", "single_amine_site", "activity_effects_ignored"))


@dataclass(frozen=True, slots=True)
class NitriteBudget:
    nitrite_mass_mg: float
    nitrite_mmol: float
    amine_mmol: float
    maximum_product_mmol: float
    product_bounds_mg: dict[str, float]
    assumptions: tuple[str, ...]


def nitrite_limited_bound(smiles: str, *, amine_mass_mg: float,
                         nitrite_mass_mg: float) -> NitriteBudget:
    """Stoichiometric upper bound for direct mononitrosation, not expected yield.

    Nitrite mass is NO2- ion equivalent, not NaNO2 salt mass. Each positional
    product is an alternative allocation of the same budget; do not sum bounds.
    Other nitrosating sources, competing reactions and impurities are excluded.
    """
    _finite(amine_mass_mg, "amine_mass_mg", nonnegative=True)
    _finite(nitrite_mass_mg, "nitrite_mass_mg", nonnegative=True)
    mol = mol_from_smiles(smiles)
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Use a single precursor molecule and its equivalent mass")
    products = enumerate_secondary_amine_nitrosation_products(smiles)
    if not products:
        raise ValueError("No supported non-amide secondary amine for direct nitrosation")
    nitrite_mmol = nitrite_mass_mg / Descriptors.MolWt(mol_from_smiles("[O-]N=O"))
    amine_mmol = amine_mass_mg / Descriptors.MolWt(mol)
    maximum = min(nitrite_mmol, amine_mmol)
    bounds = {s: maximum * Descriptors.MolWt(mol_from_smiles(s)) for s in products}
    return NitriteBudget(nitrite_mass_mg, nitrite_mmol, amine_mmol, maximum, bounds,
                         ("one_nitrite_per_mononitroso_product", "complete_conversion_upper_bound",
                          "positional_products_share_budget", "no_other_nitrosating_sources"))


def excipient_nitrite_mass(*, excipient_mass_mg: float, nitrite_ppm: float) -> float:
    """NO2- ion mass (mg) from mass/mass ppm."""
    _finite(excipient_mass_mg, "excipient_mass_mg", nonnegative=True)
    _finite(nitrite_ppm, "nitrite_ppm", nonnegative=True)
    return excipient_mass_mg * nitrite_ppm / 1_000_000
