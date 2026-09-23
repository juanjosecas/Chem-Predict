from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from chem_predict.chemistry import mol_from_smiles
from chem_predict.nitrosamines.models import NitrosationAssessment, NitrosationContext
from chem_predict.nitrosamines.sites import find_nitrosatable_centers


_SECONDARY_AMINE_NITROSATION = rdChemReactions.ReactionFromSmarts(
    "[N;H1;X3:1]([#6:2])[#6:3]>>[N:1]([#6:2])([#6:3])N=O"
)
if _SECONDARY_AMINE_NITROSATION is None:  # pragma: no cover
    raise RuntimeError("Could not compile internal nitrosation reaction SMARTS")


def assess_nitrosation_context(smiles: str, context: NitrosationContext) -> NitrosationAssessment:
    """Combine structural precursor screening with process/formulation context.

    No arbitrary probability score is generated. The result reports evidence
    flags that can later feed a validated quantitative or expert model.
    """

    centers = tuple(find_nitrosatable_centers(smiles))
    structural = bool(centers)

    if context.nitrosating_agent_present is True or context.nitrite_present is True:
        nitrosating_source = True
    elif context.nitrosating_agent_present is False and context.nitrite_present is False:
        nitrosating_source = False
    else:
        nitrosating_source = None

    acidic = None if context.ph is None else context.ph < 7.0
    flags: list[str] = []

    if structural:
        flags.append("amine_precursor_present")
    else:
        flags.append("no_screened_amine_precursor_found")

    if nitrosating_source is True:
        flags.append("nitrosating_source_present")
    elif nitrosating_source is None:
        flags.append("nitrosating_source_unknown")

    if acidic is True:
        flags.append("acidic_conditions")
    if context.storage is True:
        flags.append("storage_context")
    if context.excipient_nitrite_ppm is not None:
        flags.append("excipient_nitrite_level_provided")

    if structural and nitrosating_source is True:
        flags.append("fda_root_cause_combination_present")

    return NitrosationAssessment(
        centers=centers,
        structural_precursor_present=structural,
        nitrosating_source_supported=nitrosating_source,
        acidic_conditions=acidic,
        flags=tuple(flags),
    )


def enumerate_secondary_amine_nitrosation_products(smiles: str) -> list[str]:
    """Enumerate direct N-nitrosation products for secondary amines only.

    Tertiary/quaternary precursor pathways are intentionally not enumerated by
    this function because their formation may involve impurities, degradation,
    or dealkylation rather than a direct one-step structural transformation.
    """

    mol = mol_from_smiles(smiles)
    outcomes = _SECONDARY_AMINE_NITROSATION.RunReactants((mol,))
    products: list[str] = []
    seen: set[str] = set()

    for outcome in outcomes:
        product = Chem.Mol(outcome[0])
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
        if product_smiles not in seen:
            seen.add(product_smiles)
            products.append(product_smiles)

    return products
