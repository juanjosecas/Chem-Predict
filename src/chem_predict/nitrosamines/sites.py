from __future__ import annotations

from rdkit import Chem

from chem_predict.chemistry import mol_from_smiles
from chem_predict.nitrosamines.models import AmineKind, NitrosamineSite, NitrosatableCenter


_NITROSAMINE_QUERY = Chem.MolFromSmarts("[#7:1]([#6:2])([#6:3])-[#7:4]=[#8:5]")
if _NITROSAMINE_QUERY is None:  # pragma: no cover
    raise RuntimeError("Could not compile internal nitrosamine SMARTS")


def _carbon_double_bonded_to_heteroatom(mol: Chem.Mol, carbon_idx: int) -> bool:
    carbon = mol.GetAtomWithIdx(carbon_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() != 6:
            return True
    return False


def _amide_like_nitrogen(mol: Chem.Mol, nitrogen_idx: int) -> bool:
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    for neighbor in nitrogen.GetNeighbors():
        if neighbor.GetAtomicNum() != 6:
            continue
        if _carbon_double_bonded_to_heteroatom(mol, neighbor.GetIdx()):
            return True
    return False


def _is_methyl_carbon(atom: Chem.Atom) -> bool:
    return (
        atom.GetAtomicNum() == 6
        and atom.GetHybridization() == Chem.HybridizationType.SP3
        and atom.GetTotalNumHs() == 3
    )


def find_nitrosamine_sites(smiles: str) -> list[NitrosamineSite]:
    """Locate R2N-N=O sites and annotate FDA-CPCA scope exclusions.

    The detector intentionally finds broader N-nitrosamine-like structures than
    those eligible for the CPCA. Applicability is reported separately.
    """

    mol = mol_from_smiles(smiles)
    sites: list[NitrosamineSite] = []
    seen: set[tuple[int, int, int, int, int]] = set()

    for match in mol.GetSubstructMatches(_NITROSAMINE_QUERY, uniquify=True):
        amine_n, alpha_1, alpha_2, nitroso_n, oxygen = match
        key = (amine_n, *sorted((alpha_1, alpha_2)), nitroso_n, oxygen)
        if key in seen:
            continue
        seen.add(key)

        alpha = tuple(sorted((alpha_1, alpha_2)))
        alpha_h = tuple(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in alpha)
        exclusions: list[str] = []

        if mol.GetAtomWithIdx(amine_n).GetIsAromatic():
            exclusions.append("n_nitroso_group_in_aromatic_ring")
        if any(_carbon_double_bonded_to_heteroatom(mol, idx) for idx in alpha):
            exclusions.append("alpha_carbon_directly_double_bonded_to_heteroatom")

        sites.append(
            NitrosamineSite(
                nitrosamine_n_index=amine_n,
                nitroso_n_index=nitroso_n,
                oxygen_index=oxygen,
                alpha_carbon_indices=alpha,
                alpha_hydrogens=alpha_h,
                cpca_applicable=not exclusions,
                exclusion_reasons=tuple(exclusions),
            )
        )

    return sites


def find_nitrosatable_centers(smiles: str) -> list[NitrosatableCenter]:
    """Find amine centers relevant to nitrosamine formation screening.

    This is a structural screen, not a formation probability model. It keeps
    secondary, tertiary, dimethyl-tertiary, and quaternary centers distinct so
    downstream process/context modules can treat them differently.
    """

    mol = mol_from_smiles(smiles)
    centers: list[NitrosatableCenter] = []

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetIsAromatic():
            continue

        carbon_neighbors = tuple(
            neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6
        )
        total_h = int(atom.GetTotalNumHs())
        formal_charge = atom.GetFormalCharge()
        methyl_count = sum(_is_methyl_carbon(mol.GetAtomWithIdx(i)) for i in carbon_neighbors)

        kind: AmineKind | None = None
        if formal_charge == 0 and len(carbon_neighbors) == 2 and total_h >= 1:
            kind = AmineKind.SECONDARY
        elif formal_charge == 0 and len(carbon_neighbors) == 3 and total_h == 0:
            kind = (
                AmineKind.DIMETHYL_TERTIARY
                if methyl_count >= 2
                else AmineKind.TERTIARY
            )
        elif formal_charge > 0 and len(carbon_neighbors) == 4:
            kind = AmineKind.QUATERNARY

        if kind is None:
            continue

        centers.append(
            NitrosatableCenter(
                nitrogen_index=atom.GetIdx(),
                kind=kind,
                carbon_neighbor_indices=carbon_neighbors,
                methyl_substituent_count=methyl_count,
                total_hydrogens=total_h,
                aromatic=False,
                amide_like=_amide_like_nitrogen(mol, atom.GetIdx()),
            )
        )

    return centers
