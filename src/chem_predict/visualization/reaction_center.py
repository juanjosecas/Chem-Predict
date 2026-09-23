from __future__ import annotations

from collections.abc import Iterable

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from chem_predict.visualization.models import BondChange, ReactionCenterDepiction


def _split_reaction(reaction_smiles: str) -> tuple[str, str]:
    parts = reaction_smiles.split(">")
    if len(parts) != 3:
        raise ValueError(
            "Expected reaction SMILES in reactants>agents>products or reactants>>products form"
        )
    reactants, _, products = parts
    if not reactants or not products:
        raise ValueError("Reaction SMILES requires non-empty reactant and product sides")
    return reactants, products


def _mol_from_side(side: str, label: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(side)
    if mol is None:
        raise ValueError(f"Could not parse {label} side as SMILES: {side}")
    return mol


def _mapped_bonds(mol: Chem.Mol) -> dict[tuple[int, int], float]:
    bonds: dict[tuple[int, int], float] = {}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a <= 0 or b <= 0:
            continue
        key = tuple(sorted((a, b)))
        bonds[key] = float(bond.GetBondTypeAsDouble())
    return bonds


def find_reaction_center(reaction_smiles: str) -> tuple[BondChange, ...]:
    """Identify changed mapped bonds from atom-mapped reaction SMILES."""

    reactant_side, product_side = _split_reaction(reaction_smiles)
    reactants = _mol_from_side(reactant_side, "reactant")
    products = _mol_from_side(product_side, "product")

    before = _mapped_bonds(reactants)
    after = _mapped_bonds(products)

    changes: list[BondChange] = []
    for atom_pair in sorted(set(before) | set(after)):
        old = before.get(atom_pair, 0.0)
        new = after.get(atom_pair, 0.0)
        if old == new:
            continue
        changes.append(
            BondChange(
                atom_map_1=atom_pair[0],
                atom_map_2=atom_pair[1],
                reactant_order=old,
                product_order=new,
            )
        )
    return tuple(changes)


def _indices_for_maps(mol: Chem.Mol, atom_maps: set[int]) -> list[int]:
    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() in atom_maps
    ]


def _bond_indices_for_changes(
    mol: Chem.Mol,
    changed_pairs: set[tuple[int, int]],
) -> list[int]:
    indices: list[int] = []
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if tuple(sorted((a, b))) in changed_pairs:
            indices.append(bond.GetIdx())
    return indices


def _highlighted_svg(
    mol: Chem.Mol,
    *,
    atom_maps: set[int],
    changed_pairs: set[tuple[int, int]],
    width: int,
    height: int,
    legend: str,
) -> str:
    options = rdMolDraw2D.MolDrawOptions()
    options.addAtomIndices = False
    return Draw.MolToSVG(
        mol,
        size=(width, height),
        legend=legend,
        highlightAtoms=_indices_for_maps(mol, atom_maps),
        highlightBonds=_bond_indices_for_changes(mol, changed_pairs),
        drawOptions=options,
    )


def draw_reaction_center(
    reaction_smiles: str,
    *,
    width: int = 450,
    height: int = 300,
) -> ReactionCenterDepiction:
    """Render reactant/product sides with the mapped reaction center highlighted."""

    reactant_side, product_side = _split_reaction(reaction_smiles)
    reactants = _mol_from_side(reactant_side, "reactant")
    products = _mol_from_side(product_side, "product")
    changes = find_reaction_center(reaction_smiles)

    atom_maps = {
        atom_map
        for change in changes
        for atom_map in (change.atom_map_1, change.atom_map_2)
    }
    changed_pairs = {
        tuple(sorted((change.atom_map_1, change.atom_map_2)))
        for change in changes
    }

    return ReactionCenterDepiction(
        reaction_smiles=reaction_smiles,
        changed_bonds=changes,
        changed_atom_maps=tuple(sorted(atom_maps)),
        reactants_svg=_highlighted_svg(
            reactants,
            atom_maps=atom_maps,
            changed_pairs=changed_pairs,
            width=width,
            height=height,
            legend="Reactants",
        ),
        products_svg=_highlighted_svg(
            products,
            atom_maps=atom_maps,
            changed_pairs=changed_pairs,
            width=width,
            height=height,
            legend="Products",
        ),
    )
