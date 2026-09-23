from __future__ import annotations

from pathlib import Path
from typing import Iterable

from rdkit.Chem import Draw, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

from chem_predict.chemistry import mol_from_smiles


def molecule_svg(
    smiles: str,
    *,
    width: int = 400,
    height: int = 300,
    legend: str = "",
    atom_indices: bool = False,
    highlight_atoms: Iterable[int] = (),
    highlight_bonds: Iterable[int] = (),
) -> str:
    """Render a SMILES string to SVG using RDKit 2026.03 MolDraw2D APIs."""

    mol = mol_from_smiles(smiles)
    options = rdMolDraw2D.MolDrawOptions()
    options.addAtomIndices = atom_indices
    return Draw.MolToSVG(
        mol,
        size=(width, height),
        legend=legend,
        highlightAtoms=list(highlight_atoms),
        highlightBonds=list(highlight_bonds),
        drawOptions=options,
    )


def reaction_svg(
    reaction: str,
    *,
    width: int = 900,
    height: int = 300,
    input_format: str = "smiles",
    highlight_by_reactant: bool = False,
) -> str:
    """Render reaction SMILES or reaction SMARTS to SVG."""

    if input_format not in {"smiles", "smarts"}:
        raise ValueError("input_format must be 'smiles' or 'smarts'")

    rxn = rdChemReactions.ReactionFromSmarts(
        reaction,
        useSmiles=input_format == "smiles",
    )
    if rxn is None:
        raise ValueError(f"Could not parse reaction {input_format}: {reaction}")

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawReaction(rxn, highlightByReactant=highlight_by_reactant)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def molecule_png_bytes(
    smiles: str,
    *,
    width: int = 400,
    height: int = 300,
    legend: str = "",
) -> bytes:
    """Render a molecule to PNG bytes with the RDKit Cairo drawer."""

    mol = mol_from_smiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    drawer.DrawMolecule(mol, legend=legend)
    drawer.FinishDrawing()
    return bytes(drawer.GetDrawingText())


def save_svg(svg: str, path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(svg, encoding="utf-8")
    return destination


def save_png(data: bytes, path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)
    return destination
