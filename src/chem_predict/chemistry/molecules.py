"""Small RDKit molecule utilities used by the rest of the package."""

from __future__ import annotations

from collections.abc import Iterable

from rdkit import Chem


class InvalidSmilesError(ValueError):
    pass


def mol_from_smiles(smiles: str) -> Chem.Mol:
    """Parse and sanitize a SMILES string.

    RDKit MolFromSmiles performs sanitization by default. Invalid input is
    converted into an explicit package exception instead of propagating a
    later None failure.
    """

    if not isinstance(smiles, str) or not smiles.strip():
        raise InvalidSmilesError("SMILES must be a non-empty string")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSmilesError(f"RDKit could not parse SMILES: {smiles!r}")
    return mol


def canonicalize_smiles(smiles: str, *, isomeric: bool = True) -> str:
    mol = mol_from_smiles(smiles)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)


def deduplicate_smiles(smiles_values: Iterable[str], *, isomeric: bool = True) -> list[str]:
    """Canonicalize and deduplicate while preserving first-seen order."""

    seen: set[str] = set()
    result: list[str] = []
    for value in smiles_values:
        canonical = canonicalize_smiles(value, isomeric=isomeric)
        if canonical not in seen:
            seen.add(canonical)
            result.append(canonical)
    return result
