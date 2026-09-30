"""Transparent nearest-neighbor structural coverage, independent of ML backends.

This baseline is not an implementation of SynAD or a calibrated confidence
model. Thresholds must be assessed against held-out errors for each endpoint.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from math import isfinite

from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

from chem_predict.chemistry import mol_from_smiles


def split_reaction_smiles(reaction: str) -> tuple[str, str, str]:
    """Validate molecular (not SMARTS) reactants>agents>products."""
    if not isinstance(reaction, str) or len(parts := reaction.split(">")) != 3:
        raise ValueError("Use reactants>agents>products or reactants>>products")
    if not parts[0] or not parts[2]:
        raise ValueError("Reactants and products must be non-empty")
    for side in parts:
        if side:
            mol_from_smiles(side)
    return tuple(parts)


@dataclass(frozen=True, slots=True)
class Neighbor:
    reference_index: int
    reference: str
    similarity: float


@dataclass(frozen=True, slots=True)
class DomainAssessment:
    query: str
    in_domain: bool
    threshold: float
    nearest_similarity: float
    neighbors: tuple[Neighbor, ...]
    representation: str


class SimilarityDomain:
    """Morgan/Tanimoto coverage for molecules or role-separated reaction sides.

    Reaction fingerprints concatenate reactants, agents and products. Atom map
    numbers are stripped; salts and stereochemistry are preserved. No process
    conditions are encoded, so chemical similarity alone cannot validate yield.
    """

    def __init__(self, references: Sequence[str], *, mode: str = "molecule",
                 threshold: float = 0.5, radius: int = 2, n_bits: int = 2048):
        if isinstance(references, str) or not references:
            raise ValueError("Provide a non-empty sequence of reference SMILES")
        if mode not in {"molecule", "reaction"}:
            raise ValueError("mode must be 'molecule' or 'reaction'")
        if not isfinite(threshold) or not 0 <= threshold <= 1:
            raise ValueError("threshold must be finite and between 0 and 1")
        if isinstance(radius, bool) or not isinstance(radius, int) or radius < 0:
            raise ValueError("radius must be a non-negative integer")
        if isinstance(n_bits, bool) or not isinstance(n_bits, int) or n_bits < 1:
            raise ValueError("n_bits must be a positive integer")
        self.mode, self.threshold, self.n_bits = mode, threshold, n_bits
        self.references = tuple(references)
        self.generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=n_bits, includeChirality=True)
        self.representation = f"Morgan(radius={radius},bits={n_bits},chirality=True);{mode}"
        self.fingerprints = tuple(self._fingerprint(s) for s in references)

    def _fingerprint(self, text: str):
        sides = split_reaction_smiles(text) if self.mode == "reaction" else (text,)
        result = DataStructs.ExplicitBitVect(self.n_bits * len(sides))
        for i, side in enumerate(sides):
            if not side:
                continue
            mol = mol_from_smiles(side)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            for bit in self.generator.GetFingerprint(mol).GetOnBits():
                result.SetBit(i * self.n_bits + bit)
        return result

    def assess(self, query: str, *, k: int = 5) -> DomainAssessment:
        if isinstance(k, bool) or not isinstance(k, int) or k < 1:
            raise ValueError("k must be a positive integer")
        similarities = DataStructs.BulkTanimotoSimilarity(
            self._fingerprint(query), self.fingerprints)
        order = sorted(range(len(similarities)), key=lambda i: (-similarities[i], i))[:k]
        neighbors = tuple(Neighbor(i, self.references[i], float(similarities[i])) for i in order)
        nearest = neighbors[0].similarity
        return DomainAssessment(query, nearest >= self.threshold, self.threshold,
                                nearest, neighbors, self.representation)

    def assess_batch(self, queries: Sequence[str], *, k: int = 5) -> list[DomainAssessment]:
        if isinstance(queries, str):
            raise TypeError("Pass a sequence of query SMILES")
        return [self.assess(query, k=k) for query in queries]


__all__ = ["Neighbor", "DomainAssessment", "SimilarityDomain", "split_reaction_smiles"]
