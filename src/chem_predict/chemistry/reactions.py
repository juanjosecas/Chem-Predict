"""RDKit reaction SMARTS execution with explicit validation boundaries."""

from __future__ import annotations

from collections.abc import Sequence

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from chem_predict.chemistry.molecules import mol_from_smiles


class InvalidReactionError(ValueError):
    pass


class ReactionArityError(ValueError):
    pass


class ProductSanitizationError(ValueError):
    pass


def reaction_from_smarts(reaction_smarts: str) -> rdChemReactions.ChemicalReaction:
    if not isinstance(reaction_smarts, str) or ">>" not in reaction_smarts:
        raise InvalidReactionError("reaction_smarts must be a reaction SMARTS containing '>>'")

    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    except Exception as exc:
        raise InvalidReactionError(f"Could not parse reaction SMARTS: {reaction_smarts!r}") from exc

    if reaction is None:
        raise InvalidReactionError(f"Could not parse reaction SMARTS: {reaction_smarts!r}")
    return reaction


def apply_reaction(
    reaction_smarts: str,
    reactant_smiles: Sequence[str],
    *,
    max_products: int = 1000,
    sanitize_products: bool = True,
) -> list[tuple[str, ...]]:
    """Apply a reaction SMARTS to an ordered sequence of reactants.

    Results are canonical isomeric SMILES tuples. Duplicate product tuples are
    removed while preserving RDKit's generation order.
    """

    if max_products < 0:
        raise ValueError("max_products must be >= 0; RDKit uses 0 for unlimited generation")

    reaction = reaction_from_smarts(reaction_smarts)
    expected = reaction.GetNumReactantTemplates()
    if len(reactant_smiles) != expected:
        raise ReactionArityError(
            f"Reaction expects {expected} reactant(s), received {len(reactant_smiles)}"
        )

    reactants = tuple(mol_from_smiles(smiles) for smiles in reactant_smiles)
    raw_outcomes = reaction.RunReactants(reactants, maxProducts=max_products)

    seen: set[tuple[str, ...]] = set()
    outcomes: list[tuple[str, ...]] = []

    for raw_products in raw_outcomes:
        product_smiles: list[str] = []
        for product in raw_products:
            if sanitize_products:
                try:
                    Chem.SanitizeMol(product)
                except Exception as exc:
                    raise ProductSanitizationError(
                        f"Generated product failed RDKit sanitization for rule {reaction_smarts!r}"
                    ) from exc
            product_smiles.append(
                Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
            )

        outcome = tuple(product_smiles)
        if outcome not in seen:
            seen.add(outcome)
            outcomes.append(outcome)

    return outcomes
