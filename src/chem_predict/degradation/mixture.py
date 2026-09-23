"""Bounded one-step rule enumeration over a supplied mixture.

No kinetic, concentration, or product probability model is implied.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Any

from rdkit import Chem

from chem_predict.chemistry import apply_reaction, canonicalize_smiles
from chem_predict.chemistry.reactions import reaction_from_smarts
from chem_predict.core import Conditions
from chem_predict.rules import RuleRegistry


@dataclass(frozen=True, slots=True)
class MixtureComponent:
    smiles: str
    role: str = "unspecified"
    label: str | None = None


@dataclass(frozen=True, slots=True)
class MixturePrediction:
    reactants: tuple[str, ...]
    reactant_indices: tuple[int, ...]
    reactant_roles: tuple[str, ...]
    products: tuple[str, ...]
    rule_id: str
    rule_name: str
    priority: int
    source: str | None
    metadata: dict[str, Any] = field(default_factory=dict)


class MixtureDegradationEngine:
    def __init__(self, registry: RuleRegistry) -> None:
        self.registry = registry

    def predict(
        self,
        components: list[MixtureComponent],
        conditions: Conditions,
        *,
        max_combinations_per_rule: int = 10000,
        max_products_per_combination: int = 1000,
        allow_self_reaction: bool = False,
    ) -> list[MixturePrediction]:
        """Enumerate direct products; each component can fill one reactant slot.

        Identical component entries represent separate available species. Reusing a
        single entry across slots requires ``allow_self_reaction=True``. A rule
        exceeding the combination cap raises rather than returning partial results.
        """
        if max_combinations_per_rule < 1:
            raise ValueError("max_combinations_per_rule must be >= 1")
        if max_products_per_combination < 1:
            raise ValueError("max_products_per_combination must be >= 1")
        smiles = [canonicalize_smiles(component.smiles) for component in components]
        mols = [Chem.MolFromSmiles(value) for value in smiles]
        predictions: list[MixturePrediction] = []
        rules = sorted(
            (rule for rule in self.registry.enabled() if rule.when.matches(conditions)),
            key=lambda rule: (-rule.priority, rule.id),
        )
        for rule in rules:
            reaction = reaction_from_smarts(rule.reaction_smarts)
            slots = [
                [index for index, mol in enumerate(mols) if mol.HasSubstructMatch(template)]
                for template in reaction.GetReactants()
            ]
            if not slots or any(not candidates for candidates in slots):
                continue
            count = 0
            for indices in product(*slots):
                if not allow_self_reaction and len(set(indices)) != len(indices):
                    continue
                count += 1
                if count > max_combinations_per_rule:
                    raise ValueError(f"Rule {rule.id!r} exceeds max_combinations_per_rule")
                reactants = tuple(smiles[index] for index in indices)
                for products in apply_reaction(
                    rule.reaction_smarts, reactants, max_products=max_products_per_combination
                ):
                    predictions.append(
                        MixturePrediction(
                            reactants=reactants,
                            reactant_indices=indices,
                            reactant_roles=tuple(components[index].role for index in indices),
                            products=products,
                            rule_id=rule.id,
                            rule_name=rule.name,
                            priority=rule.priority,
                            source=rule.source,
                            metadata={"rule_metadata": dict(rule.metadata)},
                        )
                    )
        return predictions
