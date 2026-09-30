"""Bounded reaction-network enumeration over a pool of available species.

Species availability is qualitative: no concentrations, rates or depletion are
inferred. Reaction steps keep all co-reactants and co-products together.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from itertools import islice, product
from typing import Any

from chem_predict.chemistry import canonicalize_smiles, mol_from_smiles
from chem_predict.chemistry.reactions import ProductSanitizationError, reaction_from_smarts, _run_reaction
from chem_predict.core import Conditions
from chem_predict.rules import RuleRegistry


@dataclass(frozen=True, slots=True)
class ReactionStep:
    depth: int
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    rule_id: str
    rule_name: str
    priority: int
    source: str | None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ReactionNetwork:
    inputs: dict[str, str]
    species_depth: dict[str, int]
    steps: list[ReactionStep]
    warnings: list[str]
    truncated: bool = False
    conditions: dict[str, Any] = field(default_factory=dict)
    limits: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_rows(self) -> list[dict[str, Any]]:
        """Flat records suitable for pandas.DataFrame or CSV export."""
        return [
            {**asdict(step), "reactants": ".".join(step.reactants),
             "products": ".".join(step.products)}
            for step in self.steps
        ]

    def to_graph_dict(self) -> dict[str, list[dict[str, Any]]]:
        """Bipartite graph; reaction nodes preserve multi-reactant semantics."""
        ids = {smiles: f"m{i}" for i, smiles in enumerate(self.species_depth)}
        nodes = [{"id": ids[s], "kind": "molecule", "smiles": s, "depth": d}
                 for s, d in self.species_depth.items()]
        edges = []
        for i, step in enumerate(self.steps):
            rid = f"r{i}"
            nodes.append({"id": rid, "kind": "reaction", **asdict(step)})
            # Repeated edges encode repeated species in distinct template slots.
            edges.extend({"source": ids[s], "target": rid, "slot": j}
                         for j, s in enumerate(step.reactants))
            edges.extend({"source": rid, "target": ids[s], "slot": j}
                         for j, s in enumerate(step.products))
        return {"nodes": nodes, "edges": edges}


def enumerate_network(
    registry: RuleRegistry,
    components: Mapping[str, str] | Sequence[str],
    conditions: Conditions,
    *,
    max_depth: int = 2,
    max_species: int = 500,
    max_steps: int = 1000,
    max_combinations_per_rule: int = 1000,
    max_products_per_rule: int = 100,
    on_error: str = "raise",
) -> ReactionNetwork:
    """Enumerate unary and multi-reactant rules, retaining provenance and cycles.

    Each species combination is executed once per rule. Products re-enter the
    available pool on the next breadth-first round. Reuse of a species in more
    than one template slot is allowed (stoichiometry is not supplied).
    """
    for name, value in (("max_depth", max_depth), ("max_species", max_species),
                        ("max_steps", max_steps),
                        ("max_combinations_per_rule", max_combinations_per_rule),
                        ("max_products_per_rule", max_products_per_rule)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if on_error not in {"raise", "collect"}:
        raise ValueError("on_error must be 'raise' or 'collect'")
    if isinstance(components, str):
        raise TypeError("Pass a list of SMILES or a name-to-SMILES mapping")
    items = components.items() if isinstance(components, Mapping) else (
        (f"component_{i + 1}", s) for i, s in enumerate(components))
    inputs = {name: canonicalize_smiles(s) for name, s in items}
    if not inputs:
        raise ValueError("At least one component is required")
    depths = dict.fromkeys(inputs.values(), 0)
    if len(depths) > max_species:
        raise ValueError("Input species exceed max_species")
    network = ReactionNetwork(inputs, depths, [], [])
    network.conditions = {**asdict(conditions), "stresses": sorted(s.value for s in conditions.stresses)}
    network.limits = dict(max_depth=max_depth, max_species=max_species, max_steps=max_steps,
                          max_combinations_per_rule=max_combinations_per_rule,
                          max_products_per_rule=max_products_per_rule)
    rules = sorted((r for r in registry.enabled() if r.when.matches(conditions)),
                   key=lambda r: (-r.priority, r.id))
    compiled = [(r, reaction_from_smarts(r.reaction_smarts)) for r in rules]
    if any(rxn.GetNumReactantTemplates() == 0 for _, rxn in compiled):
        raise ValueError("Network rules must have at least one reactant template")
    attempted: set[tuple[str, tuple[str, ...]]] = set()
    counts = dict.fromkeys((r.id for r in rules), 0)
    molecules = {s: mol_from_smiles(s) for s in depths}

    def warn(message: str) -> None:
        if message not in network.warnings:
            network.warnings.append(message)

    frontier = set(depths)
    for depth in range(1, max_depth + 1):
        new_species: set[str] = set()
        pool = tuple(depths)
        for rule, rxn in compiled:
            candidates = [tuple(s for s in pool if molecules[s].HasSubstructMatch(template))
                          for template in rxn.GetReactants()]
            if any(not group for group in candidates):
                continue
            # Include a frontier species without scanning the full Cartesian
            # product of old species. Deduplication handles overlapping slots.
            for slot, group in enumerate(candidates):
                fresh = tuple(s for s in group if s in frontier)
                if not fresh:
                    continue
                choices = list(candidates)
                choices[slot] = fresh
                for reactants in product(*choices):
                    key = (rule.id, reactants)
                    if key in attempted:
                        continue
                    if counts[rule.id] >= max_combinations_per_rule:
                        network.truncated = True
                        warn(f"combination_limit:{rule.id}")
                        break
                    attempted.add(key)
                    counts[rule.id] += 1
                    try:
                        # One extra outcome exposes RDKit's generation cap.
                        outcomes, raw_count = _run_reaction(
                            rxn, tuple(molecules[s] for s in reactants), max_products_per_rule + 1)
                    except ProductSanitizationError as exc:
                        if on_error == "raise":
                            raise
                        warn(f"product_error:{rule.id}:{reactants}:{exc}")
                        continue
                    if raw_count > max_products_per_rule:
                        network.truncated = True
                        warn(f"product_limit:{rule.id}")
                    for products in islice(outcomes, max_products_per_rule):
                        # Skip identity transforms, preserve genuine cycles.
                        if sorted(products) == sorted(reactants):
                            continue
                        additions = set(products).difference(depths)
                        if len(depths) + len(additions) > max_species:
                            network.truncated = True
                            warn("species_limit")
                            continue
                        if len(network.steps) >= max_steps:
                            network.truncated = True
                            warn("step_limit")
                            return network
                        network.steps.append(ReactionStep(
                            depth, reactants, products, rule.id, rule.name,
                            rule.priority, rule.source, {"rule_metadata": dict(rule.metadata)}))
                        for s in products:
                            if s not in depths:
                                depths[s] = depth
                                molecules[s] = mol_from_smiles(s)
                                new_species.add(s)
        if not new_species:
            break
        frontier = new_species
    return network


def compare_scenarios(
    registry: RuleRegistry,
    components: Mapping[str, str] | Sequence[str],
    scenarios: Mapping[str, Conditions],
    **network_options: Any,
) -> dict[str, ReactionNetwork]:
    """Evaluate the same mixture under named conditions with identical limits."""
    return {name: enumerate_network(registry, components, conditions, **network_options)
            for name, conditions in scenarios.items()}
