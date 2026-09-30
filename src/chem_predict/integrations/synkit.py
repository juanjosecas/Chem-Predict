from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from chem_predict.applicability import split_reaction_smiles
from chem_predict.chemistry import mol_from_smiles


SUPPORTED_SYNKIT_MIN = "1.6.2"


class SynKitUnavailableError(ImportError):
    pass


@dataclass(frozen=True, slots=True)
class BondChange:
    source: int
    target: int
    reactant_order: float
    product_order: float


@dataclass(frozen=True, slots=True)
class SynKitITSResult:
    graph: Any
    changed_bonds: tuple[BondChange, ...]
    synkit_version: str

    def changes_to_dict(self) -> dict[str, Any]:
        """JSON-compatible change report without serializing graph internals."""
        return {"synkit_version": self.synkit_version,
                "changed_bonds": [asdict(change) for change in self.changed_bonds]}


def _load_synkit_io() -> Any:
    try:
        return import_module("synkit.IO")
    except ImportError as exc:
        raise SynKitUnavailableError(
            'SynKit is optional. Install with: python -m pip install "chem-predict[synkit]"'
        ) from exc


def installed_synkit_version() -> str | None:
    try:
        return version("synkit")
    except PackageNotFoundError:
        return None


def reaction_to_its(
    reaction_smiles: str,
    *,
    core: bool = False,
    explicit_hydrogen: bool = False,
    include_context_edges: bool = True,
) -> SynKitITSResult:
    """Convert atom-mapped reaction SMILES using SynKit 1.6.x tuple ITS API."""

    io = _load_synkit_io()
    sides = split_reaction_smiles(reaction_smiles)
    for side in (sides[0], sides[2]):
        maps = [atom.GetAtomMapNum() for atom in mol_from_smiles(side).GetAtoms()]
        if any(number <= 0 for number in maps) or len(set(maps)) != len(maps):
            raise ValueError("SynKit requires fully mapped atoms with unique positive maps on each side")
    graph = io.rsmi_to_its(
        reaction_smiles,
        core=core,
        explicit_hydrogen=explicit_hydrogen,
        format="tuple",
        include_context_edges=include_context_edges,
    )

    changes: list[BondChange] = []
    for source, target, data in graph.edges(data=True):
        order = data.get("order")
        if not isinstance(order, (tuple, list)) or len(order) != 2:
            continue
        r_order = float(order[0])
        p_order = float(order[1])
        if r_order != p_order:
            changes.append(
                BondChange(
                    source=int(source),
                    target=int(target),
                    reactant_order=r_order,
                    product_order=p_order,
                )
            )

    return SynKitITSResult(
        graph=graph,
        changed_bonds=tuple(changes),
        synkit_version=installed_synkit_version() or "unknown",
    )
