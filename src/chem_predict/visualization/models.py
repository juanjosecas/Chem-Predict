from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class BondChange:
    atom_map_1: int
    atom_map_2: int
    reactant_order: float
    product_order: float


@dataclass(frozen=True, slots=True)
class ReactionCenterDepiction:
    reaction_smiles: str
    changed_bonds: tuple[BondChange, ...]
    changed_atom_maps: tuple[int, ...]
    reactants_svg: str
    products_svg: str

    def to_html(self) -> str:
        return (
            "<div style='display:flex;gap:1rem;align-items:center'>"
            f"<div>{self.reactants_svg}</div>"
            "<div style='font-size:2rem'>&rarr;</div>"
            f"<div>{self.products_svg}</div>"
            "</div>"
        )
