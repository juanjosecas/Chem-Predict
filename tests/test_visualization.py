from chem_predict.visualization import (
    draw_reaction_center,
    find_reaction_center,
    molecule_svg,
    reaction_svg,
)


MAPPED_SUBSTITUTION = (
    "[CH3:1][Br:2].[OH-:3]"
    ">>"
    "[CH3:1][OH:3].[Br-:2]"
)


def test_molecule_svg() -> None:
    svg = molecule_svg("CCO", atom_indices=True)
    assert "<svg" in svg
    assert "</svg>" in svg


def test_reaction_svg_from_smiles() -> None:
    svg = reaction_svg(MAPPED_SUBSTITUTION)
    assert "<svg" in svg
    assert "</svg>" in svg


def test_find_reaction_center_from_atom_mapping() -> None:
    changes = find_reaction_center(MAPPED_SUBSTITUTION)
    observed = {
        (change.atom_map_1, change.atom_map_2, change.reactant_order, change.product_order)
        for change in changes
    }
    assert (1, 2, 1.0, 0.0) in observed
    assert (1, 3, 0.0, 1.0) in observed


def test_draw_reaction_center_returns_highlighted_sides() -> None:
    view = draw_reaction_center(MAPPED_SUBSTITUTION)
    assert view.changed_atom_maps == (1, 2, 3)
    assert "<svg" in view.reactants_svg
    assert "<svg" in view.products_svg
    assert "&rarr;" in view.to_html()
