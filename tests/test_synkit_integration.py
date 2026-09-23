import pytest

from chem_predict.integrations.synkit import (
    SynKitUnavailableError,
    installed_synkit_version,
    reaction_to_its,
)


REACTION = "[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]"


def test_synkit_is_optional() -> None:
    if installed_synkit_version() is not None:
        pytest.skip("SynKit is installed in this environment")
    with pytest.raises(SynKitUnavailableError):
        reaction_to_its(REACTION)


def test_synkit_tuple_its_when_installed() -> None:
    if installed_synkit_version() is None:
        pytest.skip("SynKit optional dependency is not installed")

    result = reaction_to_its(REACTION)
    changes = {
        (change.reactant_order, change.product_order)
        for change in result.changed_bonds
    }
    assert (1.0, 0.0) in changes
    assert (0.0, 1.0) in changes
    assert result.synkit_version.startswith("1.6.")
