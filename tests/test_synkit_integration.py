import pytest

from chem_predict.integrations.synkit import (
    SynKitUnavailableError,
    installed_synkit_version,
    reaction_to_its,
)


def test_synkit_is_optional() -> None:
    if installed_synkit_version() is not None:
        pytest.skip("SynKit is installed in this environment")
    with pytest.raises(SynKitUnavailableError):
        reaction_to_its("[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]")
