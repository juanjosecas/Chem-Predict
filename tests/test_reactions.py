import pytest

from chem_predict.chemistry import ReactionArityError, apply_reaction


REACTION = "[C:1]=[O:2]>>[C:1][O:2]"


def test_apply_reaction_returns_canonical_products() -> None:
    assert apply_reaction(REACTION, ["CC=O"]) == [("CCO",)]


def test_reaction_arity_is_checked() -> None:
    with pytest.raises(ReactionArityError):
        apply_reaction(REACTION, ["CC=O", "O"])
