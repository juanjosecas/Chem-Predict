import pytest

from chem_predict.chemistry import InvalidSmilesError, canonicalize_smiles, deduplicate_smiles


def test_canonicalize_smiles() -> None:
    assert canonicalize_smiles("C1=CC=CN=C1") == "c1ccncc1"


def test_deduplicate_smiles_preserves_order() -> None:
    assert deduplicate_smiles(["C1=CC=CN=C1", "c1ccncc1", "CCO"]) == ["c1ccncc1", "CCO"]


def test_invalid_smiles_raises_package_error() -> None:
    with pytest.raises(InvalidSmilesError):
        canonicalize_smiles("not-a-smiles")
