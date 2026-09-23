from chem_predict.medchem import LillyNativeRules


def test_lilly_native_rejects_nitrosamine() -> None:
    result = LillyNativeRules(min_heavy_atoms=1).assess("CCN(N=O)CC")
    assert result.rejected
    assert any(hit.rule_id == "lilly_nitrosamine" for hit in result.hits)


def test_lilly_native_ester_demerit_is_35() -> None:
    result = LillyNativeRules(min_heavy_atoms=1).assess("CC(=O)OCC")
    hit = next(hit for hit in result.hits if hit.rule_id == "lilly_ester")
    assert hit.demerit_each == 35
    assert result.total_demerits >= 35


def test_lilly_native_nitro_demerit_is_60() -> None:
    result = LillyNativeRules(min_heavy_atoms=1).assess("CC[N+](=O)[O-]")
    hit = next(hit for hit in result.hits if hit.rule_id == "lilly_nitro")
    assert hit.demerit_each == 60
    assert "lilly_nitro_is_approximate" in result.warnings


def test_lilly_native_rejects_acid_halide() -> None:
    result = LillyNativeRules(min_heavy_atoms=1).assess("CC(=O)Cl")
    assert result.rejected
    assert any(hit.rule_id == "lilly_acid_halide" for hit in result.hits)


def test_lilly_default_rejects_too_small() -> None:
    result = LillyNativeRules().assess("CCO")
    assert result.rejected
    assert any(hit.rule_id == "lilly_too_few_atoms" for hit in result.hits)
