from chem_predict.nitrosamines import (
    AmineKind,
    NitrosationContext,
    assess_cpca,
    assess_nitrosation_context,
    enumerate_secondary_amine_nitrosation_products,
    find_nitrosamine_sites,
    find_nitrosatable_centers,
)


def test_detects_ndsr_like_site_and_alpha_hydrogens() -> None:
    sites = find_nitrosamine_sites("CCN(N=O)CC")
    assert len(sites) == 1
    assert sites[0].alpha_hydrogens == (2, 2)
    assert sites[0].cpca_applicable


def test_cpca_excludes_nitroso_aromatic_ring() -> None:
    sites = find_nitrosamine_sites("O=Nn1cc2ccccc2c1")
    assert len(sites) == 1
    assert not sites[0].cpca_applicable
    assert "n_nitroso_group_in_aromatic_ring" in sites[0].exclusion_reasons


def test_cpca_category_mapping_uses_current_fda_ai() -> None:
    result = assess_cpca("CCN(N=O)CC")
    site = result.sites[0]
    assert site.alpha_hydrogen_score == 1
    assert site.potency_score == 1
    assert site.category == 1
    assert site.ai_ng_per_day == 26.5


def test_pyrrolidine_ring_is_deactivating_feature() -> None:
    result = assess_cpca("O=NN1CCCC1")
    site = result.sites[0]
    assert site.category == 4
    assert site.potency_score == 4
    assert any(
        feature.id == "pyrrolidine_ring" and feature.score == 3
        for feature in site.deactivating_features
    )


def test_morpholine_ring_feature() -> None:
    result = assess_cpca("O=NN1CCOCC1")
    site = result.sites[0]
    assert any(
        feature.id == "morpholine_ring" and feature.score == 1
        for feature in site.deactivating_features
    )


def test_secondary_and_dimethyl_tertiary_centers_are_distinct() -> None:
    secondary = find_nitrosatable_centers("CCNCC")
    assert secondary[0].kind == AmineKind.SECONDARY

    tertiary = find_nitrosatable_centers("CN(C)CC")
    assert tertiary[0].kind == AmineKind.DIMETHYL_TERTIARY


def test_context_reports_fda_root_cause_combination_without_probability() -> None:
    result = assess_nitrosation_context(
        "CCNCC",
        NitrosationContext(nitrite_present=True, ph=3.5),
    )
    assert result.structural_precursor_present
    assert result.nitrosating_source_supported is True
    assert result.acidic_conditions is True
    assert "fda_root_cause_combination_present" in result.flags


def test_secondary_amine_nitrosation_enumeration() -> None:
    products = enumerate_secondary_amine_nitrosation_products("CCNCC")
    assert products == ["CCN(CC)N=O"]


def test_tertiary_alpha_carbon_is_category_five() -> None:
    result = assess_cpca("CC(C)(C)N(N=O)CC")
    site = result.sites[0]
    assert site.category == 5
    assert site.ai_ng_per_day == 1500.0
    assert "tertiary alpha-carbon" in site.rationale[0]


def test_cpca_scope_excludes_alpha_carbonyl_structure() -> None:
    result = assess_cpca("CCN(N=O)C(C)=O")
    site = result.sites[0]
    assert site.category is None
    assert not site.site.cpca_applicable
    assert "alpha_carbon_directly_double_bonded_to_heteroatom" in site.warnings


def test_carboxylic_acid_adds_three_deactivating_points() -> None:
    result = assess_cpca("CCN(N=O)CC(=O)O")
    site = result.sites[0]
    assert site.potency_score == 4
    assert site.category == 4
    assert any(
        feature.id == "carboxylic_acid" and feature.score == 3
        for feature in site.deactivating_features
    )


def test_benzylic_aryl_feature_is_activating() -> None:
    result = assess_cpca("CCN(N=O)Cc1ccccc1")
    site = result.sites[0]
    assert any(
        feature.id == "aryl_on_alpha" and feature.score == -1
        for feature in site.activating_features
    )
