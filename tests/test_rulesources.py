from chem_predict.rulesources import (
    RetroRulesSource,
    SourceKind,
    get_source,
    list_sources,
    parse_chet_export_csv,
    parse_reaction_smiles_lines,
    parse_templates_tsv,
    to_core_rule,
)


def test_retrorules_tsv_parsing_and_core_conversion() -> None:
    text = (
        "TEMPLATE_ID\tTEMPLATE\tREACTIONS\tECS\tRADIUS_MIN\tSCORE\tVALID\tDATASETS\n"
        "RR:03-TEST\t[C:1]=[O:2]>>[C:1][O:2]\tRHEA:12345;RHEA:23456\t1.1.1.1\t4\t0.87\t1\trhea\n"
    )
    records = parse_templates_tsv(text)
    assert len(records) == 1
    record = records[0]
    assert record.id == "RR:03-TEST"
    assert record.radius == 4
    assert record.score == 0.87
    assert record.valid is True
    assert record.source_record_ids == ("RHEA:12345", "RHEA:23456")

    rule = to_core_rule(record, tags=("retrorules",))
    assert rule.reaction_smarts == record.reaction_smarts
    assert rule.metadata["radius"] == 4
    assert "retrorules" in rule.tags


def test_retrorules_uses_documented_api_shapes() -> None:
    source = RetroRulesSource()
    assert source.build_search_url(ec="1.2.1", radius=4).endswith(
        "/api/templates?ec=1.2.1&radius=4"
    )
    assert source.template_summary_url("RR:03-TEST").endswith(
        "/api/templates/RR%3A03-TEST/summary"
    )
    assert source.bulk_download_url("rhea", format="tsv") == (
        "https://retrorules.org/dl/v3.1.0/rhea/templates?format=tsv"
    )


def test_rhea_headerless_parser_finds_id_and_reaction_smiles() -> None:
    lines = [
        "12345\tCCO>>CC=O\n",
        "RHEA:23456\tCC=O>>CCO\textra\n",
    ]
    records = list(parse_reaction_smiles_lines(lines))
    assert [record.id for record in records] == ["RHEA:12345", "RHEA:23456"]
    assert records[0].reaction_smiles == "CCO>>CC=O"


def test_chet_export_keeps_provenance_and_builds_reaction_when_smiles_exist() -> None:
    text = (
        "Reaction ID,Parent SMILES,Product SMILES,Process,Reference\n"
        "rxn-1,CCO,CC=O,Abiotic Oxidation,example\n"
    )
    record = next(parse_chet_export_csv(text))
    assert record.id == "rxn-1"
    assert record.reaction_smiles == "CCO>>CC=O"
    assert record.metadata["Process"] == "Abiotic Oxidation"


def test_registry_distinguishes_direct_rules_from_primary_reactions() -> None:
    direct = list_sources(kind=SourceKind.REACTION_RULES)
    assert any(source.id == "retrorules" for source in direct)
    assert get_source("rhea").kind == SourceKind.REACTIONS
    assert get_source("envipath").redistributable_in_repo is False
