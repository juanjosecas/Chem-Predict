from chem_predict.rulesources import RetroRulesSource, to_rule_registry


source = RetroRulesSource()

# Direct API query for a small, targeted rule set.
records = source.search(ec="1.2.1", radius=4)
registry = to_rule_registry(records, tags=("retrorules", "biochemical"))

print(f"Loaded {len(tuple(registry))} targeted rules")

# For bulk work, prefer the official compressed release instead:
#
# source.download_templates(
#     "rhea",
#     "data/raw/retrorules-3.1.0-rhea.tsv.gz",
# )
#
# records = source.iter_downloaded_tsv(
#     "data/raw/retrorules-3.1.0-rhea.tsv.gz"
# )
# registry = to_rule_registry(records, tags=("retrorules", "rhea"))
