# Rule JSON format

A rule file is a JSON list. Each object maps directly to `chem_predict.core.Rule`.

```json
[
  {
    "id": "unique_rule_id",
    "name": "Human-readable name",
    "reaction_smarts": "[C:1]=[O:2]>>[C:1][O:2]",
    "when": {
      "stresses": ["reduction"],
      "ph_min": null,
      "ph_max": null,
      "temperature_min_c": null,
      "temperature_max_c": null,
      "requires_oxygen": null,
      "requires_light": null
    },
    "priority": 0,
    "enabled": true,
    "source": "citation, DOI, URL, dataset record, or curation note",
    "description": "why this transformation exists",
    "tags": ["example"],
    "metadata": {}
  }
]
```

Supported stress labels currently are:

`acid`, `base`, `neutral`, `oxidation`, `reduction`, `photolysis`, `thermal`, `humidity`, `excipient`, `nitrosation`, `custom`.

Empty `when.stresses` means the rule does not restrict stress type. Numerical limits only match when the corresponding run condition is actually supplied; unknown conditions do not silently satisfy a constrained rule.
