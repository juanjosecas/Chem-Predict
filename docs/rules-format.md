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
      "duration_min_h": null,
      "duration_max_h": null,
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

## Condition matching

Nonempty `stresses` matches **any** shared label (OR). All other specified requirements must hold together (AND). Numeric bounds are inclusive. An omitted run condition does not satisfy a numeric or boolean constraint. Duration is in hours and must be nonnegative; temperature must exceed absolute zero. Inverted windows and nonfinite numeric values are rejected.

Priority controls deterministic ordering (higher first, then rule ID); it is not a probability. A rule window filters applicability, not reaction rate. The engine does not infer stress labels from pH or add reactive species from condition flags.

See the [demonstration registry](../examples/rules_mixture_demo.json), [CLI](cli.md), and [notebook guide](notebook.md) for runnable multi-reactant examples. Record source and validation status when replacing these demonstration rules with curated rules.
