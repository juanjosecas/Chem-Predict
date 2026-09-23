# Architecture

Chem-Predict separates execution infrastructure from scientific knowledge.

## Layers

1. **Core models** contain serializable domain objects (`Rule`, `Conditions`, `Prediction`). They do not depend on RDKit.
2. **Chemistry primitives** are the narrow RDKit boundary: parsing, canonicalization, sanitization, and reaction execution.
3. **Rule registry** stores and serializes scientific transformation rules without embedding them in Python control flow.
4. **Domain engines** consume the shared core. The first is `DegradationEngine`; later modules can implement purge or nitrosamine-specific workflows without changing the chemistry layer.
5. **Scoring/property layers** should remain optional. A rule engine must still work without a trained model or external service.

## Extension rule

A new scientific capability should normally be added as one of:

- a new rule/data file;
- a new adapter/provider inside a domain namespace;
- a new engine consuming the existing models;
- a new model field only when the concept is genuinely shared across domains.

Avoid adding domain-specific switches to `chemistry/` or giant conditional blocks to `DegradationEngine`.

## Provenance

`Rule.source`, `Rule.description`, `Rule.tags`, and `Rule.metadata` are deliberately generic. A future curation schema can add DOI, patent, dataset record, curator, version, and evidence level without changing RDKit execution.

## Scoring

`Rule.priority` is deterministic ordering metadata, not a probability. Probabilistic or ML ranking should be implemented as a separate scorer so the raw generated pathway remains reproducible and auditable.
