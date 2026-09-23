# Chem-Predict

Chem-Predict is a modular Python toolkit for building auditable chemical prediction workflows around reaction rules, degradation pathways, impurity fate, physicochemical properties, purge models, and nitrosamine risk.

The repository starts deliberately small: the first layer is a stable chemical/reaction core plus registries and interfaces that can be extended without rewriting the package.

> **Status:** early research scaffold. Outputs are not validated for regulatory or safety-critical decisions.

## Design principles

- **Chemistry engine and chemical knowledge are separate.** Reaction/degradation rules are data objects, not hard-coded branches.
- **Provenance first.** Every rule can carry a source, description, tags, priority, and arbitrary metadata.
- **Fail loudly while curating.** Invalid SMILES, reaction arity errors, and invalid generated products raise explicit exceptions by default.
- **No hidden probability claims.** Rule priority is an ordering field, not a probability of degradation.
- **Dependency discipline.** External APIs are added only after checking their current official documentation.
- **Modular growth.** `degradation`, `properties`, `purge`, and `nitrosamines` can evolve independently over the shared core.

## Current architecture

```text
src/chem_predict/
├── core/          # shared domain models
├── chemistry/     # molecule and reaction primitives
├── rules/         # rule registry + JSON serialization
├── degradation/   # degradation prediction engine
├── properties/    # future property providers
├── purge/         # future impurity purge models
├── nitrosamines/  # future nitrosamine-specific logic
└── cli.py          # small command-line interface
```

## Requirements

- Python >= 3.11
- RDKit >= 2025.09.4

The initial RDKit integration was checked against the official RDKit 2026.03.6 documentation.

## Installation for development

```bash
git clone https://github.com/juanjosecas/Chem-Predict.git
cd Chem-Predict
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
pytest
```

## Quick start

Normalize a structure:

```python
from chem_predict.chemistry import canonicalize_smiles

print(canonicalize_smiles("C1=CC=CN=C1"))
# c1ccncc1
```

Apply a reaction SMARTS:

```python
from chem_predict.chemistry import apply_reaction

outcomes = apply_reaction(
    "[C:1]=[O:2]>>[C:1][O:2]",
    ["CC=O"],
)
print(outcomes)
# [('CCO',)]
```

Create a rule and run the degradation engine:

```python
from chem_predict.core import ConditionWindow, Conditions, Rule, StressType
from chem_predict.degradation import DegradationEngine
from chem_predict.rules import RuleRegistry

rule = Rule(
    id="demo_carbonyl_reduction",
    name="Demo carbonyl reduction",
    reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
    when=ConditionWindow(stresses=frozenset({StressType.REDUCTION})),
    source="Demonstration only; not a curated degradation rule",
)

registry = RuleRegistry([rule])
engine = DegradationEngine(registry)

predictions = engine.predict(
    "CC=O",
    Conditions(stresses=frozenset({StressType.REDUCTION})),
)

for prediction in predictions:
    print(prediction.rule_id, prediction.products)
```

The example transformation exists only to test the machinery. Scientific degradation rules should be separately curated and referenced.

## CLI

```bash
chem-predict normalize "C1=CC=CN=C1"
chem-predict apply "[C:1]=[O:2]>>[C:1][O:2]" "CC=O"
chem-predict predict examples/rules_demo.json "CC=O" --stress reduction
```

## Near-term roadmap

1. Curated, versioned degradation-rule library with citations and applicability constraints.
2. Multi-reactant handling for API-excipient and nitrosation chemistry.
3. Pluggable physicochemical-property providers.
4. Explicit impurity fate/purge model separating reactivity, solubility, volatility, and process operations.
5. Nitrosamine precursor/formation/persistence module.
6. Dataset adapters for reaction/degradation corpora without coupling datasets to the core engine.
7. Scoring/ranking layer kept separate from rule execution.
8. Provenance and validation reports suitable for reproducible research.

## License

MIT. See `LICENSE`.
