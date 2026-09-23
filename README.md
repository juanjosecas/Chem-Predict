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
├── medchem/        # Lilly-style quality/reactivity filters
├── rulesources/    # primary-source reaction rules/reaction evidence
├── visualization/  # molecules, reactions, reaction-center depictions
├── properties/    # future property providers
├── purge/         # future impurity purge models
├── nitrosamines/  # CPCA, precursor/context screening, reactivity evidence
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

## Nitrosamine screening

The nitrosamine module keeps potency categorization, formation context, and
reactivity/purge evidence separate:

```python
from chem_predict.nitrosamines import (
    NitrosationContext,
    assess_cpca,
    assess_nitrosation_context,
    enumerate_secondary_amine_nitrosation_products,
)

cpca = assess_cpca("CCN(N=O)CC")
print(cpca.overall_category, cpca.overall_ai_ng_per_day)

context = assess_nitrosation_context(
    "CCNCC",
    NitrosationContext(nitrite_present=True, ph=3.5),
)
print(context.flags)

print(enumerate_secondary_amine_nitrosation_products("CCNCC"))
```

See `docs/nitrosamines.md` for scope, regulatory-source hierarchy, current
limitations, and literature provenance.

## Medchem and reaction informatics

```python
from chem_predict.medchem import LillyNativeRules

assessment = LillyNativeRules().assess("CCN(N=O)CC")
print(assessment.passed, assessment.total_demerits)
```

The native Lilly port is incremental and auditable; the package also provides
an adapter to a local official Lilly-Medchem-Rules/LillyMol installation for
exact comparison.

SynKit is optional:

```bash
python -m pip install -e ".[synkit]"
```

```python
from chem_predict.integrations import reaction_to_its

result = reaction_to_its(
    "[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]"
)
print(result.changed_bonds)
```

See `docs/integrations.md` and `THIRD_PARTY_NOTICES.md`.

## Primary reaction-rule sources

Direct reaction rules and primary reaction evidence are deliberately kept
separate.

```python
from chem_predict.rulesources import RetroRulesSource, to_rule_registry

source = RetroRulesSource()
records = source.search(ec="1.2.1", radius=4)

registry = to_rule_registry(records, tags=("retrorules",))
```

Rhea reactions can be downloaded directly from the official ExPASy export:

```python
from chem_predict.rulesources import RheaSource

source = RheaSource()
source.download("data/raw/rhea-reaction-smiles.tsv")
```

CheT official CSV exports, NORMAN/Zenodo records and explicitly supplied
enviPath rule URLs also have adapters. ORD and BioTransformer are registered as
external sources with their licensing constraints rather than being copied into
the MIT repository.

See `docs/rule-sources.md`.

## Visualization

```python
from chem_predict.visualization import (
    draw_reaction_center,
    molecule_svg,
    reaction_svg,
)

mol_svg = molecule_svg("CCO", atom_indices=True)

reaction = "[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]"
rxn_svg = reaction_svg(reaction)

center = draw_reaction_center(reaction)
print(center.changed_bonds)
```

The reaction-center view highlights changed mapped atoms/bonds independently on
the reactant and product sides. See `docs/visualization.md`.

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
5. Expand nitrosamine formation/persistence models and validate the open CPCA feature catalogue.
6. Expand primary-source adapters and derive validated reaction templates from Rhea, CheT, NORMAN and ORD reaction evidence.
7. Scoring/ranking layer kept separate from rule execution.
8. Provenance and validation reports suitable for reproducible research.

## License

MIT. See `LICENSE`.

## Mixture and context enumeration

Supply the mixture components explicitly and use the existing condition windows
on rules. Multi-reactant SMARTS use reactant template order; substrate, excipient,
and reagent roles are annotations, not evidence of reactivity.

```python
from chem_predict.core import Conditions, StressType
from chem_predict.degradation import MixtureComponent, MixtureDegradationEngine
from chem_predict.rules import RuleRegistry

registry = RuleRegistry.from_json("my_curated_rules.json")
results = MixtureDegradationEngine(registry).predict(
    [MixtureComponent("CBr", role="substrate"),
     MixtureComponent("[OH-]", role="reagent")],
    Conditions(stresses=frozenset({StressType.BASE}), ph=10.0),
)
for result in results:
    print(result.reactant_indices, result.rule_id, result.products, result.source)
```

This enumerates **one-step possible products**, without kinetics, yields,
concentration effects, competing reactions, or calibrated confidence. Rule
priority only determines output order. No reaction rules are downloaded
implicitly; use the existing rule-source adapters and curate applicability
before loading a registry. The combination limit raises on overflow rather
than silently returning a partial list. Products retain RDKit isomeric SMILES,
rule source and per-component indices for auditability.
