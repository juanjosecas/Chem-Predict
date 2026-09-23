# Primary reaction-rule and transformation sources

Chem-Predict separates three kinds of knowledge:

1. **direct reaction rules** published by a source;
2. **primary reaction/transformation records** from which Chem-Predict may later
   derive rules;
3. **external/restricted knowledge bases** that must not be silently copied into
   the MIT repository.

This distinction is stored in `chem_predict.rulesources.SOURCES`.

## RetroRules

Primary URL: https://retrorules.org/

Checked: 2026-09-23. The live site reports RetroRules 3.1.0.

RetroRules is currently the preferred source for ready-to-use reaction
templates because it exposes reaction SMARTS directly and retains links to
source reactions, EC numbers, modeled radii and scores.

The adapter only uses API endpoints explicitly documented by RetroRules:

```text
GET /api/templates?ec=1.2.1&radius=4
GET /api/templates/<TEMPLATE_ID>/summary
```

Example:

```python
from chem_predict.rulesources import RetroRulesSource, to_rule_registry

source = RetroRulesSource()
records = source.search(ec="1.2.1", radius=4)

registry = to_rule_registry(records, tags=("retrorules", "biochemical"))
```

A downloaded RetroRules TSV can also be parsed without network access:

```python
from pathlib import Path
from chem_predict.rulesources import parse_templates_tsv

records = parse_templates_tsv(Path("rules.tsv").read_text())
```

RetroRules states that data are CC BY 4.0 unless stated otherwise, while also
warning that integrated source databases may carry their own terms. Preserve
the original provenance fields.

## Rhea

Primary URL: https://www.rhea-db.org/

Official directed reaction-SMILES export:

```text
https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv
```

License:

```text
https://ftp.expasy.org/databases/rhea/LICENSE.txt
```

Rhea is CC BY 4.0. It provides reactions rather than generalized reaction
rules. Chem-Predict therefore represents these as `ReactionRecord` objects.

```python
from chem_predict.rulesources import RheaSource

source = RheaSource()
source.download("data/raw/rhea-reaction-smiles.tsv")

for reaction in source.read("data/raw/rhea-reaction-smiles.tsv"):
    print(reaction.id, reaction.reaction_smiles)
```

Any SMARTS/template later extracted from these reactions must be marked as a
Chem-Predict-derived rule with the Rhea reaction IDs retained as provenance.

## US EPA CheT

Primary documentation:

https://www.epa.gov/comptox-tools/chemical-transformations-database-chet-user-guide

CheT contains observed parent/product transformations plus library/process
classification, experimental details, kinetics and references.

The official documentation describes CSV export through the CheT interface but
does not document a stable public bulk API endpoint. Chem-Predict deliberately
does not reverse-engineer one.

Export a reaction table from CheT and parse it locally:

```python
from pathlib import Path
from chem_predict.rulesources import parse_chet_export_csv

text = Path("chet_export.csv").read_text()

for record in parse_chet_export_csv(text):
    print(record.id, record.reaction_smiles)
    print(record.metadata)
```

The parser preserves every original CSV column. If SMILES columns are present,
it builds a reaction SMILES; otherwise the record can later be joined to a
chemical export through the identifiers supplied by CheT.

## NORMAN S74 REFTPS / FAIR-TPs

Concept DOI:

https://doi.org/10.5281/zenodo.4318838

S74 is continuously updated, so Chem-Predict requires an explicit Zenodo
record ID rather than silently pinning a moving dataset.

```python
from chem_predict.rulesources import NormanRefTPSSource

source = NormanRefTPSSource()
source.download_transformations(
    record_id=17895771,
    path="data/raw/S74_REFTPS_Transformations.csv",
)
```

The record ID above is an example version seen during development, not a
promise that it is the newest version. Check the concept DOI before a
production data refresh.

## Open Reaction Database

Primary documentation:

https://docs.open-reaction-database.org/

ORD provides structured single-step organic reaction data. ORD data are
CC BY-SA 4.0 and its software is Apache-2.0.

Chem-Predict currently records ORD in the source registry but does not vendor
ORD datasets or add `ord-schema` as a core dependency. A future optional
adapter should consume local ORD datasets and preserve dataset/reaction IDs.
Any derived data subject to share-alike requirements must stay outside the MIT
software source tree unless separately licensed.

## enviPath

Primary URL: https://envipath.org/

enviPath rules are particularly useful because they contain SMIRKS,
reactant/product exclusion SMARTS and likelihood metadata.

Important public knowledge packages contain CC BY-NC-SA 4.0 material.
Chem-Predict therefore exposes only an external provider for explicit rule URLs
and does not mirror those packs into the repository.

## BioTransformer

Primary repository:

https://github.com/Wishartlab-openscience/Biotransformer

BioTransformer is useful as an external metabolism/degradation predictor and
knowledge reference. Its environmental module incorporates EAWAG/enviPath
material with additional restrictions. The knowledge base is therefore not
vendored into Chem-Predict.

## Structural-alert sources

RDKit FilterCatalog and the ChEMBL-derived alert collections belong in the
medchem screening layer rather than in the reaction-rule layer. They remain
registered separately so that an alert cannot accidentally be interpreted as
a chemical transformation.

## Rule provenance policy

Every imported/generated rule should retain, at minimum:

```text
source
source version
source record/reaction IDs
retrieval date
original reaction/template string
license
whether the rule was direct or derived
derivation method/version
validation status
```

A rule generated from Rhea, CheT, NORMAN or ORD is not labelled as a native
rule from that source. The source supplies the reaction evidence; Chem-Predict
supplies the rule extraction.
