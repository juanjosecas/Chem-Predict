# Nitrosamine module

The nitrosamine namespace is split into independent layers. This is deliberate:
formation risk, structural potency categorization, reactivity/purge evidence, and
future mechanistic/ML models are different scientific questions.

## Source hierarchy

Rules and constants should be curated in this order:

1. **Current health-authority guidance and live tables.**
   - FDA, *Control of Nitrosamine Impurities in Human Drugs* (September 2024, Rev. 2):
     https://www.fda.gov/regulatory-information/search-fda-guidance-documents/control-nitrosamine-impurities-human-drugs
   - FDA, *Recommended Acceptable Intake Limits for Nitrosamine Drug Substance-Related Impurities* (August 2023):
     https://www.fda.gov/regulatory-information/search-fda-guidance-documents/recommended-acceptable-intake-limits-nitrosamine-drug-substance-related-impurities
   - FDA live acceptable-intake tables:
     https://www.fda.gov/regulatory-information/search-fda-guidance-documents/cder-nitrosamine-impurity-acceptable-intake-limits
2. **Peer-reviewed experimental/mechanistic literature.**
3. **Open implementations**, used as implementation/validation references rather
   than as normative regulatory sources.
4. Secondary reviews and vendor articles are contextual only.

Regulatory information evolves. A computed result is therefore returned with
source/version metadata and an explicit warning to verify the current FDA table.

## Public API

### Detect an existing N-nitrosamine group

\`\`\`python
from chem_predict.nitrosamines import find_nitrosamine_sites

sites = find_nitrosamine_sites("CCN(N=O)CC")
for site in sites:
    print(site.alpha_hydrogens, site.cpca_applicable)
\`\`\`

The detector separates **detection** from **CPCA applicability**. For example,
N-nitroso structures with an alpha carbon directly double-bonded to a
heteroatom can be detected while being marked outside the CPCA structural
scope.

### FDA CPCA structural assessment

\`\`\`python
from chem_predict.nitrosamines import assess_cpca

result = assess_cpca("CCN(N=O)CC")
print(result.overall_category)
print(result.overall_ai_ng_per_day)
print(result.warnings)
\`\`\`

The initial implementation covers:

- alpha-hydrogen scoring;
- the 0/2 alpha-H ethyl exception;
- tertiary alpha-carbon rule;
- carboxylic acid;
- pyrrolidine, morpholine, sulfur-containing 6-membered, 5/6-membered, and
  7-membered ring features;
- long-chain feature;
- alpha electron-withdrawing feature using an explicit finite open catalogue;
- beta-hydroxyl;
- aryl-on-alpha activating feature;
- beta-methyl activating feature;
- one or two N-nitroso sites.

The current category-to-AI mapping used by the code is:

| CPCA category | AI (ng/day) |
| ---: | ---: |
| 1 | 26.5 |
| 2 | 100 |
| 3 | 400 |
| 4 | 1500 |
| 5 | 1500 |

These are not a substitute for a nitrosamine-specific value in the current FDA
tables. Compound-specific/read-across values take precedence when FDA provides
one.

The EWG feature is intentionally flagged as needing continued curation. The FDA
table refers to a defined feature set and supporting literature; the open
implementation should expand this feature only with explicit test cases and
provenance rather than silently broadening SMARTS.

### Find nitrosatable amine centers

\`\`\`python
from chem_predict.nitrosamines import find_nitrosatable_centers

for center in find_nitrosatable_centers("CN(C)CC"):
    print(center.kind, center.methyl_substituent_count, center.amide_like)
\`\`\`

Secondary, tertiary, dimethyl-tertiary, and quaternary centers are kept
separate. Detection is not a probability of NDSRI formation.

### Add process/formulation context

\`\`\`python
from chem_predict.nitrosamines import NitrosationContext, assess_nitrosation_context

assessment = assess_nitrosation_context(
    "CCNCC",
    NitrosationContext(
        nitrite_present=True,
        ph=3.5,
        storage=True,
    ),
)
print(assessment.flags)
\`\`\`

This function reports evidence flags. It deliberately does **not** generate a
made-up risk percentage.

### Enumerate direct secondary-amine nitrosation products

\`\`\`python
from chem_predict.nitrosamines import enumerate_secondary_amine_nitrosation_products

print(enumerate_secondary_amine_nitrosation_products("CCNCC"))
# ['CCN(CC)N=O']
\`\`\`

This helper only represents a direct secondary-amine N-nitrosation
transformation. Tertiary/quaternary pathways can require impurities,
degradation, or dealkylation and are not collapsed into this one-step rule.

## Reactivity and purge evidence

Reactivity evidence is stored separately from executable transformation rules:

\`\`\`python
from chem_predict.nitrosamines import list_reactivity_evidence

for record in list_reactivity_evidence():
    print(record.reagent_or_condition, record.source_doi)
\`\`\`

The first curated record comes from Hodgin et al., *Organic Process Research &
Development* 2024, DOI 10.1021/acs.oprd.4c00217. The study combined a reaction
database survey with experiments on eight commercial nitrosamines and reported
strong consumption/reduction under sodium dithionite in 1 M aqueous NaOH at
50 °C.

That observation is **not** encoded as a universal numerical purge factor.
Substrate, concentration, reaction time, phase behavior, workup, and process
context must be represented before quantitative purge prediction is justified.

## Open implementation reference

Novartis:
https://github.com/Novartis/NDSRIs_in_silico_tool

Associated publication:
Green Chemistry (2024), DOI 10.1039/D3GC03478J.

Chem-Predict does not copy that implementation. It is used as an open
cross-check of the CPCA workflow. The repository contains mixed licensing
signals at file/repository level and it also reflects an older regulatory
snapshot, so the Chem-Predict implementation is rebuilt from current FDA
sources and independently tested.

## Future layer: mechanistic / weight-of-evidence models

The package should later accept pluggable predictors for:

- quantum-mechanical descriptors;
- CYP bioactivation/binding models;
- mutagenicity models;
- read-across/similarity evidence;
- experimental Ames and in-vivo evidence;
- calibrated uncertainty.

These should produce a weight-of-evidence object and should not overwrite the
raw CPCA result. A useful recent scientific reference is:

https://www.science.org/doi/10.1126/sciadv.aee8680

## Non-normative contextual sources

The following sources supplied during project development are useful for
workflow ideas but are not used to define regulatory rules:

- https://zamann-pharma.com/2025/01/15/the-role-of-artificial-intelligence-in-nitrosamine-risk/
- https://intuitionlabs.ai/articles/nitrosamine-impurities-fda-guidance-ai-risk-assessment

## Current limitations

- This is a screening/research implementation, not a validated regulatory tool.
- The CPCA feature catalogue requires continued validation against FDA examples
  and curated reference compounds.
- Formation currently captures structural precursors and context flags, not
  formation kinetics or concentrations.
- Direct product enumeration is intentionally limited.
- Reactivity evidence is descriptive; quantitative purge modeling belongs in
  the separate \`purge\` domain.
