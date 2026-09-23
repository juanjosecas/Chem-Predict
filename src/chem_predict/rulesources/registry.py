from __future__ import annotations

from chem_predict.rulesources.models import AccessMode, RuleSourceDescriptor, SourceKind


SOURCES: dict[str, RuleSourceDescriptor] = {
    "retrorules": RuleSourceDescriptor(
        id="retrorules",
        name="RetroRules",
        kind=SourceKind.REACTION_RULES,
        primary_url="https://retrorules.org/",
        access_mode=AccessMode.API,
        license="CC BY 4.0 by default; verify provenance-specific upstream terms",
        redistributable_in_repo=True,
        version_checked="3.1.0; checked 2026-09-23",
        notes=(
            "Primary source for ready-to-use reaction SMARTS templates. "
            "Templates retain source reaction provenance and modeled radius."
        ),
    ),
    "rhea": RuleSourceDescriptor(
        id="rhea",
        name="Rhea",
        kind=SourceKind.REACTIONS,
        primary_url="https://www.rhea-db.org/",
        access_mode=AccessMode.DOWNLOAD,
        license="CC BY 4.0",
        redistributable_in_repo=True,
        version_checked="FTP release checked 2026-09-23",
        notes=(
            "Curated biochemical reactions. Rhea supplies directed reaction SMILES; "
            "Chem-Predict treats rules derived from these reactions as derived artifacts."
        ),
    ),
    "epa_chet": RuleSourceDescriptor(
        id="epa_chet",
        name="US EPA Chemical Transformations Database (CheT)",
        kind=SourceKind.TRANSFORMATIONS,
        primary_url="https://www.epa.gov/comptox-tools/chemical-transformations-database-chet-user-guide",
        access_mode=AccessMode.EXPORT,
        license="US EPA public resource; preserve record-level literature provenance",
        redistributable_in_repo=False,
        version_checked="user guide checked 2026-09-23",
        notes=(
            "Observed parent/product transformations with process, experimental conditions, "
            "kinetics and references. Official UI documents CSV export; no undocumented bulk "
            "endpoint is hard-coded."
        ),
    ),
    "norman_reftps": RuleSourceDescriptor(
        id="norman_reftps",
        name="NORMAN S74 REFTPS / FAIR-TPs",
        kind=SourceKind.TRANSFORMATIONS,
        primary_url="https://doi.org/10.5281/zenodo.4318838",
        access_mode=AccessMode.DOWNLOAD,
        license="Open Zenodo dataset; verify the license metadata of the selected version before redistribution",
        redistributable_in_repo=False,
        version_checked="concept DOI checked 2026-09-23",
        notes=(
            "Transformation products and reactions from literature. The deposition is "
            "continuously updated; adapters therefore accept an explicit Zenodo record ID."
        ),
    ),
    "ord": RuleSourceDescriptor(
        id="ord",
        name="Open Reaction Database",
        kind=SourceKind.REACTIONS,
        primary_url="https://open-reaction-database.org/",
        access_mode=AccessMode.DOWNLOAD,
        license="Data CC BY-SA 4.0; software Apache-2.0",
        redistributable_in_repo=False,
        version_checked="documentation checked 2026-09-23",
        notes=(
            "Single-step reaction records for synthesis/process chemistry. Keep ORD datasets "
            "outside the MIT source tree and preserve share-alike provenance."
        ),
    ),
    "envipath": RuleSourceDescriptor(
        id="envipath",
        name="enviPath",
        kind=SourceKind.REACTION_RULES,
        primary_url="https://envipath.org/",
        access_mode=AccessMode.EXTERNAL,
        license="Public knowledge packages include CC BY-NC-SA 4.0 material",
        redistributable_in_repo=False,
        version_checked="API/wiki checked 2026-09-23",
        notes=(
            "Rules are represented as SMIRKS with reactant/product exclusion SMARTS and "
            "likelihood metadata. Use as an external provider; do not vendor restricted packs."
        ),
    ),
    "biotransformer": RuleSourceDescriptor(
        id="biotransformer",
        name="BioTransformer",
        kind=SourceKind.REACTION_RULES,
        primary_url="https://github.com/Wishartlab-openscience/Biotransformer",
        access_mode=AccessMode.EXTERNAL,
        license="Software LGPL-3.0; environmental module includes enviPath/EAWAG restrictions",
        redistributable_in_repo=False,
        version_checked="3.0.0 repository checked 2026-09-23",
        notes="Use as an external validation/provider layer rather than copying its knowledge base.",
    ),
    "rdkit_filters": RuleSourceDescriptor(
        id="rdkit_filters",
        name="RDKit FilterCatalog",
        kind=SourceKind.STRUCTURAL_ALERTS,
        primary_url="https://www.rdkit.org/docs/",
        access_mode=AccessMode.EXTERNAL,
        license="RDKit BSD; individual catalog provenance should be retained",
        redistributable_in_repo=False,
        version_checked="RDKit 2026.03.6 docs checked 2026-09-23",
        notes="Structural alerts, not reaction rules; useful for the medchem screening layer.",
    ),
}


def get_source(source_id: str) -> RuleSourceDescriptor:
    try:
        return SOURCES[source_id]
    except KeyError as exc:
        raise KeyError(f"Unknown rule source: {source_id}") from exc


def list_sources(*, kind: SourceKind | None = None) -> tuple[RuleSourceDescriptor, ...]:
    sources = tuple(SOURCES.values())
    if kind is None:
        return sources
    return tuple(source for source in sources if source.kind == kind)
