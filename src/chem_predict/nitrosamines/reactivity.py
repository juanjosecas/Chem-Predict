from __future__ import annotations

from chem_predict.nitrosamines.models import ReactivityEvidence


# Evidence records are intentionally descriptive. They are not executable purge
# factors or universal reaction rules.
LITERATURE_REACTIVITY_EVIDENCE: tuple[ReactivityEvidence, ...] = (
    ReactivityEvidence(
        id="hodgin_2024_dithionite_basic_aqueous_50c",
        nitrosamine_scope="Experimental screen of eight commercial N-nitrosamines",
        reagent_or_condition="Na2S2O4 in 1 M aqueous NaOH at 50 °C",
        outcome="Reported as highly effective at consuming N-nitrosamines, with reduction to parent amines discussed in the study",
        evidence_type="experimental_reactivity_screen",
        source_doi="10.1021/acs.oprd.4c00217",
        source_url="https://pubs.acs.org/doi/10.1021/acs.oprd.4c00217",
        notes="Do not convert this record into a numerical purge factor without compound- and process-specific evidence.",
    ),
)


def list_reactivity_evidence() -> tuple[ReactivityEvidence, ...]:
    return LITERATURE_REACTIVITY_EVIDENCE
