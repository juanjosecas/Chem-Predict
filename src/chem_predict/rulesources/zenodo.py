from __future__ import annotations

from pathlib import Path
from typing import Any

from chem_predict.rulesources.http import download_to_path, fetch_json


ZENODO_API = "https://zenodo.org/api/records"
NORMAN_REFTPS_CONCEPT_DOI = "10.5281/zenodo.4318838"
NORMAN_REFTPS_TRANSFORMATIONS = "S74_REFTPS_Transformations.csv"


class ZenodoRecordSource:
    def metadata(self, record_id: int | str) -> dict[str, Any]:
        payload = fetch_json(f"{ZENODO_API}/{record_id}")
        if not isinstance(payload, dict):
            raise ValueError("Unexpected Zenodo record payload")
        return payload

    def file_url(self, record_id: int | str, filename: str) -> str:
        return f"https://zenodo.org/records/{record_id}/files/{filename}?download=1"

    def download_file(
        self,
        record_id: int | str,
        filename: str,
        path: str | Path,
    ) -> Path:
        return download_to_path(self.file_url(record_id, filename), path)


class NormanRefTPSSource(ZenodoRecordSource):
    """Version-explicit downloader for NORMAN S74 transformation records."""

    concept_doi = NORMAN_REFTPS_CONCEPT_DOI

    def download_transformations(
        self,
        *,
        record_id: int | str,
        path: str | Path,
    ) -> Path:
        return self.download_file(record_id, NORMAN_REFTPS_TRANSFORMATIONS, path)
