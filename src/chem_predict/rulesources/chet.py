from __future__ import annotations

import csv
import io
from typing import Iterator

from chem_predict.rulesources.models import ReactionRecord


CHET_GUIDE_URL = (
    "https://www.epa.gov/comptox-tools/"
    "chemical-transformations-database-chet-user-guide"
)


_PARENT_SMILES_NAMES = {
    "parent smiles",
    "parent_smiles",
    "parentsmiles",
}
_PRODUCT_SMILES_NAMES = {
    "product smiles",
    "product_smiles",
    "productsmiles",
}
_ID_NAMES = {
    "reaction id",
    "reaction_id",
    "reactionid",
    "id",
}


def _normalized(row: dict[str, str]) -> dict[str, str]:
    return {
        str(key).strip().lower(): (value or "").strip()
        for key, value in row.items()
        if key is not None
    }


def parse_chet_export_csv(text: str) -> Iterator[ReactionRecord]:
    """Parse an official CheT CSV export without assuming undocumented endpoints.

    CheT exports can vary by reaction library. All original columns are kept in
    metadata. If parent/product SMILES are present, a reaction SMILES is also
    constructed; otherwise the record remains useful as provenance/condition
    evidence and can be joined to the CheT chemical export by DTXSID.
    """

    reader = csv.DictReader(io.StringIO(text))
    for row_number, raw_row in enumerate(reader, start=2):
        row = {str(k): (v or "") for k, v in raw_row.items() if k is not None}
        norm = _normalized(row)

        parent = next((norm[name] for name in _PARENT_SMILES_NAMES if norm.get(name)), None)
        product = next((norm[name] for name in _PRODUCT_SMILES_NAMES if norm.get(name)), None)
        identifier = next((norm[name] for name in _ID_NAMES if norm.get(name)), None)
        identifier = identifier or f"chet-row-{row_number}"

        reaction_smiles = f"{parent}>>{product}" if parent and product else None
        products = (product,) if product else ()

        yield ReactionRecord(
            id=identifier,
            source="epa_chet",
            reaction_smiles=reaction_smiles,
            parent_smiles=parent,
            product_smiles=products,
            metadata={"row_number": row_number, **row},
        )
