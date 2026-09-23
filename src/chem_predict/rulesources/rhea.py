from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator

from chem_predict.rulesources.http import download_to_path, fetch_text
from chem_predict.rulesources.models import ReactionRecord


RHEA_REACTION_SMILES_URL = (
    "https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv"
)
RHEA_LICENSE_URL = "https://ftp.expasy.org/databases/rhea/LICENSE.txt"


def parse_reaction_smiles_lines(lines: Iterable[str]) -> Iterator[ReactionRecord]:
    """Parse Rhea's directed reaction-SMILES TSV.

    The upstream file is intentionally treated as headerless because that is
    how the public export is documented/consumed. The parser locates the field
    containing '>>' instead of depending on a fragile fixed column count.
    """

    for line_number, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        fields = [field.strip() for field in line.split("\t")]
        reaction_smiles = next((field for field in fields if ">>" in field), None)
        if reaction_smiles is None:
            continue

        identifier = next(
            (
                field
                for field in fields
                if field != reaction_smiles
                and field.removeprefix("RHEA:").isdigit()
            ),
            None,
        )
        if identifier is None:
            identifier = f"line-{line_number}"
        if not identifier.startswith("RHEA:"):
            identifier = f"RHEA:{identifier}"

        yield ReactionRecord(
            id=identifier,
            source="rhea",
            reaction_smiles=reaction_smiles,
            metadata={"line_number": line_number},
        )


class RheaSource:
    """Primary-source adapter for Rhea directed reaction SMILES."""

    reaction_smiles_url = RHEA_REACTION_SMILES_URL

    def fetch(self, *, limit: int | None = None) -> list[ReactionRecord]:
        text = fetch_text(self.reaction_smiles_url)
        records = parse_reaction_smiles_lines(text.splitlines())
        if limit is None:
            return list(records)
        if limit < 1:
            raise ValueError("limit must be >= 1")
        out: list[ReactionRecord] = []
        for record in records:
            out.append(record)
            if len(out) >= limit:
                break
        return out

    def download(self, path: str | Path) -> Path:
        return download_to_path(self.reaction_smiles_url, path)

    def read(self, path: str | Path) -> Iterator[ReactionRecord]:
        with Path(path).open("r", encoding="utf-8") as handle:
            yield from parse_reaction_smiles_lines(handle)
