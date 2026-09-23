from __future__ import annotations

import csv
import io
from typing import Any
from urllib.parse import urlencode

from chem_predict.rulesources.http import fetch_json
from chem_predict.rulesources.models import ReactionRuleRecord, SourceFetchError


BASE_URL = "https://retrorules.org"
TEMPLATES_API = f"{BASE_URL}/api/templates"
RETRO_RULES_VERSION = "3.1.0"


def _first(row: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in row and row[key] not in (None, ""):
            return row[key]
    return None


def _split_values(value: Any) -> tuple[str, ...]:
    if value in (None, ""):
        return ()
    if isinstance(value, (list, tuple)):
        return tuple(str(item) for item in value if str(item))
    text = str(value).strip()
    separator = ";" if ";" in text else ","
    return tuple(part.strip() for part in text.split(separator) if part.strip())


def _parse_bool(value: Any) -> bool | None:
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _parse_int(value: Any) -> int | None:
    if value in (None, ""):
        return None
    return int(value)


def _parse_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    return float(value)


def parse_template_row(row: dict[str, Any]) -> ReactionRuleRecord:
    template_id = _first(row, "TEMPLATE_ID", "template_id", "id")
    template = _first(row, "TEMPLATE", "template", "smarts", "reaction_smarts")
    if not template_id or not template:
        raise ValueError("RetroRules row requires template id and reaction SMARTS")

    radius = _first(row, "RADIUS", "radius", "RADIUS_MIN", "radius_min")
    known = {
        "TEMPLATE_ID", "template_id", "id",
        "TEMPLATE", "template", "smarts", "reaction_smarts",
        "REACTIONS", "reactions",
        "ECS", "ecs", "ec_numbers",
        "RADIUS", "radius", "RADIUS_MIN", "radius_min",
        "SCORE", "score",
        "VALID", "valid",
    }

    return ReactionRuleRecord(
        id=str(template_id),
        reaction_smarts=str(template),
        source="retrorules",
        source_record_ids=_split_values(_first(row, "REACTIONS", "reactions")),
        ec_numbers=_split_values(_first(row, "ECS", "ecs", "ec_numbers")),
        radius=_parse_int(radius),
        score=_parse_float(_first(row, "SCORE", "score")),
        valid=_parse_bool(_first(row, "VALID", "valid")),
        metadata={
            "retrorules_version": RETRO_RULES_VERSION,
            **{str(k): v for k, v in row.items() if k not in known},
        },
    )


def parse_templates_tsv(text: str) -> list[ReactionRuleRecord]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    return [parse_template_row(dict(row)) for row in reader]


def _payload_rows(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [dict(item) for item in payload]
    if not isinstance(payload, dict):
        raise SourceFetchError("Unexpected RetroRules API payload")

    for key in ("results", "items", "data", "templates"):
        value = payload.get(key)
        if isinstance(value, list):
            return [dict(item) for item in value]

    if any(key in payload for key in ("TEMPLATE_ID", "template_id", "id")):
        return [dict(payload)]

    raise SourceFetchError(
        "RetroRules API response shape was not recognized; upstream API may have changed"
    )


class RetroRulesSource:
    """Primary-source adapter for the documented RetroRules template API."""

    def build_search_url(
        self,
        *,
        ec: str | None = None,
        radius: int | None = None,
    ) -> str:
        params: dict[str, str | int] = {}
        if ec:
            params["ec"] = ec
        if radius is not None:
            if not 0 <= radius <= 10:
                raise ValueError("RetroRules radius must be between 0 and 10")
            params["radius"] = radius

        query = urlencode(params)
        return TEMPLATES_API if not query else f"{TEMPLATES_API}?{query}"

    def template_summary_url(self, template_id: str) -> str:
        if not template_id:
            raise ValueError("template_id must not be empty")
        from urllib.parse import quote
        return f"{TEMPLATES_API}/{quote(template_id, safe='')}/summary"

    def search(
        self,
        *,
        ec: str | None = None,
        radius: int | None = None,
    ) -> list[ReactionRuleRecord]:
        payload = fetch_json(self.build_search_url(ec=ec, radius=radius))
        return [parse_template_row(row) for row in _payload_rows(payload)]

    def summary(self, template_id: str) -> dict[str, Any]:
        payload = fetch_json(self.template_summary_url(template_id))
        if not isinstance(payload, dict):
            raise SourceFetchError("Unexpected RetroRules template summary payload")
        return payload
