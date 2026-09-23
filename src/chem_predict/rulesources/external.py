from __future__ import annotations

from typing import Any

from chem_predict.rulesources.http import fetch_json


class EnviPathSource:
    """Thin external provider for explicitly supplied enviPath API URLs.

    Chem-Predict does not mirror public enviPath knowledge packages because
    important packages are licensed CC BY-NC-SA 4.0.
    """

    def fetch_rule(self, rule_url: str) -> dict[str, Any]:
        if not rule_url.startswith(("https://envipath.org/", "https://beta.envipath.org/")):
            raise ValueError("Expected an enviPath rule URL")
        payload = fetch_json(rule_url)
        if not isinstance(payload, dict):
            raise ValueError("Unexpected enviPath rule payload")
        return payload
