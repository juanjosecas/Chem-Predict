from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from chem_predict.rulesources.models import SourceFetchError


USER_AGENT = "Chem-Predict/0.1 (+https://github.com/juanjosecas/Chem-Predict)"


def fetch_bytes(url: str, *, timeout: float = 60.0) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read()
    except (HTTPError, URLError, TimeoutError) as exc:
        raise SourceFetchError(f"Could not fetch {url}: {exc}") from exc


def fetch_text(
    url: str,
    *,
    timeout: float = 60.0,
    encoding: str = "utf-8",
) -> str:
    return fetch_bytes(url, timeout=timeout).decode(encoding)


def fetch_json(url: str, *, timeout: float = 60.0) -> Any:
    try:
        return json.loads(fetch_text(url, timeout=timeout))
    except json.JSONDecodeError as exc:
        raise SourceFetchError(f"Source did not return valid JSON: {url}") from exc


def download_to_path(
    url: str,
    path: str | Path,
    *,
    timeout: float = 120.0,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=timeout) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)
    except (HTTPError, URLError, TimeoutError) as exc:
        destination.unlink(missing_ok=True)
        raise SourceFetchError(f"Could not download {url}: {exc}") from exc
    return destination
