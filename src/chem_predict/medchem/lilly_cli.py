from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class LillyCLIResult:
    returncode: int
    stdout: str
    stderr: str
    command: tuple[str, ...]


class LillyOfficialBackend:
    """Adapter for a local checkout of the official Lilly-Medchem-Rules 2.1.0.

    This backend is for exact upstream comparison and requires its LillyMol
    executables to be installed as documented by the upstream project.
    """

    def __init__(self, checkout: str | Path, *, python: str = "python") -> None:
        self.checkout = Path(checkout)
        self.python = python
        self.script = self.checkout / "Lilly_Medchem_Rules.py"
        if not self.script.is_file():
            raise FileNotFoundError(f"Missing Lilly_Medchem_Rules.py in {self.checkout}")

    def run_smiles(self, smiles: str, *, molecule_id: str = "query") -> LillyCLIResult:
        with tempfile.TemporaryDirectory(prefix="chem_predict_lilly_") as tmp:
            input_path = Path(tmp) / "input.smi"
            input_path.write_text(f"{smiles} {molecule_id}\n", encoding="utf-8")
            cmd = (
                self.python,
                str(self.script),
                "-nobadfiles",
                "-tabular",
                str(input_path),
            )
            completed = subprocess.run(
                cmd,
                cwd=self.checkout,
                text=True,
                capture_output=True,
                check=False,
            )
            return LillyCLIResult(
                returncode=completed.returncode,
                stdout=completed.stdout,
                stderr=completed.stderr,
                command=cmd,
            )
