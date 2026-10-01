"""Execute the example in a fresh kernel without changing tracked outputs."""
from pathlib import Path
from tempfile import TemporaryDirectory

import nbformat
from nbclient import NotebookClient


def main() -> None:
    path = Path(__file__).resolve().parents[1] / "notebooks" / "ChemPredict_reacciones_completo.ipynb"
    notebook = nbformat.read(path, as_version=4)
    nbformat.validate(notebook)
    with TemporaryDirectory(prefix="chem-predict-notebook-") as workdir:
        NotebookClient(
            notebook, timeout=180, kernel_name="python3",
            resources={"metadata": {"path": workdir}},
        ).execute()
    count = sum(cell.cell_type == "code" for cell in notebook.cells)
    print(f"Notebook executed successfully: {count} code cells")


if __name__ == "__main__":
    main()
