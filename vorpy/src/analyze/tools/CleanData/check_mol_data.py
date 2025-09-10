import os
import sys
from pathlib import Path
import tkinter as tk
from tkinter import filedialog
from typing import Optional, Union



# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[5]
sys.path.append(str(project_root))



from vorpy.src.inputs.pdb import read_pdb  # noqa: F401 (kept in case other parts rely on side effects)
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.system.system import System



"""
Checks a folder with aw, pow, prm folders and compares the number of complete cells
to the number of group atoms in the PDB.

The PDB filename must match the folder's system name (folder name minus the leading letter and underscore).
Example: folder 'A_system' -> looks for 'system.pdb' in the folder.
"""



def _extract_system_name(folder: Path) -> str:
    """Derive system name from folder per convention '<letter>_<system_name>'."""
    base = folder.name
    parts = base.split('_', 1)
    if len(parts) == 2 and parts[1]:
        return parts[1].lower()
    # Fallback: if pattern not met, use full name
    return base.lower()



def _find_pdb_file(folder: Path, system_name: str) -> Optional[Path]:
    """Return the PDB path in folder whose stem matches system_name (case-insensitive)."""
    for p in folder.iterdir():
        if p.suffix.lower() == ".pdb" and p.stem.lower() == system_name:
            return p
    return None



def _load_complete_indices(csv_path: Path, **read_kwargs) -> list[int]:
    """Load indices of atoms with Complete Cell? == True from a logs CSV."""
    logs = read_logs2(str(csv_path), **read_kwargs)
    atoms = logs["atoms"]

    # Prefer column-based 'Index' if present; otherwise use the dataframe index.
    if "Index" in atoms.columns:
        return atoms.loc[atoms["Complete Cell?"], "Index"].to_list()

    return atoms.index[atoms["Complete Cell?"]].to_list()


def _missing_in_order(population: list[int], have: set[int]) -> list[int]:
    """Return items from population missing in 'have', preserving original order."""
    return [a for a in population if a not in have]


def check_mol_data(folder: Optional[Union[str, os.PathLike]] = None) -> None:
    # Choose folder (GUI prompt if not provided)
    if folder is None:
        folder_chosen = filedialog.askdirectory()
        if not folder_chosen:
            print("No folder selected.")
            return
        folder_path = Path(folder_chosen)
    else:
        folder_path = Path(folder)

    if not folder_path.exists():
        print(f"Folder not found: {folder_path}")
        return

    # Derive system name from folder and find matching PDB
    system_name = _extract_system_name(folder_path)
    print(system_name)
    pdb_path = _find_pdb_file(folder_path, system_name)
    if pdb_path is None:
        print("No pdb file found for", system_name)
        return

    # Build the System and compute atom numbers (exclude solvent atoms) — vectorized
    system = System(str(pdb_path))
    balls = system.balls  # DataFrame with 'num' column
    if "num" not in balls.columns:
        print("Expected 'num' column in system.balls but it was not found.")
        return

    mask = ~balls["num"].isin(system.sol.atoms)
    atom_nums = balls.loc[mask, "num"].to_list()
    num_atoms = len(atom_nums)

    # Map of scheme -> (csv path, kwargs)
    files = {
        "aw":  (folder_path / "aw" / "aw_logs.csv",  dict()),
        "pow": (folder_path / "pow" / "pow_logs.csv", dict(all_=False, balls=True)),
        "prm": (folder_path / "prm" / "prm_logs.csv", dict(all_=False, balls=True)),
    }

    complete = {}
    for key, (csv_path, kwargs) in files.items():
        if not csv_path.exists():
            print(f"No {key} logs found for", system_name)
            return
        complete[key] = _load_complete_indices(csv_path, **kwargs)

    # Print summary
    print("Number of atoms:", num_atoms)

    for key in ("aw", "pow", "prm"):
        comp_list = complete[key]
        comp_set = set(comp_list)
        missing = _missing_in_order(atom_nums, comp_set)
        print(
            f"Number of complete {key} cells: {len(comp_list)}"
            f" - Missing atoms: {missing}"
        )



if __name__ == "__main__":
    check_mol_data()
