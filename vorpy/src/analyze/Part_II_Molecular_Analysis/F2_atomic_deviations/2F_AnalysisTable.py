import os
import sys
import tkinter as tk
from tkinter import filedialog
from typing import List

import pandas as pd


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


# =========================================================
# Residue / element classification helpers
# =========================================================

PROTEIN_RES = {
    'ALA', 'ARG', 'THR', 'LYS', 'GLN', 'SER', 'GLY', 'PRO', 'LEU',
    'VAL', 'HIS', 'TYR', 'GLU', 'ILE', 'PHE', 'ASP', 'MET', 'ASN', 'CYS'
}

DNA_RES = {'DA', 'DC', 'DG', 'DT'}
RNA_RES = {'U'}

TWO_LETTER_ELEMENTS = {
    "NA", "CL", "MG", "FE", "ZN", "CA", "CU", "MN", "CO",
    "SE", "BR", "LI", "KR", "XE", "CS", "BA", "SR", "AG",
    "PT", "AU", "PD", "RU", "RH", "IR", "NI"
}


def infer_element(atom_name: str) -> str:
    """
    Infer element symbol from atom Name using a PDB-like heuristic.
    """
    if not isinstance(atom_name, str):
        return "UNK"

    clean = atom_name.strip()

    # Remove leading digits
    while clean and clean[0].isdigit():
        clean = clean[1:]

    if not clean:
        return "UNK"

    c0 = clean[0].upper()
    c1 = clean[1].lower() if len(clean) > 1 else ""

    # Try two-letter element first
    if len(clean) > 1:
        candidate_two = (c0 + c1).upper()
        if candidate_two in TWO_LETTER_ELEMENTS:
            return candidate_two

    return c0


def classify_atom(res_name: str, atom_name: str) -> str:
    """
    Map (Residue, Name) to a generalized atom-type category.
    """
    res = str(res_name).strip()
    name = str(atom_name).strip()

    element = infer_element(name)

    # ---------------- protein residues ----------------
    if res in PROTEIN_RES:
        # backbone N & attached Hs
        if name == 'N':
            return 'Backbone N'
        # backbone CA
        if name == 'CA':
            return 'Backbone CA'
        # backbone C
        if name == 'C':
            return 'Backbone C'
        # backbone O / carbonyl O
        if name in {'O', 'OC1', 'OC2'}:
            return 'Backbone O'
        # Backbone H
        if name in {'H', 'H1', 'H2', 'H3'}:
            return 'Backbone H'

        # side chain atoms: grouped by element
        if element == 'C':
            return 'Sidechain C'
        if element == 'N':
            return 'Sidechain N'
        if element == 'O':
            return 'Sidechain O'
        if element == 'S':
            return 'Sidechain S'
        if element == 'H':
            return 'Sidechain H'
        return 'Sidechain Other'

    # ---------------- nucleic acids (DNA/RNA) ----------------
    if res in DNA_RES or res in RNA_RES:
        # phosphate group
        if name == 'P':
            return 'Phosphate P'
        if name in {'O1P', 'O2P', 'O2P'}:
            return 'Phosphate O'

        # sugar carbons/oxygens (contain apostrophe)
        if "'" in name:
            if name.startswith('C'):
                return 'Sugar C'
            if name.startswith('O'):
                return 'Sugar O'

        # base atoms grouped by element
        if element == 'C':
            return 'Base C'
        if element == 'N':
            return 'Base N'
        if element == 'O':
            return 'Base O'
        if element == 'H':
            return 'H'
        return 'Other'

    # ---------------- catch-all ----------------
    return 'Other'


def neighbor_count(value) -> int:
    """
    Compute neighbor count from the 'Neighbors' column.

    'Neighbors' should be a list of ints or a list of lists of ints.
    We handle both by flattening one level if needed.
    """
    if value is None:
        return 0
    if isinstance(value, list):
        if len(value) > 0 and isinstance(value[0], list):
            # Flatten one level
            return sum(len(sub) for sub in value)
        return len(value)
    # If it's not a list, treat as 0 (or you can try parsing).
    return 0


# =========================================================
# Data loading and aggregation
# =========================================================

def load_atoms_from_logs(log_files: List[str]) -> pd.DataFrame:
    """
    Load atoms tables from multiple VorPy log files and concatenate.
    """
    all_atoms = []

    for path in log_files:
        info = read_logs2(path)        # one_file=True behavior: dict for this file
        atoms = info['atoms'].copy()
        atoms['LogFile'] = os.path.basename(path)
        all_atoms.append(atoms)

    if not all_atoms:
        raise ValueError("No log files provided.")

    atoms_df = pd.concat(all_atoms, ignore_index=True)
    return atoms_df


def summarize_by_atom_type(atoms_df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify each atom, create a neighbor-count column, then aggregate
    geometry per atom-type category.

    Output columns:
      Name (category)
      Count
      Radius (mean only)
      Avg Vol (mean +/- SD)
      Avg VdW Vol (mean +/- SD)
      Avg Surf Area (mean +/- SD)
      Avg Ψ (mean +/- SD)
      Avg Nbors (mean +/- SD)
    """
    # Ensure required columns exist
    required_cols = [
        'Name', 'Residue', 'Radius', 'Volume', 'Van Der Waals Volume',
        'Surface Area', 'Sphericity', 'Neighbors'
    ]
    missing = [c for c in required_cols if c not in atoms_df.columns]
    if missing:
        raise KeyError(f"Missing required columns in atoms dataframe: {missing}")

    atoms_df = atoms_df.copy()

    # Classify atom types
    atoms_df['AtomType'] = atoms_df.apply(
        lambda row: classify_atom(row['Residue'], row['Name']),
        axis=1
    )

    # Compute neighbor counts
    atoms_df['Neighbor Count'] = atoms_df['Neighbors'].apply(neighbor_count)

    # Format "mean +/- SD" using ASCII symbols
    def fmt(series: pd.Series, decimals: int = 2) -> str:
        n = series.count()
        if n == 0:
            return "NA"
        mean = float(series.mean())
        sd = float(series.std(ddof=1)) if n > 1 else 0.0
        return f"{mean:.{decimals}f} +/- {sd:.{decimals}f}"

    grouped = atoms_df.groupby('AtomType')

    records = []
    for atom_type, g in grouped:
        n = len(g)
        if n == 0:
            continue

        record = {
            'Name': atom_type,
            'Count': n,

            # Radius: mean only (NO SD)
            'Radius': float(g['Radius'].mean()),

            # All others: mean +/- SD
            'Avg Vol':       fmt(g['Volume'], decimals=2),
            'Avg VdW Vol':   fmt(g['Van Der Waals Volume'], decimals=2),
            'Avg Surf Area': fmt(g['Surface Area'], decimals=2),
            'Avg \u03A8':    fmt(g['Sphericity'], decimals=3),
            'Avg Nbors':     fmt(g['Neighbor Count'], decimals=2),
        }
        records.append(record)

    summary = pd.DataFrame.from_records(records)
    summary = summary.sort_values(by='Name')

    return summary




# =========================================================
# GUI main: select logs, produce CSV
# =========================================================

def main():
    print("Select one or more VorPy log files (AW or Power).")
    print("You can run this script separately for AW and Power if desired.\n")

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    file_paths = []
    while True:
        file_path = filedialog.askopenfilename(
            title="Select VorPy log files",
            filetypes=[("Log files", "*.csv"), ("All files", "*.*")]
        )
        if file_path:
            file_paths.append(file_path)
        else:
            break

    if not file_paths:
        print("No files selected. Exiting.")
        return

    atoms_df = load_atoms_from_logs(list(file_paths))
    summary_df = summarize_by_atom_type(atoms_df)
    output_folder = filedialog.askdirectory(title="Select Output Folder")
    out_csv = output_folder + "/atom_type_summary.csv"
    summary_df.to_csv(out_csv, index=False)
    print(f"\nSaved atom-type summary with SD to: {out_csv}\n")


if __name__ == "__main__":
    main()
