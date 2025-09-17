import tkinter as tk
from tkinter import filedialog
import csv
import os
import sys
import numpy as np

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.system.system import System


def get_correct_SA(logs, group_atoms):
    """
    Sum 'Surface Area' for rows whose two 'Balls' are NOT both inside group_atoms.

    Speeds up by vectorizing:
      - If group_atoms are sorted & consecutive: O(n) range checks
      - Else: uses np.isin for membership (still vectorized)

    Parameters
    ----------
    logs : dict-like with key 'surfs' -> pandas.DataFrame
        DataFrame must have:
          - 'Balls' column: iterable of 2 ints (e.g., tuples)
          - 'Surface Area' column: numeric
    group_atoms : Sequence[int]
        Atom indices (may or may not be sorted/consecutive)

    Returns
    -------
    float
        Sum of surface areas where NOT (both balls in group_atoms).
    """
    surfs = logs['surfs']

    # Extract as NumPy for speed
    # balls: shape (n, 2), dtype int
    balls = np.asarray(surfs['Balls'].to_list(), dtype=int)
    sa_vals = surfs['Surface Area'].to_numpy(dtype=float, copy=False)

    # Empty group: nothing is "inside", so we sum everything
    if not group_atoms:
        return float(sa_vals.sum())

    arr = np.asarray(group_atoms, dtype=int)

    # Check strictly increasing by 1 and uniqueness
    arr_unique = np.unique(arr)
    is_consecutive = (
        arr_unique.size == arr.size and               # no duplicates
        arr_unique.size >= 1 and
        np.all(np.diff(arr_unique) == 1)              # consecutive steps of 1
    )

    if is_consecutive:
        gmin = int(arr_unique[0])
        gmax = int(arr_unique[-1])
        # Fast range check when consecutive
        both_in = (
            (balls[:, 0] >= gmin) & (balls[:, 0] <= gmax) &
            (balls[:, 1] >= gmin) & (balls[:, 1] <= gmax)
        )
    else:
        # General membership (still vectorized)
        # np.isin against the unique array is faster than Python set in a loop here
        in0 = np.isin(balls[:, 0], arr_unique, assume_unique=False)
        in1 = np.isin(balls[:, 1], arr_unique, assume_unique=False)
        both_in = in0 & in1

    # We want SA where NOT both balls are in the group
    return float(sa_vals[~both_in].sum())


import csv
from pathlib import Path


def replace_csv_cell_by_index(input_path, output_path, row_index, col_index, new_value, old_value, encoding="utf-8"):

    """
    Replace a single CSV cell (zero-based row/col) after verifying its current value.

    Parameters
    ----------
    input_path : str | Path
        Path to the input CSV file.
    output_path : str | Path
        Path to write the updated CSV file (can be the same as input_path for in-place update).
    row_index : int
        Zero-based row index of the target cell (includes header if present).
    col_index : int
        Zero-based column index of the target cell.
    new_value : Any
        The value to write into the cell.
    old_value : Any
        The value you expect currently in the cell; used as a safety check.
    encoding : str
        File encoding (default "utf-8").

    Raises
    ------
    IndexError
        If row_index or col_index is out of bounds.
    ValueError
        If the current value at (row_index, col_index) does not match old_value.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)

    # Read the whole file, preserving dialect if possible
    with input_path.open("r", newline="", encoding=encoding) as f:
        # Sniff the CSV dialect to preserve delimiter/quoting
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except csv.Error:
            dialect = csv.excel

        reader = csv.reader(f, dialect)
        rows = list(reader)

    if not rows:
        raise IndexError("CSV is empty; no rows to modify.")

    if row_index < 0 or row_index >= len(rows):
        raise IndexError(
            f"Row index {row_index} out of range (0..{len(rows)-1})."
        )

    if col_index < 0 or col_index >= len(rows[row_index]):
        raise IndexError(
            f"Column index {col_index} out of range for row {row_index} "
            f"(0..{len(rows[row_index]) - 1})."
        )

    current_value = rows[row_index][col_index]

    # Safety check to ensure we're editing the intended cell
    if str(current_value) != str(old_value):
        raise ValueError(
            f"Safety check failed at ({row_index}, {col_index}): "
            f"found '{current_value}', expected '{old_value}'. "
            f"No changes written."
        )

    # Perform replacement
    rows[row_index][col_index] = float(round(new_value, 2))

    # Write updated CSV
    with output_path.open("w", newline="", encoding=encoding) as f:
        writer = csv.writer(f, dialect)
        writer.writerows(rows)


def correct_logs_SA(folder=None):
    # Check if the folder is None
    if folder is None:
        # Get the dropbox folder
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()
    # Get the aw logs
    aw_logs_file = os.path.join(folder, 'aw/aw_logs.csv')
    # Get the pow logs
    pow_logs_file = os.path.join(folder, 'pow/pow_logs.csv')
    # Get the prm logs
    prm_logs_file = os.path.join(folder, 'prm/prm_logs.csv')
    # Read the logs
    aw_logs = read_logs2(aw_logs_file, all_=False, balls=True, surfs=True)
    pow_logs = read_logs2(pow_logs_file, all_=False, balls=True, surfs=True)
    prm_logs = read_logs2(prm_logs_file, all_=False, balls=True, surfs=True)
    # Get the subfolder
    subfolder = os.path.basename(folder)
    # Get the pdb path
    pdb_path = os.path.join(folder, subfolder[2:] + '.pdb')

    # Check if the pdb file is a real path
    if not os.path.exists(pdb_path):
        pdb_path = os.path.join(folder, subfolder + '.pdb')
        if not os.path.exists(pdb_path):
            print(f"PDB file not found for {subfolder}")
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

    # Get the correct surface area
    aw_sa = get_correct_SA(aw_logs, atom_nums)
    pow_sa = get_correct_SA(pow_logs, atom_nums)
    prm_sa = get_correct_SA(prm_logs, atom_nums)

    # Print the surface area difference
    print(f"\nSA before: aw {aw_logs['group data']['Surface Area']}, pow {pow_logs['group data']['Surface Area']}, prm {prm_logs['group data']['Surface Area']}")
    print(f"SA after: aw {aw_sa}, pow {pow_sa}, prm {prm_sa}\n")

    # Replace the logs SA
    replace_csv_cell_by_index(aw_logs_file, os.path.join(folder, 'aw/aw_logs_new.csv'), 5, 2, aw_sa, aw_logs['group data']['Surface Area'])
    replace_csv_cell_by_index(pow_logs_file, os.path.join(folder, 'pow/pow_logs_new.csv'), 5, 2, pow_sa, pow_logs['group data']['Surface Area'])
    replace_csv_cell_by_index(prm_logs_file, os.path.join(folder, 'prm/prm_logs_new.csv'), 5, 2, prm_sa, prm_logs['group data']['Surface Area'])


def correct_multiple_logs_SA(folder=None):
    # Check if the folder is None
    if folder is None:
        # Get the dropbox folder
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()
    # Loop through the folder and get the logs
    for subfolder in os.listdir(folder):
        print('\n', subfolder)
        # Correct the logs SA
        correct_logs_SA(os.path.join(folder, subfolder))

if __name__ == '__main__':
    correct_multiple_logs_SA()