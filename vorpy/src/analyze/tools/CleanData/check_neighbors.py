import os
import sys
import re
import ast
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Optional

import numpy as np
import pandas as pd



# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



@dataclass
class Config:
    output_dir: str = "neighbor_validation"
    max_examples_to_print: int = 25

    # Surface pair schema:
    # If your surfs_df has a single column with a pair (e.g. "[12, 34]"), set pair_col to that name.
    # Otherwise set i_col and j_col to the two columns that hold the indices.
    pair_col: str = "bALLS"
    i_col: str = "i"
    j_col: str = "j"

    # Atoms schema:
    index_col: str = "Index"
    neighbors_col: str = "neighbors"
    drop_self: bool = True



def choose_logs_file(initial_dir: Optional[str] = None) -> str:
    root = tk.Tk()
    root.withdraw()

    path = filedialog.askopenfilename(
        title="Select a logs CSV file",
        initialdir=initial_dir or os.getcwd(),
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )

    root.destroy()

    if not path:
        raise FileNotFoundError("No logs file selected.")

    return path



def sanitize_filename(s: str) -> str:
    keep = []
    for ch in str(s):
        if ch.isalnum() or ch in ("-", "_"):
            keep.append(ch)
        elif ch.isspace():
            keep.append("_")

    out = "".join(keep).strip("_")
    return out if out else "output"



# ------------------------------------------------------------
# 1) Add surface indices to all atom indices (dictionary of lists)
# ------------------------------------------------------------
def build_surface_neighbor_map(
    surfs_df: pd.DataFrame,
    *,
    pair_col: str,
    i_col: str,
    j_col: str,
) -> Dict[int, List[int]]:
    neigh: Dict[int, Set[int]] = {}

    # Case-insensitive column resolution
    col_map = {c.lower(): c for c in surfs_df.columns}
    pair_col_real = col_map.get(pair_col.lower())
    i_col_real = col_map.get(i_col.lower())
    j_col_real = col_map.get(j_col.lower())

    if pair_col_real is not None:
        for v in surfs_df[pair_col_real].values:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue

            if isinstance(v, (list, tuple, np.ndarray)):
                if len(v) < 2:
                    continue
                a = int(v[0])
                b = int(v[1])
            else:
                nums = re.findall(r"-?\d+", str(v))
                if len(nums) < 2:
                    continue
                a = int(nums[0])
                b = int(nums[1])

            neigh.setdefault(a, set()).add(b)
            neigh.setdefault(b, set()).add(a)

        return {k: sorted(list(v)) for k, v in neigh.items()}

    if i_col_real is None or j_col_real is None:
        raise KeyError(
            f"surfs_df must contain '{pair_col}' OR both '{i_col}' and '{j_col}'. "
            f"Columns present: {list(surfs_df.columns)}"
        )

    a_vals = pd.to_numeric(surfs_df[i_col_real], errors="coerce").to_numpy()
    b_vals = pd.to_numeric(surfs_df[j_col_real], errors="coerce").to_numpy()

    for a, b in zip(a_vals, b_vals):
        if np.isnan(a) or np.isnan(b):
            continue

        ai = int(a)
        bi = int(b)

        neigh.setdefault(ai, set()).add(bi)
        neigh.setdefault(bi, set()).add(ai)

    return {k: sorted(list(v)) for k, v in neigh.items()}





# ------------------------------------------------------------
# 2) Get the neighbors lists (dictionary of lists)
# ------------------------------------------------------------
def build_neighbors_map_from_atoms(
    atoms_df: pd.DataFrame,
    *,
    index_col: str,
    neighbors_col: str,
    drop_self: bool,
) -> Dict[int, List[int]]:
    """
    Returns: dict atom_index -> neighbor list from atoms_df[neighbors_col].

    Parsing:
      - if cell is list/tuple/np array: use it
      - else: extract all integers via regex (robust to "array([..])", "[..]", etc.)
    """
    if index_col not in atoms_df.columns:
        raise KeyError(f"atoms_df missing required column '{index_col}'")
    if neighbors_col not in atoms_df.columns:
        raise KeyError(f"atoms_df missing required column '{neighbors_col}'")

    idx = pd.to_numeric(atoms_df[index_col], errors="coerce")
    df = atoms_df.loc[~idx.isna(), [index_col, neighbors_col]].copy()
    df[index_col] = pd.to_numeric(df[index_col], errors="coerce").astype(int)

    out: Dict[int, List[int]] = {}

    for i, v in df.itertuples(index=False):
        if isinstance(v, (list, tuple, np.ndarray)):
            nbs = []
            for x in v:
                try:
                    nbs.append(int(x))
                except Exception:
                    pass
        else:
            nums = re.findall(r"-?\d+", str(v))
            nbs = [int(x) for x in nums]

        if drop_self:
            nbs = [x for x in nbs if x != int(i)]

        # De-dup preserving order
        seen = set()
        uniq = []
        for x in nbs:
            if x not in seen:
                uniq.append(x)
                seen.add(x)

        out[int(i)] = uniq

    return out


def infer_neighbors_column(atoms_df: pd.DataFrame) -> str:
    cols = list(atoms_df.columns)
    low = {c.lower().strip(): c for c in cols}

    # best guesses in priority order
    candidates = [
        "neighbors",
        "neighbor",
        "neighbours",
        "neighbors_list",
        "neighbor_list",
        "neighbor indices",
        "neighbor_indices",
        "neighbors indices",
        "neighbors_indices",
        "nbors",
        "nbors_list",
    ]

    for key in candidates:
        if key in low:
            return low[key]

    # fallback: any column containing the word 'neighbor'
    for c in cols:
        if "neighbor" in c.lower():
            return c

    raise KeyError(
        "Could not find a neighbors list column in atoms_df. "
        f"Available columns: {cols}"
    )



# ------------------------------------------------------------
# 3) Compare and print the results (including missing indices)
# ------------------------------------------------------------
def compare_and_print_neighbors(
    surf_map: Dict[int, List[int]],
    atom_map: Dict[int, List[int]],
    *,
    max_examples: int,
) -> pd.DataFrame:
    rows = []

    keys = sorted(set(surf_map.keys()) | set(atom_map.keys()))

    for i in keys:
        s = set(surf_map.get(i, []))
        a = set(atom_map.get(i, []))

        missing = sorted(list(s - a))
        extra = sorted(list(a - s))

        inter = len(s & a)
        union = len(s | a)
        jacc = float(inter / union) if union > 0 else float("nan")

        rows.append(
            {
                "Index": int(i),
                "n_surfs": int(len(s)),
                "n_atoms": int(len(a)),
                "n_missing": int(len(missing)),
                "n_extra": int(len(extra)),
                "jaccard": jacc,
                "missing": missing,
                "extra": extra,
            }
        )

    report = pd.DataFrame(rows)

    tested = report[(report["n_surfs"] > 0) | (report["n_atoms"] > 0)].copy()

    mean_j = float(tested["jaccard"].mean()) if len(tested) else float("nan")
    med_j = float(tested["jaccard"].median()) if len(tested) else float("nan")
    perfect = float((tested["jaccard"] == 1.0).mean() * 100.0) if len(tested) else float("nan")

    any_missing = float((tested["n_missing"] > 0).mean() * 100.0) if len(tested) else float("nan")
    any_extra = float((tested["n_extra"] > 0).mean() * 100.0) if len(tested) else float("nan")

    print("\nNeighbor validation summary")
    print(f"  Atoms compared : {len(tested)}")
    print(f"  Mean Jaccard   : {mean_j:.4f}")
    print(f"  Median Jaccard : {med_j:.4f}")
    print(f"  Perfect match  : {perfect:.1f}%")
    print(f"  Any missing    : {any_missing:.1f}%")
    print(f"  Any extra      : {any_extra:.1f}%")

    worst = tested.sort_values(
        by=["n_missing", "n_extra", "jaccard"],
        ascending=[False, False, True],
    ).head(max_examples)
    atom_indices = set(atom_map.keys())

    print("\nWorst mismatches (surface missing, atoms extra):")
    printed = 0
    for _, r in worst.iterrows():
        if int(r["n_missing"]) == 0 and int(r["n_extra"]) == 0:
            break
        idx = int(r["Index"])

        if idx not in atom_indices:
            continue

        print(
            f"  Index {int(r['Index'])}: "
            f"missing={r['missing']} "
            f"extra={r['extra']} "
            f"jacc={float(r['jaccard']):.3f}"
        )
        printed += 1

    if printed == 0:
        print("  None (all matched for the worst subset).")

    return report


def main() -> None:
    cfg = Config()

    logs_file = choose_logs_file()

    logs_obj = read_logs2(
        logs_file,
        return_dict=False,
        all_=True,
        balls=True,
        surfs=True,
        edges=False,
        verts=False,
    )

    atoms_df = logs_obj["atoms"]
    print("atoms columns:", list(atoms_df.columns))
    surfs_df = logs_obj["surfs"]

    if surfs_df is None or len(surfs_df) == 0:
        raise ValueError("No surfaces found in logs (surfs_df is empty). Cannot validate via surfaces.")

    print("Loaded:")
    print(f"  atoms rows: {len(atoms_df)}")
    print(f"  surfs rows: {len(surfs_df)}")
    print(f"  surfs columns: {list(surfs_df.columns)}")

    surf_map = build_surface_neighbor_map(
        surfs_df,
        pair_col=cfg.pair_col,
        i_col=cfg.i_col,
        j_col=cfg.j_col,
    )

    neighbors_col = infer_neighbors_column(atoms_df)

    atom_map = build_neighbors_map_from_atoms(
        atoms_df,
        index_col=cfg.index_col,
        neighbors_col=neighbors_col,
        drop_self=cfg.drop_self,
    )

    report = compare_and_print_neighbors(
        surf_map,
        atom_map,
        max_examples=cfg.max_examples_to_print,
    )

    mol_name = str(
        logs_obj.get("data", {}).get("name", os.path.splitext(os.path.basename(logs_file))[0])
    )
    out_prefix = sanitize_filename(mol_name)

    os.makedirs(cfg.output_dir, exist_ok=True)
    out_csv = os.path.join(cfg.output_dir, f"{out_prefix}_neighbors_vs_surfaces.csv")
    report.to_csv(out_csv, index=False)

    print("\nSaved report:")
    print(f"  {out_csv}")



if __name__ == "__main__":
    main()
