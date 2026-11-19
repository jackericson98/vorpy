import os
import sys
import ast
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Optional, Set, Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



NUC_RES_NAMES: Set[str] = {
    # DNA
    "DA", "DG", "DC", "DT",
    # RNA
    "A", "G", "C", "U",
    # Some logs may store thymine as "T"
    "T",
}



@dataclass
class Config:
    curvature_kind: str = "mean"        # "mean" or "gauss"
    use_magnitude: bool = True          # |curvature|
    weight_by_area: bool = True         # residue curvature = area-weighted across atoms
    residue_name: str = "ALA"           # residue to analyze
    output_dir: str = "alanine_dna_facing"



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



def parse_neighbors_cell(x) -> List[int]:
    """
    Robust parser for the atoms_df['neighbors'] field.
    Handles:
      - actual Python lists
      - strings like "[1, 2, 3]"
      - strings like "1,2,3"
      - NaN / empty
    """
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []

    if isinstance(x, (list, tuple, np.ndarray)):
        return [int(v) for v in x if str(v).strip() != ""]

    s = str(x).strip()
    if not s:
        return []

    # Try literal_eval first (safe for list-like strings)
    try:
        v = ast.literal_eval(s)
        if isinstance(v, (list, tuple, np.ndarray)):
            return [int(z) for z in v]
        if isinstance(v, (int, float)) and not np.isnan(v):
            return [int(v)]
    except Exception:
        pass

    # Fallback: split by comma
    s = s.strip("[](){}")
    if not s:
        return []

    parts = [p.strip() for p in s.split(",") if p.strip()]
    out = []
    for p in parts:
        try:
            out.append(int(float(p)))
        except Exception:
            continue

    return out



def build_residue_instances_from_atoms(atoms_df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """
    One row per residue instance (Chain + Residue Sequence + Residue).
    Computes area-weighted curvature across atoms in the residue.

    Required columns in atoms_df:
      Index, Residue, Residue Sequence, Chain, Surface Area,
      Average Mean Surface Curvature, Average Gaussian Surface Curvature,
      neighbors
    """
    required = [
        "Index",
        "Residue",
        "Residue Sequence",
        "Chain",
        "Surface Area",
        "Average Mean Surface Curvature",
        "Average Gaussian Surface Curvature",
        "Neighbors",
    ]
    missing = [c for c in required if c not in atoms_df.columns]
    if missing:
        raise KeyError(f"Missing required columns in atoms dataframe: {missing}")

    df = atoms_df.copy()

    df["Index"] = pd.to_numeric(df["Index"], errors="coerce")
    df["Residue Sequence"] = pd.to_numeric(df["Residue Sequence"], errors="coerce")
    df["Surface Area"] = pd.to_numeric(df["Surface Area"], errors="coerce")
    df["Average Mean Surface Curvature"] = pd.to_numeric(df["Average Mean Surface Curvature"], errors="coerce")
    df["Average Gaussian Surface Curvature"] = pd.to_numeric(df["Average Gaussian Surface Curvature"], errors="coerce")

    df = df.dropna(subset=["Index", "Residue Sequence"]).copy()
    df["Index"] = df["Index"].astype(int)
    df["Residue Sequence"] = df["Residue Sequence"].astype(int)

    if cfg.curvature_kind.lower().startswith("g"):
        curv_col = "Average Gaussian Surface Curvature"
        curv_label = "Gaussian"
    else:
        curv_col = "Average Mean Surface Curvature"
        curv_label = "Mean"

    df["_curv"] = df[curv_col].astype(float)

    if cfg.use_magnitude:
        df["_curv"] = np.abs(df["_curv"])
        curv_label = f"|{curv_label}|"

    df["_area"] = pd.to_numeric(df["Surface Area"], errors="coerce").fillna(0.0).astype(float)
    df["_area"] = np.clip(df["_area"].to_numpy(dtype=float), 0.0, None)

    # Parse neighbors into Python lists for downstream use
    df["_neighbors"] = df["Neighbors"].apply(parse_neighbors_cell)

    # Compute per-residue instance curvature (area-weighted across atoms)
    df["_wy"] = df["_area"] * df["_curv"]

    gcols = ["Chain", "Residue Sequence", "Residue"]

    agg = (
        df.groupby(gcols, as_index=False)
        .agg(
            wy_sum=("_wy", "sum"),
            area_sum=("_area", "sum"),
            curv_fallback=("_curv", "mean"),
            n_atoms=("_curv", "size"),
        )
    )

    agg["curv"] = np.where(
        agg["area_sum"] > 0.0,
        agg["wy_sum"] / agg["area_sum"],
        agg["curv_fallback"],
    )

    agg.attrs["curv_label"] = f"{curv_label} residue curvature (area-weighted)"
    agg.attrs["curv_col_used"] = curv_col

    # Keep atom-level table too (needed for neighbor-based labeling)
    df.attrs["curv_col_used"] = curv_col

    return agg, df



def label_residue_dna_facing(
    residue_instances: pd.DataFrame,
    atoms_df: pd.DataFrame,
    cfg: Config,
) -> pd.DataFrame:
    """
    Labels each residue instance as DNA-facing if any of its atoms has a neighbor atom
    whose residue name is nucleic-acid-like (DA/DG/DC/DT or A/C/G/U/T).
    """
    # Map atom index -> residue name (for neighbor lookup)
    idx_to_res = dict(zip(atoms_df["Index"].to_numpy(dtype=int), atoms_df["Residue"].astype(str).str.upper()))
    idx_to_chain = dict(zip(atoms_df["Index"].to_numpy(dtype=int), atoms_df["Chain"].astype(str)))
    idx_to_rseq = dict(zip(atoms_df["Index"].to_numpy(dtype=int), atoms_df["Residue Sequence"].to_numpy(dtype=int)))

    # Build atom membership for each residue instance (Chain, Residue Sequence, Residue)
    key_cols = ["Chain", "Residue Sequence", "Residue"]
    atoms_df["_res_key"] = list(zip(atoms_df["Chain"], atoms_df["Residue Sequence"], atoms_df["Residue"]))

    key_to_atom_indices: Dict[Tuple[str, int, str], List[int]] = {}
    for k, g in atoms_df.groupby("_res_key"):
        key_to_atom_indices[k] = g["Index"].to_numpy(dtype=int).tolist()

    # Map atom index -> parsed neighbors list
    idx_to_neighbors = dict(zip(atoms_df["Index"].to_numpy(dtype=int), atoms_df["_neighbors"]))

    # Label residues
    dna_facing = []
    dna_contact_counts = []

    for _, row in residue_instances.iterrows():
        k = (row["Chain"], int(row["Residue Sequence"]), row["Residue"])
        atom_indices = key_to_atom_indices.get(k, [])

        is_dna = False
        n_dna_contacts = 0

        for ai in atom_indices:
            for nb in idx_to_neighbors.get(ai, []):
                nb_res = idx_to_res.get(int(nb), "")
                if nb_res in NUC_RES_NAMES:
                    is_dna = True
                    n_dna_contacts += 1

        dna_facing.append(is_dna)
        dna_contact_counts.append(n_dna_contacts)

    out = residue_instances.copy()
    out["dna_facing"] = dna_facing
    out["dna_neighbor_hits"] = dna_contact_counts

    return out



def plot_violin_two_groups(
    df: pd.DataFrame,
    cfg: Config,
    title: str,
    out_prefix: str,
) -> Dict[str, str]:
    df = df.copy()
    df = df[df["Residue"].astype(str).str.upper() == cfg.residue_name.upper()].copy()
    df = df.dropna(subset=["curv"]).copy()

    grp0 = df[df["dna_facing"] == False]["curv"].to_numpy(dtype=float)  # noqa: E712
    grp1 = df[df["dna_facing"] == True]["curv"].to_numpy(dtype=float)   # noqa: E712

    if grp0.size < 2 or grp1.size < 2:
        raise ValueError(
            f"Not enough data for both groups. Non-DNA-facing n={grp0.size}, DNA-facing n={grp1.size}."
        )

    fig = plt.figure(figsize=(6.5, 5.25))
    ax = fig.add_subplot(111)

    data = [grp0, grp1]
    ax.violinplot(data, showmedians=True, showextrema=True)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [f"Not DNA-facing\n(n={grp0.size})", f"DNA-facing\n(n={grp1.size})"]
    )

    ax.set_title(title)
    ax.set_ylabel(df.attrs.get("curv_label", "|Mean residue curvature| (area-weighted)"))
    ax.grid(True, axis="y", alpha=0.3)

    # Add simple summary text
    med0 = float(np.median(grp0))
    med1 = float(np.median(grp1))
    mean0 = float(np.mean(grp0))
    mean1 = float(np.mean(grp1))

    ax.text(
        0.02,
        0.98,
        f"Median: {med0:.4f} vs {med1:.4f}\nMean: {mean0:.4f} vs {mean1:.4f}",
        transform=ax.transAxes,
        va="top",
    )

    os.makedirs(cfg.output_dir, exist_ok=True)
    png = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_dna_facing_violin.png")
    pdf = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_dna_facing_violin.pdf")
    csv = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_dna_facing_table.csv")

    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    df.to_csv(csv, index=False)

    return {"png": png, "pdf": pdf, "csv": csv}



def main() -> None:
    cfg = Config(
        curvature_kind="mean",
        use_magnitude=True,
        weight_by_area=True,
        residue_name="ALA",
        output_dir="alanine_dna_facing",
    )

    logs_file = choose_logs_file()

    logs_obj = read_logs2(
        logs_file,
        return_dict=False,
        all_=True,
        balls=True,
        surfs=False,
        edges=False,
        verts=False,
    )

    atoms_df = logs_obj["atoms"]

    residue_instances, atoms_with_neighbors = build_residue_instances_from_atoms(atoms_df, cfg)

    labeled = label_residue_dna_facing(residue_instances, atoms_with_neighbors, cfg)

    mol_name = str(
        logs_obj.get("data", {}).get(
            "name",
            os.path.splitext(os.path.basename(logs_file))[0],
        )
    )
    out_prefix = sanitize_filename(mol_name)

    result = plot_violin_two_groups(
        df=labeled,
        cfg=cfg,
        title=f"{mol_name}: {cfg.residue_name} curvature (DNA-facing vs not)",
        out_prefix=out_prefix,
    )

    # Print quick counts
    ala = labeled[labeled["Residue"].astype(str).str.upper() == cfg.residue_name.upper()]
    n0 = int((ala["dna_facing"] == False).sum())  # noqa: E712
    n1 = int((ala["dna_facing"] == True).sum())   # noqa: E712

    print(f"[{mol_name}] {cfg.residue_name} counts: Not DNA-facing={n0}, DNA-facing={n1}")
    print("Saved:")
    print(f"  PNG: {result['png']}")
    print(f"  PDF: {result['pdf']}")
    print(f"  CSV: {result['csv']}")



if __name__ == "__main__":
    main()
