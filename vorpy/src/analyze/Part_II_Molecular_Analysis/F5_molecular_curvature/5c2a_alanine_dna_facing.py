import os
import sys
import tkinter as tk
from tkinter import filedialog
from typing import Set, Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



# ----------------------------
# Configuration
# ----------------------------
ALA_NAME = "ALA"

DNA_RESIDUES: Set[str] = {
    "DA", "DG", "DC", "DT",   # DNA
    "A", "G", "C", "T", "U",  # RNA / ambiguous
}

OUTPUT_DIR = "5c2A_alanine_dna_facing"


# ----------------------------
# Helpers
# ----------------------------
def choose_logs_file() -> str:
    root = tk.Tk()
    root.withdraw()

    path = filedialog.askopenfilename(
        title="Select logs CSV file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )

    root.destroy()

    if not path:
        raise FileNotFoundError("No logs file selected.")

    return path


def build_atom_index_maps(atoms_df: pd.DataFrame):
    """
    Build quick lookup maps needed for neighbor classification.
    """
    idx_to_res = dict(
        zip(
            atoms_df["Index"].astype(int),
            atoms_df["Residue"].astype(str).str.upper(),
        )
    )

    idx_to_neighbors = dict(
        zip(
            atoms_df["Index"].astype(int),
            atoms_df["Neighbors"],
        )
    )

    return idx_to_res, idx_to_neighbors


def compute_residue_curvature(atoms_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute area-weighted |mean curvature| per residue instance.
    """
    df = atoms_df.copy()

    for col in [
        "Index",
        "Residue Sequence",
        "Surface Area",
        "Average Mean Surface Curvature",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Index", "Residue Sequence"])

    df["Index"] = df["Index"].astype(int)
    df["Residue Sequence"] = df["Residue Sequence"].astype(int)

    df["_curv"] = np.abs(df["Average Mean Surface Curvature"].astype(float))
    df["_area"] = np.clip(df["Surface Area"].astype(float), 0.0, None)
    df["_wy"] = df["_curv"] * df["_area"]

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

    return agg[gcols + ["curv", "n_atoms"]]


def label_dna_facing(
    res_df: pd.DataFrame,
    atoms_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Label each residue instance as DNA-facing if any neighbor atom
    belongs to a DNA residue.
    """
    idx_to_res, idx_to_neighbors = build_atom_index_maps(atoms_df)

    # Map residue instance -> atom indices
    atoms_df["_res_key"] = list(
        zip(
            atoms_df["Chain"],
            atoms_df["Residue Sequence"].astype(int),
            atoms_df["Residue"],
        )
    )

    res_to_atoms: Dict[tuple, List[int]] = {}
    for k, g in atoms_df.groupby("_res_key"):
        res_to_atoms[k] = g["Index"].astype(int).tolist()

    dna_flags = []

    for _, row in res_df.iterrows():
        key = (row["Chain"], row["Residue Sequence"], row["Residue"])
        atom_indices = res_to_atoms.get(key, [])

        is_dna = False

        for ai in atom_indices:
            neighs = idx_to_neighbors.get(ai, [])
            for nb in neighs:
                nb_res = idx_to_res.get(int(nb), "")
                if nb_res in DNA_RESIDUES:
                    is_dna = True
                    break
            if is_dna:
                break

        dna_flags.append(is_dna)

    out = res_df.copy()
    out["dna_facing"] = dna_flags
    return out


def plot_violin(df: pd.DataFrame, title: str):
    dna = df[df["dna_facing"]]["curv"].to_numpy()
    non = df[~df["dna_facing"]]["curv"].to_numpy()

    fig, ax = plt.subplots(figsize=(6.5, 5.0))

    ax.violinplot(
        [non, dna],
        showmedians=True,
        showextrema=True,
    )

    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [
            f"Not DNA-facing\n(n={len(non)})",
            f"DNA-facing\n(n={len(dna)})",
        ]
    )

    ax.set_ylabel("|Mean residue curvature| (area-weighted)")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)

    # Summary text
    ax.text(
        0.02,
        0.98,
        f"Median: {np.median(non):.4f} vs {np.median(dna):.4f}\n"
        f"Mean: {np.mean(non):.4f} vs {np.mean(dna):.4f}",
        transform=ax.transAxes,
        va="top",
    )

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    png = os.path.join(OUTPUT_DIR, "alanine_dna_facing_violin.png")
    pdf = os.path.join(OUTPUT_DIR, "alanine_dna_facing_violin.pdf")

    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    return png, pdf


# ----------------------------
# Main
# ----------------------------
def main():
    logs_file = choose_logs_file()

    logs = read_logs2(
        logs_file,
        return_dict=False,
        all_=True,
        balls=True,
        surfs=False,
        edges=False,
        verts=False,
    )

    atoms_df = logs["atoms"].copy()

    # Build residue curvature table (all residues), then subset to alanine
    res_df_all = compute_residue_curvature(atoms_df)
    res_df_ala = res_df_all[res_df_all["Residue"].astype(str).str.upper() == ALA_NAME].copy()

    # Label alanine residues using FULL atoms_df (so neighbors can resolve DNA residues)
    res_df_ala = label_dna_facing(res_df_ala, atoms_df)

    mol_name = logs.get("data", {}).get(
        "name", os.path.splitext(os.path.basename(logs_file))[0]
    )

    png, pdf = plot_violin(
        res_df_ala,
        title=f"{mol_name}: Alanine curvature (DNA-facing vs not)",
    )

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    res_df_ala.to_csv(
        os.path.join(OUTPUT_DIR, "alanine_dna_facing_table.csv"),
        index=False,
    )

    print(f"[{mol_name}] Alanine DNA-facing analysis complete.")
    print(f"  DNA-facing     : {(res_df_ala['dna_facing']).sum()}")
    print(f"  Not DNA-facing : {(~res_df_ala['dna_facing']).sum()}")
    print(f"  Saved: {png}")
    print(f"  Saved: {pdf}")

    # Restrict to alanine early
    atoms_df = atoms_df[atoms_df["Residue"].astype(str).str.upper() == ALA_NAME].copy()

    res_df = compute_residue_curvature(atoms_df)
    res_df = label_dna_facing(res_df, atoms_df)

    mol_name = logs.get("data", {}).get(
        "name", os.path.splitext(os.path.basename(logs_file))[0]
    )

    png, pdf = plot_violin(
        res_df,
        title=f"{mol_name}: Alanine curvature (DNA-facing vs not)",
    )

    res_df.to_csv(
        os.path.join(OUTPUT_DIR, "alanine_dna_facing_table.csv"),
        index=False,
    )

    print(f"[{mol_name}] Alanine DNA-facing analysis complete.")
    print(f"  DNA-facing     : {(res_df['dna_facing']).sum()}")
    print(f"  Not DNA-facing : {(~res_df['dna_facing']).sum()}")
    print(f"  Saved: {png}")
    print(f"  Saved: {pdf}")



if __name__ == "__main__":
    main()
