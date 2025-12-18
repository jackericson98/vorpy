import os
import sys
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



AA3_TO_AA1 = {
    # --- Protein (standard) ---
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",

    # --- Protein (alternates / ambiguous) ---
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z", "XAA": "X",

    # --- DNA ---
    "DA": "A", "DG": "G", "DC": "C", "DT": "T",
    "T": "T",

    # --- RNA ---
    "A": "A", "G": "G", "C": "C", "U": "U",
}



def residue_to_one_letter(res: str) -> str:
    r = str(res).strip().upper()

    if r in AA3_TO_AA1:
        return AA3_TO_AA1[r]

    if len(r) >= 3:
        r3 = r[:3]
        return AA3_TO_AA1.get(r3, "X")

    return "X"



@dataclass
class ViolinConfig:
    curvature_kind: str = "mean"              # "mean" or "gauss"
    use_magnitude: bool = True               # plot |curv|
    weight_by_area: bool = True              # residue curvature = area-weighted across atoms
    group_label: str = "aa1"                 # "aa1" (one-letter) or "res3" (raw Residue)
    min_count_per_group: int = 10            # filter rare residue types
    sort_by: str = "median"                  # "median" or "mean" or "count"
    output_dir: str = "curvature_panel_violin"



def choose_logs_file(initial_dir: Optional[str] = None) -> str:
    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Select a logs CSV file",
        initialdir=initial_dir or os.getcwd(),
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not file_path:
        raise FileNotFoundError("No logs file selected.")

    return file_path



def sanitize_filename(s: str) -> str:
    keep = []
    for ch in str(s):
        if ch.isalnum() or ch in ("-", "_"):
            keep.append(ch)
        elif ch.isspace():
            keep.append("_")
    out = "".join(keep).strip("_")
    return out if out else "output"



def build_residue_instances_from_atoms(atoms_df: pd.DataFrame, cfg: ViolinConfig) -> pd.DataFrame:
    """
    Returns one row per residue instance (Chain + Residue Sequence + Residue),
    with an area-weighted curvature across atoms in that residue.

    Output columns:
      Chain
      Residue Sequence
      Residue (raw)
      aa1
      curv
      area
      n_atoms
    """
    required = [
        "Residue",
        "Residue Sequence",
        "Chain",
        "Surface Area",
        "Average Mean Surface Curvature",
        "Average Gaussian Surface Curvature",
    ]
    missing = [c for c in required if c not in atoms_df.columns]
    if missing:
        raise KeyError(
            f"Missing required columns in atoms dataframe: {missing}\n"
            f"Columns present: {list(atoms_df.columns)}"
        )

    df = atoms_df.copy()

    df["Residue Sequence"] = pd.to_numeric(df["Residue Sequence"], errors="coerce")
    df["Surface Area"] = pd.to_numeric(df["Surface Area"], errors="coerce")
    df["Average Mean Surface Curvature"] = pd.to_numeric(df["Average Mean Surface Curvature"], errors="coerce")
    df["Average Gaussian Surface Curvature"] = pd.to_numeric(df["Average Gaussian Surface Curvature"], errors="coerce")

    df = df.dropna(subset=["Residue Sequence"]).copy()
    df["Residue Sequence"] = df["Residue Sequence"].astype(int)

    if cfg.curvature_kind.lower().startswith("g"):
        curv_col = "Average Gaussian Surface Curvature"
        curv_name = "Gaussian"
    else:
        curv_col = "Average Mean Surface Curvature"
        curv_name = "Mean"

    df["_curv"] = df[curv_col].astype(float)

    if cfg.use_magnitude:
        df["_curv"] = np.abs(df["_curv"])
        curv_name = f"|{curv_name}|"

    if cfg.weight_by_area:
        df["_w"] = pd.to_numeric(df["Surface Area"], errors="coerce").astype(float)
        df["_w"] = df["_w"].replace([np.inf, -np.inf], np.nan).fillna(0.0)
        df["_w"] = np.clip(df["_w"].to_numpy(dtype=float), 0.0, None)

        df["_wy"] = df["_w"] * df["_curv"]

        gcols = ["Chain", "Residue Sequence", "Residue"]

        agg = (
            df.groupby(gcols, as_index=False)
            .agg(
                wy_sum=("_wy", "sum"),
                w_sum=("_w", "sum"),
                area=("_w", "sum"),
                n_atoms=("_curv", "size"),
                curv_fallback=("_curv", "mean"),
            )
        )

        agg["curv"] = np.where(
            agg["w_sum"].to_numpy(dtype=float) > 0.0,
            agg["wy_sum"].to_numpy(dtype=float) / agg["w_sum"].to_numpy(dtype=float),
            agg["curv_fallback"].to_numpy(dtype=float),
        )

        out = agg[gcols + ["curv", "area", "n_atoms"]].copy()

    else:
        out = (
            df.groupby(["Chain", "Residue Sequence", "Residue"], as_index=False)
            .agg(
                curv=("_curv", "mean"),
                area=("Surface Area", "sum"),
                n_atoms=("_curv", "size"),
            )
        )

    out["aa1"] = out["Residue"].apply(residue_to_one_letter)

    out.attrs["curv_label"] = (
        f"{curv_name} residue curvature (area-weighted)"
        if cfg.weight_by_area
        else f"{curv_name} residue curvature"
    )
    out.attrs["curv_col_used"] = curv_col

    return out



def plot_violin_by_residue_type(
    res_inst: pd.DataFrame,
    cfg: ViolinConfig,
    title: str,
    out_prefix: str,
) -> Dict[str, str]:
    if cfg.group_label == "res3":
        group_col = "Residue"
        group_title = "Residue (raw)"
    else:
        group_col = "aa1"
        group_title = "Residue (1-letter)"

    # Filter unknowns and NaNs
    df = res_inst.copy()
    df = df.dropna(subset=["curv"]).copy()
    df[group_col] = df[group_col].astype(str)

    # Optionally drop X if you don’t want unknowns; comment out if you want to keep them
    # df = df[df[group_col] != "X"].copy()

    # Group values
    grouped: Dict[str, np.ndarray] = {}
    for k, g in df.groupby(group_col):
        vals = g["curv"].to_numpy(dtype=float)
        vals = vals[~np.isnan(vals)]
        if vals.size >= cfg.min_count_per_group:
            grouped[k] = vals

    if not grouped:
        raise ValueError(
            f"No residue types have at least {cfg.min_count_per_group} samples. "
            "Lower min_count_per_group or check input data."
        )

    labels = list(grouped.keys())

    # Sort groups for readability
    if cfg.sort_by == "count":
        labels = sorted(labels, key=lambda k: grouped[k].size, reverse=True)
    elif cfg.sort_by == "mean":
        labels = sorted(labels, key=lambda k: float(np.mean(grouped[k])), reverse=True)
    else:
        labels = sorted(labels, key=lambda k: float(np.median(grouped[k])), reverse=True)

    data = [grouped[k] for k in labels]

    fig_w = max(8.0, 0.45 * len(labels))
    fig = plt.figure(figsize=(fig_w, 5.0))
    ax = fig.add_subplot(111)

    parts = ax.violinplot(
        data,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )

    ax.set_title(title)
    ax.set_xlabel(group_title)
    ax.set_ylabel(res_inst.attrs.get("curv_label", "|Mean residue curvature| (area-weighted)"))

    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=0)

    # Add counts under each label (optional but helpful)
    counts = [grouped[k].size for k in labels]
    for i, n in enumerate(counts, start=1):
        ax.text(i, ax.get_ylim()[0], f"n={n}", ha="center", va="top", fontsize=8)

    ax.grid(True, axis="y", alpha=0.3)

    os.makedirs(cfg.output_dir, exist_ok=True)
    png_path = os.path.join(cfg.output_dir, f"{out_prefix}_violin_by_{group_col}.png")
    pdf_path = os.path.join(cfg.output_dir, f"{out_prefix}_violin_by_{group_col}.pdf")
    csv_path = os.path.join(cfg.output_dir, f"{out_prefix}_residue_instances.csv")

    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    # Save underlying data for reproducibility
    df.to_csv(csv_path, index=False)

    return {"png": png_path, "pdf": pdf_path, "csv": csv_path}



def main() -> None:
    cfg = ViolinConfig(
        curvature_kind="mean",
        use_magnitude=True,
        weight_by_area=True,
        group_label="aa1",           # "aa1" or "res3"
        min_count_per_group=10,
        sort_by="median",            # "median" / "mean" / "count"
        output_dir="curvature_panel_violin",
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

    res_inst = build_residue_instances_from_atoms(atoms_df, cfg)

    molecule_name = str(
        logs_obj.get("data", {}).get("name", os.path.splitext(os.path.basename(logs_file))[0])
    )
    out_prefix = sanitize_filename(molecule_name)

    result = plot_violin_by_residue_type(
        res_inst=res_inst,
        cfg=cfg,
        title=f"{molecule_name}: residue curvature distribution by residue type",
        out_prefix=out_prefix,
    )

    print(f"[{molecule_name}] Using curvature column: {res_inst.attrs.get('curv_col_used', '')}")
    print("Saved:")
    print(f"  PNG: {result['png']}")
    print(f"  PDF: {result['pdf']}")
    print(f"  CSV: {result['csv']}")



if __name__ == "__main__":
    main()
