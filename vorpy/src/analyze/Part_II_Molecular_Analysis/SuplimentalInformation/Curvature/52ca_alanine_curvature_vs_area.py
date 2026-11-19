import os
import sys
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# ----------------------------
# Path + imports
# ----------------------------
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



# ----------------------------
# Config
# ----------------------------
@dataclass
class Config:
    curvature_kind: str = "mean"        # "mean" or "gauss"
    use_magnitude: bool = True          # |curvature|
    weight_by_area: bool = True         # area-weighted residue curvature
    residue_name: str = "ALA"           # residue to analyze
    output_dir: str = "alanine_curvature_analysis"



# ----------------------------
# File chooser
# ----------------------------
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



# ----------------------------
# Build residue-level table
# ----------------------------
def build_residue_instances(atoms_df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
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
        raise KeyError(f"Missing columns: {missing}")

    df = atoms_df.copy()

    df["Residue Sequence"] = pd.to_numeric(df["Residue Sequence"], errors="coerce")
    df["Surface Area"] = pd.to_numeric(df["Surface Area"], errors="coerce")
    df["Average Mean Surface Curvature"] = pd.to_numeric(
        df["Average Mean Surface Curvature"], errors="coerce"
    )
    df["Average Gaussian Surface Curvature"] = pd.to_numeric(
        df["Average Gaussian Surface Curvature"], errors="coerce"
    )

    df = df.dropna(subset=["Residue Sequence"]).copy()
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

    # Area-weighted aggregation per residue instance
    df["_w"] = pd.to_numeric(df["Surface Area"], errors="coerce").fillna(0.0)
    df["_w"] = np.clip(df["_w"].to_numpy(dtype=float), 0.0, None)
    df["_wy"] = df["_w"] * df["_curv"]

    gcols = ["Chain", "Residue Sequence", "Residue"]

    agg = (
        df.groupby(gcols, as_index=False)
        .agg(
            wy_sum=("_wy", "sum"),
            w_sum=("_w", "sum"),
            area=("Surface Area", "sum"),
            n_atoms=("_curv", "size"),
            curv_fallback=("_curv", "mean"),
        )
    )

    agg["curv"] = np.where(
        agg["w_sum"] > 0.0,
        agg["wy_sum"] / agg["w_sum"],
        agg["curv_fallback"],
    )

    out = agg[gcols + ["curv", "area", "n_atoms"]].copy()
    out.attrs["curv_label"] = f"{curv_label} residue curvature (area-weighted)"

    return out



# ----------------------------
# Alanine-only scatter
# ----------------------------
def plot_residue_scatter(res_df: pd.DataFrame, cfg: Config, title: str) -> None:
    df = res_df[res_df["Residue"] == cfg.residue_name].copy()

    if df.empty:
        raise ValueError(f"No residues found for {cfg.residue_name}")

    # Clean
    df = df[(df["area"] > 0.0) & (~df["curv"].isna())].copy()

    x = df["area"].to_numpy(dtype=float)
    y = df["curv"].to_numpy(dtype=float)

    # Correlation (rank-based is most meaningful here)
    spear = np.corrcoef(
        np.argsort(np.argsort(x)),
        np.argsort(np.argsort(y)),
    )[0, 1]

    fig = plt.figure(figsize=(6.5, 5.5))
    ax = fig.add_subplot(111)

    ax.scatter(x, y, s=40, alpha=0.8)

    ax.set_xlabel("Alanine residue surface area")
    ax.set_ylabel(res_df.attrs.get("curv_label", "|Mean residue curvature|"))
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    ax.text(
        0.04,
        0.96,
        f"n = {len(df)}\nSpearman ρ ≈ {spear:.3f}",
        transform=ax.transAxes,
        va="top",
    )

    os.makedirs(cfg.output_dir, exist_ok=True)
    png = os.path.join(cfg.output_dir, "alanine_curvature_vs_area.png")
    pdf = os.path.join(cfg.output_dir, "alanine_curvature_vs_area.pdf")
    csv = os.path.join(cfg.output_dir, "alanine_curvature_vs_area.csv")

    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    df.to_csv(csv, index=False)

    print("Saved:")
    print(f"  PNG: {png}")
    print(f"  PDF: {pdf}")
    print(f"  CSV: {csv}")



# ----------------------------
# Main
# ----------------------------
def main() -> None:
    cfg = Config()

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

    res_df = build_residue_instances(atoms_df, cfg)

    mol_name = str(
        logs_obj.get("data", {}).get(
            "name",
            os.path.splitext(os.path.basename(logs_file))[0],
        )
    )

    plot_residue_scatter(
        res_df,
        cfg,
        title=f"{mol_name}: alanine curvature vs surface area",
    )



if __name__ == "__main__":
    main()
