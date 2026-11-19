import os
import sys
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Optional, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



@dataclass
class Config:
    curvature_kind: str = "mean"        # "mean" or "gauss"
    use_magnitude: bool = True          # |curvature|
    weight_by_area: bool = True         # area-weighted residue curvature across atoms
    residue_name: str = "ALA"           # residue to analyze
    output_dir: str = "alanine_curvature_by_chain"



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



def build_residue_instances(atoms_df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """
    Returns one row per residue instance (Chain + Residue Sequence + Residue),
    with area-weighted curvature across atoms in that residue.

    Output columns:
      Chain
      Residue Sequence
      Residue
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
        raise KeyError(f"Missing columns in atoms dataframe: {missing}")

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

    df["_w"] = pd.to_numeric(df["Surface Area"], errors="coerce").fillna(0.0).astype(float)
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
    out.attrs["curv_col_used"] = curv_col

    return out



def per_chain_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("Chain", as_index=False)
        .agg(
            n=("curv", "size"),
            curv_mean=("curv", "mean"),
            curv_median=("curv", "median"),
            curv_std=("curv", "std"),
            area_mean=("area", "mean"),
            area_median=("area", "median"),
        )
        .sort_values("n", ascending=False)
    )



def plot_residue_scatter_by_chain(
    res_df: pd.DataFrame,
    cfg: Config,
    title: str,
    out_prefix: str,
    use_log_area: bool = False,
) -> Dict[str, str]:
    df = res_df[res_df["Residue"].astype(str).str.upper() == cfg.residue_name.upper()].copy()

    if df.empty:
        raise ValueError(f"No residues found for {cfg.residue_name}")

    df = df[(df["area"] > 0.0) & (~df["curv"].isna())].copy()

    if use_log_area:
        df["area_plot"] = np.log10(df["area"].to_numpy(dtype=float))
        x_label = "log10(residue surface area)"
    else:
        df["area_plot"] = df["area"].to_numpy(dtype=float)
        x_label = "Residue surface area"

    chains = sorted(df["Chain"].astype(str).unique().tolist())

    fig = plt.figure(figsize=(7.25, 5.5))
    ax = fig.add_subplot(111)

    for ch in chains:
        sub = df[df["Chain"].astype(str) == ch]
        ax.scatter(
            sub["area_plot"].to_numpy(dtype=float),
            sub["curv"].to_numpy(dtype=float),
            s=42,
            alpha=0.80,
            label=f"Chain {ch} (n={len(sub)})",
        )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(res_df.attrs.get("curv_label", "|Mean residue curvature| (area-weighted)"))
    ax.grid(True, alpha=0.3)

    ax.legend(loc="best", fontsize=9)

    os.makedirs(cfg.output_dir, exist_ok=True)
    png = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_curv_vs_area_by_chain.png")
    pdf = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_curv_vs_area_by_chain.pdf")
    csv = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_points.csv")
    summary_csv = os.path.join(cfg.output_dir, f"{out_prefix}_{cfg.residue_name}_chain_summary.csv")

    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    df.to_csv(csv, index=False)

    summ = per_chain_summary(df)
    summ.to_csv(summary_csv, index=False)

    return {"png": png, "pdf": pdf, "csv": csv, "summary_csv": summary_csv}



def main() -> None:
    cfg = Config(
        curvature_kind="mean",
        use_magnitude=True,
        weight_by_area=True,
        residue_name="ALA",
        output_dir="alanine_curvature_by_chain",
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

    res_df = build_residue_instances(atoms_df, cfg)

    mol_name = str(
        logs_obj.get("data", {}).get(
            "name",
            os.path.splitext(os.path.basename(logs_file))[0],
        )
    )
    out_prefix = sanitize_filename(mol_name)

    result = plot_residue_scatter_by_chain(
        res_df=res_df,
        cfg=cfg,
        title=f"{mol_name}: {cfg.residue_name} curvature vs surface area (by chain)",
        out_prefix=out_prefix,
        use_log_area=False,
    )

    print(f"[{mol_name}] Using curvature column: {res_df.attrs.get('curv_col_used', '')}")
    print("Saved:")
    print(f"  PNG: {result['png']}")
    print(f"  PDF: {result['pdf']}")
    print(f"  Points CSV: {result['csv']}")
    print(f"  Chain summary CSV: {result['summary_csv']}")



if __name__ == "__main__":
    main()
