import os
import sys
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

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
class CurvConfig:
    curvature_kind: str = "mean"                 # "mean" or "gauss"
    use_magnitude: bool = True                  # True -> |curvature|
    weight_by_area: bool = True                 # residue aggregation uses atom Surface Area as weight
    output_dir: str = "curvature_panel_c2"

    do_area_coupling: bool = True               # 5c2A
    do_lorenz: bool = False                     # optional legacy panel



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



def _rankdata(a: np.ndarray) -> np.ndarray:
    """
    Minimal rankdata implementation (average ranks for ties).
    Avoids scipy dependency.
    """
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)

    sorted_a = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sorted_a[j + 1] == sorted_a[i]:
            j += 1
        if j > i:
            avg = 0.5 * (i + 1 + j + 1)
            ranks[order[i:j + 1]] = avg
        i = j + 1

    return ranks



def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """
    Spearman correlation computed as Pearson correlation on ranks.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    rx = _rankdata(x)
    ry = _rankdata(y)

    rx = rx - float(np.mean(rx))
    ry = ry - float(np.mean(ry))

    denom = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
    if denom <= 0.0:
        return float("nan")

    return float(np.sum(rx * ry) / denom)



def build_residue_df_from_atoms(atoms_df: pd.DataFrame, cfg: CurvConfig) -> pd.DataFrame:
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

        grouped = agg[gcols + ["curv", "area", "n_atoms"]].copy()

    else:
        grouped = (
            df.groupby(["Chain", "Residue Sequence", "Residue"], as_index=False)
            .agg(curv=("_curv", "mean"), area=("Surface Area", "sum"), n_atoms=("_curv", "size"))
        )

    grouped["aa1"] = grouped["Residue"].apply(residue_to_one_letter)
    grouped["res_label"] = grouped.apply(
        lambda r: f"{r['aa1']}{int(r['Residue Sequence'])}",
        axis=1,
    )

    grouped.attrs["curv_label"] = (
        f"{curv_name} residue curvature (area-weighted)"
        if cfg.weight_by_area
        else f"{curv_name} residue curvature"
    )
    grouped.attrs["curv_col_used"] = curv_col

    return grouped



def aggregate_across_chains(res_df_chain: pd.DataFrame) -> pd.DataFrame:
    gcols = ["Residue Sequence", "Residue", "aa1"]

    out = (
        res_df_chain
        .groupby(gcols, as_index=False)
        .agg(
            curv_mean=("curv", "mean"),
            curv_sd=("curv", "std"),
            n_repeats=("curv", "count"),
            area_sum=("area", "sum"),
        )
    )

    out["curv_sd"] = out["curv_sd"].fillna(0.0)
    out["res_label"] = out.apply(lambda r: f"{r['aa1']}{int(r['Residue Sequence'])}", axis=1)

    out.attrs["curv_label"] = res_df_chain.attrs.get("curv_label", "Residue curvature")
    out.attrs["curv_col_used"] = res_df_chain.attrs.get("curv_col_used", "")

    return out



def plot_panel_c2_area_coupling(
    res_df: pd.DataFrame,
    cfg: CurvConfig,
    title: str,
    out_prefix: str,
    use_log_area: bool = True,
    fit_line: bool = True,
) -> Dict[str, object]:
    """
    5c2A: Scatter of curvature magnitude vs residue surface area.

    Accepts either:
      - per-chain residue df: columns include ['area','curv','res_label','Chain',...]
      - across-chain aggregated df: columns include ['area_sum','curv_mean','curv_sd','n_repeats','res_label',...]
    """
    if "curv_mean" in res_df.columns:
        y = res_df["curv_mean"].to_numpy(dtype=float)
        x = res_df["area_sum"].to_numpy(dtype=float) if "area_sum" in res_df.columns else np.full_like(y, np.nan)
        labels = res_df["res_label"].to_numpy(dtype=str) if "res_label" in res_df.columns else np.array([""] * len(y))
        mode = "across-chain mean"
    else:
        y = res_df["curv"].to_numpy(dtype=float)
        x = res_df["area"].to_numpy(dtype=float)
        labels = res_df["res_label"].to_numpy(dtype=str) if "res_label" in res_df.columns else np.array([""] * len(y))
        mode = "per-chain"

    keep = (~np.isnan(x)) & (~np.isnan(y)) & (x > 0.0) & (y >= 0.0)
    x = x[keep]
    y = y[keep]
    labels = labels[keep]

    if x.size < 3:
        raise ValueError("Not enough valid residue points to plot area-coupling (need at least 3).")

    x_plot = np.log10(x) if use_log_area else x
    x_label = "log10(residue surface area)" if use_log_area else "Residue surface area"

    pearson = float(np.corrcoef(x_plot, y)[0, 1]) if x_plot.size > 1 else float("nan")
    spear = spearman_rho(x_plot, y)

    m = b = float("nan")
    if fit_line and x_plot.size >= 2:
        m, b = np.polyfit(x_plot, y, deg=1)

    fig = plt.figure(figsize=(6.75, 5.25))
    ax = fig.add_subplot(111)

    ax.scatter(x_plot, y, s=28, alpha=0.85)

    if fit_line and np.isfinite(m) and np.isfinite(b):
        xs = np.linspace(float(np.min(x_plot)), float(np.max(x_plot)), 200)
        ax.plot(xs, m * xs + b, linewidth=2.0, label="Linear trend")

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("|Mean residue curvature| (area-weighted)")
    ax.grid(True, alpha=0.3)

    stat_txt = f"Mode: {mode}\nSpearman ρ = {spear:.3f}\nPearson r = {pearson:.3f}"
    ax.text(0.04, 0.96, stat_txt, transform=ax.transAxes, va="top")

    if fit_line and np.isfinite(m) and np.isfinite(b):
        ax.legend()

    os.makedirs(cfg.output_dir, exist_ok=True)
    png_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c2A_area_coupling.png")
    pdf_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c2A_area_coupling.pdf")
    svg_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c2A_area_coupling.svg")
    csv_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c2A_area_coupling.csv")

    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)

    out_table = pd.DataFrame(
        {
            "res_label": labels,
            "area": x,
            "area_plotted": x_plot,
            "curv": y,
        }
    )
    out_table.to_csv(csv_path, index=False)

    return {
        "png": png_path,
        "pdf": pdf_path,
        "svg": svg_path,
        "csv": csv_path,
        "spearman_rho": spear,
        "pearson_r": pearson,
        "fit_m": m,
        "fit_b": b,
        "n_points": int(x.size),
        "mode": mode,
    }



def lorenz_curve(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    v = np.clip(v, 0.0, None)

    if v.size == 0 or float(np.sum(v)) <= 0.0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])

    v_sorted = np.sort(v)
    cum_v = np.cumsum(v_sorted)
    total = float(cum_v[-1])

    y = np.concatenate(([0.0], cum_v / total))
    x = np.linspace(0.0, 1.0, num=y.size)

    return x, y



def gini_from_lorenz(x: np.ndarray, y: np.ndarray) -> float:
    auc = float(np.trapezoid(y, x))
    return float(1.0 - 2.0 * auc)



def top_k_contrib(values: np.ndarray, percents: Tuple[float, ...]) -> Dict[float, float]:
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    v = np.clip(v, 0.0, None)

    total = float(np.sum(v))
    if total <= 0.0 or v.size == 0:
        return {p: float("nan") for p in percents}

    v_desc = np.sort(v)[::-1]
    n = v_desc.size

    out = {}
    for p in percents:
        k = max(1, int(np.ceil((p / 100.0) * n)))
        out[p] = float(np.sum(v_desc[:k]) / total)

    return out



def plot_panel_c2_lorenz(
    res_df: pd.DataFrame,
    cfg: CurvConfig,
    title: str,
    out_prefix: str,
) -> Dict[str, object]:
    y = res_df["curv_mean"].to_numpy(dtype=float) if "curv_mean" in res_df.columns else res_df["curv"].to_numpy(dtype=float)

    x, y_lor = lorenz_curve(y)
    gini = gini_from_lorenz(x, y_lor)
    top = top_k_contrib(y, (5.0, 10.0, 20.0))

    fig = plt.figure(figsize=(6.25, 5.25))
    ax = fig.add_subplot(111)

    ax.plot(x, y_lor, linewidth=2.0, label="Lorenz curve")
    ax.plot([0, 1], [0, 1], linewidth=1.0, linestyle="--", label="Uniform (reference)")

    ax.set_title(title)
    ax.set_xlabel("Cumulative fraction of residues")
    ax.set_ylabel("Cumulative fraction of curvature mass")
    ax.grid(True, alpha=0.3)
    ax.legend()

    annotation = (
        f"Gini = {gini:.3f}\n"
        f"Top 5% → {100.0 * top[5.0]:.1f}%\n"
        f"Top 10% → {100.0 * top[10.0]:.1f}%\n"
        f"Top 20% → {100.0 * top[20.0]:.1f}%"
    )
    ax.text(0.05, 0.80, annotation, transform=ax.transAxes, fontsize=10, va="top")

    os.makedirs(cfg.output_dir, exist_ok=True)
    png_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c2_lorenz.png")
    pdf_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c2_lorenz.pdf")

    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    return {
        "png": png_path,
        "pdf": pdf_path,
        "gini": gini,
        "top_5": top[5.0],
        "top_10": top[10.0],
        "top_20": top[20.0],
    }



def main() -> None:
    cfg = CurvConfig(
        curvature_kind="mean",
        use_magnitude=True,
        weight_by_area=True,
        output_dir="curvature_panel_c2",
        do_area_coupling=True,
        do_lorenz=False,
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

    res_df_chain = build_residue_df_from_atoms(atoms_df, cfg)
    res_df = aggregate_across_chains(res_df_chain)

    molecule_name = str(
        logs_obj.get("data", {}).get("name", os.path.splitext(os.path.basename(logs_file))[0])
    )
    out_prefix = sanitize_filename(molecule_name)

    if cfg.do_area_coupling:
        area_result = plot_panel_c2_area_coupling(
            res_df=res_df,
            cfg=cfg,
            title=f"{molecule_name}: curvature vs residue surface area",
            out_prefix=out_prefix,
            use_log_area=True,
            fit_line=True,
        )

        print(
            f"[{molecule_name}] 5c2A area-coupling: n={area_result['n_points']} | "
            f"Spearman ρ={area_result['spearman_rho']:.3f} | Pearson r={area_result['pearson_r']:.3f}"
        )
        print(f"Saved 5c2A PNG: {area_result['png']}")
        print(f"Saved 5c2A CSV: {area_result['csv']}")

    if cfg.do_lorenz:
        lorenz_result = plot_panel_c2_lorenz(
            res_df=res_df,
            cfg=cfg,
            title=f"{molecule_name}: curvature concentration",
            out_prefix=out_prefix,
        )

        print(
            f"[{molecule_name}] Lorenz: Gini={lorenz_result['gini']:.3f} | "
            f"Top 5%={100.0 * lorenz_result['top_5']:.1f}% | "
            f"Top 10%={100.0 * lorenz_result['top_10']:.1f}% | "
            f"Top 20%={100.0 * lorenz_result['top_20']:.1f}%"
        )
        print(f"Saved Lorenz PNG: {lorenz_result['png']}")



if __name__ == "__main__":
    main()
