import os
import sys
import tkinter as tk
from tkinter import filedialog
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



@dataclass
class CurvConfig:
    curvature_kind: str = "mean"                 # "mean" or "gauss"
    use_magnitude: bool = True                  # True -> |curvature|
    top_percent: float = 5.0                    # mark top X% residues
    smooth_window: int = 9                      # odd integer recommended
    clip_percentiles: Optional[Tuple[float, float]] = (2.0, 98.0)  # None to disable
    weight_by_area: bool = True                 # residue aggregation uses atom Surface Area as weight
    output_dir: str = "curvature_panel_c1"


AA3_TO_AA1 = {
    # --- Protein (standard) ---
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",

    # --- Protein (alternates / ambiguous) ---
    "SEC": "U",
    "PYL": "O",
    "ASX": "B",
    "GLX": "Z",
    "XAA": "X",

    "DA": "A",
    "A": "A",
    "DG": "G",
    "G": "G",
    "DC": "C",
    "C": "C",
    "DT": "T",
    "T":  "T",
    "DU": "U",
    "U":  "U",
}

# Global or module-level collector
UNKNOWN_RESIDUES = defaultdict(int)


def residue_to_one_letter(res: str) -> str:
    r = str(res).strip().upper()

    # Direct hit
    if r in AA3_TO_AA1:
        return AA3_TO_AA1[r]

    # Try first 3 characters (protein-style)
    if len(r) >= 3:
        r3 = r[:3]
        if r3 in AA3_TO_AA1:
            return AA3_TO_AA1[r3]

    # Unknown residue → log it
    UNKNOWN_RESIDUES[r] += 1
    return "X"



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


def moving_average(y: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return y.copy()

    if window % 2 == 0:
        window += 1

    kernel = np.ones(window, dtype=float) / float(window)
    pad = window // 2
    y_pad = np.pad(y, (pad, pad), mode="reflect")
    return np.convolve(y_pad, kernel, mode="valid")


def percentile_clip(y: np.ndarray, lo: float, hi: float) -> np.ndarray:
    y = np.asarray(y, dtype=float)

    lo_v = np.nanpercentile(y, lo)
    hi_v = np.nanpercentile(y, hi)

    return np.clip(y, lo_v, hi_v)


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

        # weighted curvature contribution
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

        # weighted mean; fallback to simple mean when weights are zero
        agg["curv"] = np.where(
            agg["w_sum"].to_numpy(dtype=float) > 0.0,
            agg["wy_sum"].to_numpy(dtype=float) / agg["w_sum"].to_numpy(dtype=float),
            agg["curv_fallback"].to_numpy(dtype=float),
        )

        grouped = agg[gcols + ["curv", "area", "n_atoms"]].copy()


    else:
        grouped = (
            df.groupby(["Chain", "Residue Sequence", "Residue"], as_index=False)
            .agg(curv=("__curv__", "mean"))
        )

        grouped["area"] = df.groupby(["Chain", "Residue Sequence", "Residue"])["Surface Area"].sum().values
        grouped["n_atoms"] = df.groupby(["Chain", "Residue Sequence", "Residue"]).size().values

    grouped["res_label"] = grouped.apply(
        lambda r: f"{r['Residue']}{int(r['Residue Sequence'])}" if str(r["Residue"]).strip() else str(int(r["Residue Sequence"])),
        axis=1,
    )

    grouped.attrs["curv_label"] = f"{curv_name} residue curvature (area-weighted)" if cfg.weight_by_area else f"{curv_name} residue curvature"
    grouped.attrs["curv_col_used"] = curv_col

    return grouped


def aggregate_across_chains(res_df_chain: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate curvature across repeated chains at each aligned residue index.
    Requires columns: Residue Sequence, Residue, curv (and optionally area, n_atoms).
    """
    required = ["Residue Sequence", "Residue", "curv"]
    missing = [c for c in required if c not in res_df_chain.columns]
    if missing:
        raise KeyError(f"aggregate_across_chains missing columns: {missing}")

    gcols = ["Residue Sequence", "Residue"]

    out = (
        res_df_chain
        .groupby(gcols, as_index=False)
        .agg(
            curv_mean=("curv", "mean"),
            curv_sd=("curv", "std"),
            n_repeats=("curv", "count"),
        )
    )

    out["curv_sd"] = out["curv_sd"].fillna(0.0)
    out["aa1"] = out["Residue"].apply(residue_to_one_letter)
    out["res_label"] = out.apply(
        lambda r: f"{r['aa1']}{int(r['Residue Sequence'])}",
        axis=1,
    )

    out.attrs["curv_label"] = res_df_chain.attrs.get("curv_label", "Residue curvature")
    out.attrs["curv_col_used"] = res_df_chain.attrs.get("curv_col_used", "")

    return out


def expand_to_full_residue_range(res_df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands residue dataframe to include missing residue indices as NaNs.
    This allows matplotlib to show gaps rather than compressing the x-axis.
    """
    min_res = int(res_df["Residue Sequence"].min())
    max_res = int(res_df["Residue Sequence"].max())

    full_index = pd.DataFrame(
        {"Residue Sequence": np.arange(min_res, max_res + 1, dtype=int)}
    )

    expanded = full_index.merge(
        res_df,
        on="Residue Sequence",
        how="left",
        sort=True,
    )

    return expanded


def aggregate_across_repeats(res_df_chain: pd.DataFrame) -> pd.DataFrame:
    """
    Take per-chain residue curvature and aggregate across repeated sequences
    (e.g., tetramer chains) at each residue position.

    Output:
      Residue Sequence
      Residue
      curv_mean (mean across chains)
      curv_sd
      n_repeats (number of chains contributing)
      res_label (e.g., GLU330)
    """
    gcols = ["Residue Sequence", "Residue"]

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

    out["res_label"] = out.apply(
        lambda r: f"{r['Residue']}{int(r['Residue Sequence'])}" if str(r["Residue"]).strip() else str(int(r["Residue Sequence"])),
        axis=1,
    )

    # Preserve label metadata for plotting
    out.attrs["curv_label"] = res_df_chain.attrs.get("curv_label", "Residue curvature")
    out.attrs["curv_col_used"] = res_df_chain.attrs.get("curv_col_used", "")

    return out


def plot_panel_c1(res_df: pd.DataFrame, cfg: CurvConfig, title: str, out_prefix: str) -> Dict[str, object]:

    # res_df now represents aggregates across repeats
    # res_df is now across-chain aggregate
    y = res_df["curv_mean"].to_numpy(dtype=float)
    y_sd = res_df["curv_sd"].to_numpy(dtype=float)
    aa1 = res_df["aa1"].to_numpy(dtype=str)
    labels = res_df["res_label"].to_numpy(dtype=str)

    x_idx = res_df["Residue Sequence"].to_numpy(dtype=int)

    order = np.argsort(x_idx)
    x_idx = x_idx[order]
    y = y[order]
    y_sd = y_sd[order]
    aa1 = aa1[order]
    labels = labels[order]

    keep = ~np.isnan(y)
    x_idx = x_idx[keep]
    y = y[keep]
    y_sd = y_sd[keep]
    aa1 = aa1[keep]
    labels = labels[keep]

    # Plot positions 1..N but label them with residue letters
    x = np.arange(1, len(y) + 1, dtype=int)


    order = np.argsort(x)
    x = x[order]
    y = y[order]
    labels = labels[order]

    keep = ~np.isnan(y)
    x = x[keep]
    y = y[keep]
    labels = labels[keep]

    y_plot = y.copy()
    if cfg.clip_percentiles is not None:
        y_plot = percentile_clip(y_plot, cfg.clip_percentiles[0], cfg.clip_percentiles[1])

    y_smooth = moving_average(y_plot, cfg.smooth_window)


    n = len(y)
    k = max(1, int(np.ceil((cfg.top_percent / 100.0) * n)))
    top_idx = np.argsort(y)[-k:]

    total_mass = float(np.sum(y))
    top_mass = float(np.sum(y[top_idx]))
    frac_mass = (top_mass / total_mass) if total_mass > 0 else np.nan

    fig = plt.figure(figsize=(11, 4.75))
    ax = fig.add_subplot(111)

    # If you prefer SEM instead of SD:
    # y_err = y_sd / np.sqrt(np.maximum(n_rep, 1))
    y_err = y_sd

    lower = y_plot - y_err
    upper = y_plot + y_err
    ax.fill_between(x, lower, upper, alpha=0.2, linewidth=0.0, label="Across-chain variability")

    ax.plot(
        x,
        y_plot,
        linestyle="None",
        marker="o",
        markersize=4,
        alpha=0.8,
        label="Residue curvature",
    )

    ax.plot(
        x,
        y_smooth,
        linewidth=2.0,
        alpha=0.9,
        label=f"Smoothed (window={cfg.smooth_window})",
    )

    ax.scatter(
        x[top_idx],
        y_plot[top_idx],
        s=45,
        zorder=5,
        label=f"Top {cfg.top_percent:.1f}% residues",
    )

    ax.set_title(title)
    ax.set_xlabel("Aligned residue position (aggregated across repeated chains)")
    ax.set_ylabel(res_df.attrs.get("curv_label", "Residue curvature") + " (mean across chains)")

    ax.grid(True, alpha=0.3)
    ax.legend()

    png_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c1.png")
    pdf_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c1.pdf")
    csv_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c1_residue_table.csv")
    svg_path = os.path.join(cfg.output_dir, f"{out_prefix}_panel_c1.svg")
    ax.set_xticks(x)
    ax.set_xticklabels(aa1, fontsize=9)
    ax.set_xlabel("Aligned residue position (sequence letters; averaged across chains)")

    # Optional: print the sequence string in the console for easy copy/paste
    seq_str = "".join(list(aa1))
    print(f"Sequence window (x-axis): {seq_str}")

    # Optional: show mapping back to absolute residue indices
    ax.text(
        0.01,
        -0.22,
        f"Residue Sequence range: {int(x_idx.min())}-{int(x_idx.max())}",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
    )

    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.show()
    plt.close(fig)

    out_table = pd.DataFrame(
        {
            "res_seq": x,
            "res_label": labels,
            "curv_raw": y,
            "curv_plot": y_plot,
            "curv_smooth": y_smooth,
            "is_top": np.isin(np.arange(len(x)), top_idx),
        }
    )
    out_table.to_csv(csv_path, index=False)

    return {
        "png": png_path,
        "pdf": pdf_path,
        "csv": csv_path,
        "n_residues": n,
        "k_top": k,
        "top_percent": cfg.top_percent,
        "frac_mass": frac_mass,
        "curv_col_used": res_df.attrs.get("curv_col_used", ""),
    }


def main() -> None:
    cfg = CurvConfig(
        curvature_kind="mean",
        use_magnitude=True,
        top_percent=5.0,
        smooth_window=9,
        clip_percentiles=(2.0, 98.0),
        weight_by_area=True,
        output_dir="curvature_panel_c1",
    )

    logs_file = choose_logs_file()

    logs_file_prefix = os.path.dirname(logs_file)

    logs_obj = read_logs2(logs_file, return_dict=False, all_=True, balls=True, surfs=False, edges=False, verts=False)
    atoms_df = logs_obj["atoms"]

    res_df_chain = build_residue_df_from_atoms(atoms_df, cfg)
    res_df = aggregate_across_chains(res_df_chain)

    # Sanity check: confirm repeats are actually 4 most of the time
    print("Repeat counts (n_repeats) value counts:")
    print(res_df["n_repeats"].value_counts().sort_index())

    molecule_name = str(logs_obj.get("data", {}).get("name", os.path.splitext(os.path.basename(logs_file))[0]))
    title = f"{molecule_name}: residue curvature localization"
    out_prefix = molecule_name.replace(" ", "_")

    result = plot_panel_c1(res_df, cfg, title=title, out_prefix=os.path.join(logs_file_prefix, out_prefix))

    if np.isfinite(result["frac_mass"]):
        print(
            f"[{molecule_name}] Using column '{result['curv_col_used']}'. "
            f"Top {result['top_percent']:.1f}% of residues (n={result['k_top']}/{result['n_residues']}) "
            f"contribute {100.0 * result['frac_mass']:.1f}% of total curvature mass."
        )
    else:
        print(f"[{molecule_name}] Curvature mass fraction could not be computed (total mass was zero or NaN).")

    print("Saved:")
    print(f"  PNG: {result['png']}")
    print(f"  PDF: {result['pdf']}")
    print(f"  CSV: {result['csv']}")
    if UNKNOWN_RESIDUES:
        print("\n[WARNING] Unmapped residue names encountered (mapped to 'X'):")
        for res, count in sorted(UNKNOWN_RESIDUES.items(), key=lambda x: -x[1]):
            print(f"  {res:>6s} : {count}")
    else:
        print("\nNo unmapped residue names encountered.")


if __name__ == "__main__":
    main()
