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
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z", "XAA": "X",
}

DNA_3_TO_BASE = {"DA": "A", "DG": "G", "DC": "C", "DT": "T"}
RNA_1_TO_BASE = {"A": "A", "G": "G", "C": "C", "U": "U"}

# Some files may store thymine as just "T"
DNA_ALT = {"T": "T"}  # treated as DNA thymine by default


def residue_to_token(res: str) -> str:
    """
    Returns an unambiguous label token:
      - Protein AAs: 'A', 'R', ...
      - DNA bases:   'dA', 'dC', 'dG', 'dT'
      - RNA bases:   'rA', 'rC', 'rG', 'rU'
      - Unknown:     'X'
    """
    r = str(res).strip().upper()

    # DNA (explicit)
    if r in DNA_3_TO_BASE:
        return "d" + DNA_3_TO_BASE[r]

    # DNA (alternate thymine encoding)
    if r in DNA_ALT:
        return "d" + DNA_ALT[r]

    # RNA (single-letter bases)
    if r in RNA_1_TO_BASE:
        return "r" + RNA_1_TO_BASE[r]

    # Protein (3-letter)
    if r in AA3_TO_AA1:
        return AA3_TO_AA1[r]

    # Fallback: try first 3 letters for protein
    if len(r) >= 3 and r[:3] in AA3_TO_AA1:
        return AA3_TO_AA1[r[:3]]

    return "X"


def token_class(token: str) -> str:
    """
    protein vs nucleic for plotting order + separator.
    """
    t = str(token)
    if t.startswith("d") or t.startswith("r"):
        return "nucleic"
    return "protein"


def print_token_outliers_iqr(
    res_inst: pd.DataFrame,
    token: str = "A",
    curv_col: str = "curv",
    k: float = 1.5,
    max_rows: int = 50,
    save_csv_path: str = None,
) -> None:
    """
    Print IQR-based outliers for a given residue token (e.g., 'A' for alanine).

    Outlier rule:
        curv < Q1 - k*IQR  OR  curv > Q3 + k*IQR
    """

    if curv_col not in res_inst.columns:
        raise KeyError(f"Expected curvature column '{curv_col}' not found. Columns: {list(res_inst.columns)}")

    if "token" not in res_inst.columns:
        raise KeyError("Expected column 'token' not found in res_inst. Did build_residue_instances_from_atoms add it?")

    df = res_inst[res_inst["token"] == token].copy()

    if df.empty:
        print(f"[OUTLIERS] No rows found for token='{token}'.")
        return

    df[curv_col] = pd.to_numeric(df[curv_col], errors="coerce")
    df = df.dropna(subset=[curv_col])

    if df.empty:
        print(f"[OUTLIERS] Token='{token}' has no finite curvature values in '{curv_col}'.")
        return

    q1 = float(df[curv_col].quantile(0.25))
    q3 = float(df[curv_col].quantile(0.75))
    iqr = q3 - q1

    lo = q1 - k * iqr
    hi = q3 + k * iqr

    out = df[(df[curv_col] < lo) | (df[curv_col] > hi)].copy()
    out = out.sort_values(curv_col)

    print(f"\n[OUTLIERS] token='{token}' using '{curv_col}' (k={k})")
    print(f"  Q1={q1:.6f}, Q3={q3:.6f}, IQR={iqr:.6f}")
    print(f"  Bounds: [{lo:.6f}, {hi:.6f}]")
    print(f"  Found {len(out)} outliers (of {len(df)} total '{token}' residues)\n")

    if out.empty:
        return

    # Prefer a rich but safe set of columns if they exist
    preferred_cols = [
        "Chain",
        "Residue Sequence",
        "Residue",
        "token",
        curv_col,
        "n_atoms",
        "area_sum",
        "area",
        "Surface Area",
        "Complete Cell?",
        "Number of Neighbors",
    ]
    cols = [c for c in preferred_cols if c in out.columns]

    if not cols:
        cols = [curv_col, "token"]

    print(out[cols].head(max_rows).to_string(index=False))

    if len(out) > max_rows:
        print(f"... ({len(out) - max_rows} more)")

    if save_csv_path is not None:
        out.to_csv(save_csv_path, index=False)
        print(f"\n[OUTLIERS] Saved outliers CSV: {save_csv_path}")



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
    curvature_kind: str = "mean"             # "mean" or "gauss"
    use_magnitude: bool = True               # plot |curv|
    weight_by_area: bool = True              # residue curvature = area-weighted across atoms
    group_label: str = "token"               # "token" (disambiguated) or "aa1" or "res3"
    min_count_per_group: int = 5             # filter rare residue types
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

    out["token"] = out["Residue"].apply(residue_to_token)
    out["token_class"] = out["token"].apply(token_class)

    out.attrs["curv_label"] = (
        f"Res Curv Dist"
        if cfg.weight_by_area
        else f"{curv_name} residue curvature"
    )
    out.attrs["curv_col_used"] = curv_col

    return out


def apply_percentile_cutoff_per_group(
    df: pd.DataFrame,
    group_col: str,
    value_col: str = "curv",
    low: float = 2.5,
    high: float = 97.5,
    min_n: int = 20,
) -> pd.DataFrame:
    """
    Apply a percentile cutoff (e.g. 2.5–97.5%) per group.
    Groups with fewer than min_n samples are left untouched.
    """

    keep_frames = []

    for g, sub in df.groupby(group_col):
        vals = sub[value_col].to_numpy(dtype=float)
        vals = vals[~np.isnan(vals)]

        if vals.size < min_n:
            keep_frames.append(sub)
            continue

        lo = np.percentile(vals, low)
        hi = np.percentile(vals, high)

        trimmed = sub[(sub[value_col] >= lo) & (sub[value_col] <= hi)]
        keep_frames.append(trimmed)

    return pd.concat(keep_frames, ignore_index=True)


def apply_classwise_y_clipping(
    df: pd.DataFrame,
    value_col: str = "curv",
    token_col: str = "token",
    protein_limits: tuple[float, float] = (0.20, 0.37),
    nucleic_limits: tuple[float, float] = (0.22, 0.30),
) -> pd.DataFrame:
    """
    Clip curvature values by residue class (protein vs nucleic).
    """

    df = df.copy()

    is_protein = df[token_col].apply(token_class) == "protein"
    is_nucleic = df[token_col].apply(token_class) == "nucleic"

    df.loc[is_protein, value_col] = df.loc[is_protein, value_col].clip(
        lower=protein_limits[0],
        upper=protein_limits[1],
    )

    df.loc[is_nucleic, value_col] = df.loc[is_nucleic, value_col].clip(
        lower=nucleic_limits[0],
        upper=nucleic_limits[1],
    )

    return df


def plot_violin_by_residue_type(
    res_inst: pd.DataFrame,
    cfg: ViolinConfig,
    title: str,
    out_prefix: str,
) -> Dict[str, str]:
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.size": 24,
        "axes.titlesize": 28,
        "axes.labelsize": 26,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "legend.fontsize": 22,
        "figure.titlesize": 30,
    })

    if cfg.group_label == "res3":
        group_col = "Residue"
        group_title = "Residue (raw)"
    elif cfg.group_label == "token":
        group_col = "token"
        group_title = "Residue (protein vs DNA/RNA)"
    else:
        group_col = "aa1"
        group_title = "Residue"


    # Filter unknowns and NaNs
    df = res_inst.copy()
    df = df.dropna(subset=["curv"]).copy()
    df[group_col] = df[group_col].astype(str)

    df = res_inst.copy()
    df = df.dropna(subset=["curv"]).copy()
    df[group_col] = df[group_col].astype(str)

    # NEW: 95% cutoff per residue group
    df = apply_percentile_cutoff_per_group(
        df,
        group_col=group_col,
        value_col="curv",
        low=2.5,
        high=97.5,
        min_n=20,
    )


    # NEW: classwise y-range enforcement
    df = apply_classwise_y_clipping(
        df,
        value_col="curv",
        token_col=group_col,
        protein_limits=(0.20, 0.37),
        nucleic_limits=(0.22, 0.30),
    )

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

    # Enforce protein-left / nucleic-right ordering
    protein_labels = [k for k in labels if token_class(k) == "protein"]
    nucleic_labels = [k for k in labels if token_class(k) == "nucleic"]

    labels = protein_labels + nucleic_labels
    boundary = len(protein_labels)  # divider goes after last protein violin


    data = [grouped[k] for k in labels]
    global_min = min(np.min(vals) for vals in data)
    ymin = 0.95 * global_min

    fig_w = max(8.0, 0.45 * len(labels))
    fig = plt.figure(figsize=(fig_w, 5.0))
    ax = fig.add_subplot(111)

    parts = ax.violinplot(
        data,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )
    ax.set_ylim(bottom=ymin)

    ax.set_title(title)
    ax.set_xlabel(group_title)
    ax.set_ylabel(res_inst.attrs.get("curv_label", "Mean res curv"))

    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=0)

    ax.set_ylim(0.20, 0.37)

    if boundary > 0 and boundary < len(labels):
        ax.axvline(boundary + 0.5, linewidth=2.5, linestyle="--")
        ax.text(
            boundary + 0.5,
            ax.get_ylim()[1],
            "protein | nucleic",
            ha="center",
            va="bottom",
            fontsize=16,
        )


    ax.tick_params(axis="x", labelsize=mpl.rcParams["xtick.labelsize"])
    ax.tick_params(axis="y", labelsize=mpl.rcParams["ytick.labelsize"])
    ax.tick_params(
        axis="both",
        which="major",
        length=10,
        width=2.5,
    )
    ax.tick_params(top=False, right=False)

    for spine in ax.spines.values():
        spine.set_linewidth(2.5)

    # Add counts under each label (optional but helpful)
    counts = [grouped[k].size for k in labels]

    for i, (vals, n) in enumerate(zip(data, counts), start=1):
        y_min = min(vals) * 0.99
        ax.text(
            i,
            y_min,  # slightly below the median
            f"{n}",
            ha="center",
            va="top",
            fontsize=10,
            color="black",
        )

    ax.grid(True, axis="y", alpha=0.3)

    os.makedirs(cfg.output_dir, exist_ok=True)
    png_path = os.path.join(cfg.output_dir, f"{out_prefix}_violin_by_{group_col}.png")
    pdf_path = os.path.join(cfg.output_dir, f"{out_prefix}_violin_by_{group_col}.pdf")
    svg_path = os.path.join(cfg.output_dir, f"{out_prefix}_violin_by_{group_col}.svg")
    csv_path = os.path.join(cfg.output_dir, f"{out_prefix}_residue_instances.csv")

    fig.tight_layout()
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.savefig(svg_path, bbox_inches="tight")
    plt.show()
    plt.close(fig)

    # Save underlying data for reproducibility
    df.to_csv(csv_path, index=False)

    return {"png": png_path, "pdf": pdf_path, "svg": svg_path, "csv": csv_path}


def main() -> None:
    cfg = ViolinConfig(
        curvature_kind="mean",
        use_magnitude=True,
        weight_by_area=True,
        group_label="token",
        min_count_per_group=5,
        sort_by="median",
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

    # Print alanine outliers (protein alanine token is 'A')
    outlier_csv = os.path.join(cfg.output_dir, "outliers_token_A.csv")
    print_token_outliers_iqr(
        res_inst=res_inst,
        token="A",
        curv_col="curv",
        k=1.5,
        max_rows=50,
        save_csv_path=outlier_csv,
    )


    molecule_name = str(
        logs_obj.get("data", {}).get("name", os.path.splitext(os.path.basename(logs_file))[0])
    )
    out_prefix = sanitize_filename(molecule_name)

    result = plot_violin_by_residue_type(
        res_inst=res_inst,
        cfg=cfg,
        title=f"{molecule_name}",
        out_prefix=out_prefix,
    )

    print(f"[{molecule_name}] Using curvature column: {res_inst.attrs.get('curv_col_used', '')}")
    print("Saved:")
    print(f"  PNG: {result['png']}")
    print(f"  PDF: {result['pdf']}")
    print(f"  SVG: {result['svg']}")
    print(f"  CSV: {result['csv']}")


if __name__ == "__main__":
    main()
