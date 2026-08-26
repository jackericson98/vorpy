"""Figure 3C — Residue-type AW vs Power deviations.

Creates a 2 x 3 residue-level comparison:

    Protein       : Volume | Surface Area | Contacts
    Nucleic acids : Volume | Surface Area | Contacts

Each faint point is one residue instance. Large markers show the mean signed
Power-vs-AW percent difference for each residue type, with SEM error bars.

The residue metrics come from F3_common:
    - Volume: sum of constituent atom volumes.
    - Surface Area: only Voronoi surfaces crossing the residue boundary.
    - Contacts: number of Voronoi surfaces crossing the residue boundary.

DNA and RNA residue names are kept distinct. For example, DA and A are not
silently merged.

Outputs
-------
F3C_residue_type_deviations.png
F3C_residue_type_deviations.svg
F3C_residue_type_summary.csv
"""
from __future__ import annotations

import os
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# VorPy imports
# ---------------------------------------------------------------------------

vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
if vorpy_root not in sys.path:
    sys.path.append(vorpy_root)

from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common import (
    add_deviation_columns,
    build_atomic_metrics,
    build_residue_metrics,
    discover_systems,
    read_pair,
)


# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

DATA_ROOT = None
FIGURE_DIR = None

EXCLUDE_KEYS = ["A", "B", "C"]

FIGSIZE = (19, 11)
DPI = 300
SHOW = True
SAVE_PNG = True
SAVE_SVG = True

# Raw residue points
POINT_SIZE = 12
POINT_ALPHA = 0.14
JITTER_WIDTH = 0.13
RANDOM_SEED = 17

# Mean ± SEM
MEAN_SIZE = 52
MEAN_EDGE_WIDTH = 0.9
ERROR_CAPSIZE = 3

# If a residue type has fewer observations than this, it is omitted.
MIN_RESIDUE_COUNT = 3

# Shared y-axis behavior.
# Percentile trimming only controls plotting limits; it does NOT remove data
# from means/SEMs or from the saved CSV.
Y_PERCENTILE = 99.5
MIN_ABS_Y_RANGE = {
    "Volume": 5.0,
    "Surface Area": 5.0,
    "Contacts": 10.0,
}

AA_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
]

DNA_ORDER = ["DA", "DC", "DG", "DT", "DI"]
RNA_ORDER = ["A", "C", "G", "U", "I", "RA", "RC", "RG", "RU"]
NUCLEIC_ORDER = DNA_ORDER + RNA_ORDER

PROTEIN_RESIDUES = set(AA_ORDER)
DNA_RESIDUES = set(DNA_ORDER)
RNA_RESIDUES = set(RNA_ORDER)

# Chemical-category colors for residue means.
# Raw points use the same color with low alpha.
AA_CATEGORY = {
    # hydrophobic / nonpolar
    "ALA": "hydrophobic", "VAL": "hydrophobic", "ILE": "hydrophobic",
    "LEU": "hydrophobic", "MET": "hydrophobic", "PHE": "hydrophobic",
    "TRP": "hydrophobic", "PRO": "hydrophobic",

    # polar uncharged
    "SER": "polar", "THR": "polar", "ASN": "polar",
    "GLN": "polar", "TYR": "polar", "CYS": "polar",

    # acidic
    "ASP": "acidic", "GLU": "acidic",

    # basic
    "LYS": "basic", "ARG": "basic", "HIS": "basic",

    # special
    "GLY": "special",
}

CATEGORY_COLORS = {
    "hydrophobic": "#7f7f7f",
    "polar": "#2ca02c",
    "acidic": "#d62728",
    "basic": "#1f77b4",
    "special": "#9467bd",
    "purine": "#ff7f0e",
    "pyrimidine": "#17becf",
    "other": "#8c564b",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def select_directory(title: str) -> str:
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    path = filedialog.askdirectory(title=title)
    root.destroy()
    return path


def classify_system(atoms: pd.DataFrame) -> Optional[str]:
    """Classify system by the dominant recognized biomolecular residue class."""
    if "Residue" not in atoms.columns:
        return None

    residues = atoms["Residue"].astype(str).str.strip().str.upper()

    counts = {
        "protein": int(residues.isin(PROTEIN_RESIDUES).sum()),
        "dna": int(residues.isin(DNA_RESIDUES).sum()),
        "rna": int(residues.isin(RNA_RESIDUES).sum()),
    }

    best = max(counts, key=counts.get)
    return best if counts[best] > 0 else None


def residue_category(residue: str, broad_class: str) -> str:
    residue = str(residue).upper()

    if broad_class == "Protein":
        return AA_CATEGORY.get(residue, "other")

    # Purines: adenine / guanine
    if residue in {"DA", "DG", "DI", "A", "G", "I", "RA", "RG"}:
        return "purine"

    # Pyrimidines: cytosine / thymine / uracil
    if residue in {"DC", "DT", "C", "U", "RC", "RU"}:
        return "pyrimidine"

    return "other"


def broad_class(system_class: str) -> str:
    if system_class == "protein":
        return "Protein"
    if system_class in {"dna", "rna"}:
        return "Nucleic acid"
    raise ValueError(system_class)


def retain_class_residues(
    residue_df: pd.DataFrame,
    system_class: str,
) -> pd.DataFrame:
    if system_class == "protein":
        allowed = PROTEIN_RESIDUES
    elif system_class == "dna":
        allowed = DNA_RESIDUES
    elif system_class == "rna":
        allowed = RNA_RESIDUES
    else:
        return residue_df.iloc[0:0].copy()

    return residue_df[
        residue_df["Residue"].astype(str).str.upper().isin(allowed)
    ].copy()


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_residue_data(
    data_root: str,
    exclude_keys: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    systems = discover_systems(
        data_root,
        exclude_keys=exclude_keys,
    )

    print(f"\nFound {len(systems)} valid Figure 3 system folders.")

    frames: List[pd.DataFrame] = []

    for system in systems:
        print(f"\nProcessing {system.name} ...")

        aw_logs, power_logs = read_pair(
            system,
            need_surfs=True,
        )

        system_class = classify_system(aw_logs["atoms"])
        if system_class is None:
            print("  skipped: no recognized protein/DNA/RNA residue class")
            continue

        atom_df = build_atomic_metrics(
            aw_logs,
            power_logs,
        )

        residue_df = build_residue_metrics(
            atom_df,
            aw_logs,
            power_logs,
        )

        if residue_df.empty:
            print("  skipped: no matched residue data")
            continue

        residue_df = retain_class_residues(
            residue_df,
            system_class=system_class,
        )

        if residue_df.empty:
            print("  skipped: no residues retained after class filtering")
            continue

        residue_df = add_deviation_columns(residue_df)
        residue_df["System"] = system.name
        residue_df["SystemClass"] = system_class
        residue_df["BroadClass"] = broad_class(system_class)
        residue_df["Category"] = [
            residue_category(res, broad_class(system_class))
            for res in residue_df["Residue"]
        ]

        print(
            f"  classified as {system_class}; "
            f"retained residues: {len(residue_df):,}"
        )

        frames.append(residue_df)

    if not frames:
        return pd.DataFrame()

    result = pd.concat(frames, ignore_index=True)

    print("\nCollected residue data:")
    for name, subdf in result.groupby("BroadClass"):
        print(f"  {name}: {len(subdf):,} residues")

    return result


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

METRIC_DIFF_COLUMNS = {
    "Volume": "Volume Signed % Diff",
    "Surface Area": "Surface Area Signed % Diff",
    "Contacts": "Contacts Signed % Diff",
}


def sem(values) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]

    if len(vals) <= 1:
        return 0.0

    return float(np.std(vals, ddof=1) / np.sqrt(len(vals)))


def summarize_residue_types(
    residue_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for broad_name in ["Protein", "Nucleic acid"]:
        class_df = residue_df[
            residue_df["BroadClass"] == broad_name
        ]

        for metric, diff_col in METRIC_DIFF_COLUMNS.items():
            for residue_name, group_df in class_df.groupby("Residue"):
                vals = pd.to_numeric(
                    group_df[diff_col],
                    errors="coerce",
                ).dropna().to_numpy(float)

                if len(vals) < MIN_RESIDUE_COUNT:
                    continue

                rows.append(
                    {
                        "BroadClass": broad_name,
                        "Metric": metric,
                        "Residue": residue_name,
                        "Category": residue_category(
                            residue_name,
                            broad_name,
                        ),
                        "Count": len(vals),
                        "Mean Signed % Diff": float(np.mean(vals)),
                        "Median Signed % Diff": float(np.median(vals)),
                        "SD": (
                            float(np.std(vals, ddof=1))
                            if len(vals) > 1
                            else 0.0
                        ),
                        "SEM": sem(vals),
                        "Mean Abs % Diff": float(np.mean(np.abs(vals))),
                    }
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def residue_order_for_class(
    residue_df: pd.DataFrame,
    broad_class_name: str,
) -> List[str]:
    present = set(
        residue_df.loc[
            residue_df["BroadClass"] == broad_class_name,
            "Residue",
        ].astype(str)
    )

    if broad_class_name == "Protein":
        return [x for x in AA_ORDER if x in present]

    ordered = [x for x in NUCLEIC_ORDER if x in present]

    # Preserve any unexpected nucleic labels at the end rather than dropping.
    extras = sorted(present - set(ordered))
    return ordered + extras


def get_shared_y_limits(
    residue_df: pd.DataFrame,
    metric: str,
) -> Tuple[float, float]:
    diff_col = METRIC_DIFF_COLUMNS[metric]

    vals = pd.to_numeric(
        residue_df[diff_col],
        errors="coerce",
    ).dropna().to_numpy(float)

    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return -1.0, 1.0

    tail = (100.0 - Y_PERCENTILE) / 2.0
    lo = float(np.percentile(vals, tail))
    hi = float(np.percentile(vals, 100.0 - tail))

    max_abs = max(
        abs(lo),
        abs(hi),
        float(MIN_ABS_Y_RANGE.get(metric, 1.0)),
    )

    pad = max_abs * 0.08
    return -(max_abs + pad), max_abs + pad


def plot_residue_metric(
    ax,
    residue_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    broad_class_name: str,
    metric: str,
    y_limits: Tuple[float, float],
    rng: np.random.Generator,
):
    diff_col = METRIC_DIFF_COLUMNS[metric]

    class_df = residue_df[
        residue_df["BroadClass"] == broad_class_name
    ].copy()

    order = residue_order_for_class(
        residue_df,
        broad_class_name,
    )

    if not order:
        ax.text(
            0.5,
            0.5,
            f"No {broad_class_name} data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return

    positions = {
        residue: i
        for i, residue in enumerate(order)
    }

    # ------------------------------------------------------------------
    # Individual residue observations.
    # ------------------------------------------------------------------
    for residue in order:
        group_df = class_df[
            class_df["Residue"] == residue
        ]

        vals = pd.to_numeric(
            group_df[diff_col],
            errors="coerce",
        ).dropna().to_numpy(float)

        if len(vals) == 0:
            continue

        category = residue_category(
            residue,
            broad_class_name,
        )
        color = CATEGORY_COLORS.get(
            category,
            CATEGORY_COLORS["other"],
        )

        x0 = positions[residue]
        xs = x0 + rng.uniform(
            -JITTER_WIDTH,
            JITTER_WIDTH,
            len(vals),
        )

        ax.scatter(
            xs,
            vals,
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            color=color,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

    # ------------------------------------------------------------------
    # Mean ± SEM.
    # ------------------------------------------------------------------
    panel_summary = summary_df[
        (summary_df["BroadClass"] == broad_class_name) &
        (summary_df["Metric"] == metric)
    ].copy()

    for _, row in panel_summary.iterrows():
        residue = row["Residue"]

        if residue not in positions:
            continue

        category = row["Category"]
        color = CATEGORY_COLORS.get(
            category,
            CATEGORY_COLORS["other"],
        )

        x = positions[residue]
        y = float(row["Mean Signed % Diff"])
        error = float(row["SEM"])

        ax.errorbar(
            x,
            y,
            yerr=error,
            fmt="o",
            markersize=np.sqrt(MEAN_SIZE),
            color=color,
            markeredgecolor="white",
            markeredgewidth=MEAN_EDGE_WIDTH,
            ecolor=color,
            elinewidth=1.4,
            capsize=ERROR_CAPSIZE,
            zorder=5,
        )

    ax.axhline(
        0.0,
        linestyle="--",
        linewidth=1.4,
        color="black",
        alpha=0.8,
        zorder=0,
    )

    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(*y_limits)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(
        order,
        rotation=45 if broad_class_name == "Protein" else 0,
        ha="right" if broad_class_name == "Protein" else "center",
        fontsize=8.5 if broad_class_name == "Protein" else 10,
    )

    ax.set_ylabel(
        "Power vs AW signed difference (%)",
        fontsize=10,
    )
    ax.set_title(
        metric,
        fontsize=13,
        fontweight="bold",
    )

    ax.tick_params(
        axis="y",
        labelsize=9,
        width=1.2,
        length=5,
    )
    ax.tick_params(
        axis="x",
        width=1.2,
        length=4,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def print_summary_table(summary_df: pd.DataFrame):
    if summary_df.empty:
        return

    for broad_name in ["Protein", "Nucleic acid"]:
        for metric in ["Volume", "Surface Area", "Contacts"]:
            panel = summary_df[
                (summary_df["BroadClass"] == broad_name) &
                (summary_df["Metric"] == metric)
            ].copy()

            if panel.empty:
                continue

            panel = panel.sort_values(
                "Mean Signed % Diff"
            )

            print()
            print("=" * 100)
            print(
                f"F3C RESIDUE TYPES — "
                f"{broad_name.upper()} | {metric.upper()}"
            )
            print("=" * 100)

            print(
                panel[
                    [
                        "Residue",
                        "Category",
                        "Count",
                        "Mean Signed % Diff",
                        "Median Signed % Diff",
                        "SD",
                        "SEM",
                        "Mean Abs % Diff",
                    ]
                ].to_string(
                    index=False,
                    formatters={
                        "Mean Signed % Diff": lambda x: f"{x:8.3f}",
                        "Median Signed % Diff": lambda x: f"{x:8.3f}",
                        "SD": lambda x: f"{x:8.3f}",
                        "SEM": lambda x: f"{x:8.3f}",
                        "Mean Abs % Diff": lambda x: f"{x:8.3f}",
                    },
                )
            )


def make_figure(
    residue_df: pd.DataFrame,
    figure_dir: str,
):
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_df = summarize_residue_types(
        residue_df
    )

    print_summary_table(summary_df)

    metrics = [
        "Volume",
        "Surface Area",
        "Contacts",
    ]
    rows = [
        "Protein",
        "Nucleic acid",
    ]

    shared_y = {
        metric: get_shared_y_limits(
            residue_df,
            metric,
        )
        for metric in metrics
    }

    print()
    print("F3C SHARED Y LIMITS")
    for metric, limits in shared_y.items():
        print(
            f"  {metric}: "
            f"{limits[0]:.3f} -> {limits[1]:.3f}"
        )

    fig, axes = plt.subplots(
        2,
        3,
        figsize=FIGSIZE,
        constrained_layout=False,
    )

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    for row_i, broad_name in enumerate(rows):
        for col_i, metric in enumerate(metrics):
            plot_residue_metric(
                axes[row_i, col_i],
                residue_df=residue_df,
                summary_df=summary_df,
                broad_class_name=broad_name,
                metric=metric,
                y_limits=shared_y[metric],
                rng=rng,
            )

        axes[row_i, 0].text(
            -0.22,
            0.5,
            broad_name,
            transform=axes[row_i, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=15,
            fontweight="bold",
        )

    fig.suptitle(
        "Figure 3C — Residue-level differences between AW and Power partitions",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )

    fig.subplots_adjust(
        left=0.075,
        right=0.99,
        top=0.92,
        bottom=0.10,
        wspace=0.25,
        hspace=0.35,
    )

    if SAVE_PNG:
        path = figure_dir / "F3C_residue_type_deviations.png"
        fig.savefig(
            path,
            dpi=DPI,
            bbox_inches="tight",
        )
        print(f"Saved: {path}")

    if SAVE_SVG:
        path = figure_dir / "F3C_residue_type_deviations.svg"
        fig.savefig(
            path,
            bbox_inches="tight",
        )
        print(f"Saved: {path}")

    summary_path = (
        figure_dir /
        "F3C_residue_type_summary.csv"
    )
    summary_df.to_csv(
        summary_path,
        index=False,
    )
    print(f"Saved: {summary_path}")

    if SHOW:
        plt.show()

    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(
    data_root: Optional[str] = DATA_ROOT,
    figure_dir: Optional[str] = FIGURE_DIR,
):
    if data_root is None:
        data_root = select_directory(
            "Select Figure 3 data folder"
        )

    if not data_root:
        print("No data folder selected.")
        return

    if figure_dir is None:
        figure_dir = select_directory(
            "Select figures/Figure_3 folder"
        )

    if not figure_dir:
        print("No figure output folder selected.")
        return

    residue_df = collect_residue_data(
        data_root=data_root,
        exclude_keys=EXCLUDE_KEYS,
    )

    if residue_df.empty:
        print("No residue data were collected.")
        return

    make_figure(
        residue_df,
        figure_dir=figure_dir,
    )


if __name__ == "__main__":
    main()