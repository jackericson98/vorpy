"""VorPy atom-level correlation analysis."""

from __future__ import annotations

import re
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# =============================================================================
# USER SETTINGS
# =============================================================================

# Turn metrics on/off for the correlation matrix.
METRICS = {
    "Radius": True,
    "Volume": True,
    "Surface Area": True,
    "Maximum Mean Curvature": True,
    "Average Mean Surface Curvature": True,
    "Maximum Gaussian Curvature": True,
    "Average Gaussian Surface Curvature": True,
    "Integrated Mean Curvature": True,
    "Integrated Mean Curvature Squared": True,
    "Integrated Gaussian Curvature": True,
    "Representative Surface Energy": True,
}

# Pair to show in the detailed scatter plot.
X_METRIC = "Surface Area"
Y_METRIC = "Integrated Mean Curvature"

SAVE_PNG = True
SAVE_SVG = True
SHOW = True
DPI = 300
SHOW_LINEAR_FIT = True


# =============================================================================
# VORPY IMPORT
# =============================================================================

def find_vorpy_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "vorpy").is_dir():
            return candidate
    raise RuntimeError("Could not locate the VorPy repository root.")


SCRIPT_DIR = Path(__file__).resolve().parent
VORPY_ROOT = find_vorpy_root(SCRIPT_DIR)

if str(VORPY_ROOT) not in sys.path:
    sys.path.insert(0, str(VORPY_ROOT))

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


# =============================================================================
# FILE SELECTION
# =============================================================================

def select_log_files() -> list[Path]:
    """
    Choose VorPy aw_logs.csv files from one or more folders.

    The picker reopens after every successful selection. Select one or several
    logs from the current folder, then choose logs from another folder on the
    next popup. Press Cancel (or close the picker) when finished.

    Duplicate paths are removed while preserving selection order.
    """

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)

    selected_files: list[Path] = []
    seen: set[Path] = set()
    selection_number = 1

    while True:
        files = filedialog.askopenfilenames(
            parent=root,
            title=(
                f"Select aw_logs.csv file(s) - round {selection_number} "
                f"({len(selected_files)} selected). Cancel when finished."
            ),
            filetypes=[
                ("VorPy logs", "aw_logs.csv"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )

        if not files:
            break

        added = 0

        for file_path in files:
            path = Path(file_path).resolve()

            if path in seen:
                print(f"Already selected; skipping duplicate: {path}")
                continue

            seen.add(path)
            selected_files.append(path)
            added += 1

        print(
            f"Selection round {selection_number}: "
            f"added {added} file(s); "
            f"{len(selected_files)} total."
        )

        selection_number += 1

    root.destroy()
    return selected_files


# =============================================================================
# DATA READING
# =============================================================================

def read_atoms(path: Path) -> pd.DataFrame:
    """Read the atom table from one VorPy log."""
    result = read_logs2(str(path))
    atoms = result.get("atoms")

    if atoms is None:
        raise ValueError("read_logs2 did not return an 'atoms' table.")

    if not isinstance(atoms, pd.DataFrame):
        atoms = pd.DataFrame(atoms)

    atoms = atoms.copy()

    if atoms.empty:
        raise ValueError("Atom table is empty.")

    atoms.insert(0, "Source Log", str(path))
    atoms.insert(1, "Source File", path.name)
    return atoms


def combine_logs(paths: list[Path]) -> pd.DataFrame:
    tables = []

    print("\n" + "=" * 80)
    print("READING ATOMIC DATA")
    print("=" * 80)

    for i, path in enumerate(paths, 1):
        print(f"[{i:>3}/{len(paths)}] {path}")
        try:
            table = read_atoms(path)
            print(f"    {len(table):,} atoms")
            tables.append(table)
        except Exception as exc:
            print(f"    ERROR: {exc}")

    if not tables:
        return pd.DataFrame()

    return pd.concat(tables, ignore_index=True, sort=False)


# =============================================================================
# METRICS / CORRELATION
# =============================================================================

def enabled_metrics() -> list[str]:
    return [name for name, enabled in METRICS.items() if enabled]


def prepare_metrics(
    atoms: pd.DataFrame,
    requested: list[str],
) -> tuple[pd.DataFrame, list[str]]:

    available = []

    print("\nMetric availability:")

    for metric in requested:
        if metric not in atoms.columns:
            print(f"  MISSING: {metric}")
            continue

        atoms[metric] = pd.to_numeric(atoms[metric], errors="coerce")
        n = int(atoms[metric].notna().sum())

        if n < 2:
            print(f"  INSUFFICIENT: {metric} ({n} values)")
            continue

        available.append(metric)
        print(f"  OK: {metric} ({n:,} values)")

    return atoms, available


def correlation_matrix(
    atoms: pd.DataFrame,
    metrics: list[str],
) -> pd.DataFrame:
    return atoms[metrics].corr(method="pearson", min_periods=2)


# =============================================================================
# CORRELATION MATRIX PLOT
# =============================================================================

def plot_matrix(corr: pd.DataFrame, output_dir: Path):
    n = len(corr)
    size = max(8.0, 0.9 * n + 3.0)

    fig, ax = plt.subplots(figsize=(size, size))

    image = ax.imshow(
        corr.to_numpy(float),
        vmin=-1,
        vmax=1,
        aspect="equal",
    )

    labels = list(corr.columns)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    ax.set_title(
        "Atomic Metric Pearson Correlation Matrix",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )

    for row in range(n):
        for col in range(n):
            value = corr.iloc[row, col]
            label = f"{value:.2f}" if np.isfinite(value) else "NA"
            ax.text(
                col,
                row,
                label,
                ha="center",
                va="center",
                fontsize=8,
            )

    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r", rotation=270, labelpad=18)

    fig.tight_layout()

    if SAVE_PNG:
        path = output_dir / "atomic_correlation_matrix.png"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {path}")

    if SAVE_SVG:
        path = output_dir / "atomic_correlation_matrix.svg"
        fig.savefig(path, bbox_inches="tight")
        print(f"Saved: {path}")

    if SHOW:
        plt.show()
    else:
        plt.close(fig)


# =============================================================================
# SELECTED PAIR SCATTER PLOT
# =============================================================================

def safe_name(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def plot_pair(
    atoms: pd.DataFrame,
    x_metric: str,
    y_metric: str,
    output_dir: Path,
) -> dict | None:

    pair = atoms[[x_metric, y_metric]].copy()
    pair[x_metric] = pd.to_numeric(pair[x_metric], errors="coerce")
    pair[y_metric] = pd.to_numeric(pair[y_metric], errors="coerce")
    pair = pair.replace([np.inf, -np.inf], np.nan).dropna()

    if len(pair) < 2:
        print(f"Not enough paired values for {x_metric} vs {y_metric}.")
        return None

    x = pair[x_metric].to_numpy(float)
    y = pair[y_metric].to_numpy(float)

    if np.std(x) == 0 or np.std(y) == 0:
        r = np.nan
    else:
        r = float(np.corrcoef(x, y)[0, 1])

    if len(np.unique(x)) >= 2:
        slope, intercept = np.polyfit(x, y, 1)
    else:
        slope = intercept = np.nan

    r2 = r * r if np.isfinite(r) else np.nan

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.scatter(x, y, s=16, alpha=0.45)

    if SHOW_LINEAR_FIT and np.isfinite(slope):
        x_fit = np.linspace(x.min(), x.max(), 200)
        ax.plot(
            x_fit,
            slope * x_fit + intercept,
            linewidth=1.5,
            label="Linear fit",
        )
        ax.legend(frameon=False)

    ax.set_title(
        f"Atomic {x_metric} vs {y_metric}",
        fontsize=14,
        fontweight="bold",
    )
    ax.set_xlabel(x_metric)
    ax.set_ylabel(y_metric)
    ax.grid(alpha=0.25)

    r_text = f"{r:.4f}" if np.isfinite(r) else "NA"
    r2_text = f"{r2:.4f}" if np.isfinite(r2) else "NA"

    ax.text(
        0.03,
        0.97,
        f"N = {len(pair):,}\nPearson r = {r_text}\nR² = {r2_text}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )

    fig.tight_layout()

    base = f"atomic_scatter_{safe_name(x_metric)}_vs_{safe_name(y_metric)}"

    if SAVE_PNG:
        path = output_dir / f"{base}.png"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {path}")

    if SAVE_SVG:
        path = output_dir / f"{base}.svg"
        fig.savefig(path, bbox_inches="tight")
        print(f"Saved: {path}")

    if SHOW:
        plt.show()
    else:
        plt.close(fig)

    return {
        "X Metric": x_metric,
        "Y Metric": y_metric,
        "N": len(pair),
        "Pearson r": r,
        "R Squared": r2,
        "Slope": slope,
        "Intercept": intercept,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "=" * 80)
    print("VORPY ATOMIC CORRELATION ANALYSIS")
    print("=" * 80)

    paths = select_log_files()

    if not paths:
        print("No files selected.")
        return

    atoms = combine_logs(paths)

    if atoms.empty:
        print("No atomic data could be read.")
        return

    atoms, metrics = prepare_metrics(atoms, enabled_metrics())

    if len(metrics) < 2:
        print("Fewer than two usable metrics were found.")
        return

    output_dir = paths[0].parent / "atomic_correlation_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    combined_path = output_dir / "atomic_correlation_combined_data.csv"
    atoms.to_csv(combined_path, index=False)
    print(f"\nSaved: {combined_path}")

    corr = correlation_matrix(atoms, metrics)
    corr_path = output_dir / "atomic_correlation_matrix.csv"
    corr.to_csv(corr_path)
    print(f"Saved: {corr_path}")

    plot_matrix(corr, output_dir)

    print("\n" + "=" * 80)
    print("CORRELATION MATRIX")
    print("=" * 80)
    print(corr.round(4).to_string())

    if X_METRIC not in metrics:
        print(f"\nConfigured X_METRIC is unavailable: {X_METRIC}")
        return

    if Y_METRIC not in metrics:
        print(f"\nConfigured Y_METRIC is unavailable: {Y_METRIC}")
        return

    stats = plot_pair(
        atoms,
        X_METRIC,
        Y_METRIC,
        output_dir,
    )

    if stats is not None:
        stats_path = output_dir / "atomic_selected_pair_statistics.csv"
        pd.DataFrame([stats]).to_csv(stats_path, index=False)

        print("\n" + "=" * 80)
        print("SELECTED CORRELATION")
        print("=" * 80)
        print(f"X:         {X_METRIC}")
        print(f"Y:         {Y_METRIC}")
        print(f"N:         {stats['N']:,}")
        print(f"Pearson r: {stats['Pearson r']:.8f}")
        print(f"R squared: {stats['R Squared']:.8f}")
        print(f"Saved:     {stats_path}")

    print("\nResults:")
    print(f"  {output_dir}")


if __name__ == "__main__":
    main()
