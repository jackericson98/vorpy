"""
VorPy Frame Variation Analysis
==============================

Plots group-level geometric quantities across simulation frames.

Expected directory structure
----------------------------

Selected root/
    frames_new/
        f1/
            aw/
                aw/
                    aw_logs.csv
        f2/
            aw/
                aw/
                    aw_logs.csv
        f3/
            aw/
                aw/
                    aw_logs.csv
        ...

The script recursively searches beneath the selected directory, so the
selected directory can be frames_new itself or a directory above it.

Default enabled metrics
-----------------------
Volume
Surface Area
Integrated Mean Curvature
Integrated Mean Curvature Squared
Integrated Gaussian Curvature

Outputs
-------
frame_variation_data.csv
frame_variation_summary.csv
one PNG and SVG per enabled metric
"""

from __future__ import annotations

import os
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

# Metrics from the "Group Information" section of the VorPy log.
#
# Turn individual metrics on/off here. Every enabled metric is:
#   1. retained as a raw value,
#   2. converted to signed % difference from its first valid frame,
#   3. saved to the output CSV, and
#   4. plotted in its own figure.
#
# Percentage difference:
#
#   100 * (frame_value - first_frame_value) / abs(first_frame_value)
#
METRICS = {
    "Volume": True,
    "Surface Area": True,
    "Integrated Mean Curvature": True,
    "Integrated Mean Curvature Squared": True,
    "Integrated Gaussian Curvature": True,

    # Topology / Gauss-Bonnet diagnostics
    "Euler Characteristic": False,
    "Gauss-Bonnet Expected": False,
    "Gauss-Bonnet Error": False,
    "Gauss-Bonnet Relative Error": False,

    # Enable if desired/present in the group data
    "Representative Surface Energy": False,
}

ENABLED_METRICS = [
    metric
    for metric, enabled in METRICS.items()
    if enabled
]

# Save figures.
SAVE_PNG = True
SAVE_SVG = True

# Show plots interactively.
SHOW = True

# Figure quality.
DPI = 300
FIGSIZE = (11, 6)

# Add a horizontal mean line.
SHOW_MEAN = False

# Add +/- 1 standard deviation around the mean.
SHOW_STD = False


# =============================================================================
# VORPY IMPORT
# =============================================================================

# This assumes the script lives somewhere inside the VorPy repository.
#
# We search upward until we find the directory containing the "vorpy" package.
# This makes the script less dependent on exactly how deep it sits in
# src/analyze/.

def find_vorpy_root(start: Path) -> Path:
    """Find the repository root containing the vorpy package."""

    start = start.resolve()

    for candidate in [start, *start.parents]:
        if (candidate / "vorpy").is_dir():
            return candidate

    raise RuntimeError(
        "Could not locate the VorPy repository root.\n"
        "Make sure this script is located somewhere inside the VorPy repository."
    )


SCRIPT_DIR = Path(__file__).resolve().parent
VORPY_ROOT = find_vorpy_root(SCRIPT_DIR)

if str(VORPY_ROOT) not in sys.path:
    sys.path.insert(0, str(VORPY_ROOT))


# -------------------------------------------------------------------------
# IMPORTANT
#
# Change this import only if read_logs2.py lives somewhere different in
# your current repository.
#
# Based on your existing analysis structure, this is the expected import.
# -------------------------------------------------------------------------

try:
    from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2

except ImportError:
    # Fallback: attempt common historical location.
    try:
        from vorpy.src.analyze.read_logs2 import read_logs2

    except ImportError as exc:
        raise ImportError(
            "\nCould not import read_logs2.\n\n"
            "Update the read_logs2 import near the top of this script to match "
            "its actual location in your VorPy repository.\n"
        ) from exc


# =============================================================================
# DIRECTORY SELECTION
# =============================================================================

def select_directory(title: str) -> str:
    """Open a Tkinter directory-selection dialog."""

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)

    path = filedialog.askdirectory(title=title)

    root.destroy()

    return path


# =============================================================================
# FRAME / PATH HELPERS
# =============================================================================

def frame_number_from_path(path: Path) -> int | None:
    """
    Extract the numeric frame number from a directory named f<number>.

    Examples
    --------
    f1      -> 1
    f20     -> 20
    f1000   -> 1000
    """

    for part in reversed(path.parts):
        match = re.fullmatch(r"f(\d+)", part, flags=re.IGNORECASE)

        if match:
            return int(match.group(1))

    return None


def frame_name_from_path(path: Path) -> str | None:
    """Return the f<number> directory name associated with a log."""

    for part in reversed(path.parts):
        if re.fullmatch(r"f\d+", part, flags=re.IGNORECASE):
            return part

    return None


def find_frame_logs(root: Path) -> list[Path]:
    """
    Recursively locate frame AW logs.

    Required ending:

        f#/aw/aw/aw_logs.csv

    Files elsewhere with the same filename are ignored.
    """

    logs = []

    for path in root.rglob("aw_logs.csv"):

        parts_lower = [part.lower() for part in path.parts]

        # Need:
        #
        # f#/aw/aw/aw_logs.csv

        if len(parts_lower) < 4:
            continue

        if parts_lower[-1] != "aw_logs.csv":
            continue

        if parts_lower[-2] != "aw":
            continue

        if parts_lower[-3] != "aw":
            continue

        frame_dir = path.parts[-4]

        if not re.fullmatch(r"f\d+", frame_dir, flags=re.IGNORECASE):
            continue

        logs.append(path)

    # Natural numeric frame ordering.
    logs.sort(
        key=lambda p: (
            frame_number_from_path(p)
            if frame_number_from_path(p) is not None
            else float("inf")
        )
    )

    return logs


# =============================================================================
# DATA READING
# =============================================================================

def read_frame_log(log_path: Path) -> dict:
    """
    Read one VorPy frame log using read_logs2.

    Returns the complete read_logs2 result.
    """

    return read_logs2(str(log_path))


def collect_frame_data(
    log_files: list[Path],
    metrics: list[str],
) -> pd.DataFrame:
    """
    Read requested group-level metrics from every frame.
    """

    rows = []

    print()
    print("=" * 80)
    print("READING FRAME LOGS")
    print("=" * 80)

    for i, log_path in enumerate(log_files, start=1):

        frame_number = frame_number_from_path(log_path)
        frame_name = frame_name_from_path(log_path)

        print(
            f"[{i:>4}/{len(log_files)}] "
            f"{frame_name or 'unknown'}"
        )

        try:
            logs = read_frame_log(log_path)

        except Exception as exc:
            print(f"    ERROR reading log: {exc}")
            continue

        group_data = logs.get("group data", {})

        row = {
            "Frame": frame_name,
            "Frame Number": frame_number,
            "Log Path": str(log_path),
        }

        missing = []

        for metric in metrics:

            if metric in group_data:
                row[metric] = group_data[metric]

            else:
                row[metric] = np.nan
                missing.append(metric)

        if missing:
            print(
                "    Missing metric(s): "
                + ", ".join(missing)
            )

        rows.append(row)

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # Make absolutely sure frames are ordered numerically.
    df = df.sort_values(
        "Frame Number",
        na_position="last",
    ).reset_index(drop=True)

    return df


# =============================================================================
# PERCENT DIFFERENCE FROM FIRST FRAME
# =============================================================================

def add_percent_difference_columns(
    df: pd.DataFrame,
    metrics: list[str],
) -> pd.DataFrame:
    """
    Add signed percentage-difference columns relative to the first valid
    frame for each metric.

    The first valid frame for a metric is assigned 0%. Raw values are kept.
    """

    df = df.copy()

    for metric in metrics:

        if metric not in df.columns:
            continue

        values = pd.to_numeric(
            df[metric],
            errors="coerce",
        )

        valid_values = values.dropna()

        if valid_values.empty:
            print(f"WARNING: No valid values available for {metric}.")
            continue

        baseline = float(valid_values.iloc[0])
        percent_column = f"{metric} % Difference From First Frame"

        if baseline == 0:
            df[percent_column] = np.nan
            print(
                f"WARNING: Cannot calculate percentage difference for "
                f"{metric}; first valid frame is zero."
            )
            continue

        df[percent_column] = (
            100.0
            * (values - baseline)
            / abs(baseline)
        )

    return df


# =============================================================================
# STATISTICS
# =============================================================================

def calculate_metric_statistics(
    df: pd.DataFrame,
    metric: str,
) -> dict:
    """Calculate frame-to-frame variation statistics."""

    values = pd.to_numeric(
        df[metric],
        errors="coerce",
    ).dropna()

    n = len(values)

    if n == 0:
        return {
            "Metric": metric,
            "N": 0,
            "Mean": np.nan,
            "Standard Deviation": np.nan,
            "Minimum": np.nan,
            "Maximum": np.nan,
            "Range": np.nan,
            "CV (%)": np.nan,
            "Range / Mean (%)": np.nan,
        }

    mean = float(values.mean())

    if n > 1:
        std = float(values.std(ddof=1))
    else:
        std = 0.0

    minimum = float(values.min())
    maximum = float(values.max())
    value_range = maximum - minimum

    if mean != 0:
        cv = 100.0 * std / abs(mean)
        range_percent = 100.0 * value_range / abs(mean)
    else:
        cv = np.nan
        range_percent = np.nan

    return {
        "Metric": metric,
        "N": n,
        "Mean": mean,
        "Standard Deviation": std,
        "Minimum": minimum,
        "Maximum": maximum,
        "Range": value_range,
        "CV (%)": cv,
        "Range / Mean (%)": range_percent,
    }


def build_summary(
    df: pd.DataFrame,
    metrics: list[str],
) -> pd.DataFrame:
    """Create summary statistics for all requested metrics."""

    rows = []

    for metric in metrics:

        if metric not in df.columns:
            continue

        rows.append(
            calculate_metric_statistics(
                df=df,
                metric=metric,
            )
        )

    return pd.DataFrame(rows)


# =============================================================================
# PLOTTING
# =============================================================================

def safe_filename(text: str) -> str:
    """Convert a metric name into a filesystem-friendly name."""

    text = text.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def plot_metric(
    df: pd.DataFrame,
    metric: str,
    output_dir: Path,
):
    """
    Plot one metric across frames.
    """

    plot_df = df[
        ["Frame Number", metric]
    ].copy()

    plot_df["Frame Number"] = pd.to_numeric(
        plot_df["Frame Number"],
        errors="coerce",
    )

    plot_df[metric] = pd.to_numeric(
        plot_df[metric],
        errors="coerce",
    )

    plot_df = plot_df.dropna()

    if plot_df.empty:
        print(f"No valid values found for: {metric}")
        return

    plot_df = plot_df.sort_values("Frame Number")

    x = plot_df["Frame Number"].to_numpy(float)
    raw_y = plot_df[metric].to_numpy(float)

    # ---------------------------------------------------------------------
    # Convert values to percentage difference from the first frame
    # ---------------------------------------------------------------------

    baseline = raw_y[0]

    if baseline == 0:
        print(
            f"Cannot calculate percentage difference for {metric}: "
            "the first frame has a value of zero."
        )
        return

    y = 100.0 * (raw_y - baseline) / abs(baseline)

    mean = float(np.mean(y))

    if len(y) > 1:
        std = float(np.std(y, ddof=1))
    else:
        std = 0.0

    # ---------------------------------------------------------------------
    # Plot
    # ---------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=FIGSIZE,
    )

    # Main frame-to-frame series.
    ax.plot(
        x,
        y,
        marker="o",
        markersize=4,
        linewidth=1.25,
        label=metric,
    )

    # Mean.
    if SHOW_MEAN:

        ax.axhline(
            mean,
            linestyle="--",
            linewidth=1.25,
            label=f"Mean = {mean:.6g}",
        )

    # Standard deviation envelope.
    if SHOW_STD and len(y) > 1:

        ax.axhspan(
            mean - std,
            mean + std,
            alpha=0.15,
            label=f"Mean ± 1 SD ({std:.6g})",
        )

    # ---------------------------------------------------------------------
    # Labels
    # ---------------------------------------------------------------------

    ax.set_title(
        f"{metric} Variation Across Frames",
        fontsize=14,
        fontweight="bold",
    )

    ax.set_xlabel(
        "Frame",
        fontsize=11,
    )

    ax.set_ylabel(
        "Difference from First Frame (%)",
        fontsize=11,
    )

    ax.axhline(
        0.0,
        linewidth=1.0,
        linestyle="--",
        alpha=0.6,
        label="First frame baseline",
    )

    ax.grid(
        alpha=0.25,
    )

    ax.legend(
        frameon=False,
    )

    # Avoid excessive margins around first/last frame.
    if len(x) > 1:
        x_span = max(x) - min(x)
        pad = max(0.5, 0.02 * x_span)

        ax.set_xlim(
            min(x) - pad,
            max(x) + pad,
        )

    fig.tight_layout()

    # ---------------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------------

    filename = (
        "frame_variation_"
        + safe_filename(metric)
    )

    if SAVE_PNG:

        png_path = output_dir / f"{filename}.png"

        fig.savefig(
            png_path,
            dpi=DPI,
            bbox_inches="tight",
        )

        print(f"Saved: {png_path}")

    if SAVE_SVG:

        svg_path = output_dir / f"{filename}.svg"

        fig.savefig(
            svg_path,
            bbox_inches="tight",
        )

        print(f"Saved: {svg_path}")

    if SHOW:
        plt.show()
    else:
        plt.close(fig)


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

def print_summary(summary: pd.DataFrame):
    """Print readable variation statistics."""

    print()
    print("=" * 80)
    print("FRAME VARIATION SUMMARY")
    print("=" * 80)

    if summary.empty:
        print("No statistics available.")
        return

    for _, row in summary.iterrows():

        print()
        print(row["Metric"])
        print("-" * 80)

        print(f"Frames:                 {int(row['N']):,}")
        print(f"Mean:                   {row['Mean']:.8g}")
        print(f"Standard deviation:     {row['Standard Deviation']:.8g}")
        print(f"Minimum:                {row['Minimum']:.8g}")
        print(f"Maximum:                {row['Maximum']:.8g}")
        print(f"Range:                  {row['Range']:.8g}")

        if np.isfinite(row["CV (%)"]):
            print(f"Coefficient variation:  {row['CV (%)']:.4f} %")
        else:
            print("Coefficient variation:  N/A")

        if np.isfinite(row["Range / Mean (%)"]):
            print(f"Range / mean:           {row['Range / Mean (%)']:.4f} %")
        else:
            print("Range / mean:           N/A")


# =============================================================================
# MAIN
# =============================================================================

def main():

    print()
    print("=" * 80)
    print("VORPY FRAME VARIATION ANALYSIS")
    print("=" * 80)

    # ---------------------------------------------------------------------
    # Select data directory
    # ---------------------------------------------------------------------

    selected = select_directory(
        "Select frames_new or a directory containing frames_new"
    )

    if not selected:
        print("No directory selected.")
        return

    data_root = Path(selected).resolve()

    print()
    print(f"Selected directory:")
    print(f"  {data_root}")

    # ---------------------------------------------------------------------
    # Find logs
    # ---------------------------------------------------------------------

    log_files = find_frame_logs(data_root)

    print()
    print(f"Found {len(log_files):,} frame log(s).")

    if not log_files:
        print()
        print("No logs matching:")
        print()
        print("    f#/aw/aw/aw_logs.csv")
        print()
        print("were found beneath the selected directory.")
        return

    print()
    print("First log:")
    print(f"  {log_files[0]}")

    print()
    print("Last log:")
    print(f"  {log_files[-1]}")

    # ---------------------------------------------------------------------
    # Read data
    # ---------------------------------------------------------------------

    frame_df = collect_frame_data(
        log_files=log_files,
        metrics=ENABLED_METRICS,
    )

    if frame_df.empty:
        print("No frame data could be read.")
        return

    # ---------------------------------------------------------------------
    # Add signed percentage difference from the first valid frame
    # ---------------------------------------------------------------------

    frame_df = add_percent_difference_columns(
        df=frame_df,
        metrics=ENABLED_METRICS,
    )

    # ---------------------------------------------------------------------
    # Output directory
    # ---------------------------------------------------------------------

    output_dir = data_root / "frame_variation_analysis"
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------------------------------------------------------------------
    # Save raw frame values
    # ---------------------------------------------------------------------

    data_csv = output_dir / "frame_variation_data.csv"

    frame_df.to_csv(
        data_csv,
        index=False,
    )

    print()
    print(f"Saved frame data:")
    print(f"  {data_csv}")

    # ---------------------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------------------

    summary = build_summary(
        df=frame_df,
        metrics=ENABLED_METRICS,
    )

    summary_csv = output_dir / "frame_variation_summary.csv"

    summary.to_csv(
        summary_csv,
        index=False,
    )

    print(f"Saved summary:")
    print(f"  {summary_csv}")

    print_summary(summary)

    # ---------------------------------------------------------------------
    # Plot metrics
    # ---------------------------------------------------------------------

    print()
    print("=" * 80)
    print("GENERATING PLOTS")
    print("=" * 80)
    print()

    for metric in ENABLED_METRICS:

        if metric not in frame_df.columns:
            print(f"Skipping missing metric: {metric}")
            continue

        plot_metric(
            df=frame_df,
            metric=metric,
            output_dir=output_dir,
        )

    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)

    print()
    print(f"Results:")
    print(f"  {output_dir}")
    print()


if __name__ == "__main__":
    main()