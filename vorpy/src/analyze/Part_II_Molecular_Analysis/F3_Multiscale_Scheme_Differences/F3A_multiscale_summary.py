"""Figure 3A — multiscale AW vs Power deviation summary.

Produces three compact panels (Volume, Surface Area, Contacts) comparing the
magnitude of Power-vs-AW deviations at atomic, residue, and molecular scales.
"""
from __future__ import annotations

import os
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from F3_common import (
    POWER_COLOR,
    add_deviation_columns,
    build_atomic_metrics,
    build_molecule_metrics,
    build_residue_metrics,
    discover_systems,
    read_pair,
    summarize_scale,
)


EXCLUDE_KEYS = ["A", "B", "C"]
FIGSIZE = (12.0, 4.2)
DPI = 300


def choose_directory(title: str) -> str:
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    path = filedialog.askdirectory(title=title)
    root.destroy()
    return path


def collect_multiscale_summary(data_root: str, exclude_keys=None):
    """Collect per-system summaries plus pooled atom/residue observations."""
    systems = discover_systems(data_root, exclude_keys=exclude_keys)
    if not systems:
        raise FileNotFoundError("No matched aw/aw_logs.csv + pow/pow_logs.csv systems were found.")

    atom_frames = []
    residue_frames = []
    molecule_rows = []

    for paths in systems:
        print(f"Processing {paths.name}...")
        aw_logs, power_logs = read_pair(paths, need_surfs=True)

        atoms = add_deviation_columns(build_atomic_metrics(aw_logs, power_logs))
        atoms["System"] = paths.name
        atom_frames.append(atoms)

        residues = add_deviation_columns(build_residue_metrics(atoms, aw_logs, power_logs))
        residues["System"] = paths.name
        residue_frames.append(residues)

        molecule = build_molecule_metrics(aw_logs, power_logs)
        molecule["System"] = paths.name
        molecule_rows.append(molecule)

    atom_df = pd.concat(atom_frames, ignore_index=True)
    residue_df = pd.concat(residue_frames, ignore_index=True)
    molecule_df = add_deviation_columns(pd.DataFrame(molecule_rows))

    summary = pd.concat(
        [
            summarize_scale(atom_df, "Atom"),
            summarize_scale(residue_df, "Residue"),
            summarize_scale(molecule_df, "Molecule"),
        ],
        ignore_index=True,
    )

    return summary, atom_df, residue_df, molecule_df


def plot_summary(summary: pd.DataFrame, save_paths: list[str] | None = None, show: bool = True):
    metrics = ["Volume", "Surface Area", "Contacts"]
    scales = ["Atom", "Residue", "Molecule"]

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE)

    for ax, metric in zip(axes, metrics):
        sub = summary[summary["Metric"] == metric].set_index("Scale")
        available = [scale for scale in scales if scale in sub.index and np.isfinite(sub.loc[scale, "Mean Abs % Diff"])]

        if not available:
            ax.text(0.5, 0.5, "Metric unavailable", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(metric)
            ax.set_xticks([])
            continue

        y = [sub.loc[scale, "Mean Abs % Diff"] for scale in available]
        err = [sub.loc[scale, "SEM"] for scale in available]
        x = np.arange(len(available))

        ax.bar(x, y, yerr=err, capsize=4, color=POWER_COLOR, alpha=0.82, linewidth=0)
        ax.set_xticks(x)
        ax.set_xticklabels(available)
        ax.set_title(metric, fontsize=13, fontweight="bold")
        ax.axhline(0, color="0.25", linewidth=1)
        ax.tick_params(axis="both", labelsize=10, width=1.2, length=5)

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    axes[0].set_ylabel("Mean absolute difference from AW (%)", fontsize=11)
    fig.tight_layout(w_pad=2.0)

    if save_paths:
        for save_path in save_paths:
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            fig.savefig(save_path, dpi=DPI, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def main():
    data_root = choose_directory("Select Figure 3 data folder")
    if not data_root:
        return

    figure_dir = choose_directory("Select figures/Figure_3 output folder")
    if not figure_dir:
        return

    summary, atom_df, residue_df, molecule_df = collect_multiscale_summary(
        data_root,
        exclude_keys=EXCLUDE_KEYS,
    )

    # Save lightweight analysis tables beside the figure during development.
    # Remove these three lines if you want the figure directory to contain only assets.
    summary.to_csv(os.path.join(figure_dir, "F3A_multiscale_summary.csv"), index=False)

    print("\nFigure 3A summary:")
    print(summary.to_string(index=False))

    plot_summary(
        summary,
        save_paths=[os.path.join(figure_dir, "F3A_multiscale_summary.png"),
                    os.path.join(figure_dir, "F3A_multiscale_summary.svg"),
                    os.path.join(figure_dir, "F3A_multiscale_summary.pdf")],
        show=True,
    )


if __name__ == "__main__":
    main()
