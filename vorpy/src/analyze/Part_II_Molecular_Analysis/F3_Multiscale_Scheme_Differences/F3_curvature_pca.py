#!/usr/bin/env python3
"""Figure 3 atomic geometry/curvature correlations and PCA.

Run this script directly from ``F3_old``. Unless paths are supplied on the
command line, it opens Tkinter dialogs for:

1. the Figure 3 molecular data root containing folders such as D_Hairpin and
   E_Cambrin; and
2. the output folder for the PCA tables and figures.

The script uses the same ``discover_systems`` and ``read_pair`` functions as
the other Figure 3 scripts. It matches AW and Power atom rows within each
system, pools those matched atoms across systems, calculates Pearson and
Spearman correlations, standardizes the variables, and runs PCA.

Recognized descriptors include radius; AW, Power, van der Waals, and named SA
volumes; AW/Power surface areas and contacts; local maximum/average mean and
Gaussian curvatures; and the integrated curvature terms. It also derives the
isolated-sphere references ``4*pi*r^2`` and ``4*pi*r^3/3`` from atomic radius.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# VorPy imports
# ---------------------------------------------------------------------------

vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
if vorpy_root not in sys.path:
    sys.path.append(vorpy_root)

from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common import (
    discover_systems,
    read_pair,
)


# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

DATA_ROOT = None
OUTPUT_ROOT = None
EXCLUDE_KEYS = ["A", "B", "C"]
RESULTS_FOLDER_NAME = "F3_atomic_curvature_pca"


ALIASES = {
    "atom_index": ("Index", "Atom Index", "atom_index", "index"),
    "atom_name": ("Name", "Atom Name", "atom_name"),
    "residue": ("Residue", "Residue Name", "residue"),
    "residue_sequence": ("Residue Sequence", "Residue Number", "residue_sequence"),
    "chain": ("Chain", "chain"),
    "radius": ("Radius", "radius"),
    "volume": ("Volume", "volume"),
    "vdw_volume": ("Van Der Waals Volume", "VDW Volume", "vdw_volume"),
    "sa_volume": (
        "SA Volume",
        "Solvent Accessible Volume",
        "solvent_accessible_volume",
        "sa_volume",
    ),
    "surface_area": ("Surface Area", "surface_area", "SA"),
    "contacts": ("Number of Neighbors", "Contacts", "contacts"),
    "maximum_mean_curvature": (
        "Maximum Mean Curvature",
        "maximum_mean_curvature",
    ),
    "mean_curvature": (
        "Average Mean Surface Curvature",
        "Area-Weighted Mean Curvature",
        "Mean Curvature",
        "mean_curvature",
    ),
    "integrated_mean_curvature": (
        "Integrated Mean Curvature",
        "int_mean_curv",
        "integrated_mean_curvature",
    ),
    "integrated_mean_curvature_squared": (
        "Integrated Mean Curvature Squared",
        "int_mean_curv_sq",
        "integrated_mean_curvature_squared",
    ),
    "gaussian_curvature": (
        "Average Gaussian Surface Curvature",
        "Area-Weighted Gaussian Curvature",
        "Gaussian Curvature",
        "gaussian_curvature",
    ),
    "maximum_gaussian_curvature": (
        "Maximum Gaussian Curvature",
        "maximum_gaussian_curvature",
    ),
    "integrated_gaussian_curvature": (
        "Integrated Gaussian Curvature",
        "int_gauss_curv",
        "integrated_gaussian_curvature",
    ),
    "gaussian_face": ("Integrated Gaussian Curvature Face", "int_gauss_curv_face"),
    "gaussian_edge": ("Integrated Gaussian Curvature Edge", "int_gauss_curv_edge"),
    "gaussian_vertex": ("Integrated Gaussian Curvature Vertex", "int_gauss_curv_vertex"),
    "mean_face": ("Integrated Mean Curvature Face", "int_mean_curv_face"),
    "mean_edge": ("Integrated Mean Curvature Edge", "int_mean_curv_edge"),
}


def select_directory(title: str) -> str:
    """Open a topmost Tkinter directory picker and return the selected path."""
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    path = filedialog.askdirectory(title=title)
    root.destroy()
    return path


def canonical(text: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).strip().lower())


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        frame = pd.read_excel(path)
    else:
        frame = pd.read_csv(path, sep=None, engine="python")
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def find_column(frame: pd.DataFrame, names: tuple[str, ...] | list[str]) -> str | None:
    lookup = {canonical(column): column for column in frame.columns}
    for name in names:
        if canonical(name) in lookup:
            return lookup[canonical(name)]
    return None


def numeric_column(frame: pd.DataFrame, alias: str) -> pd.Series | None:
    column = find_column(frame, ALIASES[alias])
    if column is None:
        return None
    return pd.to_numeric(frame[column], errors="coerce")


def choose_key(frame: pd.DataFrame, requested: str | None) -> str:
    if requested:
        column = find_column(frame, [requested])
        if column is None:
            raise ValueError(f"Join key {requested!r} was not found. Available columns: {list(frame.columns)}")
        return column
    column = find_column(frame, ALIASES["atom_index"])
    if column is None:
        raise ValueError("Could not find an atom index column. Pass --key with the shared join column.")
    return column


def normalized_atom_ids(values: pd.Series) -> pd.Series:
    """Normalize integer-like atom indexes so AW and Power keys match."""
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().all() and np.allclose(numeric, np.round(numeric)):
        return numeric.round().astype("int64").astype(str)
    return values.astype(str).str.strip()


def prepare_source(frame: pd.DataFrame, scheme: str, requested_key: str | None) -> pd.DataFrame:
    key = choose_key(frame, requested_key)
    output = pd.DataFrame({"atom_id": normalized_atom_ids(frame[key])})

    if scheme == "aw":
        for metadata in ("atom_name", "residue", "residue_sequence", "chain"):
            source = find_column(frame, ALIASES[metadata])
            if source is not None:
                output[metadata] = frame[source]

    aw_base = (
        "radius",
        "volume",
        "vdw_volume",
        "sa_volume",
        "surface_area",
        "contacts",
    )
    comparison_base = ("volume", "surface_area", "contacts")
    curvature = (
        "maximum_mean_curvature",
        "mean_curvature",
        "integrated_mean_curvature",
        "integrated_mean_curvature_squared",
        "maximum_gaussian_curvature",
        "gaussian_curvature",
        "integrated_gaussian_curvature",
        "gaussian_face",
        "gaussian_edge",
        "gaussian_vertex",
        "mean_face",
        "mean_edge",
    )
    requested = aw_base + curvature if scheme == "aw" else comparison_base
    for alias in requested:
        values = numeric_column(frame, alias)
        if values is not None:
            output[f"{scheme}_{alias}"] = values

    if scheme == "aw" and "aw_radius" in output:
        # Isolated-sphere reference values. These are intentionally named
        # "sphere", not "solvent accessible", because no probe radius is used.
        radius = output["aw_radius"]
        output["sphere_volume"] = (4.0 / 3.0) * np.pi * radius**3
        output["sphere_surface_area"] = 4.0 * np.pi * radius**2

    if output["atom_id"].duplicated().any():
        duplicated = output.loc[output["atom_id"].duplicated(), "atom_id"].head().tolist()
        raise ValueError(f"Duplicate {scheme} atom IDs prevent a one-to-one join: {duplicated}")
    return output


def collect_figure3_atoms(
    data_root: Path,
    requested_key: str | None,
    exclude_keys: list[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Read and pool matched AW/Power atoms from valid Figure 3 systems."""
    systems = discover_systems(str(data_root), exclude_keys=exclude_keys)
    print(f"\nFound {len(systems)} valid Figure 3 system folders.")

    frames: list[pd.DataFrame] = []
    system_reports: list[dict[str, object]] = []

    for system in systems:
        print(f"\nProcessing {system.name} ...")
        try:
            aw_logs, power_logs = read_pair(system, need_surfs=False)
            aw_atoms = aw_logs["atoms"]
            power_atoms = power_logs["atoms"]

            aw = prepare_source(aw_atoms, "aw", requested_key)
            power = prepare_source(power_atoms, "pow", requested_key)
            matched = aw.merge(power, on="atom_id", how="inner", validate="one_to_one")

            if matched.empty:
                raise ValueError("no matching AW/Power atom indexes")

            matched.insert(0, "system", system.name)
            frames.append(matched)
            system_reports.append(
                {
                    "system": system.name,
                    "aw_atoms": int(len(aw)),
                    "power_atoms": int(len(power)),
                    "matched_atoms": int(len(matched)),
                    "aw_columns": list(aw_atoms.columns),
                    "power_columns": list(power_atoms.columns),
                }
            )
            print(
                f"  AW atoms: {len(aw):,}; Power atoms: {len(power):,}; "
                f"matched: {len(matched):,}"
            )
        except Exception as error:
            system_reports.append(
                {"system": system.name, "status": "skipped", "reason": str(error)}
            )
            print(f"  skipped: {error}")

    if not frames:
        raise ValueError("No systems produced matched AW/Power atom data.")

    pooled = pd.concat(frames, ignore_index=True)
    pooled.insert(1, "atom_uid", pooled["system"] + ":" + pooled["atom_id"])
    return pooled, {
        "data_root": str(data_root),
        "systems_discovered": int(len(systems)),
        "systems_used": int(len(frames)),
        "systems": system_reports,
    }


def plot_heatmap(correlation: pd.DataFrame, path: Path, title: str) -> None:
    size = max(8.0, 0.72 * len(correlation.columns))
    fig, ax = plt.subplots(figsize=(size, size))
    sns.heatmap(
        correlation,
        vmin=-1,
        vmax=1,
        center=0,
        cmap="vlag",
        square=True,
        annot=len(correlation.columns) <= 12,
        fmt=".2f",
        ax=ax,
    )
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_scree(explained: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    numbers = np.arange(1, len(explained) + 1)
    ax.bar(numbers, explained["explained_variance_ratio"], color="#2b6f8a")
    ax.plot(numbers, explained["cumulative_variance_ratio"], marker="o", color="#b5473c")
    ax.set(xlabel="Principal component", ylabel="Fraction of variance", title="PCA explained variance")
    ax.set_xticks(numbers)
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_scores(scores: pd.DataFrame, explained: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    if "system" in scores and scores["system"].nunique() <= 15:
        hue = "system"
    elif "atom_name" in scores and scores["atom_name"].nunique() <= 15:
        hue = "atom_name"
    else:
        hue = None
    sns.scatterplot(data=scores, x="PC1", y="PC2", hue=hue, s=55, alpha=0.8, ax=ax)
    x_pct = 100 * explained.loc[0, "explained_variance_ratio"]
    y_pct = 100 * explained.loc[1, "explained_variance_ratio"]
    ax.set(xlabel=f"PC1 ({x_pct:.1f}%)", ylabel=f"PC2 ({y_pct:.1f}%)", title="Atomic geometry PCA")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_loadings(loadings: pd.DataFrame, path: Path) -> None:
    shown = loadings.iloc[:, : min(4, loadings.shape[1])]
    fig, ax = plt.subplots(figsize=(7.5, max(5.0, 0.45 * len(shown))))
    sns.heatmap(shown, center=0, cmap="vlag", annot=True, fmt=".2f", ax=ax)
    ax.set_title("PCA loadings")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DATA_ROOT,
        help="Figure 3 molecular data root; opens a Tkinter picker when omitted",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="Parent output folder; opens a Tkinter picker when omitted",
    )
    parser.add_argument("--key", help="Shared atom join column; defaults to Index/Atom Index")
    parser.add_argument("--components", type=int, default=6)
    parser.add_argument(
        "--exclude-key",
        action="append",
        dest="exclude_keys",
        help=(
            "Leading Figure 3 system letter to exclude. Repeat for multiple keys. "
            "Defaults to A, B, and C."
        ),
    )
    parser.add_argument(
        "--missing",
        choices=("median", "complete"),
        default="median",
        help="Median-impute missing values or retain complete rows only",
    )
    parser.add_argument(
        "--variance-tolerance",
        type=float,
        default=1e-10,
        help="Drop PCA features with variance at or below this threshold",
    )
    parser.add_argument(
        "--include-topological-gaussian",
        action="store_true",
        help=(
            "Force integrated Gaussian curvature into PCA even when every atom is "
            "within the Gauss-Bonnet tolerance of 4*pi. It is always retained in "
            "the merged data and correlation tables."
        ),
    )
    parser.add_argument(
        "--gauss-bonnet-tolerance",
        type=float,
        default=0.05,
        help=(
            "Absolute tolerance for treating atomic integrated Gaussian curvature "
            "as the topological constant 4*pi (default: 0.05 in dimensionless "
            "angle measure)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    data_root = args.data_root
    if data_root is None:
        selected = select_directory(
            "STEP 1/2 — Select MOLECULAR DATA ROOT containing D_Hairpin, E_Cambrin, etc."
        )
        if not selected:
            print("No molecular data root selected.")
            return
        data_root = Path(selected)

    output_root = args.output_root
    if output_root is None:
        selected = select_directory(
            "STEP 2/2 — Select OUTPUT folder for the atomic curvature PCA"
        )
        if not selected:
            print("No output folder selected.")
            return
        output_root = Path(selected)

    output = output_root / RESULTS_FOLDER_NAME
    output.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    exclude_keys = args.exclude_keys if args.exclude_keys is not None else EXCLUDE_KEYS
    print()
    print(f"Molecular data root: {data_root}")
    print(f"PCA output folder:   {output}")
    print(f"Excluded system keys: {', '.join(exclude_keys) if exclude_keys else '(none)'}")

    merged, collection_report = collect_figure3_atoms(
        data_root=data_root,
        requested_key=args.key,
        exclude_keys=exclude_keys,
    )

    feature_columns = [
        column
        for column in merged.columns
        if column.startswith(("aw_", "pow_", "sa_", "sphere_"))
    ]
    if len(feature_columns) < 2:
        raise ValueError(f"PCA needs at least two recognized numeric features; found {feature_columns}")

    merged.to_csv(output / "atomic_descriptors_merged.csv", index=False)
    raw = merged[feature_columns].replace([np.inf, -np.inf], np.nan)
    missing_counts = raw.isna().sum()
    all_missing = raw.columns[raw.isna().all()].tolist()
    raw = raw.drop(columns=all_missing)

    if args.missing == "complete":
        kept_rows = raw.notna().all(axis=1)
        analysis = raw.loc[kept_rows].copy()
        analysis_metadata = merged.loc[kept_rows].reset_index(drop=True)
    else:
        imputer = SimpleImputer(strategy="median")
        analysis = pd.DataFrame(imputer.fit_transform(raw), columns=raw.columns, index=raw.index)
        analysis_metadata = merged.reset_index(drop=True)

    # Keep the full usable descriptor set for correlation analysis. PCA applies
    # stricter filtering below because constant/topological variables cannot
    # carry meaningful standardized variance.
    correlation_analysis = analysis.copy()

    topological_columns = []
    gaussian_column = "aw_integrated_gaussian_curvature"
    if gaussian_column in analysis.columns and not args.include_topological_gaussian:
        gaussian_values = analysis[gaussian_column].to_numpy(dtype=float)
        if (
            gaussian_values.size
            and np.all(np.isfinite(gaussian_values))
            and np.max(np.abs(gaussian_values - 4.0 * np.pi))
            <= args.gauss_bonnet_tolerance
        ):
            topological_columns.append(gaussian_column)
            analysis = analysis.drop(columns=gaussian_column)

    variances = analysis.var(axis=0, ddof=0)
    constant_columns = variances[variances <= args.variance_tolerance].index.tolist()
    analysis = analysis.drop(columns=constant_columns)
    if analysis.shape[1] < 2:
        raise ValueError(f"Fewer than two varying features remain after filtering: {list(analysis.columns)}")
    if len(analysis) < 3:
        raise ValueError("PCA needs at least three usable atom rows.")

    pearson = correlation_analysis.corr(method="pearson")
    spearman = correlation_analysis.corr(method="spearman")
    pearson.to_csv(output / "correlation_pearson.csv")
    spearman.to_csv(output / "correlation_spearman.csv")
    plot_heatmap(pearson, output / "correlation_pearson.png", "Atomic descriptor Pearson correlations")
    plot_heatmap(spearman, output / "correlation_spearman.png", "Atomic descriptor Spearman correlations")

    scaler = StandardScaler()
    scaled_values = scaler.fit_transform(analysis)
    scaled = pd.DataFrame(scaled_values, columns=analysis.columns)
    scaled.to_csv(output / "atomic_descriptors_standardized.csv", index=False)

    component_count = min(args.components, analysis.shape[0], analysis.shape[1])
    pca = PCA(n_components=component_count)
    score_values = pca.fit_transform(scaled_values)
    pc_names = [f"PC{number}" for number in range(1, component_count + 1)]
    scores = pd.DataFrame(score_values, columns=pc_names)
    for metadata in (
        "system",
        "atom_uid",
        "atom_id",
        "atom_name",
        "residue",
        "residue_sequence",
        "chain",
    ):
        if metadata in analysis_metadata:
            scores.insert(len(scores.columns), metadata, analysis_metadata[metadata].to_numpy())
    scores.to_csv(output / "pca_scores.csv", index=False)

    loadings = pd.DataFrame(pca.components_.T, index=analysis.columns, columns=pc_names)
    loadings.to_csv(output / "pca_loadings.csv")
    explained = pd.DataFrame(
        {
            "component": pc_names,
            "explained_variance": pca.explained_variance_,
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
        }
    )
    explained.to_csv(output / "pca_explained_variance.csv", index=False)
    plot_scree(explained, output / "pca_scree.png")
    if component_count >= 2:
        plot_scores(scores, explained, output / "pca_scores_pc1_pc2.png")
    plot_loadings(loadings, output / "pca_loadings.png")

    gaussian_column = "aw_integrated_gaussian_curvature"
    gauss_bonnet = None
    if gaussian_column in merged:
        values = pd.to_numeric(merged[gaussian_column], errors="coerce")
        deviations = values - 4 * np.pi
        gauss_bonnet = {
            "count": int(values.notna().sum()),
            "mean": float(values.mean()),
            "standard_deviation": float(values.std(ddof=0)),
            "mean_absolute_error_from_4pi": float(deviations.abs().mean()),
            "maximum_absolute_error_from_4pi": float(deviations.abs().max()),
        }

    report = {
        "collection": collection_report,
        "output_folder": str(output),
        "excluded_system_keys": exclude_keys,
        "merged_atom_count": int(len(merged)),
        "pca_atom_count": int(len(analysis)),
        "recognized_features": feature_columns,
        "pca_features": list(analysis.columns),
        "constant_or_near_constant_features_excluded": constant_columns,
        "topologically_constrained_features_excluded": topological_columns,
        "all_missing_features_excluded": all_missing,
        "missing_values_by_feature": {key: int(value) for key, value in missing_counts.items()},
        "gauss_bonnet_atom_check": gauss_bonnet,
    }
    (output / "analysis_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"Analyzed {len(analysis)} atoms with {analysis.shape[1]} varying features.")
    if constant_columns:
        print("Excluded constant/near-constant PCA features:", ", ".join(constant_columns))
    if topological_columns:
        print(
            "Excluded topologically constrained PCA features:",
            ", ".join(topological_columns),
            "(retained in raw and correlation outputs)",
        )
    print(f"Results written to: {output.resolve()}")

    try:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes("-topmost", 1)
        messagebox.showinfo(
            "F3 atomic curvature PCA complete",
            f"Analyzed {len(analysis):,} matched atoms.\n\nResults:\n{output.resolve()}",
            parent=root,
        )
        root.destroy()
    except tk.TclError:
        # Command-line/headless runs still complete normally.
        pass


if __name__ == "__main__":
    main()
