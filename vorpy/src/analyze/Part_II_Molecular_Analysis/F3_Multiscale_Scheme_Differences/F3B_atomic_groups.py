"""Figure 3B — Atomic AW vs Power comparisons.

Creates a 2 x 3 Figure 3B panel matrix:

    Protein       : Volume | Surface Area | Contacts
    Nucleic acids : Volume | Surface Area | Contacts

Only systems discovered by F3_common.discover_systems() are considered, so the
data root is restricted to folders named <letter>_<molecule name> and the
leading letter is still filtered by EXCLUDE_KEYS.

The script classifies each system from its residue content, then pools matched
atoms across systems. Atom groups use the same canonicalization rules already
developed for F2C_All_Data_Sets_Clustered.py.

Outputs
-------
F3B_atomic_groups.png
F3B_atomic_groups.svg
F3B_protein_group_summary.csv
F3B_nucleic_group_summary.csv
"""
from __future__ import annotations

import os
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse


# ---------------------------------------------------------------------------
# VorPy imports
# ---------------------------------------------------------------------------

vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
if vorpy_root not in sys.path:
    sys.path.append(vorpy_root)

from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common import (
    build_atomic_metrics,
    discover_systems,
    read_pair,
)

from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_All_Data_Sets_Clustered import (
    canonicalize_atom_name,
)


# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

DATA_ROOT = None
FIGURE_DIR = None

EXCLUDE_KEYS = ["A", "B", "C"]

MIN_GROUP_COUNT = 25

# Figure 3B should show a small, representative set of chemically meaningful
# groups rather than simply the mathematically tightest clusters.
REPRESENTATIVE_GROUPS = 6

# Raw scatter only. Full data are still used for group statistics/selection.
PROTEIN_PLOT_FRACTION = 0.20
PLOT_RANDOM_SEED = 17

# Candidate filtering: discard only the most diffuse 30% before scoring.
COMPACTNESS_QUANTILE = 0.70

# Automatic representative-selection weights.
WEIGHT_TIGHTNESS = 0.35
WEIGHT_SEPARATION = 0.30
WEIGHT_ATOM_DIVERSITY = 0.20
WEIGHT_COUNT = 0.15

# Optional manual overrides after inspecting the printed grouping tables.
# Leave a list empty to use automatic selection for that panel.
# Final Figure 3B representative groups.
#
# Use the SAME atom groups for Volume, Surface Area, and Contacts so the
# three columns can be compared directly.
#
# Protein:
#   CA_* = generic alpha carbon (CA)
#   CB_* = generic beta carbon (CB)
#   O    = backbone oxygen
#   NZ   = lysine terminal nitrogen
#   HN   = backbone/amide hydrogen
#   SG   = cysteine sulfur
#
# Nucleic acids:
#   C5'
#   C_sugar
#   N9
#   H_base_exocyclic
#   O_backbone
#   P
PROTEIN_REPRESENTATIVE_GROUPS = [
    "CA_*",
    "CB_*",
    "O",
    "NZ",
    "HN",
    "SG",
]

NUCLEIC_REPRESENTATIVE_GROUPS = [
    "C5'",
    "C_sugar",
    "N9",
    "H_base_exocyclic",
    "O_backbone",
    "P",
]

HARD_CODED_GROUPS = {
    ("Protein", "Volume"): PROTEIN_REPRESENTATIVE_GROUPS,
    ("Protein", "Surface Area"): PROTEIN_REPRESENTATIVE_GROUPS,
    ("Protein", "Contacts"): PROTEIN_REPRESENTATIVE_GROUPS,

    ("Nucleic acid", "Volume"): NUCLEIC_REPRESENTATIVE_GROUPS,
    ("Nucleic acid", "Surface Area"): NUCLEIC_REPRESENTATIVE_GROUPS,
    ("Nucleic acid", "Contacts"): NUCLEIC_REPRESENTATIVE_GROUPS,
}

# Require one representative from every major element present in each class.
# Protein residues contain H/C/N/O/S; nucleic acids contain H/C/N/O/P.
REQUIRED_ELEMENTS = {
    "Protein": ["H", "C", "N", "O", "S"],
    "Nucleic acid": ["H", "C", "N", "O", "P"],
}

LABEL_SELECTED_GROUPS = True
ELLIPSE_N_STD = 1.25

# Minimum shared upper limits. These prevent selected representative groups
# from being clipped by percentile-based axis trimming.
MIN_SHARED_UPPER_LIMIT = {
    "Volume": 32.0,
    "Surface Area": 60.0,
    "Contacts": None,
}

POINT_SIZE = 8
POINT_ALPHA = 0.10
GROUP_POINT_SIZE = 55

FIGSIZE = (18, 11)
DPI = 300

SAVE_PNG = True
SAVE_SVG = True
SHOW = True

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DI"}
RNA_RESIDUES = {"A", "C", "G", "U", "I", "RA", "RC", "RG", "RU"}

ELEMENT_COLORS = {
    "H": "#1f77b4",
    "C": "#ff7f0e",
    "N": "#2ca02c",
    "O": "#d62728",
    "P": "#9467bd",
    "S": "#8c564b",
    "SE": "#8c564b",
}


# ---------------------------------------------------------------------------
# Folder selection
# ---------------------------------------------------------------------------

def select_directory(title: str) -> str:
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    path = filedialog.askdirectory(title=title)
    root.destroy()
    return path


# ---------------------------------------------------------------------------
# System classification
# ---------------------------------------------------------------------------

def classify_system(atoms: pd.DataFrame) -> Optional[str]:
    """Classify a molecular system as protein, dna, or rna from residues."""
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


def broad_class(system_class: str) -> str:
    if system_class == "protein":
        return "Protein"
    if system_class in {"dna", "rna"}:
        return "Nucleic acid"
    raise ValueError(system_class)


# ---------------------------------------------------------------------------
# Atom grouping
# ---------------------------------------------------------------------------

def element_from_group(group_name: str) -> str:
    name = str(group_name).strip().upper()
    base = name.split("_")[0].replace("*", "")

    if base.startswith("SE"):
        return "SE"

    for element in ("H", "C", "N", "O", "P", "S"):
        if base.startswith(element):
            return element

    return "Other"


def canonicalize_atomic_dataframe(
    atom_df: pd.DataFrame,
    system_class: str,
    system_name: str,
) -> pd.DataFrame:
    df = atom_df.copy()

    df["System"] = system_name
    df["SystemClass"] = system_class
    df["BroadClass"] = broad_class(system_class)

    df["CanonicalName"] = [
        canonicalize_atom_name(
            atom_name=str(atom),
            molecule_class=system_class,
            residue_name=str(residue),
        )
        for atom, residue in zip(df["Atom"], df["Residue"])
    ]

    df["ElementClass"] = df["CanonicalName"].apply(element_from_group)
    return df


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_atomic_data(
    data_root: str,
    exclude_keys: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    systems = discover_systems(data_root, exclude_keys=exclude_keys)

    print(f"\nFound {len(systems)} valid Figure 3 system folders.")

    frames: List[pd.DataFrame] = []

    for system in systems:
        print(f"\nProcessing {system.name} ...")

        aw_logs, power_logs = read_pair(system, need_surfs=True)
        system_class = classify_system(aw_logs["atoms"])

        if system_class is None:
            print("  skipped: no dominant protein/DNA/RNA residue class")
            continue

        print(f"  classified as: {system_class}")

        atom_df = build_atomic_metrics(aw_logs, power_logs)
        if atom_df.empty:
            print("  skipped: no matched AW/Power atoms")
            continue

        atom_df = canonicalize_atomic_dataframe(
            atom_df=atom_df,
            system_class=system_class,
            system_name=system.name,
        )

        # Remove solvent, ions, ligands, etc. from the class-specific analysis.
        if system_class == "protein":
            keep = atom_df["Residue"].isin(PROTEIN_RESIDUES)
        elif system_class == "dna":
            keep = atom_df["Residue"].isin(DNA_RESIDUES)
        else:
            keep = atom_df["Residue"].isin(RNA_RESIDUES)

        atom_df = atom_df.loc[keep].copy()

        print(f"  retained atoms: {len(atom_df):,}")
        frames.append(atom_df)

    if not frames:
        return pd.DataFrame()

    result = pd.concat(frames, ignore_index=True)

    print("\nCollected atomic data:")
    for name, subdf in result.groupby("BroadClass"):
        print(f"  {name}: {len(subdf):,} atoms")

    return result


# ---------------------------------------------------------------------------
# Group statistics
# ---------------------------------------------------------------------------

METRIC_COLUMNS = {
    "Volume": ("AW Volume", "Power Volume"),
    "Surface Area": ("AW Surface Area", "Power Surface Area"),
    "Contacts": ("AW Contacts", "Power Contacts"),
}


def compute_group_stats(
    df: pd.DataFrame,
    metric: str,
    min_count: int = MIN_GROUP_COUNT,
) -> pd.DataFrame:
    aw_col, power_col = METRIC_COLUMNS[metric]

    rows = []

    for group_name, group_df in df.groupby("CanonicalName"):
        values = (
            group_df[[aw_col, power_col]]
            .apply(pd.to_numeric, errors="coerce")
            .dropna()
        )

        if len(values) < min_count:
            continue

        aw = values[aw_col].to_numpy(float)
        power = values[power_col].to_numpy(float)

        mean_aw = float(np.mean(aw))
        mean_power = float(np.mean(power))
        sd_aw = float(np.std(aw, ddof=1)) if len(aw) > 1 else 0.0
        sd_power = float(np.std(power, ddof=1)) if len(power) > 1 else 0.0

        spread_2d = float(np.sqrt(sd_aw ** 2 + sd_power ** 2))
        center_mag = float(np.sqrt(mean_aw ** 2 + mean_power ** 2))
        normalized_spread = (
            spread_2d / center_mag
            if center_mag > 0.0
            else np.inf
        )

        delta = mean_power - mean_aw
        pct = np.nan if mean_aw == 0 else 100.0 * delta / mean_aw

        rows.append(
            {
                "Metric": metric,
                "GroupName": group_name,
                "AtomType": str(group_name).split("_")[0],
                "Element": element_from_group(group_name),
                "Count": len(values),
                "Mean_AW": mean_aw,
                "Mean_Power": mean_power,
                "SD_AW": sd_aw,
                "SD_Power": sd_power,
                "Spread_2D": spread_2d,
                "Normalized_Spread": normalized_spread,
                "Mean_Delta": delta,
                "Mean_Signed_%_Diff": pct,
                "Abs_Mean_Signed_%_Diff": abs(pct) if np.isfinite(pct) else np.nan,
            }
        )

    stats = pd.DataFrame(rows)
    if stats.empty:
        return stats

    return stats.sort_values(
        ["Abs_Mean_Signed_%_Diff", "Count"],
        ascending=[False, False],
    ).reset_index(drop=True)


def _normalize_series(values: pd.Series, invert: bool = False) -> pd.Series:
    """Scale finite values to [0, 1]."""
    vals = pd.to_numeric(values, errors="coerce").astype(float)
    finite = np.isfinite(vals)
    out = pd.Series(np.zeros(len(vals), dtype=float), index=vals.index)

    if not finite.any():
        return out

    finite_vals = vals[finite]
    lo = float(finite_vals.min())
    hi = float(finite_vals.max())

    if hi <= lo:
        out.loc[finite] = 1.0
    else:
        out.loc[finite] = (finite_vals - lo) / (hi - lo)

    if invert:
        out.loc[finite] = 1.0 - out.loc[finite]

    return out


def _normalized_group_coordinates(candidates: pd.DataFrame) -> np.ndarray:
    """Normalize AW/Power group centers independently to [0, 1]."""
    coords = candidates[["Mean_AW", "Mean_Power"]].astype(float).to_numpy()
    mins = np.nanmin(coords, axis=0)
    maxs = np.nanmax(coords, axis=0)
    spans = maxs - mins
    spans[spans == 0.0] = 1.0
    return (coords - mins) / spans


def _manual_group_selection(
    stats_df: pd.DataFrame,
    broad_class_name: str,
    metric: str,
) -> pd.DataFrame:
    requested = HARD_CODED_GROUPS.get((broad_class_name, metric), [])

    if not requested:
        return pd.DataFrame()

    lookup = stats_df.set_index("GroupName", drop=False)
    rows = []

    for name in requested:
        if name not in lookup.index:
            print(
                f"WARNING: manual F3B group '{name}' not found for "
                f"{broad_class_name} / {metric}"
            )
            continue

        row = lookup.loc[name]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        rows.append(row.copy())

    if not rows:
        return pd.DataFrame()

    selected = pd.DataFrame(rows).reset_index(drop=True)
    selected["RepresentativeRank"] = np.arange(1, len(selected) + 1)
    selected["SelectionMode"] = "manual"
    return selected


def _base_group_score(candidates: pd.DataFrame) -> pd.Series:
    """
    Score individual groups before spatial-diversity selection.

    Higher is better:
      - compact groups,
      - large groups,
      - groups with a clear AW/Power shift.
    """
    tightness = _normalize_series(
        candidates["Normalized_Spread"],
        invert=True,
    )
    count = _normalize_series(
        np.log1p(candidates["Count"].astype(float))
    )
    deviation = _normalize_series(
        np.abs(candidates["Mean_Signed_%_Diff"].astype(float))
    )

    # Deviation is useful, but compactness and count remain primary because
    # these examples need to be visually coherent and representative.
    return (
        0.50 * tightness
        + 0.30 * count
        + 0.20 * deviation
    )


def select_representative_groups(
    stats_df: pd.DataFrame,
    broad_class_name: str,
    metric: str,
    n_groups: int = REPRESENTATIVE_GROUPS,
    compactness_quantile: float = COMPACTNESS_QUANTILE,
) -> pd.DataFrame:
    """
    Select Figure 3B groups with explicit elemental representation.

    Stage 1
    -------
    Pick the best available H/C/N/O/S group for proteins or H/C/N/O/P group
    for nucleic acids.

    Stage 2
    -------
    Fill any remaining slots with groups that are compact, abundant, and
    spatially separated from those already selected.

    HARD_CODED_GROUPS still overrides the automatic choice.
    """
    manual = _manual_group_selection(
        stats_df,
        broad_class_name=broad_class_name,
        metric=metric,
    )
    if not manual.empty:
        return manual

    if stats_df.empty:
        return stats_df.copy()

    candidates = stats_df[
        np.isfinite(stats_df["Normalized_Spread"]) &
        np.isfinite(stats_df["Mean_AW"]) &
        np.isfinite(stats_df["Mean_Power"])
    ].copy()

    if candidates.empty:
        return candidates

    # Keep this fairly permissive so O/H groups are not automatically removed
    # simply because they are broader than tightly packed carbon classes.
    spread_cutoff = float(
        candidates["Normalized_Spread"].quantile(compactness_quantile)
    )

    compact_candidates = candidates[
        candidates["Normalized_Spread"] <= spread_cutoff
    ].copy()

    # Element representation takes precedence over the compactness cutoff.
    # If an element has no group inside the cutoff, use the best group for that
    # element from the complete eligible set instead.
    candidates["BaseScore"] = _base_group_score(candidates)
    compact_candidates["BaseScore"] = _base_group_score(compact_candidates)

    selected_rows = []
    selected_names = set()

    required = REQUIRED_ELEMENTS.get(broad_class_name, [])

    # ------------------------------------------------------------------
    # Stage 1: one best group for every required element.
    # ------------------------------------------------------------------
    for element in required:
        pool = compact_candidates[
            compact_candidates["Element"] == element
        ].copy()

        if pool.empty:
            pool = candidates[
                candidates["Element"] == element
            ].copy()

        if pool.empty:
            print(
                f"WARNING: no {element} grouping found for "
                f"{broad_class_name} / {metric}"
            )
            continue

        chosen = pool.sort_values(
            ["BaseScore", "Normalized_Spread", "Count"],
            ascending=[False, True, False],
        ).iloc[0].copy()

        selected_rows.append(chosen)
        selected_names.add(chosen["GroupName"])

    # ------------------------------------------------------------------
    # Stage 2: fill extra slots using spatial separation + base quality.
    # ------------------------------------------------------------------
    remaining_slots = max(0, int(n_groups) - len(selected_rows))

    if remaining_slots > 0:
        remaining = compact_candidates[
            ~compact_candidates["GroupName"].isin(selected_names)
        ].copy()

        if not remaining.empty:
            # Normalize all available centers on the combined candidate range.
            coord_source = candidates.reset_index(drop=True)
            all_coords = _normalized_group_coordinates(coord_source)
            coord_map = {
                coord_source.loc[i, "GroupName"]: all_coords[i]
                for i in coord_source.index
            }

            for _ in range(remaining_slots):
                if remaining.empty:
                    break

                best_idx = None
                best_score = -np.inf

                selected_coords = [
                    coord_map[row["GroupName"]]
                    for row in selected_rows
                    if row["GroupName"] in coord_map
                ]

                for idx, row in remaining.iterrows():
                    coord = coord_map.get(row["GroupName"])

                    if coord is None or not selected_coords:
                        separation = 1.0
                    else:
                        separation = min(
                            float(np.linalg.norm(coord - other))
                            for other in selected_coords
                        ) / np.sqrt(2.0)
                        separation = min(max(separation, 0.0), 1.0)

                    score = (
                        0.70 * float(row["BaseScore"])
                        + 0.30 * separation
                    )

                    if score > best_score:
                        best_score = score
                        best_idx = idx

                if best_idx is None:
                    break

                chosen = remaining.loc[best_idx].copy()
                selected_rows.append(chosen)
                selected_names.add(chosen["GroupName"])
                remaining = remaining.drop(index=best_idx)

    if not selected_rows:
        return pd.DataFrame(columns=stats_df.columns)

    selected = pd.DataFrame(selected_rows).reset_index(drop=True)
    selected["RepresentativeRank"] = np.arange(1, len(selected) + 1)
    selected["SelectionMode"] = "automatic_element_balanced"

    return selected


def print_grouping_table(
    stats_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    broad_class_name: str,
    metric: str,
):
    """Print all groupings so manual representative choices are easy."""
    print()
    print("=" * 120)
    print(f"F3B GROUPINGS — {broad_class_name.upper()} | {metric.upper()}")
    print("=" * 120)

    if stats_df.empty:
        print("No eligible groups.")
        return

    selected_names = (
        set(selected_df["GroupName"])
        if not selected_df.empty
        else set()
    )

    table = stats_df.copy()
    table["Selected"] = table["GroupName"].isin(selected_names)

    display_cols = [
        "Selected",
        "GroupName",
        "AtomType",
        "Element",
        "Count",
        "Mean_AW",
        "Mean_Power",
        "Mean_Signed_%_Diff",
        "SD_AW",
        "SD_Power",
        "Normalized_Spread",
    ]

    table = table[display_cols].sort_values(
        ["Selected", "Mean_Signed_%_Diff", "Normalized_Spread"],
        ascending=[False, True, True],
    )

    print(
        table.to_string(
            index=False,
            formatters={
                "Mean_AW": lambda x: f"{x:8.3f}",
                "Mean_Power": lambda x: f"{x:8.3f}",
                "Mean_Signed_%_Diff": lambda x: f"{x:8.3f}",
                "SD_AW": lambda x: f"{x:7.3f}",
                "SD_Power": lambda x: f"{x:7.3f}",
                "Normalized_Spread": lambda x: f"{x:8.4f}",
            },
        )
    )

    if not selected_df.empty:
        print()
        print("SELECTED REPRESENTATIVES:")
        for _, row in selected_df.iterrows():
            print(
                f"  {int(row['RepresentativeRank'])}. "
                f"{row['GroupName']} "
                f"(type={row['AtomType']}, element={row['Element']}, "
                f"n={int(row['Count'])}, "
                f"AW={row['Mean_AW']:.3f}, "
                f"Pow={row['Mean_Power']:.3f}, "
                f"diff={row['Mean_Signed_%_Diff']:.3f}%, "
                f"spread={row['Normalized_Spread']:.4f})"
            )


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def add_covariance_ellipse(
    ax,
    x,
    y,
    color,
    n_std: float = ELLIPSE_N_STD,
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    if len(x) < 3:
        return

    cov = np.cov(x, y)
    if np.any(~np.isfinite(cov)):
        return

    vals, vecs = np.linalg.eigh(cov)
    if np.any(vals < 0):
        return

    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width = 2.0 * n_std * np.sqrt(vals[0])
    height = 2.0 * n_std * np.sqrt(vals[1])

    fill = Ellipse(
        (np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        edgecolor=color,
        alpha=0.08,
        linewidth=1.2,
        zorder=2,
    )
    ax.add_patch(fill)

    edge = Ellipse(
        (np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=angle,
        facecolor="none",
        edgecolor=color,
        alpha=0.65,
        linewidth=1.3,
        zorder=2.1,
    )
    edge.set_path_effects(
        [pe.Stroke(linewidth=3.0, foreground="white"), pe.Normal()]
    )
    ax.add_patch(edge)


def get_axis_limits(
    x: np.ndarray,
    y: np.ndarray,
    metric: str,
) -> Tuple[float, float]:
    values = np.concatenate([x, y])
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return 0.0, 1.0

    if metric == "Contacts":
        return max(0.0, float(np.min(values)) - 0.75), float(np.max(values)) + 0.75

    low = float(np.percentile(values, 0.25))
    high = float(np.percentile(values, 99.75))
    span = high - low
    pad = 0.05 * span if span > 0 else 1.0

    return max(0.0, low - pad), high + pad


def label_representative_groups(
    ax,
    selected_df: pd.DataFrame,
):
    if selected_df.empty or not LABEL_SELECTED_GROUPS:
        return

    for _, row in selected_df.iterrows():
        color = ELEMENT_COLORS.get(row["Element"], "0.35")
        ax.annotate(
            str(row["GroupName"]),
            xy=(row["Mean_AW"], row["Mean_Power"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8.5,
            color=color,
            fontweight="bold",
            zorder=7,
        )


def get_plot_sample(
    plot_df: pd.DataFrame,
    broad_class_name: str,
) -> pd.DataFrame:
    """Downsample only the raw protein scatter layer."""
    if broad_class_name != "Protein":
        return plot_df

    if PROTEIN_PLOT_FRACTION >= 1.0:
        return plot_df

    if PROTEIN_PLOT_FRACTION <= 0.0:
        return plot_df.iloc[0:0].copy()

    return plot_df.sample(
        frac=PROTEIN_PLOT_FRACTION,
        random_state=PLOT_RANDOM_SEED,
    )


def plot_atomic_metric(
    ax,
    df: pd.DataFrame,
    metric: str,
    broad_class_name: str,
    shared_limits: Optional[Tuple[float, float]] = None,
):
    aw_col, power_col = METRIC_COLUMNS[metric]

    plot_df = df[df["BroadClass"] == broad_class_name].copy()
    plot_df[aw_col] = pd.to_numeric(plot_df[aw_col], errors="coerce")
    plot_df[power_col] = pd.to_numeric(plot_df[power_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[aw_col, power_col])

    if plot_df.empty:
        ax.text(
            0.5,
            0.5,
            f"No {broad_class_name} data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return pd.DataFrame()

    # Raw atoms. Protein is downsampled for plotting only.
    raw_plot_df = get_plot_sample(
        plot_df,
        broad_class_name=broad_class_name,
    )

    for element, element_df in raw_plot_df.groupby("ElementClass"):
        color = ELEMENT_COLORS.get(element, "0.45")
        ax.scatter(
            element_df[aw_col],
            element_df[power_col],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            color=color,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

    stats = compute_group_stats(plot_df, metric=metric)

    representative_stats = select_representative_groups(
        stats,
        broad_class_name=broad_class_name,
        metric=metric,
    )

    print_grouping_table(
        stats_df=stats,
        selected_df=representative_stats,
        broad_class_name=broad_class_name,
        metric=metric,
    )

    # Continuous metrics get covariance envelopes. Contacts are discrete
    # integer counts, so ellipses are intentionally omitted there.
    for _, row in representative_stats.iterrows():
        group_df = plot_df[
            plot_df["CanonicalName"] == row["GroupName"]
        ]
        color = ELEMENT_COLORS.get(row["Element"], "0.35")

        if metric != "Contacts":
            add_covariance_ellipse(
                ax,
                group_df[aw_col].to_numpy(float),
                group_df[power_col].to_numpy(float),
                color=color,
            )

        # The center remains useful for all three metrics. For Contacts it
        # summarizes a discrete group without implying a continuous Gaussian
        # population.
        ax.scatter(
            row["Mean_AW"],
            row["Mean_Power"],
            s=GROUP_POINT_SIZE + (15 if metric == "Contacts" else 0),
            color=color,
            edgecolor="white",
            linewidth=0.9,
            zorder=5,
        )

    x = plot_df[aw_col].to_numpy(float)
    y = plot_df[power_col].to_numpy(float)

    if shared_limits is None:
        lo, hi = get_axis_limits(x, y, metric)
    else:
        lo, hi = shared_limits

    ax.plot(
        [lo, hi],
        [lo, hi],
        linestyle="--",
        color="black",
        linewidth=1.5,
        alpha=0.8,
        zorder=0,
    )

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    if metric == "Volume":
        unit = r" ($\AA^3$)"
    elif metric == "Surface Area":
        unit = r" ($\AA^2$)"
    else:
        unit = ""

    ax.set_xlabel(f"AW {metric}{unit}", fontsize=11)
    ax.set_ylabel(f"Power {metric}{unit}", fontsize=11)
    ax.set_title(metric, fontsize=13, fontweight="bold")

    ax.tick_params(axis="both", labelsize=9, width=1.2, length=5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    label_representative_groups(ax, representative_stats)

    stats = stats.copy()

    representative_names = set(
        representative_stats["GroupName"]
    ) if not representative_stats.empty else set()

    stats["Is_Representative"] = stats["GroupName"].isin(
        representative_names
    )

    if not representative_stats.empty:
        rank_map = dict(
            zip(
                representative_stats["GroupName"],
                representative_stats["RepresentativeRank"],
            )
        )

    return stats


def add_element_legend(fig):
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=7,
            label=element,
        )
        for element, color in ELEMENT_COLORS.items()
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        fontsize=10,
        bbox_to_anchor=(0.5, 0.01),
    )


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def compute_shared_metric_limits(
    atom_df: pd.DataFrame,
    metric: str,
) -> Tuple[float, float]:
    """Use one AW/Power range for protein and nucleic panels of a metric."""
    aw_col, power_col = METRIC_COLUMNS[metric]

    vals = atom_df[[aw_col, power_col]].apply(
        pd.to_numeric,
        errors="coerce",
    )

    vals = vals.to_numpy(float).ravel()
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return 0.0, 1.0

    if metric == "Contacts":
        return (
            max(0.0, float(np.min(vals)) - 0.75),
            float(np.max(vals)) + 0.75,
        )

    low = float(np.percentile(vals, 0.25))
    high = float(np.percentile(vals, 99.75))
    span = high - low
    pad = 0.05 * span if span > 0 else 1.0

    lower = max(0.0, low - pad)
    upper = high + pad

    min_upper = MIN_SHARED_UPPER_LIMIT.get(metric)
    if min_upper is not None:
        upper = max(upper, float(min_upper))

    return lower, upper


def make_figure(atom_df: pd.DataFrame, figure_dir: str):
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=FIGSIZE,
        constrained_layout=False,
    )

    row_classes = ["Protein", "Nucleic acid"]
    metrics = ["Volume", "Surface Area", "Contacts"]

    shared_limits = {
        metric: compute_shared_metric_limits(atom_df, metric)
        for metric in metrics
    }

    print()
    print("F3B SHARED AXIS LIMITS")
    for metric, limits in shared_limits.items():
        print(f"  {metric}: {limits[0]:.3f} -> {limits[1]:.3f}")

    all_stats: Dict[str, List[pd.DataFrame]] = {
        "Protein": [],
        "Nucleic acid": [],
    }

    for row_i, broad_name in enumerate(row_classes):
        for col_i, metric in enumerate(metrics):
            stats = plot_atomic_metric(
                axes[row_i, col_i],
                df=atom_df,
                metric=metric,
                broad_class_name=broad_name,
                shared_limits=shared_limits[metric],
            )

            if not stats.empty:
                stats = stats.copy()
                stats.insert(0, "BroadClass", broad_name)
                all_stats[broad_name].append(stats)

        axes[row_i, 0].text(
            -0.28,
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
        "Figure 3B — Atomic differences between AW and Power partitions",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )

    add_element_legend(fig)

    fig.subplots_adjust(
        left=0.08,
        right=0.985,
        top=0.92,
        bottom=0.09,
        wspace=0.27,
        hspace=0.28,
    )

    if SAVE_PNG:
        path = figure_dir / "F3B_atomic_groups.png"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {path}")

    if SAVE_SVG:
        path = figure_dir / "F3B_atomic_groups.svg"
        fig.savefig(path, bbox_inches="tight")
        print(f"Saved: {path}")

    for broad_name, pieces in all_stats.items():
        if not pieces:
            continue

        summary = pd.concat(pieces, ignore_index=True)

        filename = (
            "F3B_protein_group_summary.csv"
            if broad_name == "Protein"
            else "F3B_nucleic_group_summary.csv"
        )

        path = figure_dir / filename
        summary.to_csv(path, index=False)
        print(f"Saved: {path}")

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
        data_root = select_directory("Select Figure 3 data folder")

    if not data_root:
        print("No data folder selected.")
        return

    if figure_dir is None:
        figure_dir = select_directory("Select figures/Figure_3 folder")

    if not figure_dir:
        print("No figure output folder selected.")
        return

    atom_df = collect_atomic_data(
        data_root=data_root,
        exclude_keys=EXCLUDE_KEYS,
    )

    if atom_df.empty:
        print("No protein or nucleic-acid atomic data were collected.")
        return

    make_figure(atom_df, figure_dir=figure_dir)


if __name__ == "__main__":
    main()