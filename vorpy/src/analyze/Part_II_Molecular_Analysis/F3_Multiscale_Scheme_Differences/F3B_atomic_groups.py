"""Figure 3B — Atomic AW vs Power comparisons.

Creates 12 separate Figure 3B plots (PNG and SVG):

    Protein       : Volume | Surface Area | Contacts
    DNA           : Volume | Surface Area | Contacts
    RNA           : Volume | Surface Area | Contacts
    Small molecule: Volume | Surface Area | Contacts

Only systems discovered by F3_common.discover_systems() are considered, so the
data root is restricted to folders named <letter>_<molecule name> and the
leading letter is still filtered by EXCLUDE_KEYS.

The script classifies each system from its residue content, then pools matched
atoms across systems. Atom groups use the same canonicalization rules already
developed for F2C_All_Data_Sets_Clustered.py.

Outputs
-------
F3B_<class>_<metric>.png and .svg (24 image files)
F3B_<class>_group_summary.csv (four CSV files)
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
    add_deviation_columns,
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

EXCLUDE_KEYS = ["A"]

# ---------------------------------------------------------------------------
# Instance cache + Figure 3D candidate exports
# ---------------------------------------------------------------------------
# Set True only when the underlying logs, parsing, or metric construction changes.
REBUILD_INSTANCE_CACHE = True
INSTANCE_CACHE_NAME = "F3B_atomic_instance_cache.csv"

# Number of examples retained in each candidate class / molecular class.
OUTLIER_TOP_N = 25

# Minimum observations required for a canonical atom type to receive a stable
# within-type outlier score.
MIN_WITHIN_TYPE_COUNT = 5

MIN_GROUP_COUNT = 25
MIN_RANK_SYSTEMS = 2

SMALL_MOLECULE_MIN_GROUP_COUNT = 2
SMALL_MOLECULE_MIN_RANK_SYSTEMS = 1
TOP_RANK_N = 5

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

LABEL_SELECTED_GROUPS = True
ELLIPSE_N_STD = 1.25

# Raw atomic scatter
POINT_SIZE = 14
POINT_ALPHA = 0.30

# Representative group centers
GROUP_POINT_SIZE = 75

# Ellipse appearance
ELLIPSE_FILL_ALPHA = 0.06
ELLIPSE_EDGE_ALPHA = 0.65
ELLIPSE_LINEWIDTH = 0.65
ELLIPSE_HALO_LINEWIDTH = 0.0

# Typography
TICK_LABEL_SIZE = 20
AXIS_LABEL_SIZE = 20
PLOT_TITLE_SIZE = 18
GROUP_LABEL_SIZE = 13

# Figure dimensions
FIGSIZE = (7.2, 6.4)
DPI = 300

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
    "DNA": ["H", "C", "N", "O", "P"],
    "RNA": ["H", "C", "N", "O", "P"],
    "Small molecule": ["H", "C", "N", "O", "S", "P"],
}

# Minimum shared upper limits. These prevent selected representative groups
# from being clipped by percentile-based axis trimming.
MIN_SHARED_UPPER_LIMIT = {
    "Volume": 32.0,
    "Surface Area": 60.0,
    "Contacts": None,
}

SAVE_PNG = True
SAVE_SVG = True
SHOW = False  # Set True to display all 12 figures interactively

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DI"}
RNA_RESIDUES = {"A", "C", "G", "U", "I", "RA", "RC", "RG", "RU"}

# Non-molecular components excluded from the small-molecule class.
EXCLUDED_SMALL_RESIDUES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE",
    "NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN", "CU", "CO", "NI",
    "SO4", "PO4", "NO3", "NH4"}
MOLECULAR_CLASSES = ["Protein", "DNA", "RNA", "Small molecule"]

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
    """Classify systems by their dominant polymer residue type, or as small molecules."""
    if "Residue" not in atoms.columns:
        return None
    residues = atoms["Residue"].astype(str).str.strip().str.upper()
    counts = {
        "protein": int(residues.isin(PROTEIN_RESIDUES).sum()),
        "dna": int(residues.isin(DNA_RESIDUES).sum()),
        "rna": int(residues.isin(RNA_RESIDUES).sum()),
    }
    best = max(counts, key=counts.get)
    if counts[best] > 0:
        return best
    return "small_molecule" if (~residues.isin(EXCLUDED_SMALL_RESIDUES)).any() else None


def broad_class(system_class: str) -> str:
    return {"protein": "Protein", "dna": "DNA", "rna": "RNA",
            "small_molecule": "Small molecule"}[system_class]


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

    if system_class == "small_molecule":
        # Preserve chemical atom labels rather than forcing a protein/nucleotide
        # canonicalization rule onto arbitrary ligand names.
        df["CanonicalName"] = df["Atom"].astype(str).str.strip().str.upper()
    else:
        df["CanonicalName"] = [
            canonicalize_atom_name(atom_name=str(atom), molecule_class=system_class,
                                   residue_name=str(residue))
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

        # Keep only the dominant polymer class, or non-solvent small molecules.
        atom_df["Residue"] = atom_df["Residue"].astype(str).str.strip().str.upper()
        if system_class == "protein":
            keep = atom_df["Residue"].isin(PROTEIN_RESIDUES)
        elif system_class == "dna":
            keep = atom_df["Residue"].isin(DNA_RESIDUES)
        elif system_class == "rna":
            keep = atom_df["Residue"].isin(RNA_RESIDUES)
        else:
            keep = ~atom_df["Residue"].isin(EXCLUDED_SMALL_RESIDUES)

        atom_df = atom_df.loc[keep].copy()
        if atom_df.empty:
            print("  skipped: no retained atoms after residue filtering")
            continue

        # Preserve instance-level signed and absolute deviations for Figure 3D
        # candidate selection. This does not alter Figure 3B plotting.
        atom_df = add_deviation_columns(atom_df)

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

        rms_difference = float(np.sqrt(np.mean((power - aw) ** 2)))
        n_systems = int(group_df.loc[values.index, "System"].nunique()) if "System" in group_df else 0
        delta = mean_power - mean_aw
        pct = np.nan if mean_aw == 0 else 100.0 * delta / mean_aw

        rows.append(
            {
                "Metric": metric,
                "GroupName": group_name,
                "AtomType": str(group_name).split("_")[0],
                "Element": element_from_group(group_name),
                "Count": len(values),
                "Systems": n_systems,
                "RMS_AW_Power_Difference": rms_difference,
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


MAX_GROUPS_PER_ELEMENT = 3
MIN_SYSTEMS_FOR_REPRESENTATIVES = 2


def select_representative_groups(
    stats_df: pd.DataFrame,
    broad_class_name: str,
    metric: str,
    n_groups: int = REPRESENTATIVE_GROUPS,
    compactness_quantile: float = COMPACTNESS_QUANTILE,
) -> pd.DataFrame:
    """Select compact, well-sampled, chemically varied groups on both sides of y=x.

    Constraints (when feasible): at most three representatives per element,
    at least one Power > AW and one Power < AW. Prefer groups seen in two or
    more systems, but allow single-system exploratory choices for RNA.
    """
    if stats_df.empty:
        return stats_df.copy()

    candidates = stats_df.loc[
        np.isfinite(stats_df["Normalized_Spread"])
        & np.isfinite(stats_df["Mean_AW"])
        & np.isfinite(stats_df["Mean_Power"])
    ].copy()
    if candidates.empty:
        return candidates

    # Some datasets (RNA, for example) contain only one independent system.
    required_systems = (SMALL_MOLECULE_MIN_RANK_SYSTEMS
                        if broad_class_name == "Small molecule" else MIN_SYSTEMS_FOR_REPRESENTATIVES)
    if "Systems" in candidates.columns and (candidates["Systems"] >= required_systems).any():
        candidates = candidates[candidates["Systems"] >= required_systems].copy()
    else:
        print(f"NOTE: {broad_class_name} / {metric}: exploratory single-system representative selection")

    candidates["BaseScore"] = _base_group_score(candidates)
    candidates["Side"] = np.where(
        candidates["Mean_Power"] > candidates["Mean_AW"], 1,
        np.where(candidates["Mean_Power"] < candidates["Mean_AW"], -1, 0),
    )
    candidates = candidates.sort_values(
        ["BaseScore", "Normalized_Spread", "Count", "GroupName"],
        ascending=[False, True, False, True],
    ).reset_index(drop=True)
    available_elements = set(candidates["Element"].astype(str))
    desired_elements = [e for e in REQUIRED_ELEMENTS.get(broad_class_name, []) if e in available_elements]
    if not desired_elements:
        desired_elements = list(dict.fromkeys(candidates["Element"].astype(str)))
    # Preserve coverage only when there are enough slots; do not exceed the cap.
    desired_elements = desired_elements[:n_groups]

    selected = []
    selected_names = set()
    element_counts = {}
    coord_source = candidates.reset_index(drop=True)
    coord_map = dict(zip(coord_source["GroupName"], _normalized_group_coordinates(coord_source)))

    def add(row):
        name = str(row["GroupName"])
        element = str(row["Element"])
        selected.append(row.copy())
        selected_names.add(name)
        element_counts[element] = element_counts.get(element, 0) + 1

    def eligible(pool):
        return pool.loc[
            (~pool["GroupName"].isin(selected_names))
            & (pool["Element"].astype(str).map(lambda e: element_counts.get(e, 0) < MAX_GROUPS_PER_ELEMENT))
        ]

    def best(pool):
        pool = eligible(pool)
        if pool.empty:
            return None
        scores = []
        for _, row in pool.iterrows():
            coord = coord_map.get(row["GroupName"])
            if coord is None or not selected:
                separation = 1.0
            else:
                separation = min(float(np.linalg.norm(coord - coord_map[r["GroupName"]]))
                                 for r in selected if r["GroupName"] in coord_map) / np.sqrt(2.0)
            scores.append(0.85 * float(row["BaseScore"]) + 0.15 * separation)
        return pool.iloc[int(np.argmax(scores))]

    # First cover the required chemical elements, leaving room to fix sign balance.
    for element in desired_elements:
        if len(selected) >= n_groups:
            break
        row = best(candidates[candidates["Element"].astype(str) == element])
        if row is not None:
            add(row)

    # Guarantee representatives on both sides whenever the eligible data allow it.
    for side in (1, -1):
        if any(int(r["Side"]) == side for r in selected):
            continue
        row = best(candidates[candidates["Side"] == side])
        if row is not None and len(selected) < n_groups:
            add(row)
        elif row is not None and selected:
            # Replace the weakest group without losing required element coverage,
            # the opposite side, or the per-element limit.
            choices = []
            for i, old in enumerate(selected):
                old_el = str(old["Element"])
                new_el = str(row["Element"])
                if old_el != new_el and element_counts.get(new_el, 0) >= MAX_GROUPS_PER_ELEMENT:
                    continue
                if old_el in desired_elements and element_counts[old_el] == 1 and old_el != new_el:
                    continue
                remaining = [r for j, r in enumerate(selected) if j != i]
                if any(int(r["Side"]) == -side for r in selected) and not any(int(r["Side"]) == -side for r in remaining) and int(row["Side"]) != -side:
                    continue
                choices.append(i)
            if choices:
                i = min(choices, key=lambda j: float(selected[j]["BaseScore"]))
                old = selected.pop(i)
                element_counts[str(old["Element"])] -= 1
                selected_names.remove(str(old["GroupName"]))
                add(row)

    while len(selected) < n_groups:
        row = best(candidates)
        if row is None:
            break
        add(row)

    sides = {int(r["Side"]) for r in selected}
    for side, label in ((1, "above"), (-1, "below")):
        if side not in sides:
            print(f"WARNING: {broad_class_name} / {metric}: no representative {label} y=x under selection constraints")
    if len(selected) < n_groups:
        print(f"WARNING: {broad_class_name} / {metric}: selected only {len(selected)}/{n_groups} groups")

    result = pd.DataFrame(selected).reset_index(drop=True)
    result["RepresentativeRank"] = np.arange(1, len(result) + 1)
    result["SelectionMode"] = "automatic_element_and_sign_balanced"
    return result.drop(columns=["BaseScore", "Side"], errors="ignore")


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
        alpha=ELLIPSE_FILL_ALPHA,
        linewidth=0,
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
        alpha=ELLIPSE_EDGE_ALPHA,
        linewidth=ELLIPSE_LINEWIDTH,
        zorder=2.1,
    )
    if ELLIPSE_HALO_LINEWIDTH > 0:
        edge.set_path_effects([
            pe.Stroke(
                linewidth=ELLIPSE_HALO_LINEWIDTH,
                foreground="white"
            ),
            pe.Normal()
        ])
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
            fontsize=GROUP_LABEL_SIZE,
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

    min_count = (SMALL_MOLECULE_MIN_GROUP_COUNT
                 if broad_class_name == "Small molecule" else MIN_GROUP_COUNT)
    stats = compute_group_stats(plot_df, metric=metric, min_count=min_count)

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

    ax.set_xlabel(f"AW {metric}{unit}", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(f"Power {metric}{unit}", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(metric, fontsize=PLOT_TITLE_SIZE, fontweight="bold")

    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE, width=1.2, length=5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    label_representative_groups(ax, representative_stats)

    stats = stats.copy()

    if stats.empty:
        stats["Is_Representative"] = pd.Series(dtype=bool)
    else:
        representative_names = (
            set(representative_stats["GroupName"])
            if not representative_stats.empty
            else set()
        )
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



# ---------------------------------------------------------------------------
# Instance cache + outlier analysis for Figure 3D
# ---------------------------------------------------------------------------

def _robust_high_z(values: pd.Series) -> pd.Series:
    """MAD-based one-sided robust z: positive means unusually high."""
    vals = pd.to_numeric(values, errors="coerce").astype(float)
    finite = np.isfinite(vals)

    result = pd.Series(np.nan, index=vals.index, dtype=float)
    if not finite.any():
        return result

    valid = vals[finite]
    median = float(np.median(valid))
    mad = float(np.median(np.abs(valid - median)))

    if mad <= 1e-12:
        result.loc[finite] = 0.0
        return result

    result.loc[finite] = 0.67448975 * (valid - median) / mad
    return result


def _within_group_percentile(
    df: pd.DataFrame,
    value_col: str,
    group_cols: List[str],
    min_count: int,
) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype=float)

    for _, group in df.groupby(group_cols, dropna=False):
        vals = pd.to_numeric(group[value_col], errors="coerce")
        valid = vals.notna() & np.isfinite(vals)

        if int(valid.sum()) < min_count:
            continue

        ranks = vals.loc[valid].rank(method="average", pct=True)
        out.loc[ranks.index] = ranks.astype(float)

    return out


def _within_group_robust_z(
    df: pd.DataFrame,
    value_col: str,
    group_cols: List[str],
    min_count: int,
) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype=float)

    for _, group in df.groupby(group_cols, dropna=False):
        vals = pd.to_numeric(group[value_col], errors="coerce")
        valid = vals.notna() & np.isfinite(vals)

        if int(valid.sum()) < min_count:
            continue

        z = _robust_high_z(vals.loc[valid])
        out.loc[z.index] = z

    return out


def build_atomic_outlier_table(atom_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add global and within-canonical-type outlier scores to every atom instance.

    Joint scores use Volume and Surface Area only because these are the
    continuous geometric quantities intended for Figure 3D.
    """
    df = atom_df.copy()

    for metric in ["Volume", "Surface Area", "Contacts"]:
        abs_col = f"{metric} Abs % Diff"
        if abs_col not in df.columns:
            signed_col = f"{metric} Signed % Diff"
            if signed_col in df.columns:
                df[abs_col] = np.abs(pd.to_numeric(df[signed_col], errors="coerce"))

    metric_cols = {
        "Volume": "Volume Abs % Diff",
        "Surface Area": "Surface Area Abs % Diff",
    }

    for metric, value_col in metric_cols.items():
        vals = pd.to_numeric(df[value_col], errors="coerce")
        df[f"{metric} Percentile"] = vals.rank(method="average", pct=True)

        df[f"{metric} Within-Type Percentile"] = _within_group_percentile(
            df,
            value_col=value_col,
            group_cols=["BroadClass", "CanonicalName"],
            min_count=MIN_WITHIN_TYPE_COUNT,
        )

        df[f"{metric} Within-Type Robust Z"] = _within_group_robust_z(
            df,
            value_col=value_col,
            group_cols=["BroadClass", "CanonicalName"],
            min_count=MIN_WITHIN_TYPE_COUNT,
        )

    df["Joint Outlier Score"] = df[
        ["Volume Percentile", "Surface Area Percentile"]
    ].mean(axis=1, skipna=False)

    df["Within-Type Joint Outlier Score"] = df[
        ["Volume Within-Type Percentile", "Surface Area Within-Type Percentile"]
    ].mean(axis=1, skipna=False)

    df["Within-Type Mean Robust Z"] = df[
        ["Volume Within-Type Robust Z", "Surface Area Within-Type Robust Z"]
    ].mean(axis=1, skipna=False)

    volume_abs = pd.to_numeric(df["Volume Abs % Diff"], errors="coerce")
    sa_abs = pd.to_numeric(df["Surface Area Abs % Diff"], errors="coerce")
    df["Max Metric"] = np.where(volume_abs >= sa_abs, "Volume", "Surface Area")

    # Explicit ranks make sorting/filtering straightforward outside Python.
    df["Absolute Joint Rank"] = (
        df.groupby("BroadClass")["Joint Outlier Score"]
        .rank(method="min", ascending=False)
    )
    df["Within-Type Joint Rank"] = (
        df.groupby("BroadClass")["Within-Type Joint Outlier Score"]
        .rank(method="min", ascending=False)
    )
    df["Low-Deviation Rank"] = (
        df.groupby("BroadClass")["Joint Outlier Score"]
        .rank(method="min", ascending=True)
    )

    return df


def summarize_atomic_outlier_groups(outlier_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for (broad, group_name), group in outlier_df.groupby(
        ["BroadClass", "CanonicalName"],
        dropna=False,
    ):
        rows.append({
            "BroadClass": broad,
            "CanonicalName": group_name,
            "Element": group["ElementClass"].iloc[0] if "ElementClass" in group else "",
            "Count": len(group),
            "Systems": group["System"].nunique() if "System" in group else np.nan,
            "Median Volume Abs % Diff": pd.to_numeric(
                group["Volume Abs % Diff"], errors="coerce"
            ).median(),
            "Median Surface Area Abs % Diff": pd.to_numeric(
                group["Surface Area Abs % Diff"], errors="coerce"
            ).median(),
            "Mean Joint Outlier Score": pd.to_numeric(
                group["Joint Outlier Score"], errors="coerce"
            ).mean(),
            "95th % Joint Outlier Score": pd.to_numeric(
                group["Joint Outlier Score"], errors="coerce"
            ).quantile(0.95),
            "Maximum Joint Outlier Score": pd.to_numeric(
                group["Joint Outlier Score"], errors="coerce"
            ).max(),
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values(
            ["BroadClass", "Median Volume Abs % Diff", "Median Surface Area Abs % Diff"],
            ascending=[True, False, False],
        )
    return result


def build_atomic_top_candidates(outlier_df: pd.DataFrame) -> pd.DataFrame:
    pieces = []

    for broad_name in MOLECULAR_CLASSES:
        class_df = outlier_df[outlier_df["BroadClass"] == broad_name].copy()
        if class_df.empty:
            continue

        selections = [
            (
                "Absolute high",
                class_df.sort_values("Joint Outlier Score", ascending=False).head(OUTLIER_TOP_N),
            ),
            (
                "Within-type high",
                class_df.sort_values(
                    ["Within-Type Joint Outlier Score", "Within-Type Mean Robust Z"],
                    ascending=[False, False],
                ).head(OUTLIER_TOP_N),
            ),
            (
                "Low-deviation control",
                class_df.sort_values("Joint Outlier Score", ascending=True).head(OUTLIER_TOP_N),
            ),
        ]

        for candidate_set, selected in selections:
            selected = selected.copy()
            selected.insert(0, "Candidate Set", candidate_set)
            pieces.append(selected)

    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def select_signed_atomic_examples(atom_df: pd.DataFrame) -> pd.DataFrame:
    """Choose two AW>Power and two AW<Power volume cells per class.

    Rank by absolute volume difference (percentage), using surface-area
    disagreement as a tiebreaker. Prefer distinct systems and residues.
    Never manufacture a missing direction: return fewer than two if necessary.
    """
    df = atom_df.copy()
    for col in ("AW Volume", "Power Volume", "AW Surface Area", "Power Surface Area"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["AW Volume", "Power Volume"])
    df = df.loc[(df["AW Volume"] > 0) & (df["Power Volume"] > 0)].copy()
    df["AW Minus Power Volume"] = df["AW Volume"] - df["Power Volume"]
    df["Volume Signed % (AW relative to Power)"] = (
        100 * df["AW Minus Power Volume"] / df["Power Volume"]
    )
    df["AW Minus Power Surface Area"] = df["AW Surface Area"] - df["Power Surface Area"]
    df["Surface Area Signed % (AW relative to Power)"] = (
        100 * df["AW Minus Power Surface Area"] / df["Power Surface Area"].where(df["Power Surface Area"] > 0)
    )
    df["Volume Magnitude %"] = df["Volume Signed % (AW relative to Power)"].abs()
    df["Surface Area Magnitude %"] = df["Surface Area Signed % (AW relative to Power)"].abs()
    selected = []
    for broad in MOLECULAR_CLASSES:
        class_df = df.loc[df["BroadClass"] == broad]
        for label, sign in (("AW larger", 1), ("AW smaller", -1)):
            pool = class_df.loc[class_df["AW Minus Power Volume"] * sign > 0].copy()
            pool = pool.sort_values(
                ["Volume Magnitude %", "Surface Area Magnitude %"],
                ascending=[False, False], kind="stable"
            )
            picks = []
            for _, row in pool.iterrows():
                if len(picks) == 2:
                    break
                # Favor different systems when available, then different residues.
                if picks and row.get("System") == picks[0].get("System") and pool["System"].nunique() > 1:
                    continue
                picks.append(row)
            if len(picks) < 2:
                for _, row in pool.iterrows():
                    if len(picks) == 2:
                        break
                    if any(row.name == prev.name for prev in picks):
                        continue
                    picks.append(row)
            if len(picks) < 2:
                print(f"WARNING: {broad}: only {len(picks)} candidate(s) with {label} volume")
            for rank, row in enumerate(picks, 1):
                record = row.copy()
                record["Selection Direction"] = label
                record["Selection Rank"] = rank
                selected.append(record)
    return pd.DataFrame(selected)


def atomic_cache_path(figure_dir: str) -> Path:
    return Path(figure_dir) / INSTANCE_CACHE_NAME


def print_atomic_cache_status(figure_dir: str):
    path = atomic_cache_path(figure_dir)
    print()
    print("=" * 80)
    print("F3B INSTANCE CACHE")
    print("=" * 80)
    print(f"Path:    {path}")
    print(f"Exists:  {path.exists()}")
    print(f"Rebuild: {REBUILD_INSTANCE_CACHE}")
    print()


def load_or_collect_atomic_data(
    figure_dir: str,
    data_root: Optional[str] = None,
    exclude_keys=None,
) -> pd.DataFrame:
    path = atomic_cache_path(figure_dir)

    if path.exists() and not REBUILD_INSTANCE_CACHE:
        print(f"Loading cached F3B atomic instances: {path}")
        cached = pd.read_csv(path)
        if not {"SystemClass", "BroadClass", "CanonicalName"}.issubset(cached.columns):
            raise ValueError("Outdated instance cache: set REBUILD_INSTANCE_CACHE=True")
        cached["BroadClass"] = cached["SystemClass"].map(broad_class)
        print(f"  loaded {len(cached):,} atom instances")
        print("  Note: caches built before small-molecule support may omit those systems; rebuild to include them.")
        return cached

    if not data_root:
        data_root = select_directory(
            "Select MOLECULAR DATA ROOT containing D_Hairpin, E_Cambrin, etc."
        )

    if not data_root:
        print("No Figure 3 data folder selected; cannot rebuild atomic cache.")
        return pd.DataFrame()

    print("Building F3B atomic instance cache from AW/Power logs ...")

    atom_df = collect_atomic_data(
        data_root=data_root,
        exclude_keys=exclude_keys,
    )

    if not atom_df.empty:
        Path(figure_dir).mkdir(parents=True, exist_ok=True)
        atom_df.to_csv(path, index=False)
        print(f"Saved atomic instance cache: {path}")

    return atom_df


def export_atomic_outliers(atom_df: pd.DataFrame, figure_dir: str):
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    outliers = build_atomic_outlier_table(atom_df)
    group_summary = summarize_atomic_outlier_groups(outliers)
    candidates = build_atomic_top_candidates(outliers)
    signed_examples = select_signed_atomic_examples(outliers)

    instance_path = figure_dir / "F3B_atomic_instance_outliers.csv"
    group_path = figure_dir / "F3B_atomic_group_outliers.csv"
    candidate_path = figure_dir / "F3B_atomic_top_candidates.csv"

    outliers.to_csv(instance_path, index=False)
    group_summary.to_csv(group_path, index=False)
    candidates.to_csv(candidate_path, index=False)
    signed_path = figure_dir / "F3B_signed_atomic_examples.csv"
    signed_examples.to_csv(signed_path, index=False)

    print(f"Saved: {instance_path}")
    print(f"Saved: {group_path}")
    print(f"Saved: {candidate_path}")
    print(f"Saved: {signed_path}")
    if not signed_examples.empty:
        cols = [c for c in ("BroadClass", "Selection Direction", "Selection Rank", "System", "Index", "Atom", "Residue", "Residue Sequence", "Chain", "AW Volume", "Power Volume", "Volume Signed % (AW relative to Power)", "AW Surface Area", "Power Surface Area") if c in signed_examples.columns]
        print("\nFOUR SIGNED VOLUME CANDIDATES PER MOLECULAR CLASS:")
        print(signed_examples[cols].to_string(index=False))

    preview_cols = [
        c for c in [
            "Candidate Set", "BroadClass", "System", "Index", "Atom",
            "CanonicalName", "Residue", "Residue Sequence", "Chain",
            "Volume Abs % Diff", "Surface Area Abs % Diff",
            "Joint Outlier Score", "Within-Type Joint Outlier Score",
            "Within-Type Mean Robust Z", "Max Metric",
        ]
        if c in candidates.columns
    ]

    if not candidates.empty and preview_cols:
        print()
        print("=" * 120)
        print("F3B TOP ATOMIC CANDIDATES FOR FIGURE 3D")
        print("=" * 120)
        print(candidates[preview_cols].head(60).to_string(index=False))


def print_top_five_rankings(stats: pd.DataFrame, class_name: str, metric: str):
    """Three independent group rankings, using complete (not sampled) data."""
    print("\n" + "=" * 90)
    print(f"TOP-FIVE RANKINGS — {class_name.upper()} | {metric.upper()}")
    print("=" * 90)
    if stats.empty:
        print("No eligible canonical groups.")
        return

    min_count = SMALL_MOLECULE_MIN_GROUP_COUNT if class_name == "Small molecule" else MIN_GROUP_COUNT
    min_systems = SMALL_MOLECULE_MIN_RANK_SYSTEMS if class_name == "Small molecule" else MIN_RANK_SYSTEMS
    eligible = stats.loc[(stats["Count"] >= min_count) &
                         (stats["Systems"] >= min_systems)].copy()
    specs = [
        ("MOST NUMEROUS (all groups)", stats, "Count", False),
        (f"TIGHTEST (>= {min_count} atoms, >= {min_systems} systems)",
         eligible, "Normalized_Spread", True),
        (f"FURTHEST FROM y=x (>= {min_count} atoms, >= {min_systems} systems)",
         eligible, "RMS_AW_Power_Difference", False),
    ]
    cols = ["GroupName", "Count", "Systems", "Normalized_Spread",
            "RMS_AW_Power_Difference", "Mean_AW", "Mean_Power"]
    for title, frame, sort_col, ascending in specs:
        print(f"\n{title}:")
        if frame.empty:
            print("  No groups satisfy the selection criteria.")
            continue
        ranked = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=[sort_col])
        ranked = ranked.sort_values([sort_col, "Count"], ascending=[ascending, False]).head(TOP_RANK_N)
        print(ranked[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))


def make_figure(atom_df: pd.DataFrame, figure_dir: str):
    """Save one independent PNG/SVG scatter plot per molecular class and metric."""
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    metrics = ["Volume", "Surface Area", "Contacts"]

    for class_name in MOLECULAR_CLASSES:
        class_df = atom_df.loc[atom_df["BroadClass"] == class_name]
        summaries = []
        for metric in metrics:
            fig, ax = plt.subplots(figsize=FIGSIZE)
            # Keep matching x/y ranges within each panel, without letting
            # small-molecule outliers compress all polymer plots.
            limits = compute_shared_metric_limits(class_df, metric)
            stats = plot_atomic_metric(ax, atom_df, metric, class_name,
                                       shared_limits=limits)
            ax.set_title(f"{class_name} — {metric}", fontsize=PLOT_TITLE_SIZE, fontweight="bold")
            fig.tight_layout()
            stem = f"F3B_{class_name.lower().replace(' ', '_')}_{metric.lower().replace(' ', '_')}"
            if SAVE_PNG:
                path = figure_dir / f"{stem}.png"
                fig.savefig(path, dpi=DPI)
                print(f"Saved: {path}")
            if SAVE_SVG:
                path = figure_dir / f"{stem}.svg"
                fig.savefig(path)
                print(f"Saved: {path}")
            if SHOW:
                plt.show()
            plt.close(fig)
            print_top_five_rankings(stats, class_name, metric)
            if not stats.empty:
                stats = stats.copy()
                stats.insert(0, "BroadClass", class_name)
                summaries.append(stats)

        if summaries:
            path = figure_dir / f"F3B_{class_name.lower().replace(' ', '_')}_group_summary.csv"
            pd.concat(summaries, ignore_index=True).to_csv(path, index=False)
            print(f"Saved: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(
    data_root: Optional[str] = DATA_ROOT,
    figure_dir: Optional[str] = FIGURE_DIR,
):
    # Rebuild workflow:
    #   1) choose the molecular DATA ROOT first
    #   2) choose the Figure 3 OUTPUT folder second
    #
    # Cached workflow:
    #   only the OUTPUT folder is needed.
    if REBUILD_INSTANCE_CACHE:
        if data_root is None:
            data_root = select_directory(
                "STEP 1/2 — Select MOLECULAR DATA ROOT containing D_Hairpin, E_Cambrin, etc."
            )

        if not data_root:
            print("No molecular data root selected.")
            return

        if figure_dir is None:
            figure_dir = select_directory(
                "STEP 2/2 — Select FIGURE OUTPUT folder for F3B files/cache"
            )
    else:
        if figure_dir is None:
            figure_dir = select_directory(
                "Select FIGURE OUTPUT folder containing the F3B cache"
            )

    if not figure_dir:
        print("No figure output folder selected.")
        return

    print()
    print(f"Molecular data root: {data_root if data_root else '(not needed; using cache)'}")
    print(f"Figure output dir:   {figure_dir}")

    print_atomic_cache_status(figure_dir)

    atom_df = load_or_collect_atomic_data(
        figure_dir=figure_dir,
        data_root=data_root,
        exclude_keys=EXCLUDE_KEYS,
    )

    if atom_df.empty:
        print("No molecular atomic data were collected.")
        return

    export_atomic_outliers(
        atom_df,
        figure_dir=figure_dir,
    )

    make_figure(atom_df, figure_dir=figure_dir)


if __name__ == "__main__":
    main()