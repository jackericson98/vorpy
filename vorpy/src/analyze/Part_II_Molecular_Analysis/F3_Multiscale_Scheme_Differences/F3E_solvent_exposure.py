"""Figure 3E — Continuous solvent-exposure dependence of AW/Power deviations.

This version replaces pooled low/intermediate/high exposure violins with a
continuous exposure analysis.

Layout
------
                    Volume | Surface Area | Contacts
    Atom
    Residue

Atom row
--------
Only the same representative atom groups used in Figure 3B are shown:

Protein:
    CA_*, CB_*, O, NZ, HN, SG

Nucleic acid:
    C5', C_sugar, N9, H_base_exocyclic, O_backbone, P

For each atom group:
    x = AW solvent-facing surface area (%)
    y = absolute Power-vs-AW percent deviation

Faint points show individual atoms.
A binned mean trend is drawn for groups with enough observations.

Residue row
-----------
All residue types are retained, but residue identity is controlled explicitly.

For each residue type, subtract its own mean exposure and its own mean
absolute deviation. This is a fixed-effect / within-type transformation:

    x = AW SolFacingPct - mean(AW SolFacingPct | residue type)
    y = Abs % Diff      - mean(Abs % Diff      | residue type)

The plotted residue relationship therefore asks:

    Within the SAME residue identity, does being more solvent exposed than
    usual correspond to a larger or smaller AW/Power deviation?

Protein and nucleic-acid residuals are plotted separately.

Methodological choice
---------------------
AW defines solvent exposure. Power is never allowed to change the exposure
variable.

The full-system PDB is used only to recover metadata for surface partner
indices omitted from the logs, especially solvent atoms. AW/Power geometric
metrics still come from the logs.

Outputs
-------
F3E_continuous_solvent_exposure.png
F3E_continuous_solvent_exposure.svg
F3E_atom_group_trends.csv
F3E_residue_fixed_effect_trends.csv
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
from matplotlib.lines import Line2D


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
    surface_pairs,
)

from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_All_Data_Sets_Clustered import (
    canonicalize_atom_name,
)


# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

DATA_ROOT = None
FIGURE_DIR = None

EXCLUDE_KEYS = ["A", "B", "C", "K", "L"]

FIGSIZE = (19, 10.8)
DPI = 300
SHOW = True
SAVE_PNG = True
SAVE_SVG = True

SOLVENT_RESIDUES = {
    "SOL", "HOH", "WAT", "H2O",
    "TIP3", "TIP3P", "TIP4P", "SPC", "SPCE",
}

AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DI"}
RNA_RESIDUES = {"A", "C", "G", "U", "I", "RA", "RC", "RG", "RU"}

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

# Plot only a fraction of raw atom points; statistics use all observations.
ATOM_POINT_FRACTION = 0.35
ATOM_POINT_ALPHA = 0.12
ATOM_POINT_SIZE = 9
RANDOM_SEED = 17

# Binned trend requirements.
TREND_BINS = 8
MIN_GROUP_COUNT = 25
MIN_BIN_COUNT = 5
TREND_LINE_WIDTH = 2.0
TREND_MARKER_SIZE = 4.5

# Residue fixed-effect scatter.
RESIDUE_POINT_ALPHA = 0.12
RESIDUE_POINT_SIZE = 10

# Display trimming only; regression/statistics use all finite data.
DISPLAY_PERCENTILE = 99.0

# Atom-group colors. Chosen to stay stable across all three metrics.
ATOM_GROUP_COLORS = {
    "CA_*": "#1f77b4",
    "CB_*": "#ff7f0e",
    "O": "#d62728",
    "NZ": "#2ca02c",
    "HN": "#9467bd",
    "SG": "#8c564b",
    "C5'": "#17becf",
    "C_sugar": "#bcbd22",
    "N9": "#2ca02c",
    "H_base_exocyclic": "#9467bd",
    "O_backbone": "#d62728",
    "P": "#e377c2",
}

CLASS_COLORS = {
    "Protein": "#4c78a8",
    "Nucleic acid": "#f58518",
}

METRICS = ["Volume", "Surface Area", "Contacts"]


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def select_directory(title: str) -> str:
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    path = filedialog.askdirectory(title=title)
    root.destroy()
    return path


def classify_system(atoms: pd.DataFrame) -> Optional[str]:
    if "Residue" not in atoms.columns:
        return None

    residues = atoms["Residue"].astype(str).str.strip().str.upper()

    counts = {
        "protein": int(residues.isin(AA_RESIDUES).sum()),
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


def allowed_residues(system_class: str):
    if system_class == "protein":
        return AA_RESIDUES
    if system_class == "dna":
        return DNA_RESIDUES
    if system_class == "rna":
        return RNA_RESIDUES
    return set()


def expected_pdb_path(system) -> Path:
    # SchemePaths.aw points to .../<system>/aw/aw_logs.csv
    system_folder = Path(system.aw).resolve().parent.parent

    molecule_name = getattr(system, "molecule_name", None)
    if not molecule_name:
        folder_name = system_folder.name
        molecule_name = (
            folder_name.split("_", 1)[1]
            if "_" in folder_name
            else folder_name
        )

    expected = system_folder / f"{molecule_name}.pdb"

    if expected.exists():
        return expected

    target = expected.name.lower()
    for candidate in system_folder.iterdir():
        if candidate.is_file() and candidate.name.lower() == target:
            return candidate

    return expected


def parse_full_pdb_metadata(pdb_path: Path) -> Dict[int, Dict]:
    """
    Zero-based ATOM/HETATM file order = original VorPy index.
    """
    metadata = {}
    atom_index = 0

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            rec = line[0:6].strip().upper()
            if rec not in {"ATOM", "HETATM"}:
                continue

            atom_name = line[12:16].strip()
            residue = line[17:20].strip().upper()
            chain = line[21:22].strip()

            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                resseq = -1

            metadata[atom_index] = {
                "Atom": atom_name,
                "Residue": residue,
                "Chain": chain,
                "Residue Sequence": resseq,
            }
            atom_index += 1

    return metadata


def validate_pdb_mapping(
    metadata: Dict[int, Dict],
    log_atoms: pd.DataFrame,
    system_name: str,
):
    checked = 0
    matched = 0

    for _, atom in log_atoms.iterrows():
        try:
            idx = int(atom["Index"])
        except (TypeError, ValueError):
            continue

        checked += 1
        pdb_atom = metadata.get(idx)
        if pdb_atom is None:
            continue

        if (
            str(atom.get("Residue", "")).strip().upper()
            == str(pdb_atom.get("Residue", "")).strip().upper()
        ):
            matched += 1

    fraction = matched / checked if checked else 0.0

    print(
        f"  PDB/log residue index agreement: "
        f"{matched:,}/{checked:,} ({100.0 * fraction:.2f}%)"
    )

    if fraction < 0.95:
        raise RuntimeError(
            f"PDB index mapping failed for {system_name}: "
            f"{100.0 * fraction:.2f}% residue agreement"
        )


# ---------------------------------------------------------------------------
# Solvent-exposure calculations
# ---------------------------------------------------------------------------

def attach_atom_exposure(
    atom_df: pd.DataFrame,
    aw_logs: Dict,
    metadata: Dict[int, Dict],
) -> pd.DataFrame:
    df = atom_df.copy()

    total_area = {int(i): 0.0 for i in df["Index"]}
    solvent_area = {int(i): 0.0 for i in df["Index"]}

    for b1, b2, area in surface_pairs(aw_logs["surfs"]):
        if b1 in total_area:
            total_area[b1] += area
            if metadata.get(b2, {}).get("Residue", "") in SOLVENT_RESIDUES:
                solvent_area[b1] += area

        if b2 in total_area:
            total_area[b2] += area
            if metadata.get(b1, {}).get("Residue", "") in SOLVENT_RESIDUES:
                solvent_area[b2] += area

    df["AW SolFacingPct"] = [
        (
            100.0 * solvent_area[int(idx)] / total_area[int(idx)]
            if total_area[int(idx)] > 0.0
            else np.nan
        )
        for idx in df["Index"]
    ]

    return df


def attach_residue_exposure(
    residue_df: pd.DataFrame,
    atom_df: pd.DataFrame,
    aw_logs: Dict,
    metadata: Dict[int, Dict],
) -> pd.DataFrame:
    df = residue_df.copy()

    idx_to_residue = {
        int(row["Index"]): (
            str(row["Chain"]),
            str(row["Residue"]),
            int(row["Residue Sequence"]),
        )
        for _, row in atom_df.iterrows()
    }

    residue_keys = {
        (
            str(row["Chain"]),
            str(row["Residue"]),
            int(row["Residue Sequence"]),
        )
        for _, row in df.iterrows()
    }

    total_area = {k: 0.0 for k in residue_keys}
    solvent_area = {k: 0.0 for k in residue_keys}

    for b1, b2, area in surface_pairs(aw_logs["surfs"]):
        r1 = idx_to_residue.get(b1)
        r2 = idx_to_residue.get(b2)

        if r1 is not None and r1 == r2:
            continue

        if r1 in total_area:
            total_area[r1] += area
            if metadata.get(b2, {}).get("Residue", "") in SOLVENT_RESIDUES:
                solvent_area[r1] += area

        if r2 in total_area:
            total_area[r2] += area
            if metadata.get(b1, {}).get("Residue", "") in SOLVENT_RESIDUES:
                solvent_area[r2] += area

    exposure_lookup = {
        key: (
            100.0 * solvent_area[key] / total_area[key]
            if total_area[key] > 0.0
            else np.nan
        )
        for key in residue_keys
    }

    df["AW SolFacingPct"] = [
        exposure_lookup.get(
            (
                str(row["Chain"]),
                str(row["Residue"]),
                int(row["Residue Sequence"]),
            ),
            np.nan,
        )
        for _, row in df.iterrows()
    ]

    return df


# ---------------------------------------------------------------------------
# Figure 3B atom group reuse
# ---------------------------------------------------------------------------

def add_atom_groups(
    atom_df: pd.DataFrame,
    system_class: str,
) -> pd.DataFrame:
    df = atom_df.copy()

    df["CanonicalName"] = [
        canonicalize_atom_name(
            atom_name=str(atom),
            molecule_class=system_class,
            residue_name=str(residue),
        )
        for atom, residue in zip(df["Atom"], df["Residue"])
    ]

    wanted = (
        PROTEIN_REPRESENTATIVE_GROUPS
        if system_class == "protein"
        else NUCLEIC_REPRESENTATIVE_GROUPS
    )

    return df[df["CanonicalName"].isin(wanted)].copy()


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_data(
    data_root: str,
    exclude_keys: Optional[Iterable[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    systems = discover_systems(
        data_root,
        exclude_keys=exclude_keys,
    )

    atom_frames = []
    residue_frames = []

    print(f"\nFound {len(systems)} valid Figure 3 systems.")

    for system in systems:
        print(f"\nProcessing {system.name} ...")

        aw_logs, power_logs = read_pair(
            system,
            need_surfs=True,
        )

        system_class = classify_system(aw_logs["atoms"])
        if system_class is None:
            print("  skipped: unclassified")
            continue

        pdb_path = expected_pdb_path(system)
        if not pdb_path.exists():
            raise FileNotFoundError(
                f"Full-system PDB not found: {pdb_path}"
            )

        metadata = parse_full_pdb_metadata(pdb_path)
        validate_pdb_mapping(
            metadata,
            aw_logs["atoms"],
            system.name,
        )

        allowed = allowed_residues(system_class)
        class_name = broad_class(system_class)

        # --------------------------
        # Atom data
        # --------------------------
        atom_df = build_atomic_metrics(
            aw_logs,
            power_logs,
        )

        atom_df = atom_df[
            atom_df["Residue"].isin(allowed)
        ].copy()

        atom_df = add_deviation_columns(atom_df)
        atom_df = attach_atom_exposure(
            atom_df,
            aw_logs,
            metadata,
        )
        atom_df = add_atom_groups(
            atom_df,
            system_class=system_class,
        )

        atom_df["System"] = system.name
        atom_df["BroadClass"] = class_name

        atom_frames.append(atom_df)

        # --------------------------
        # Residue data
        # --------------------------
        # Build residue metrics from all matched biomolecular atoms, not only
        # the Figure 3B representative subset.
        all_atom_df = build_atomic_metrics(
            aw_logs,
            power_logs,
        )
        all_atom_df = all_atom_df[
            all_atom_df["Residue"].isin(allowed)
        ].copy()

        residue_df = build_residue_metrics(
            all_atom_df,
            aw_logs,
            power_logs,
        )
        residue_df = residue_df[
            residue_df["Residue"].isin(allowed)
        ].copy()

        residue_df = add_deviation_columns(residue_df)
        residue_df = attach_residue_exposure(
            residue_df,
            all_atom_df,
            aw_logs,
            metadata,
        )

        residue_df["System"] = system.name
        residue_df["BroadClass"] = class_name

        residue_frames.append(residue_df)

        print(
            f"  retained selected atoms: {len(atom_df):,}; "
            f"residues: {len(residue_df):,}"
        )

    atoms = (
        pd.concat(atom_frames, ignore_index=True)
        if atom_frames
        else pd.DataFrame()
    )
    residues = (
        pd.concat(residue_frames, ignore_index=True)
        if residue_frames
        else pd.DataFrame()
    )

    return atoms, residues


# ---------------------------------------------------------------------------
# Trend statistics
# ---------------------------------------------------------------------------

def linear_stats(x, y) -> Tuple[float, float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    if len(x) < 3 or np.ptp(x) <= 0.0:
        return np.nan, np.nan, np.nan, len(x)

    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept

    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = (
        1.0 - ss_res / ss_tot
        if ss_tot > 0.0
        else np.nan
    )

    return float(slope), float(intercept), float(r2), len(x)


def quantile_binned_means(
    x,
    y,
    n_bins: int = TREND_BINS,
) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "x": pd.to_numeric(x, errors="coerce"),
            "y": pd.to_numeric(y, errors="coerce"),
        }
    ).dropna()

    if len(df) < MIN_GROUP_COUNT:
        return pd.DataFrame()

    # qcut can collapse repeated exposure values; duplicates='drop' handles it.
    try:
        df["bin"] = pd.qcut(
            df["x"],
            q=n_bins,
            duplicates="drop",
        )
    except ValueError:
        return pd.DataFrame()

    out = (
        df.groupby("bin", observed=True)
        .agg(
            x_mean=("x", "mean"),
            y_mean=("y", "mean"),
            y_sem=(
                "y",
                lambda s: (
                    s.std(ddof=1) / np.sqrt(len(s))
                    if len(s) > 1
                    else 0.0
                ),
            ),
            count=("y", "size"),
        )
        .reset_index(drop=True)
    )

    return out[out["count"] >= MIN_BIN_COUNT].copy()


def build_atom_trend_table(atom_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for metric in METRICS:
        y_col = f"{metric} Abs % Diff"

        for (class_name, group_name), group_df in atom_df.groupby(
            ["BroadClass", "CanonicalName"]
        ):
            slope, intercept, r2, n = linear_stats(
                group_df["AW SolFacingPct"],
                group_df[y_col],
            )

            rows.append(
                {
                    "Metric": metric,
                    "BroadClass": class_name,
                    "GroupName": group_name,
                    "Count": n,
                    "Mean Exposure %": float(
                        pd.to_numeric(
                            group_df["AW SolFacingPct"],
                            errors="coerce",
                        ).mean()
                    ),
                    "Mean Abs % Diff": float(
                        pd.to_numeric(
                            group_df[y_col],
                            errors="coerce",
                        ).mean()
                    ),
                    "Slope (%diff per exposure %)": slope,
                    "Intercept": intercept,
                    "R2": r2,
                }
            )

    return pd.DataFrame(rows)


def add_residue_fixed_effects(
    residue_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Residualize exposure and each metric by residue identity within broad class.
    """
    df = residue_df.copy()

    group_cols = ["BroadClass", "Residue"]

    df["Exposure Within Type"] = (
        pd.to_numeric(
            df["AW SolFacingPct"],
            errors="coerce",
        )
        - df.groupby(group_cols)["AW SolFacingPct"].transform("mean")
    )

    for metric in METRICS:
        y_col = f"{metric} Abs % Diff"

        numeric_y = pd.to_numeric(
            df[y_col],
            errors="coerce",
        )

        df[f"{metric} Within Type"] = (
            numeric_y
            - df.assign(_y=numeric_y)
                .groupby(group_cols)["_y"]
                .transform("mean")
        )

    return df


def build_residue_trend_table(
    residue_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for metric in METRICS:
        y_col = f"{metric} Within Type"

        for class_name, class_df in residue_df.groupby("BroadClass"):
            slope, intercept, r2, n = linear_stats(
                class_df["Exposure Within Type"],
                class_df[y_col],
            )

            rows.append(
                {
                    "Metric": metric,
                    "BroadClass": class_name,
                    "Count": n,
                    "Slope": slope,
                    "Intercept": intercept,
                    "R2": r2,
                }
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def display_upper(values, percentile=DISPLAY_PERCENTILE):
    vals = pd.to_numeric(
        pd.Series(values),
        errors="coerce",
    ).dropna().to_numpy(float)

    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return 1.0

    return max(
        float(np.percentile(vals, percentile)) * 1.10,
        1.0,
    )


def symmetric_limit(values, percentile=DISPLAY_PERCENTILE):
    vals = pd.to_numeric(
        pd.Series(values),
        errors="coerce",
    ).dropna().to_numpy(float)

    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return 1.0

    abs_lim = float(
        np.percentile(np.abs(vals), percentile)
    )

    return max(abs_lim * 1.10, 0.5)


def plot_atom_panel(
    ax,
    atom_df: pd.DataFrame,
    metric: str,
    rng: np.random.Generator,
):
    y_col = f"{metric} Abs % Diff"

    panel = atom_df[
        ["AW SolFacingPct", y_col, "BroadClass", "CanonicalName"]
    ].copy()

    panel["AW SolFacingPct"] = pd.to_numeric(
        panel["AW SolFacingPct"],
        errors="coerce",
    )
    panel[y_col] = pd.to_numeric(
        panel[y_col],
        errors="coerce",
    )
    panel = panel.dropna()

    groups_in_order = (
        PROTEIN_REPRESENTATIVE_GROUPS
        + NUCLEIC_REPRESENTATIVE_GROUPS
    )

    for group_name in groups_in_order:
        group_df = panel[
            panel["CanonicalName"] == group_name
        ]

        if len(group_df) < MIN_GROUP_COUNT:
            continue

        color = ATOM_GROUP_COLORS.get(group_name, "0.4")

        # Plot sample only.
        if ATOM_POINT_FRACTION < 1.0:
            n_sample = max(
                1,
                int(round(len(group_df) * ATOM_POINT_FRACTION)),
            )
            sample_idx = rng.choice(
                group_df.index.to_numpy(),
                size=min(n_sample, len(group_df)),
                replace=False,
            )
            sample = group_df.loc[sample_idx]
        else:
            sample = group_df

        ax.scatter(
            sample["AW SolFacingPct"],
            sample[y_col],
            s=ATOM_POINT_SIZE,
            alpha=ATOM_POINT_ALPHA,
            color=color,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

        trend = quantile_binned_means(
            group_df["AW SolFacingPct"],
            group_df[y_col],
        )

        if not trend.empty:
            ax.plot(
                trend["x_mean"],
                trend["y_mean"],
                marker="o",
                markersize=TREND_MARKER_SIZE,
                linewidth=TREND_LINE_WIDTH,
                color=color,
                label=group_name,
                zorder=4,
            )

    ax.set_xlabel(
        "AW solvent-facing surface area (%)",
        fontsize=10,
    )
    ax.set_ylabel(
        "Absolute Power vs AW difference (%)",
        fontsize=10,
    )
    ax.set_title(
        metric,
        fontsize=13,
        fontweight="bold",
    )

    ax.set_ylim(
        0.0,
        display_upper(panel[y_col]),
    )

    x_max = display_upper(
        panel["AW SolFacingPct"],
        percentile=99.5,
    )
    ax.set_xlim(0.0, x_max)

    ax.tick_params(
        axis="both",
        labelsize=9,
        width=1.2,
        length=5,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def plot_residue_panel(
    ax,
    residue_df: pd.DataFrame,
    metric: str,
):
    x_col = "Exposure Within Type"
    y_col = f"{metric} Within Type"

    panel = residue_df[
        [x_col, y_col, "BroadClass"]
    ].copy()

    panel[x_col] = pd.to_numeric(
        panel[x_col],
        errors="coerce",
    )
    panel[y_col] = pd.to_numeric(
        panel[y_col],
        errors="coerce",
    )
    panel = panel.dropna()

    for class_name in ["Protein", "Nucleic acid"]:
        class_df = panel[
            panel["BroadClass"] == class_name
        ]

        if class_df.empty:
            continue

        color = CLASS_COLORS[class_name]

        ax.scatter(
            class_df[x_col],
            class_df[y_col],
            s=RESIDUE_POINT_SIZE,
            alpha=RESIDUE_POINT_ALPHA,
            color=color,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

        trend = quantile_binned_means(
            class_df[x_col],
            class_df[y_col],
            n_bins=TREND_BINS,
        )

        if not trend.empty:
            ax.plot(
                trend["x_mean"],
                trend["y_mean"],
                marker="o",
                markersize=TREND_MARKER_SIZE + 0.5,
                linewidth=TREND_LINE_WIDTH,
                color=color,
                label=class_name,
                zorder=4,
            )

    ax.axhline(
        0.0,
        linestyle="--",
        linewidth=1.0,
        color="black",
        alpha=0.7,
        zorder=0,
    )
    ax.axvline(
        0.0,
        linestyle="--",
        linewidth=1.0,
        color="black",
        alpha=0.7,
        zorder=0,
    )

    xlim = symmetric_limit(panel[x_col], percentile=99.0)
    ylim = symmetric_limit(panel[y_col], percentile=99.0)

    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-ylim, ylim)

    ax.set_xlabel(
        "Exposure relative to residue-type mean (%)",
        fontsize=10,
    )
    ax.set_ylabel(
        "Deviation relative to residue-type mean (%)",
        fontsize=10,
    )
    ax.set_title(
        metric,
        fontsize=13,
        fontweight="bold",
    )

    ax.tick_params(
        axis="both",
        labelsize=9,
        width=1.2,
        length=5,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def add_legends(fig):
    atom_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linewidth=2,
            color=ATOM_GROUP_COLORS.get(name, "0.4"),
            label=name,
            markersize=5,
        )
        for name in (
            PROTEIN_REPRESENTATIVE_GROUPS
            + NUCLEIC_REPRESENTATIVE_GROUPS
        )
    ]

    residue_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linewidth=2,
            color=color,
            label=name,
            markersize=5,
        )
        for name, color in CLASS_COLORS.items()
    ]

    atom_legend = fig.legend(
        handles=atom_handles,
        title="Atom groups",
        loc="lower left",
        bbox_to_anchor=(0.08, 0.005),
        ncol=6,
        frameon=False,
        fontsize=8.5,
        title_fontsize=9,
    )
    fig.add_artist(atom_legend)

    fig.legend(
        handles=residue_handles,
        title="Residue class",
        loc="lower right",
        bbox_to_anchor=(0.92, 0.005),
        ncol=2,
        frameon=False,
        fontsize=9,
        title_fontsize=9,
    )


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(
    atom_df: pd.DataFrame,
    residue_df: pd.DataFrame,
    figure_dir: str,
):
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    residue_fe = add_residue_fixed_effects(
        residue_df
    )

    atom_trends = build_atom_trend_table(
        atom_df
    )
    residue_trends = build_residue_trend_table(
        residue_fe
    )

    print()
    print("=" * 115)
    print("F3E ATOM GROUP EXPOSURE TRENDS")
    print("=" * 115)
    print(
        atom_trends.sort_values(
            ["Metric", "BroadClass", "GroupName"]
        ).to_string(
            index=False,
            formatters={
                "Mean Exposure %": lambda x: f"{x:8.3f}",
                "Mean Abs % Diff": lambda x: f"{x:8.3f}",
                "Slope (%diff per exposure %)": lambda x: f"{x:10.5f}",
                "Intercept": lambda x: f"{x:9.4f}",
                "R2": lambda x: f"{x:8.4f}",
            },
        )
    )

    print()
    print("=" * 90)
    print("F3E RESIDUE WITHIN-TYPE EXPOSURE TRENDS")
    print("=" * 90)
    print(
        residue_trends.to_string(
            index=False,
            formatters={
                "Slope": lambda x: f"{x:10.5f}",
                "Intercept": lambda x: f"{x:9.4f}",
                "R2": lambda x: f"{x:8.4f}",
            },
        )
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

    for col_i, metric in enumerate(METRICS):
        plot_atom_panel(
            axes[0, col_i],
            atom_df=atom_df,
            metric=metric,
            rng=rng,
        )

        plot_residue_panel(
            axes[1, col_i],
            residue_df=residue_fe,
            metric=metric,
        )

    axes[0, 0].text(
        -0.22,
        0.5,
        "Atom",
        transform=axes[0, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=15,
        fontweight="bold",
    )

    axes[1, 0].text(
        -0.22,
        0.5,
        "Residue\n(type-controlled)",
        transform=axes[1, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=14,
        fontweight="bold",
    )

    fig.suptitle(
        "Figure 3E — AW/Power sensitivity as a function of solvent exposure",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )

    add_legends(fig)

    fig.subplots_adjust(
        left=0.08,
        right=0.99,
        top=0.92,
        bottom=0.16,
        wspace=0.26,
        hspace=0.30,
    )

    png_path = (
        figure_dir /
        "F3E_continuous_solvent_exposure.png"
    )
    svg_path = (
        figure_dir /
        "F3E_continuous_solvent_exposure.svg"
    )

    if SAVE_PNG:
        fig.savefig(
            png_path,
            dpi=DPI,
            bbox_inches="tight",
        )
        print(f"Saved: {png_path}")

    if SAVE_SVG:
        fig.savefig(
            svg_path,
            bbox_inches="tight",
        )
        print(f"Saved: {svg_path}")

    atom_csv = (
        figure_dir /
        "F3E_atom_group_trends.csv"
    )
    residue_csv = (
        figure_dir /
        "F3E_residue_fixed_effect_trends.csv"
    )

    atom_trends.to_csv(
        atom_csv,
        index=False,
    )
    residue_trends.to_csv(
        residue_csv,
        index=False,
    )

    print(f"Saved: {atom_csv}")
    print(f"Saved: {residue_csv}")

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

    atom_df, residue_df = collect_data(
        data_root=data_root,
        exclude_keys=EXCLUDE_KEYS,
    )

    if atom_df.empty or residue_df.empty:
        print(
            "Continuous 3E requires both atom and residue data."
        )
        return

    make_figure(
        atom_df=atom_df,
        residue_df=residue_df,
        figure_dir=figure_dir,
    )


if __name__ == "__main__":
    main()