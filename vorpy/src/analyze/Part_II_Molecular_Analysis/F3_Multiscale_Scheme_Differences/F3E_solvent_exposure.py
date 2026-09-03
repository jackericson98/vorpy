"""Figure 3E — Direct solvent-exposure dependence of AW/Power deviations.

Layout:
                    Volume | Surface Area | Contacts
    Atom
    Residue

Each panel shows:
    x = AW solvent-facing surface area (%)
    y = absolute Power-vs-AW percent deviation

Protein and nucleic-acid observations are shown separately with faint points
and binned mean trends. Chemistry is controlled in exported statistics rather
than encoded into the visualization.
"""
from __future__ import annotations

import os
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog
from typing import Dict, Iterable, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

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

DATA_ROOT = None
FIGURE_DIR = None
EXCLUDE_KEYS = ["A", "B", "C"]

# ---------------------------------------------------------------------------
# Cache settings
# ---------------------------------------------------------------------------
# Building the exposure tables is expensive. Once built, all later figure
# experiments should read these CSVs instead of re-reading logs/PDBs.
REBUILD_CACHE = False
ATOM_CACHE_NAME = "F3E_atom_exposure_cache.csv"
RESIDUE_CACHE_NAME = "F3E_residue_exposure_cache.csv"


FIGSIZE = (18.5, 10.5)
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

METRICS = ["Volume", "Surface Area"]

CLASS_COLORS = {
    "Protein": "#4c78a8",
    "Nucleic acid": "#f58518",
}

POINT_SIZE = 8
POINT_ALPHA = 0.055
ATOM_PLOT_FRACTION = 1.0
RESIDUE_PLOT_FRACTION = 1.0
RANDOM_SEED = 17

# Hard caps keep the raw visual layer informative without allowing the
# very large protein populations to dominate the figure.
ATOM_MAX_POINTS_PER_GROUP = 500
RESIDUE_MAX_POINTS_PER_GROUP = 500

# Exposure bins are data-adaptive within each scale (Atom vs Residue).
# Protein and nucleic acids share the same cutoffs within a scale.
EXPOSURE_QUANTILES = (1.0 / 3.0, 2.0 / 3.0)
EXPOSURE_ORDER = ["Low", "Medium", "High"]

# Compact summary styling.
CLASS_OFFSET = 0.18
RAW_JITTER = 0.07
SUMMARY_MARKER_SIZE = 82
IQR_LINE_WIDTH = 3.8
MEDIAN_LINE_WIDTH = 2.0

DISPLAY_PERCENTILE = 99.0

# Figure 3E focuses on continuous geometric quantities.
# Contacts are intentionally excluded from this exposure panel.


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
    system_folder = Path(system.aw).resolve().parent.parent
    molecule_name = getattr(system, "molecule_name", None)
    if not molecule_name:
        folder_name = system_folder.name
        molecule_name = folder_name.split("_", 1)[1] if "_" in folder_name else folder_name

    expected = system_folder / f"{molecule_name}.pdb"
    if expected.exists():
        return expected

    target = expected.name.lower()
    for candidate in system_folder.iterdir():
        if candidate.is_file() and candidate.name.lower() == target:
            return candidate
    return expected


def parse_full_pdb_metadata(pdb_path: Path) -> Dict[int, Dict]:
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


def validate_pdb_mapping(metadata: Dict[int, Dict], log_atoms: pd.DataFrame, system_name: str):
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


def attach_atom_exposure(atom_df: pd.DataFrame, aw_logs: Dict, metadata: Dict[int, Dict]) -> pd.DataFrame:
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
        100.0 * solvent_area[int(idx)] / total_area[int(idx)]
        if total_area[int(idx)] > 0.0 else np.nan
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
        key: 100.0 * solvent_area[key] / total_area[key]
        if total_area[key] > 0.0 else np.nan
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


def collect_data(data_root: str, exclude_keys: Optional[Iterable[str]] = None):
    systems = discover_systems(data_root, exclude_keys=exclude_keys)
    atom_frames = []
    residue_frames = []

    print(f"\nFound {len(systems)} valid Figure 3 systems.")

    for system in systems:
        print(f"\nProcessing {system.name} ...")

        aw_logs, power_logs = read_pair(system, need_surfs=True)
        system_class = classify_system(aw_logs["atoms"])

        if system_class is None:
            print("  skipped: unclassified")
            continue

        pdb_path = expected_pdb_path(system)
        if not pdb_path.exists():
            raise FileNotFoundError(f"Full-system PDB not found: {pdb_path}")

        metadata = parse_full_pdb_metadata(pdb_path)
        validate_pdb_mapping(metadata, aw_logs["atoms"], system.name)

        allowed = allowed_residues(system_class)
        class_name = broad_class(system_class)

        atom_df = build_atomic_metrics(aw_logs, power_logs)
        atom_df = atom_df[atom_df["Residue"].isin(allowed)].copy()
        atom_df = add_deviation_columns(atom_df)
        atom_df = attach_atom_exposure(atom_df, aw_logs, metadata)
        atom_df["System"] = system.name
        atom_df["BroadClass"] = class_name
        atom_df["CanonicalName"] = [
            canonicalize_atom_name(
                atom_name=str(atom),
                molecule_class=system_class,
                residue_name=str(residue),
            )
            for atom, residue in zip(atom_df["Atom"], atom_df["Residue"])
        ]
        atom_frames.append(atom_df)

        residue_df = build_residue_metrics(atom_df, aw_logs, power_logs)
        residue_df = residue_df[residue_df["Residue"].isin(allowed)].copy()
        residue_df = add_deviation_columns(residue_df)
        residue_df = attach_residue_exposure(residue_df, atom_df, aw_logs, metadata)
        residue_df["System"] = system.name
        residue_df["BroadClass"] = class_name
        residue_frames.append(residue_df)

        print(f"  retained atoms: {len(atom_df):,}; residues: {len(residue_df):,}")

    atoms = pd.concat(atom_frames, ignore_index=True) if atom_frames else pd.DataFrame()
    residues = pd.concat(residue_frames, ignore_index=True) if residue_frames else pd.DataFrame()
    return atoms, residues


def class_exposure_cutoffs(df: pd.DataFrame) -> Dict[str, tuple]:
    """Return Low/Medium/High tertile cutoffs separately for each molecular class."""
    cutoffs = {}

    for broad_name in ["Protein", "Nucleic acid"]:
        vals = pd.to_numeric(
            df.loc[df["BroadClass"] == broad_name, "AW SolFacingPct"],
            errors="coerce",
        ).dropna().to_numpy(float)

        vals = vals[np.isfinite(vals)]

        if len(vals) == 0:
            cutoffs[broad_name] = (np.nan, np.nan)
            continue

        low, high = np.quantile(vals, EXPOSURE_QUANTILES)
        cutoffs[broad_name] = (float(low), float(high))

    return cutoffs


def assign_exposure_categories(
    df: pd.DataFrame,
    cutoffs: Dict[str, tuple],
) -> pd.DataFrame:
    """
    Assign exposure tertiles within molecular class.

    Thus Low/Medium/High means relatively buried/intermediate/exposed within
    Protein or Nucleic acid, rather than using protein-dominated global cutoffs.
    """
    out = df.copy()
    out["Exposure Group"] = None

    exposure = pd.to_numeric(
        out["AW SolFacingPct"],
        errors="coerce",
    )

    for broad_name, (low_cut, high_cut) in cutoffs.items():
        if not np.isfinite(low_cut) or not np.isfinite(high_cut):
            continue

        class_mask = out["BroadClass"] == broad_name

        out.loc[
            class_mask & (exposure <= low_cut),
            "Exposure Group",
        ] = "Low"

        out.loc[
            class_mask
            & (exposure > low_cut)
            & (exposure <= high_cut),
            "Exposure Group",
        ] = "Medium"

        out.loc[
            class_mask & (exposure > high_cut),
            "Exposure Group",
        ] = "High"

    return out


def summarize_exposure_groups(
    df: pd.DataFrame,
    scale: str,
):
    cutoffs = class_exposure_cutoffs(df)
    work = assign_exposure_categories(df, cutoffs=cutoffs)

    rows = []

    for metric in METRICS:
        y_col = f"{metric} Abs % Diff"

        for broad_name in ["Protein", "Nucleic acid"]:
            low_cut, high_cut = cutoffs[broad_name]

            for exposure_group in EXPOSURE_ORDER:
                vals = pd.to_numeric(
                    work.loc[
                        (work["BroadClass"] == broad_name)
                        & (work["Exposure Group"] == exposure_group),
                        y_col,
                    ],
                    errors="coerce",
                ).dropna().to_numpy(float)

                vals = vals[np.isfinite(vals)]

                if len(vals) == 0:
                    continue

                rows.append({
                    "Scale": scale,
                    "Metric": metric,
                    "BroadClass": broad_name,
                    "Exposure Group": exposure_group,
                    "Count": len(vals),
                    "Exposure Low Cutoff": low_cut,
                    "Exposure High Cutoff": high_cut,
                    "Mean": float(np.mean(vals)),
                    "Median": float(np.median(vals)),
                    "Q1": float(np.percentile(vals, 25)),
                    "Q3": float(np.percentile(vals, 75)),
                    "SD": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                })

    return work, pd.DataFrame(rows), cutoffs


def linear_stats(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(float)

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    if len(x) < 3 or np.ptp(x) <= 0.0:
        return np.nan, np.nan, np.nan, len(x)

    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return float(slope), float(intercept), float(r2), len(x)


def controlled_slopes_atom(atom_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        y_col = f"{metric} Abs % Diff"
        for (broad_name, atom_type), group_df in atom_df.groupby(["BroadClass", "CanonicalName"]):
            slope, intercept, r2, n = linear_stats(
                group_df["AW SolFacingPct"],
                group_df[y_col],
            )
            rows.append({
                "Metric": metric,
                "BroadClass": broad_name,
                "AtomType": atom_type,
                "Count": n,
                "Slope": slope,
                "Intercept": intercept,
                "R2": r2,
            })
    return pd.DataFrame(rows)


def controlled_slopes_residue(residue_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        y_col = f"{metric} Abs % Diff"
        for (broad_name, residue_name), group_df in residue_df.groupby(["BroadClass", "Residue"]):
            slope, intercept, r2, n = linear_stats(
                group_df["AW SolFacingPct"],
                group_df[y_col],
            )
            rows.append({
                "Metric": metric,
                "BroadClass": broad_name,
                "Residue": residue_name,
                "Count": n,
                "Slope": slope,
                "Intercept": intercept,
                "R2": r2,
            })
    return pd.DataFrame(rows)


def overall_summary(atom_df: pd.DataFrame, residue_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scale, df in [("Atom", atom_df), ("Residue", residue_df)]:
        for metric in METRICS:
            y_col = f"{metric} Abs % Diff"
            for broad_name, group_df in df.groupby("BroadClass"):
                slope, intercept, r2, n = linear_stats(
                    group_df["AW SolFacingPct"],
                    group_df[y_col],
                )
                rows.append({
                    "Scale": scale,
                    "Metric": metric,
                    "BroadClass": broad_name,
                    "Count": n,
                    "Mean Exposure %": pd.to_numeric(
                        group_df["AW SolFacingPct"], errors="coerce"
                    ).mean(),
                    "Mean Abs % Diff": pd.to_numeric(
                        group_df[y_col], errors="coerce"
                    ).mean(),
                    "Slope": slope,
                    "Intercept": intercept,
                    "R2": r2,
                })
    return pd.DataFrame(rows)


def sample_for_plot(
    df: pd.DataFrame,
    fraction: float,
    rng: np.random.Generator,
    max_points: Optional[int] = None,
) -> pd.DataFrame:
    """Downsample only the faint raw-point layer; summaries always use all data."""
    if len(df) <= 1:
        return df

    n = len(df)

    if fraction < 1.0:
        n = max(1, int(round(n * fraction)))

    if max_points is not None:
        n = min(n, max_points)

    if n >= len(df):
        return df

    idx = rng.choice(
        df.index.to_numpy(),
        size=n,
        replace=False,
    )
    return df.loc[idx]


def axis_upper(values, percentile=DISPLAY_PERCENTILE):
    vals = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return 1.0

    return max(float(np.percentile(vals, percentile)) * 1.08, 1.0)


def plot_panel(
    ax,
    df,
    scale,
    metric,
    rng,
    cutoffs,
):
    y_col = f"{metric} Abs % Diff"
    plot_fraction = (
        ATOM_PLOT_FRACTION
        if scale == "Atom"
        else RESIDUE_PLOT_FRACTION
    )

    work = assign_exposure_categories(
        df,
        cutoffs=cutoffs,
    )

    base_x = {
        "Low": 0.0,
        "Medium": 1.0,
        "High": 2.0,
    }

    offsets = {
        "Protein": -CLASS_OFFSET,
        "Nucleic acid": CLASS_OFFSET,
    }

    for broad_name in ["Protein", "Nucleic acid"]:
        color = CLASS_COLORS[broad_name]

        for exposure_group in EXPOSURE_ORDER:
            group_df = work[
                (work["BroadClass"] == broad_name)
                & (work["Exposure Group"] == exposure_group)
            ].copy()

            group_df[y_col] = pd.to_numeric(
                group_df[y_col],
                errors="coerce",
            )
            group_df = group_df.dropna(subset=[y_col])

            if group_df.empty:
                continue

            # Downsample only the raw visual layer.
            max_points = (
                ATOM_MAX_POINTS_PER_GROUP
                if scale == "Atom"
                else RESIDUE_MAX_POINTS_PER_GROUP
            )

            sample = sample_for_plot(
                group_df,
                plot_fraction,
                rng,
                max_points=max_points,
            )

            x_center = base_x[exposure_group] + offsets[broad_name]
            jitter = rng.uniform(
                -RAW_JITTER,
                RAW_JITTER,
                len(sample),
            )

            ax.scatter(
                np.full(len(sample), x_center) + jitter,
                sample[y_col],
                s=POINT_SIZE,
                alpha=POINT_ALPHA,
                color=color,
                linewidths=0,
                rasterized=True,
                zorder=1,
            )

            vals = group_df[y_col].to_numpy(float)

            median = float(np.median(vals))
            q1 = float(np.percentile(vals, 25))
            q3 = float(np.percentile(vals, 75))

            # Thick IQR bar.
            ax.vlines(
                x_center,
                q1,
                q3,
                color=color,
                linewidth=IQR_LINE_WIDTH,
                zorder=4,
            )

            # Median marker.
            ax.scatter(
                [x_center],
                [median],
                s=SUMMARY_MARKER_SIZE,
                color=color,
                edgecolor="white",
                linewidth=1.2,
                zorder=5,
            )

            # Short median tick reinforces the statistic.
            ax.hlines(
                median,
                x_center - 0.055,
                x_center + 0.055,
                color="white",
                linewidth=MEDIAN_LINE_WIDTH,
                zorder=6,
            )

            # Count above the panel.
            ax.text(
                x_center,
                0.98,
                f"n={len(vals)}",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=8,
                color=color,
            )

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(EXPOSURE_ORDER, fontsize=10)

    ax.set_xlim(-0.55, 2.55)

    panel_y = pd.to_numeric(
        work[y_col],
        errors="coerce",
    )
    ax.set_ylim(
        0.0,
        axis_upper(panel_y, percentile=99.0),
    )

    ax.set_xlabel(
        "Relative AW solvent exposure",
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

    ax.tick_params(
        axis="both",
        labelsize=9,
        width=1.2,
        length=5,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def add_legend(fig):
    handles = [
        Line2D(
            [0], [0],
            marker="o",
            linewidth=2.2,
            color=color,
            label=name,
            markersize=6,
        )
        for name, color in CLASS_COLORS.items()
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=10,
        bbox_to_anchor=(0.5, 0.02),
    )


def cache_paths(figure_dir):
    figure_dir = Path(figure_dir)
    return (
        figure_dir / ATOM_CACHE_NAME,
        figure_dir / RESIDUE_CACHE_NAME,
    )


def save_exposure_cache(atom_df, residue_df, figure_dir):
    """Save only the scalar fields needed for future Figure 3E experiments."""
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    atom_cache, residue_cache = cache_paths(figure_dir)

    atom_columns = [
        "System",
        "BroadClass",
        "Index",
        "Name",
        "Atom",
        "Residue",
        "Residue Sequence",
        "Chain",
        "CanonicalName",
        "AW SolFacingPct",
        "Volume Abs % Diff",
        "Surface Area Abs % Diff",
        "Contacts Abs % Diff",
        "Volume Signed % Diff",
        "Surface Area Signed % Diff",
        "Contacts Signed % Diff",
    ]
    atom_columns = [c for c in atom_columns if c in atom_df.columns]

    residue_columns = [
        "System",
        "BroadClass",
        "Residue",
        "Residue Sequence",
        "Chain",
        "AW SolFacingPct",
        "Volume Abs % Diff",
        "Surface Area Abs % Diff",
        "Contacts Abs % Diff",
        "Volume Signed % Diff",
        "Surface Area Signed % Diff",
        "Contacts Signed % Diff",
    ]
    residue_columns = [c for c in residue_columns if c in residue_df.columns]

    atom_df[atom_columns].to_csv(atom_cache, index=False)
    residue_df[residue_columns].to_csv(residue_cache, index=False)

    print(f"Saved atom exposure cache: {atom_cache}")
    print(f"Saved residue exposure cache: {residue_cache}")


def load_exposure_cache(figure_dir):
    atom_cache, residue_cache = cache_paths(figure_dir)

    if not atom_cache.exists() or not residue_cache.exists():
        return None, None

    print(f"Loading cached atom data: {atom_cache}")
    print(f"Loading cached residue data: {residue_cache}")

    return (
        pd.read_csv(atom_cache),
        pd.read_csv(residue_cache),
    )


def get_or_build_exposure_data(data_root, figure_dir, exclude_keys=None):
    """
    Load cached exposure/deviation tables when available.

    Set REBUILD_CACHE = True only when the underlying AW/Power logs, PDBs,
    parsing logic, or solvent-exposure calculation has changed.
    """
    if not REBUILD_CACHE:
        atom_df, residue_df = load_exposure_cache(figure_dir)

        if atom_df is not None and residue_df is not None:
            print(
                f"Using cached Figure 3E data: "
                f"{len(atom_df):,} atoms, {len(residue_df):,} residues"
            )
            return atom_df, residue_df

    print("Building Figure 3E exposure cache from logs/PDBs ...")

    atom_df, residue_df = collect_data(
        data_root=data_root,
        exclude_keys=exclude_keys,
    )

    save_exposure_cache(
        atom_df,
        residue_df,
        figure_dir=figure_dir,
    )

    return atom_df, residue_df


def build_cutoff_table(atom_cutoffs, residue_cutoffs):
    rows = []

    for scale, cutoffs in [
        ("Atom", atom_cutoffs),
        ("Residue", residue_cutoffs),
    ]:
        for broad_name in ["Protein", "Nucleic acid"]:
            low_cut, high_cut = cutoffs[broad_name]

            rows.extend([
                {
                    "Scale": scale,
                    "BroadClass": broad_name,
                    "Exposure Group": "Low",
                    "Range": f"<= {low_cut:.3f}%",
                    "Lower Bound %": 0.0,
                    "Upper Bound %": low_cut,
                },
                {
                    "Scale": scale,
                    "BroadClass": broad_name,
                    "Exposure Group": "Medium",
                    "Range": f"> {low_cut:.3f}% to <= {high_cut:.3f}%",
                    "Lower Bound %": low_cut,
                    "Upper Bound %": high_cut,
                },
                {
                    "Scale": scale,
                    "BroadClass": broad_name,
                    "Exposure Group": "High",
                    "Range": f"> {high_cut:.3f}%",
                    "Lower Bound %": high_cut,
                    "Upper Bound %": np.nan,
                },
            ])

    return pd.DataFrame(rows)


def format_cutoff_line(scale_name: str, cutoffs: Dict[str, tuple]) -> str:
    parts = []

    for broad_name in ["Protein", "Nucleic acid"]:
        low_cut, high_cut = cutoffs[broad_name]

        parts.append(
            f"{broad_name}: "
            f"Low <= {low_cut:.1f}% | "
            f"Medium {low_cut:.1f}-{high_cut:.1f}% | "
            f"High > {high_cut:.1f}%"
        )

    return f"{scale_name}: " + "     ".join(parts)


def make_figure(atom_df, residue_df, figure_dir):
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    atom_work, atom_summary, atom_cutoffs = summarize_exposure_groups(
        atom_df,
        scale="Atom",
    )
    residue_work, residue_summary, residue_cutoffs = summarize_exposure_groups(
        residue_df,
        scale="Residue",
    )

    summary = pd.concat(
        [atom_summary, residue_summary],
        ignore_index=True,
    )

    atom_controlled = controlled_slopes_atom(atom_df)
    residue_controlled = controlled_slopes_residue(residue_df)
    cutoff_table = build_cutoff_table(
        atom_cutoffs,
        residue_cutoffs,
    )

    print()
    print("=" * 110)
    print("F3E CLASS-SPECIFIC LOW / MEDIUM / HIGH EXPOSURE SUMMARY")
    print("=" * 110)

    for scale_name, cutoffs in [
        ("Atom", atom_cutoffs),
        ("Residue", residue_cutoffs),
    ]:
        for broad_name in ["Protein", "Nucleic acid"]:
            low_cut, high_cut = cutoffs[broad_name]
            print(
                f"{scale_name:7s} | {broad_name:12s}: "
                f"Low <= {low_cut:.3f}% | "
                f"Medium <= {high_cut:.3f}% | "
                f"High > {high_cut:.3f}%"
            )

    print()
    print(summary.to_string(index=False))

    fig, axes = plt.subplots(
        2,
        len(METRICS),
        figsize=(12.5 if len(METRICS) == 2 else 18.5, 10.5),
        constrained_layout=False,
        squeeze=False,
    )

    rng = np.random.default_rng(RANDOM_SEED)

    for col_i, metric in enumerate(METRICS):
        plot_panel(
            axes[0, col_i],
            atom_work,
            "Atom",
            metric,
            rng,
            cutoffs=atom_cutoffs,
        )
        plot_panel(
            axes[1, col_i],
            residue_work,
            "Residue",
            metric,
            rng,
            cutoffs=residue_cutoffs,
        )

    axes[0, 0].text(
        -0.20, 0.5, "Atom",
        transform=axes[0, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=15,
        fontweight="bold",
    )
    axes[1, 0].text(
        -0.20, 0.5, "Residue",
        transform=axes[1, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=15,
        fontweight="bold",
    )

    fig.suptitle(
        "Figure 3E — Scheme sensitivity across relative solvent exposure",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )

    fig.text(
        0.5,
        0.948,
        "Low / Medium / High are class-specific AW solvent-exposure tertiles; "
        "markers show medians and bars show IQR",
        ha="center",
        va="center",
        fontsize=9.5,
    )

    fig.text(
        0.5,
        0.925,
        format_cutoff_line("Atom", atom_cutoffs),
        ha="center",
        va="center",
        fontsize=8.7,
    )

    fig.text(
        0.5,
        0.905,
        format_cutoff_line("Residue", residue_cutoffs),
        ha="center",
        va="center",
        fontsize=8.7,
    )

    add_legend(fig)

    fig.subplots_adjust(
        left=0.08,
        right=0.99,
        top=0.865,
        bottom=0.10,
        wspace=0.26,
        hspace=0.30,
    )

    png_path = figure_dir / "F3E_exposure_multiscale_ranges.png"
    svg_path = figure_dir / "F3E_exposure_multiscale_ranges.svg"

    if SAVE_PNG:
        fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {png_path}")

    if SAVE_SVG:
        fig.savefig(svg_path, bbox_inches="tight")
        print(f"Saved: {svg_path}")

    summary.to_csv(
        figure_dir / "F3E_exposure_multiscale_summary.csv",
        index=False,
    )
    atom_controlled.to_csv(
        figure_dir / "F3E_atom_type_controlled_slopes.csv",
        index=False,
    )
    residue_controlled.to_csv(
        figure_dir / "F3E_residue_type_controlled_slopes.csv",
        index=False,
    )
    cutoff_table.to_csv(
        figure_dir / "F3E_exposure_group_cutoffs.csv",
        index=False,
    )

    if SHOW:
        plt.show()

    plt.close(fig)


def main(data_root: Optional[str] = DATA_ROOT, figure_dir: Optional[str] = FIGURE_DIR):
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

    atom_df, residue_df = get_or_build_exposure_data(
        data_root=data_root,
        figure_dir=figure_dir,
        exclude_keys=EXCLUDE_KEYS,
    )

    if atom_df.empty or residue_df.empty:
        print("3E requires both atom and residue data.")
        return

    make_figure(atom_df, residue_df, figure_dir)


if __name__ == "__main__":
    main()
