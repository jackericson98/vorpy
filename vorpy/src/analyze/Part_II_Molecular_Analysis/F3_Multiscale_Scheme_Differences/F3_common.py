"""Shared data utilities for Figure 3.

Intended VorPy location
----------------------
vorpy/src/analyze/Part_II_Molecular_Analysis/F3_multiscale_differences/

The Figure 3 scripts compare additively weighted (AW) and Power Voronoi
partitions at atomic, residue, and molecular scales.  This module owns the
metric definitions so every panel uses the same matching and deviation rules.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


AW_COLOR = "#1f77b4"
POWER_COLOR = "#d62728"
ZERO_COLOR = "0.25"


@dataclass(frozen=True)
class SchemePaths:
    """Paths and identifiers for one Figure 3 molecular system."""

    name: str
    key: str
    molecule_name: str
    aw: str
    power: str
    primitive: Optional[str] = None


def signed_percent_difference(power: float, aw: float) -> float:
    """Return 100 * (Power - AW) / AW."""
    if not np.isfinite(aw) or aw == 0:
        return np.nan
    return 100.0 * (float(power) - float(aw)) / float(aw)


def absolute_percent_difference(power: float, aw: float) -> float:
    """Absolute magnitude of :func:`signed_percent_difference`."""
    value = signed_percent_difference(power, aw)
    return abs(value) if np.isfinite(value) else np.nan


def mean_absolute_percent_difference(power: Sequence[float], aw: Sequence[float]) -> float:
    """Mean absolute percent difference over matched observations."""
    power = np.asarray(power, dtype=float)
    aw = np.asarray(aw, dtype=float)
    valid = np.isfinite(power) & np.isfinite(aw) & (aw != 0)
    if not np.any(valid):
        return np.nan
    return float(np.mean(np.abs(100.0 * (power[valid] - aw[valid]) / aw[valid])))


def sem(values: Sequence[float]) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) <= 1:
        return 0.0
    return float(np.std(vals, ddof=1) / np.sqrt(len(vals)))


SYSTEM_FOLDER_RE = re.compile(r"^([A-Za-z])_(.+)$")


def discover_systems(
    data_root: str,
    exclude_keys: Optional[Iterable[str]] = None,
    verbose: bool = True,
) -> List[SchemePaths]:
    """Discover valid Figure 3 molecular-system folders.

    A valid system is an immediate child of ``data_root`` named
    ``<single letter>_<molecule name>`` and containing:

        aw/aw_logs.csv
        pow/pow_logs.csv

    Primitive is optional. When ``verbose`` is True, candidate directories are
    reported with an explicit rejection reason so a wrong root or folder layout
    is immediately visible.
    """
    root = os.path.abspath(os.path.expanduser(str(data_root)))

    if not os.path.isdir(root):
        if verbose:
            print(f"Figure 3 data root does not exist or is not a directory: {root}")
        return []

    exclude = {str(x).strip().upper() for x in (exclude_keys or [])}
    systems: List[SchemePaths] = []

    if verbose:
        print(f"\nDiscovering Figure 3 systems in: {root}")
        print(f"Excluded leading keys: {sorted(exclude) if exclude else 'none'}")

    entries = sorted(os.listdir(root))
    directory_count = 0
    pattern_count = 0

    for entry in entries:
        folder = os.path.join(root, entry)
        if not os.path.isdir(folder):
            continue

        directory_count += 1
        match = SYSTEM_FOLDER_RE.fullmatch(entry)

        if match is None:
            if verbose:
                print(f"  SKIP {entry}: name does not match <letter>_<molecule>")
            continue

        pattern_count += 1
        key = match.group(1).upper()
        molecule_name = match.group(2).strip()

        if not molecule_name:
            if verbose:
                print(f"  SKIP {entry}: molecule name is empty")
            continue

        if key in exclude:
            if verbose:
                print(f"  SKIP {entry}: leading key {key} is excluded")
            continue

        aw = os.path.join(folder, "aw", "aw_logs.csv")
        power = os.path.join(folder, "pow", "pow_logs.csv")
        primitive = os.path.join(folder, "prm", "prm_logs.csv")

        missing = []
        if not os.path.isfile(aw):
            missing.append(r"aw\aw_logs.csv")
        if not os.path.isfile(power):
            missing.append(r"pow\pow_logs.csv")

        if missing:
            if verbose:
                print(f"  SKIP {entry}: missing {', '.join(missing)}")
            continue

        if verbose:
            print(f"  OK   {entry}")

        systems.append(
            SchemePaths(
                name=entry,
                key=key,
                molecule_name=molecule_name,
                aw=aw,
                power=power,
                primitive=primitive if os.path.isfile(primitive) else None,
            )
        )

    if verbose:
        print(
            f"Discovery summary: {directory_count} immediate directories, "
            f"{pattern_count} matched the naming pattern, "
            f"{len(systems)} valid systems."
        )
        if directory_count == 0:
            print("  The selected root contains no subdirectories.")
        elif pattern_count == 0:
            print(
                "  None of the immediate subdirectories match the required "
                "<letter>_<molecule> naming convention."
            )

    return systems


def read_pair(paths: SchemePaths, need_surfs: bool = True):
    """Read matched AW and Power log files."""
    aw = read_logs2(paths.aw, all_=False, balls=True, surfs=need_surfs)
    power = read_logs2(paths.power, all_=False, balls=True, surfs=need_surfs)
    return aw, power


def _atom_name_column(df: pd.DataFrame) -> Optional[str]:
    for col in ("Atom", "Name", "Atom Name"):
        if col in df.columns:
            return col
    return None


def _chain_value(row: pd.Series) -> str:
    for col in ("Chain", "Chain ID", "ChainID", "Subunit"):
        if col in row.index:
            value = row[col]
            return "" if pd.isna(value) else str(value).strip()
    return ""


def atom_key(row: pd.Series) -> Tuple:
    """Stable atom identity used to match AW and Power rows."""
    name_col = _atom_name_column(pd.DataFrame([row]))
    atom_name = str(row[name_col]).strip().upper() if name_col else ""
    return (
        int(row["Index"]),
        _chain_value(row),
        str(row.get("Residue", "")).strip().upper(),
        int(row.get("Residue Sequence", -1)),
        atom_name,
    )


def residue_key(row: pd.Series) -> Tuple[str, str, int]:
    return (
        _chain_value(row),
        str(row.get("Residue", "")).strip().upper(),
        int(row.get("Residue Sequence", -1)),
    )


def match_atoms(aw_atoms: pd.DataFrame, power_atoms: pd.DataFrame) -> pd.DataFrame:
    """Return one row per matched atom with AW/Power volume and metadata."""
    power_lookup = {atom_key(row): row for _, row in power_atoms.iterrows()}
    records = []

    name_col = _atom_name_column(aw_atoms)
    for _, aw_row in aw_atoms.iterrows():
        key = atom_key(aw_row)
        power_row = power_lookup.get(key)
        if power_row is None:
            continue

        records.append(
            {
                "Index": int(aw_row["Index"]),
                "Chain": _chain_value(aw_row),
                "Residue": str(aw_row.get("Residue", "")).strip().upper(),
                "Residue Sequence": int(aw_row.get("Residue Sequence", -1)),
                "Atom": str(aw_row[name_col]).strip().upper() if name_col else "",
                "AW Volume": float(aw_row["Volume"]),
                "Power Volume": float(power_row["Volume"]),
            }
        )

    return pd.DataFrame.from_records(records)


def surface_pairs(surfs: pd.DataFrame) -> List[Tuple[int, int, float]]:
    """Normalize surface rows to ``(ball1, ball2, area)`` tuples."""
    rows: List[Tuple[int, int, float]] = []
    for _, surf in surfs.iterrows():
        balls = surf.get("Balls")
        if not isinstance(balls, (list, tuple, np.ndarray)) or len(balls) != 2:
            continue
        try:
            b1, b2 = int(balls[0]), int(balls[1])
            area = float(surf["Surface Area"])
        except (TypeError, ValueError, KeyError):
            continue
        rows.append((b1, b2, area))
    return rows


def atom_surface_metrics(atoms: pd.DataFrame, surfs: pd.DataFrame) -> pd.DataFrame:
    """Compute per-atom total cell surface area and number of Voronoi contacts."""
    area = {int(i): 0.0 for i in atoms["Index"]}
    contacts = {int(i): 0 for i in atoms["Index"]}

    for b1, b2, sa in surface_pairs(surfs):
        if b1 in area:
            area[b1] += sa
            contacts[b1] += 1
        if b2 in area:
            area[b2] += sa
            contacts[b2] += 1

    return pd.DataFrame(
        {
            "Index": list(area.keys()),
            "Surface Area": list(area.values()),
            "Contacts": [contacts[i] for i in area],
        }
    )


def build_atomic_metrics(aw_logs: Dict, power_logs: Dict) -> pd.DataFrame:
    """Matched atom-level AW/Power volume, SA, and contact metrics."""
    matched = match_atoms(aw_logs["atoms"], power_logs["atoms"])
    if matched.empty:
        return matched

    aw_surf = atom_surface_metrics(aw_logs["atoms"], aw_logs["surfs"]).rename(
        columns={"Surface Area": "AW Surface Area", "Contacts": "AW Contacts"}
    )
    pow_surf = atom_surface_metrics(power_logs["atoms"], power_logs["surfs"]).rename(
        columns={"Surface Area": "Power Surface Area", "Contacts": "Power Contacts"}
    )

    matched = matched.merge(aw_surf, on="Index", how="left").merge(pow_surf, on="Index", how="left")
    return matched


def build_residue_metrics(atom_metrics: pd.DataFrame, aw_logs: Dict, power_logs: Dict) -> pd.DataFrame:
    """Aggregate matched metrics to residues.

    Volume is summed over atoms. Surface area and contacts include only surfaces
    crossing the residue boundary, so internal atom-atom surfaces are excluded.
    """
    if atom_metrics.empty:
        return pd.DataFrame()

    idx_to_res = {
        int(row["Index"]): (str(row["Chain"]), str(row["Residue"]), int(row["Residue Sequence"]))
        for _, row in atom_metrics.iterrows()
    }
    residues = sorted(set(idx_to_res.values()))

    out = {
        key: {
            "Chain": key[0],
            "Residue": key[1],
            "Residue Sequence": key[2],
            "AW Volume": 0.0,
            "Power Volume": 0.0,
            "AW Surface Area": 0.0,
            "Power Surface Area": 0.0,
            "AW Contacts": 0,
            "Power Contacts": 0,
        }
        for key in residues
    }

    for _, atom in atom_metrics.iterrows():
        key = (str(atom["Chain"]), str(atom["Residue"]), int(atom["Residue Sequence"]))
        out[key]["AW Volume"] += float(atom["AW Volume"])
        out[key]["Power Volume"] += float(atom["Power Volume"])

    for scheme, logs in (("AW", aw_logs), ("Power", power_logs)):
        for b1, b2, sa in surface_pairs(logs["surfs"]):
            r1 = idx_to_res.get(b1)
            r2 = idx_to_res.get(b2)
            if r1 == r2:
                continue
            if r1 in out:
                out[r1][f"{scheme} Surface Area"] += sa
                out[r1][f"{scheme} Contacts"] += 1
            if r2 in out:
                out[r2][f"{scheme} Surface Area"] += sa
                out[r2][f"{scheme} Contacts"] += 1

    return pd.DataFrame(out.values())


def _group_metric(group_data, candidates: Sequence[str]) -> float:
    """Read a scalar from ``logs['group data']`` using known candidate names."""
    for key in candidates:
        try:
            if key in group_data:
                return float(group_data[key])
        except TypeError:
            pass
    return np.nan


def build_molecule_metrics(aw_logs: Dict, power_logs: Dict) -> Dict[str, float]:
    """Read whole-system/group totals.

    Volume and SA come directly from group data.  Contact counts are read when
    the log format exposes a group-level contact field; otherwise they remain
    NaN rather than silently substituting an inequivalent definition.
    """
    aw_group = aw_logs["group data"]
    pow_group = power_logs["group data"]

    return {
        "AW Volume": _group_metric(aw_group, ("Volume",)),
        "Power Volume": _group_metric(pow_group, ("Volume",)),
        "AW Surface Area": _group_metric(aw_group, ("Surface Area", "SurfaceArea")),
        "Power Surface Area": _group_metric(pow_group, ("Surface Area", "SurfaceArea")),
        "AW Contacts": _group_metric(aw_group, ("Contacts", "Number of Contacts", "Contact Count")),
        "Power Contacts": _group_metric(pow_group, ("Contacts", "Number of Contacts", "Contact Count")),
    }


def add_deviation_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add signed and absolute AW→Power deviations for all available metrics."""
    df = df.copy()
    for metric in ("Volume", "Surface Area", "Contacts"):
        aw_col = f"AW {metric}"
        pow_col = f"Power {metric}"
        if aw_col not in df.columns or pow_col not in df.columns:
            continue
        aw = pd.to_numeric(df[aw_col], errors="coerce")
        power = pd.to_numeric(df[pow_col], errors="coerce")
        signed = 100.0 * (power - aw) / aw.replace(0, np.nan)
        df[f"{metric} Signed % Diff"] = signed
        df[f"{metric} Abs % Diff"] = signed.abs()
    return df


def summarize_scale(df: pd.DataFrame, scale: str) -> pd.DataFrame:
    """Summarize mean absolute deviations for one atom/residue dataframe."""
    records = []
    for metric in ("Volume", "Surface Area", "Contacts"):
        col = f"{metric} Abs % Diff"
        if col not in df:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").dropna().to_numpy(float)
        records.append(
            {
                "Scale": scale,
                "Metric": metric,
                "Mean Abs % Diff": float(np.mean(vals)) if len(vals) else np.nan,
                "SEM": sem(vals),
                "N": int(len(vals)),
            }
        )
    return pd.DataFrame(records)
