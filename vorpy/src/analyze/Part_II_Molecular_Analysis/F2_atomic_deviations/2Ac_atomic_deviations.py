import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Set
import sys

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.plot_templates.scatter import scatter


@dataclass
class MoleculeInputs:
    name: str
    pow_logs: str
    aw_logs: str


def _safe_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _extract_element(atom_name: str) -> str:
    """
    Best-effort element extraction from your 'Name' field.

    Handles cases like:
      - "C", "CA", "CL", "NA"
      - "C-HN", "C-HB" (returns "C")
      - "OW", "HW1" (returns "O" or "H" heuristically)
    """
    if atom_name is None:
        return "UNK"

    s = str(atom_name).strip().upper()

    # Common water/ion heuristics used in your codebase
    if s.startswith("OW"):
        return "O"
    if s.startswith("HW") or s.startswith("H0") or s.startswith("H"):
        return "H"
    if s.startswith("CL"):
        return "CL"
    if s.startswith("NA"):
        return "NA"
    if s.startswith("MG"):
        return "MG"
    if s.startswith("K"):
        return "K"

    # If you have atom-class names like "C-HN", take prefix before '-'
    if "-" in s:
        s = s.split("-")[0]

    # Take leading letters until first non-letter
    letters = []
    for ch in s:
        if ch.isalpha():
            letters.append(ch)
        else:
            break

    if not letters:
        return "UNK"

    # Prefer 2-letter elements when appropriate
    candidate = "".join(letters)
    if len(candidate) >= 2 and candidate[:2] in {"CL", "NA", "MG"}:
        return candidate[:2]

    return candidate[0]


def _key_cols() -> List[str]:
    return ["Index", "Name", "Residue", "Residue Sequence", "Chain"]


def _normalize_atoms_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize column names / types used downstream.
    """
    out = df.copy()

    # Ensure expected fields exist
    for c in _key_cols():
        if c not in out.columns:
            out[c] = np.nan

    if "Volume" in out.columns:
        out["Volume"] = _safe_float_series(out["Volume"])
    if "Surface Area" in out.columns:
        out["Surface Area"] = _safe_float_series(out["Surface Area"])

    if "Number of Neighbors" in out.columns:
        out["Number of Neighbors"] = pd.to_numeric(out["Number of Neighbors"], errors="coerce")

    # Curvature: you have both mean and gauss average columns in newer logs
    for c in [
        "Average Mean Surface Curvature",
        "Average Gaussian Surface Curvature",
        "Maximum Mean Curvature",
        "Maximum Gaussian Curvature",
    ]:
        if c in out.columns:
            out[c] = _safe_float_series(out[c])

    out["Element"] = out["Name"].apply(_extract_element)

    return out


def _merge_pow_aw(pow_atoms: pd.DataFrame, aw_atoms: pd.DataFrame) -> pd.DataFrame:
    """
    Merge on atom identity. Primary key = Index.
    """
    pow_df = _normalize_atoms_df(pow_atoms)
    aw_df = _normalize_atoms_df(aw_atoms)

    key = ["Index"]

    merged = pow_df.merge(
        aw_df,
        on=key,
        suffixes=("_pow", "_aw"),
        how="inner",
    )

    return merged


def _abs_percent_diff(a: pd.Series, b: pd.Series, denom: pd.Series) -> pd.Series:
    eps = 1e-12
    d = (a - b).abs()
    return (d / (denom.abs() + eps)) * 100.0


def _choose_curvature_col(merged: pd.DataFrame, curvature_kind: str) -> str:
    """
    Use AW curvature by default (you can swap), but keep it consistent.
    """
    curvature_kind = curvature_kind.strip().lower()
    if curvature_kind in {"mean", "avg_mean", "mean_avg"}:
        return "Average Mean Surface Curvature_aw" if "Average Mean Surface Curvature_aw" in merged.columns else "Average Mean Surface Curvature_pow"
    if curvature_kind in {"gauss", "gaussian", "avg_gauss"}:
        return "Average Gaussian Surface Curvature_aw" if "Average Gaussian Surface Curvature_aw" in merged.columns else "Average Gaussian Surface Curvature_pow"

    raise ValueError(f"Unknown curvature_kind='{curvature_kind}'. Use 'mean' or 'gauss'.")


def _apply_quality_filters(
    df: pd.DataFrame,
    require_complete_cells: bool,
    min_neighbors: Optional[int],
    max_neighbors: Optional[int],
) -> pd.DataFrame:
    out = df.copy()

    # Drop non-finite geometry
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["Volume_pow", "Volume_aw", "Surface Area_pow", "Surface Area_aw"], how="any")

    # Optional: require complete cells if logs provide it
    if require_complete_cells:
        c_pow = "Complete Cell?_pow"
        c_aw = "Complete Cell?_aw"
        if c_pow in out.columns and c_aw in out.columns:
            out = out[(out[c_pow] == True) & (out[c_aw] == True)]

    # Optional: neighbor window
    if min_neighbors is not None:
        n_pow = "Number of Neighbors_pow"
        n_aw = "Number of Neighbors_aw"
        if n_pow in out.columns:
            out = out[out[n_pow] >= min_neighbors]
        if n_aw in out.columns:
            out = out[out[n_aw] >= min_neighbors]

    if max_neighbors is not None:
        n_pow = "Number of Neighbors_pow"
        n_aw = "Number of Neighbors_aw"
        if n_pow in out.columns:
            out = out[out[n_pow] <= max_neighbors]
        if n_aw in out.columns:
            out = out[out[n_aw] <= max_neighbors]

    return out.reset_index(drop=True)


def select_exemplar_candidates(
    merged: pd.DataFrame,
    curvature_kind: str,
    n_candidates_per_category: int,
    exclude_indices: Optional[Set[int]] = None,
    require_complete_cells: bool = False,
    min_neighbors: Optional[int] = None,
    max_neighbors: Optional[int] = None,
    add_outlier: bool = True,
) -> pd.DataFrame:
    """
    Returns multiple candidates per exemplar category so you can choose
    alternatives when a mesh rendering looks odd.

    Categories:
      - high_curv_high_dV: highest |%ΔV| among top curvature decile
      - low_curv_high_dV : highest |%ΔV| among bottom curvature decile
      - median_dV_control: atoms near median |%ΔV|
      - global_max_dV_outlier: optional global max |%ΔV| (1 row)

    Output includes:
      - exemplar_category
      - candidate_rank (1..k within each category)
    """
    df = merged.copy()

    # Core metrics
    df["abs_pct_dV_pow"] = _abs_percent_diff(df["Volume_pow"], df["Volume_aw"], denom=df["Volume_pow"])
    df["abs_pct_dSA_pow"] = _abs_percent_diff(df["Surface Area_pow"], df["Surface Area_aw"], denom=df["Surface Area_pow"])

    curv_col = _choose_curvature_col(df, curvature_kind=curvature_kind)
    df["curv_used"] = df[curv_col]

    df = df.dropna(subset=["abs_pct_dV_pow", "curv_used"]).reset_index(drop=True)

    # Exclude known-bad indices (e.g., weird render)
    if exclude_indices:
        df = df[~df["Index"].astype(int).isin(set(int(x) for x in exclude_indices))].reset_index(drop=True)

    # Quality filters
    df = _apply_quality_filters(
        df=df,
        require_complete_cells=require_complete_cells,
        min_neighbors=min_neighbors,
        max_neighbors=max_neighbors,
    )

    # Curvature deciles
    q_lo = df["curv_used"].quantile(0.10)
    q_hi = df["curv_used"].quantile(0.90)

    high_curv = df[df["curv_used"] >= q_hi].sort_values("abs_pct_dV_pow", ascending=False)
    low_curv = df[df["curv_used"] <= q_lo].sort_values("abs_pct_dV_pow", ascending=False)

    picks: List[pd.DataFrame] = []

    if len(high_curv) > 0:
        cand = high_curv.head(n_candidates_per_category).copy()
        cand["exemplar_category"] = "high_curv_high_dV"
        cand["candidate_rank"] = np.arange(1, len(cand) + 1)
        picks.append(cand)

    if len(low_curv) > 0:
        cand = low_curv.head(n_candidates_per_category).copy()
        cand["exemplar_category"] = "low_curv_high_dV"
        cand["candidate_rank"] = np.arange(1, len(cand) + 1)
        picks.append(cand)

    # Median control: take a small window around the median and then pick k rows
    df_sorted = df.sort_values("abs_pct_dV_pow")
    if len(df_sorted) > 0:
        med_i = len(df_sorted) // 2
        half_window = max(25, n_candidates_per_category * 10)
        lo = max(0, med_i - half_window)
        hi = min(len(df_sorted), med_i + half_window + 1)

        window = df_sorted.iloc[lo:hi].copy()

        # Prefer "typical" neighbors if available: closest to median neighbor count
        if "Number of Neighbors_pow" in window.columns:
            med_n = window["Number of Neighbors_pow"].median()
            window["median_neighbor_dist"] = (window["Number of Neighbors_pow"] - med_n).abs()
            window = window.sort_values(["abs_pct_dV_pow", "median_neighbor_dist"], ascending=[True, True])
            window = window.drop(columns=["median_neighbor_dist"], errors="ignore")
        else:
            window = window.sort_values("abs_pct_dV_pow", ascending=True)

        cand = window.head(n_candidates_per_category).copy()
        cand["exemplar_category"] = "median_dV_control"
        cand["candidate_rank"] = np.arange(1, len(cand) + 1)
        picks.append(cand)

    if add_outlier and len(df_sorted) > 0:
        out_row = df_sorted.tail(1).copy()
        out_row["exemplar_category"] = "global_max_dV_outlier"
        out_row["candidate_rank"] = 1
        picks.append(out_row)

    out = pd.concat(picks, ignore_index=True)

    keep = [
        "Index",
        "Name_pow",
        "Element_pow",
        "Residue_pow",
        "Residue Sequence_pow",
        "Chain_pow",
        "Volume_pow",
        "Volume_aw",
        "Surface Area_pow",
        "Surface Area_aw",
        "Number of Neighbors_pow",
        "Number of Neighbors_aw",
        "curv_used",
        "abs_pct_dV_pow",
        "abs_pct_dSA_pow",
        "exemplar_category",
        "candidate_rank",
    ]

    existing = [c for c in keep if c in out.columns]
    out = out[existing].copy()

    out = out.sort_values(
        ["exemplar_category", "candidate_rank", "abs_pct_dV_pow"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    return out


def run_for_molecule(
    mol: MoleculeInputs,
    out_dir: str,
    curvature_kind: str,
    n_candidates_per_category: int,
    exclude_indices: Optional[Set[int]] = None,
    require_complete_cells: bool = False,
    min_neighbors: Optional[int] = None,
    max_neighbors: Optional[int] = None,
) -> Dict:
    pow_data = read_logs2(mol.pow_logs, return_dict=False, all_=True)
    aw_data = read_logs2(mol.aw_logs, return_dict=False, all_=True)

    pow_atoms = pow_data["atoms"]
    aw_atoms = aw_data["atoms"]

    merged = _merge_pow_aw(pow_atoms=pow_atoms, aw_atoms=aw_atoms)

    merged["abs_pct_dV_pow"] = _abs_percent_diff(merged["Volume_pow"], merged["Volume_aw"], denom=merged["Volume_pow"])
    merged["abs_pct_dSA_pow"] = _abs_percent_diff(merged["Surface Area_pow"], merged["Surface Area_aw"], denom=merged["Surface Area_pow"])

    curv_col = _choose_curvature_col(merged, curvature_kind=curvature_kind)
    merged["curv_used"] = merged[curv_col]

    os.makedirs(out_dir, exist_ok=True)

    merged_csv = os.path.join(out_dir, f"{mol.name}_pow_vs_aw_atoms_merged.csv")
    merged.to_csv(merged_csv, index=False)

    exemplars = select_exemplar_candidates(
        merged=merged,
        curvature_kind=curvature_kind,
        n_candidates_per_category=n_candidates_per_category,
        exclude_indices=exclude_indices,
        require_complete_cells=require_complete_cells,
        min_neighbors=min_neighbors,
        max_neighbors=max_neighbors,
        add_outlier=True,
    )

    exemplars_csv = os.path.join(out_dir, f"{mol.name}_exemplar_candidates.csv")
    exemplars.to_csv(exemplars_csv, index=False)

    summary = {
        "molecule": mol.name,
        "pow_logs": mol.pow_logs,
        "aw_logs": mol.aw_logs,
        "n_pow_atoms": int(len(pow_atoms)),
        "n_aw_atoms": int(len(aw_atoms)),
        "n_merged_atoms": int(len(merged)),
        "curvature_kind": curvature_kind,
        "n_candidates_per_category": int(n_candidates_per_category),
        "exclude_indices": sorted(list(exclude_indices)) if exclude_indices else [],
        "require_complete_cells": bool(require_complete_cells),
        "min_neighbors": min_neighbors,
        "max_neighbors": max_neighbors,
        "merged_csv": merged_csv,
        "exemplars_csv": exemplars_csv,
    }

    summary_json = os.path.join(out_dir, f"{mol.name}_summary.json")
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    """
    HARD-CODED INPUTS (edit these paths to your local machine):
    """

    root_dir = filedialog.askdirectory(title="Select the molecules directory")
    out_dir = filedialog.askdirectory(title="Select output directory")

    molecules = [
        MoleculeInputs(
            name="hammerhead",
            pow_logs=os.path.join(root_dir, "G_Hammerhead", "pow", "hammerhead_pow_logs.csv"),
            aw_logs=os.path.join(root_dir, "G_Hammerhead", "aw", "hammerhead_aw_logs.csv"),
        ),
        MoleculeInputs(
            name="p53tet",
            pow_logs=os.path.join(root_dir, "H_p53tet", "pow", "p53tet_pow_logs.csv"),
            aw_logs=os.path.join(root_dir, "H_p53tet", "aw", "p53tet_aw_logs.csv"),
        ),
    ]

    curvature_kind = "mean"  # or "gauss"

    # How many alternatives you want per category (recommended: 3 or 5)
    n_candidates_per_category = 10

    # Known-bad / weird renderings you want to exclude from exemplar selection
    exclude_indices_by_molecule = {
        "p53tet": {639, 656, 657, 1335, 1334, 640, 654},
        "hammerhead": set(),
    }

    # Optional quality gates (safe to leave False/None)
    require_complete_cells = False
    min_neighbors = None
    max_neighbors = None

    summaries = []
    for mol in molecules:
        summaries.append(
            run_for_molecule(
                mol=mol,
                out_dir=out_dir,
                curvature_kind=curvature_kind,
                n_candidates_per_category=n_candidates_per_category,
                exclude_indices=exclude_indices_by_molecule.get(mol.name, set()),
                require_complete_cells=require_complete_cells,
                min_neighbors=min_neighbors,
                max_neighbors=max_neighbors,
            )
        )

    with open(os.path.join(out_dir, "all_summaries.json"), "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)


if __name__ == "__main__":


    main()
