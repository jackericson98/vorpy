"""
combine_curvature_peaks_folder.py

Select a single root folder (like 4A_Peaks.py), crawl its subfolders, and collect
curvature bin files named:

    <MODEL_NAME>_curvature_bins.csv

Then detect peaks per (SourceFile, Pair) and write ONE combined CSV.

Assumptions about each curvature_bins.csv (from your Hairpin example / 4A export):
  required columns: Bin, BinLow, BinHigh, Pair, Count
  optional columns (ignored): Residue1, Residue2, CurvMean, CurvStd, SurfAreaMean, ...

By default:
  - Aggregates Count across residues for each (Pair, Bin) within each file.
  - Detects peaks on the aggregated histogram per Pair.
  - Outputs combined_peaks.csv in the selected root folder (unless overridden).

Usage:
  python combine_curvature_peaks_folder.py
  python combine_curvature_peaks_folder.py --exclude A B --smooth_window 7 --top_k 5
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog




@dataclass
class Peak:
    source_file: str
    source_path: str
    model_name: str
    model_folder: str
    pair: str
    peak_rank: int
    peak_index: int
    bin_id: int
    bin_low: float
    bin_high: float
    bin_center: float
    height: float
    prominence: float
    left_base_center: float
    right_base_center: float




def _moving_average(y: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return y.astype(float)

    if window % 2 == 0:
        window += 1

    window = max(1, int(window))
    kernel = np.ones(window, dtype=float) / float(window)

    return np.convolve(y.astype(float), kernel, mode="same")




def _local_maxima_indices(y: np.ndarray) -> np.ndarray:
    if y.size < 3:
        return np.array([], dtype=int)

    left = y[:-2]
    mid = y[1:-1]
    right = y[2:]

    mask = (mid > left) & (mid >= right)

    return np.where(mask)[0] + 1




def _peak_prominence(y: np.ndarray, i: int) -> Tuple[float, int, int]:
    """
    Robust, scipy-free prominence estimate.

    Returns:
      (prominence, left_base_index, right_base_index)

    Guarantees:
      0 <= left_base_index <= i <= right_base_index <= n-1
    """
    n = int(y.size)
    if n == 0:
        return 0.0, 0, 0

    i = int(np.clip(i, 0, n - 1))
    peak_h = float(y[i])

    # --- Find a left boundary where the signal rises above the peak (or edge) ---
    left = i
    while left > 0 and float(y[left - 1]) <= peak_h:
        left -= 1

    # --- Find a right boundary where the signal rises above the peak (or edge) ---
    right = i
    while right < n - 1 and float(y[right + 1]) <= peak_h:
        right += 1

    # Now compute the minima on both sides within these safe ranges
    # Left range includes i
    left_slice = y[left : i + 1]
    left_min_rel = int(np.argmin(left_slice))
    left_min_idx = left + left_min_rel
    left_min = float(y[left_min_idx])

    # Right range includes i
    right_slice = y[i : right + 1]
    right_min_rel = int(np.argmin(right_slice))
    right_min_idx = i + right_min_rel
    right_min = float(y[right_min_idx])

    ref = max(left_min, right_min)
    prom = peak_h - ref

    # Clamp indices defensively
    left_min_idx = int(np.clip(left_min_idx, 0, n - 1))
    right_min_idx = int(np.clip(right_min_idx, 0, n - 1))

    return float(prom), left_min_idx, right_min_idx





def _load_and_aggregate(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath)

    required = {"Bin", "BinLow", "BinHigh", "Pair", "Count"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{filepath}: missing required columns: {missing}")

    df = df.copy()
    df["BinCenter"] = 0.5 * (df["BinLow"].astype(float) + df["BinHigh"].astype(float))

    # Aggregate across residue pairs if present (or any other grouping columns)
    agg = (
        df.groupby(["Pair", "Bin", "BinLow", "BinHigh", "BinCenter"], as_index=False)["Count"]
        .sum()
        .rename(columns={"Count": "AggCount"})
    )

    return agg




def _detect_peaks_for_pair(
    pair_df: pd.DataFrame,
    *,
    source_file: str,
    source_path: str,
    model_name: str,
    model_folder: str,
    smooth_window: int,
    min_height: float | None,
    min_prominence: float | None,
    top_k: int | None,
) -> List[Peak]:
    pair = str(pair_df["Pair"].iloc[0])

    pair_df = pair_df.sort_values(["BinLow", "BinHigh", "Bin"], kind="mergesort").reset_index(drop=True)

    y_raw = pair_df["AggCount"].to_numpy(dtype=float)
    y_smooth = _moving_average(y_raw, smooth_window)

    peak_idxs = _local_maxima_indices(y_smooth)

    peaks: List[Tuple[int, float, float, int, int]] = []
    for idx in peak_idxs:
        height = float(y_smooth[idx])
        if min_height is not None and height < float(min_height):
            continue

        prom, left_base_idx, right_base_idx = _peak_prominence(y_smooth, int(idx))
        if min_prominence is not None and prom < float(min_prominence):
            continue

        peaks.append((int(idx), height, prom, left_base_idx, right_base_idx))

    peaks.sort(key=lambda t: (t[2], t[1]), reverse=True)

    if top_k is not None and int(top_k) > 0:
        peaks = peaks[: int(top_k)]

    out: List[Peak] = []
    for rank, (idx, height, prom, lbi, rbi) in enumerate(peaks, start=1):
        row = pair_df.iloc[int(idx)]
        n = len(pair_df)
        lbi = int(np.clip(int(lbi), 0, n - 1))
        rbi = int(np.clip(int(rbi), 0, n - 1))

        lb = pair_df.iloc[lbi]
        rb = pair_df.iloc[rbi]

        out.append(
            Peak(
                source_file=source_file,
                source_path=source_path,
                model_name=model_name,
                model_folder=model_folder,
                pair=pair,
                peak_rank=int(rank),
                peak_index=int(idx),
                bin_id=int(row["Bin"]),
                bin_low=float(row["BinLow"]),
                bin_high=float(row["BinHigh"]),
                bin_center=float(row["BinCenter"]),
                height=float(height),
                prominence=float(prom),
                left_base_center=float(lb["BinCenter"]),
                right_base_center=float(rb["BinCenter"]),
            )
        )

    return out


def _find_curvature_files(root_folder: str, excluded_prefixes: List[str] | None = None) -> List[Dict[str, str]]:
    """
    Mimics the 4A_Peaks.py folder crawl style:
      - root_folder contains multiple model subfolders
      - each model subfolder contains <MODEL>_curvature_bins.csv somewhere (often directly in that folder)
    This implementation looks for curvature bin CSVs directly inside each immediate subfolder.
    If you also want it to recurse deeper, flip RECURSE=True below.
    """
    if excluded_prefixes is None:
        excluded_prefixes = []

    RECURSE = False

    hits: List[Dict[str, str]] = []

    for subfolder in os.listdir(root_folder):
        if not subfolder:
            continue

        if subfolder[0] in excluded_prefixes:
            continue

        model_folder = os.path.join(root_folder, subfolder)
        if not os.path.isdir(model_folder):
            continue

        if RECURSE:
            for dirpath, _, filenames in os.walk(model_folder):
                for fn in filenames:
                    if fn.endswith("_curvature_bins.csv"):
                        model_name = fn[: -len("_curvature_bins.csv")]
                        hits.append(
                            {
                                "model_name": model_name,
                                "model_folder": model_folder,
                                "source_path": os.path.join(dirpath, fn),
                                "source_file": fn,
                            }
                        )
        else:
            for fn in os.listdir(model_folder):
                if fn.endswith("_curvature_bins.csv"):
                    model_name = fn[: -len("_curvature_bins.csv")]
                    hits.append(
                        {
                            "model_name": model_name,
                            "model_folder": model_folder,
                            "source_path": os.path.join(model_folder, fn),
                            "source_file": fn,
                        }
                    )

    return hits


def combine_peaks_from_folder(
    root_folder: str,
    *,
    excluded_prefixes: List[str] | None = None,
    smooth_window: int = 5,
    min_height: float | None = None,
    min_prominence: float | None = None,
    top_k: int | None = 5,
) -> pd.DataFrame:
    entries = _find_curvature_files(root_folder, excluded_prefixes=excluded_prefixes)

    if not entries:
        return pd.DataFrame(
            columns=[
                "SourceFile",
                "SourcePath",
                "Model",
                "ModelFolder",
                "Pair",
                "PeakRank",
                "PeakIndex",
                "Bin",
                "BinLow",
                "BinHigh",
                "BinCenter",
                "Height",
                "Prominence",
                "LeftBaseCenter",
                "RightBaseCenter",
            ]
        )

    all_peaks: List[Peak] = []

    for e in entries:
        source_path = e["source_path"]
        source_file = e["source_file"]
        model_name = e["model_name"]
        model_folder = e["model_folder"]

        try:
            agg = _load_and_aggregate(source_path)
        except Exception as exc:
            print(f"Skipping {source_path}: {exc}")
            continue

        for _, pair_df in agg.groupby("Pair", sort=False):
            peaks = _detect_peaks_for_pair(
                pair_df,
                source_file=source_file,
                source_path=source_path,
                model_name=model_name,
                model_folder=model_folder,
                smooth_window=smooth_window,
                min_height=min_height,
                min_prominence=min_prominence,
                top_k=top_k,
            )
            all_peaks.extend(peaks)

    rows = []
    for p in all_peaks:
        rows.append(
            {
                "SourceFile": p.source_file,
                "SourcePath": p.source_path,
                "Model": p.model_name,
                "ModelFolder": p.model_folder,
                "Pair": p.pair,
                "PeakRank": p.peak_rank,
                "PeakIndex": p.peak_index,
                "Bin": p.bin_id,
                "BinLow": p.bin_low,
                "BinHigh": p.bin_high,
                "BinCenter": p.bin_center,
                "Height": p.height,
                "Prominence": p.prominence,
                "LeftBaseCenter": p.left_base_center,
                "RightBaseCenter": p.right_base_center,
            }
        )

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        return out_df

    out_df = out_df.sort_values(
        ["Pair", "Model", "SourceFile", "PeakRank", "Prominence", "Height"],
        ascending=[True, True, True, True, False, False],
        kind="mergesort",
    ).reset_index(drop=True)

    return out_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Folder-based combine of curvature peaks (Figure 4A).")
    parser.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Exclude model subfolders by FIRST letter (same idea as 4A_Peaks.py excluded list).",
    )
    parser.add_argument(
        "--smooth_window",
        type=int,
        default=5,
        help="Odd moving-average window for smoothing (default 5).",
    )
    parser.add_argument(
        "--min_height",
        type=float,
        default=None,
        help="Minimum peak height after smoothing (default None).",
    )
    parser.add_argument(
        "--min_prominence",
        type=float,
        default=None,
        help="Minimum peak prominence (default None).",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        default=5,
        help="Keep top K peaks per Pair per file. Use <=0 for all peaks.",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output CSV path. Default is <selected_root>/combined_peaks.csv",
    )

    args = parser.parse_args()

    # Select root folder (same UX style as 4A_Peaks.py)
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)

    root_folder = filedialog.askdirectory(title="Choose Molecular Data Directory")
    if not root_folder:
        print("No folder selected.")
        return

    top_k = None
    if args.top_k is not None and int(args.top_k) > 0:
        top_k = int(args.top_k)

    peaks_df = combine_peaks_from_folder(
        root_folder,
        excluded_prefixes=list(args.exclude) if args.exclude else [],
        smooth_window=int(args.smooth_window),
        min_height=args.min_height,
        min_prominence=args.min_prominence,
        top_k=top_k,
    )

    out_path = args.out
    if out_path is None:
        out_path = os.path.join(root_folder, "combined_peaks.csv")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    peaks_df.to_csv(out_path, index=False)

    print(f"Wrote combined peaks CSV: {os.path.abspath(out_path)}")
    print(f"Peak rows: {len(peaks_df)}")


if __name__ == "__main__":


    main()
