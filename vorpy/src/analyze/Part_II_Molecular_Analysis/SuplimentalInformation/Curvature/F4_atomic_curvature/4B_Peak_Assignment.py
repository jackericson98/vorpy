import os
import sys
import re

import tkinter as tk
from tkinter import filedialog

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Get the path to the root vorpy folder (reuse your pattern)
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)



from vorpy.src.analyze.tools.plot_templates.line import line_plot


LETTER_CODES = {
    "NA5": "A",
    "EDTA": "B",
    "DB1976": "C",
    "HAIRPIN": "D",
    "CAMBRIN": "E",
    "BDNA": "F",
    "HAMMERHEAD": "G",
    "P53TET": "H",
    "T4LP": "I",
    "STREPTAVIDIN": "J",
    "NCP": "K",
    "BSA": "L",
    "BSA_20": "L",
}

COLOR_MAP = {
    "NA5": "#000000",
    "EDTA": "#0077FF",
    "DB1976": "#FF7F00",
    "HAIRPIN": "#008F00",
    "CAMBRIN": "#00C8FF",
    "BDNA": "#B300FF",
    "HAMMERHEAD": "#FF0000",
    "P53TET": "#FF00AA",
    "T4LP": "#00FF00",
    "STREPTAVIDIN": "#A52A2A",
    "NCP": "#8C8C8C",
    "BSA": "#FFD700",
    "BSA_20": "#FFD700",
}


import numpy as np
import pandas as pd



def _peak_signature(peak_df: pd.DataFrame, topk: int = 10) -> tuple[list[str], dict[str, float]]:
    """
    Returns:
      - topk NormPair list by PairCount
      - normalized weight dict over NormPair (PairCount / sum)
    """
    d = (
        peak_df.groupby("NormPair", sort=False)["PairCount"]
        .sum()
        .sort_values(ascending=False)
    )

    top = d.head(topk).index.tolist()

    total = float(d.sum())
    weights = {}
    if total > 0:
        weights = {k: float(v) / total for k, v in d.items()}

    return top, weights


def _cosine_similarity(wA: dict[str, float], wB: dict[str, float]) -> float:
    keys = sorted(set(wA) | set(wB))
    a = np.array([wA.get(k, 0.0) for k in keys], dtype=float)
    b = np.array([wB.get(k, 0.0) for k in keys], dtype=float)

    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        return 0.0

    return float(np.dot(a, b) / (na * nb))


def merge_consecutive_similar_peaks(
    peaks_df: pd.DataFrame,
    bin_group: pd.DataFrame,
    percents: np.ndarray,
    *,
    topk: int = 10,
    jaccard_thresh: float = 0.8,
    cosine_thresh: float = 0.95,
) -> pd.DataFrame:
    """
    Merge consecutive peaks if:
      1) No zero valley between them (bins between have nonzero percent)
      2) Similar composition (Jaccard on top-k + cosine on weights)

    Assumes one model at a time in peaks_df.
    """
    if peaks_df.empty:
        return peaks_df

    # Work peak-by-peak in PeakIndex order
    peak_ids = sorted(peaks_df["PeakIndex"].unique().tolist())
    merged_groups = []

    i = 0
    while i < len(peak_ids):
        pid = peak_ids[i]
        cur = peaks_df[peaks_df["PeakIndex"] == pid].copy()

        # Try to merge forward repeatedly
        j = i + 1
        while j < len(peak_ids):
            pid2 = peak_ids[j]
            nxt = peaks_df[peaks_df["PeakIndex"] == pid2].copy()

            # --- condition (1): no zero valley between ranges ---
            R_A = int(cur.iloc[0]["RangeBinRight"])
            L_B = int(nxt.iloc[0]["RangeBinLeft"])

            if L_B <= R_A + 1:
                valley_ok = True
            else:
                # Map Bin -> index in bin_group (Bin is an int column)
                # Faster: create mapping once
                valley_bins = bin_group[(bin_group["Bin"] > R_A) & (bin_group["Bin"] < L_B)]
                if valley_bins.empty:
                    valley_ok = True
                else:
                    # if any in-between bin is exactly 0, treat as separated
                    valley_ok = bool((valley_bins["TotalCount"] > 0).all())

            if not valley_ok:
                break

            # --- condition (2): composition similarity ---
            topA, wA = _peak_signature(cur, topk=topk)
            topB, wB = _peak_signature(nxt, topk=topk)

            setA = set(topA)
            setB = set(topB)
            jacc = len(setA & setB) / float(len(setA | setB) or 1)

            cos = _cosine_similarity(wA, wB)

            if jacc < jaccard_thresh or cos < cosine_thresh:
                break

            # ✅ Merge nxt into cur
            cur = pd.concat([cur, nxt], ignore_index=True)

            # Expand the merged range bounds
            cur.loc[:, "RangeBinLeft"] = min(R_A, L_B, int(cur["RangeBinLeft"].min()))
            cur.loc[:, "RangeBinRight"] = max(int(cur["RangeBinRight"].max()), int(nxt["RangeBinRight"].max()))
            cur.loc[:, "RangeLow"] = float(cur["RangeLow"].min())
            cur.loc[:, "RangeHigh"] = float(cur["RangeHigh"].max())

            # Recompute peak-range percent after merging bins
            left_bin = int(cur.iloc[0]["RangeBinLeft"])
            right_bin = int(cur.iloc[0]["RangeBinRight"])

            # Sum TotalCount over those bins from bin_group (more stable than summing PairCount)
            sub = bin_group[(bin_group["Bin"] >= left_bin) & (bin_group["Bin"] <= right_bin)]
            range_total = float(sub["TotalCount"].sum())
            total_surfs = float(bin_group["TotalCount"].sum())
            range_percent_all = 100.0 * range_total / total_surfs if total_surfs > 0 else 0.0

            cur.loc[:, "RangePercentAllSurfs"] = range_percent_all

            # Recompute per-pair PeakShare + SurfaceType% within merged range
            # PairCount sums
            summed = (
                cur.groupby("NormPair", sort=False)
                .agg(
                    PairCount=("PairCount", "sum"),
                    CurvMean=("CurvMean", "mean"),
                    CurvStd=("CurvStd", "mean"),
                    SurfaceTypePercent=("SurfaceTypePercent", "max"),  # keep existing global % per type
                )
                .reset_index()
            )

            range_pair_total = float(summed["PairCount"].sum())
            summed["PeakShare"] = summed["PairCount"] / range_pair_total if range_pair_total > 0 else 0.0

            # Carry over shared columns from cur (use first row)
            shared_cols = [c for c in cur.columns if c not in {"NormPair", "PairCount", "PeakShare", "CurvMean", "CurvStd"}]
            base = cur.iloc[0][shared_cols].to_dict()

            # Build new cur table in same schema
            new_cur = []
            for _, r in summed.iterrows():
                row = dict(base)
                row["NormPair"] = r["NormPair"]
                row["PairCount"] = int(r["PairCount"])
                row["PeakShare"] = float(r["PeakShare"])
                row["CurvMean"] = float(r["CurvMean"])
                row["CurvStd"] = float(r["CurvStd"])
                row["SurfaceTypePercent"] = float(r.get("SurfaceTypePercent", 0.0))
                new_cur.append(row)

            cur = pd.DataFrame(new_cur)

            # consume pid2
            j += 1

        # assign a new PeakIndex (compressed later if you want)
        merged_groups.append(cur)
        i = j

    out = pd.concat(merged_groups, ignore_index=True)

    # Re-number PeakIndex sequentially (nice for the txt)
    new_ids = {old: new for new, old in enumerate(sorted(out["PeakIndex"].unique()), start=1)}
    out["PeakIndex"] = out["PeakIndex"].map(new_ids)

    return out


def merge_consecutive_almost_identical_peaks(
    peaks_df: pd.DataFrame,
    bin_group: pd.DataFrame,
    *,
    topk: int = 12,
    jaccard_thresh: float = 0.85,
    count_rel_tol: float = 0.05,   # 5% per-pair relative tolerance
    min_pair_count: int = 5,       # ignore tiny pairs when comparing
) -> pd.DataFrame:
    """
    Merge consecutive peaks when:
      1) There is no zero-count bin between their ranges, AND
      2) Their dominant NormPair composition is essentially the same:
         - Jaccard(topK NormPairs) >= threshold
         - For shared dominant pairs, PairCount agrees within relative tolerance
    """
    if peaks_df.empty:
        return peaks_df

    peak_ids = sorted(peaks_df["PeakIndex"].unique().tolist())
    blocks: list[pd.DataFrame] = []

    def sig(df_peak: pd.DataFrame):
        s = (
            df_peak.groupby("NormPair")["PairCount"]
            .sum()
            .sort_values(ascending=False)
        )
        # drop tiny pairs to avoid noise-driven non-merges
        s = s[s >= min_pair_count]
        top = s.head(topk)
        return top  # Series: NormPair -> PairCount

    i = 0
    while i < len(peak_ids):
        pid = peak_ids[i]
        cur = peaks_df[peaks_df["PeakIndex"] == pid].copy()

        j = i + 1
        while j < len(peak_ids):
            pid2 = peak_ids[j]
            nxt = peaks_df[peaks_df["PeakIndex"] == pid2].copy()

            # ---- (1) valley check: no zero-count bin between ranges ----
            R_A = int(cur.iloc[0]["RangeBinRight"])
            L_B = int(nxt.iloc[0]["RangeBinLeft"])

            if L_B > R_A + 1:
                mid = bin_group[(bin_group["Bin"] > R_A) & (bin_group["Bin"] < L_B)]
                if not mid.empty and (mid["TotalCount"] == 0).any():
                    break

            # ---- (2) composition similarity ----
            A = sig(cur)
            B = sig(nxt)

            setA = set(A.index.tolist())
            setB = set(B.index.tolist())
            if len(setA | setB) == 0:
                break

            jacc = len(setA & setB) / float(len(setA | setB))
            if jacc < jaccard_thresh:
                break

            # For overlapping dominant pairs, require close counts (relative)
            common = list(setA & setB)
            ok = True
            for p in common:
                a = float(A[p])
                b = float(B[p])
                denom = max(a, b, 1.0)
                if abs(a - b) / denom > count_rel_tol:
                    ok = False
                    break

            if not ok:
                break

            # ✅ Merge nxt into cur: concatenate rows (we'll recompute shares)
            cur = pd.concat([cur, nxt], ignore_index=True)

            # Expand merged range
            cur_left = min(int(cur["RangeBinLeft"].min()), int(nxt["RangeBinLeft"].min()))
            cur_right = max(int(cur["RangeBinRight"].max()), int(nxt["RangeBinRight"].max()))
            cur_low = float(min(cur["RangeLow"].min(), nxt["RangeLow"].min()))
            cur_high = float(max(cur["RangeHigh"].max(), nxt["RangeHigh"].max()))

            cur.loc[:, "RangeBinLeft"] = cur_left
            cur.loc[:, "RangeBinRight"] = cur_right
            cur.loc[:, "RangeLow"] = cur_low
            cur.loc[:, "RangeHigh"] = cur_high

            # Update range percent
            total_surfs = float(bin_group["TotalCount"].sum())
            sub = bin_group[(bin_group["Bin"] >= cur_left) & (bin_group["Bin"] <= cur_right)]
            range_total = float(sub["TotalCount"].sum())
            cur.loc[:, "RangePercentAllSurfs"] = 100.0 * range_total / total_surfs if total_surfs > 0 else 0.0

            # Recompute per-pair totals and PeakShare within merged peak
            summed = (
                cur.groupby("NormPair", sort=False)
                .agg(
                    PairCount=("PairCount", "sum"),
                    CurvMean=("CurvMean", "mean"),
                    CurvStd=("CurvStd", "mean"),
                    PairFracOfAll=("PairFracOfAll", "max"),  # keep global “surface type %”
                )
                .reset_index()
            )
            denom = float(summed["PairCount"].sum())
            summed["PairFracInRange"] = summed["PairCount"] / denom if denom > 0 else 0.0

            # carry shared peak fields from first row
            base = cur.iloc[0].to_dict()
            keep_cols = { "Model","PeakIndex","PeakBin","PeakCenter",
                          "RangeBinLeft","RangeBinRight","RangeLow","RangeHigh","RangePercentAllSurfs" }

            rows = []
            for _, r in summed.iterrows():
                row = {k: base[k] for k in keep_cols}
                row.update({
                    "NormPair": r["NormPair"],
                    "PairCount": int(r["PairCount"]),
                    "PairFracInRange": float(r["PairFracInRange"]),
                    "PairFracOfAll": float(r["PairFracOfAll"]),
                    "CurvMean": float(r["CurvMean"]),
                    "CurvStd": float(r["CurvStd"]),
                })
                rows.append(row)

            cur = pd.DataFrame(rows)

            j += 1

        blocks.append(cur)
        i = j

    out = pd.concat(blocks, ignore_index=True)

    # Renumber PeakIndex sequentially
    old_ids = sorted(out["PeakIndex"].unique().tolist())
    id_map = {old: new for new, old in enumerate(old_ids, start=1)}
    out["PeakIndex"] = out["PeakIndex"].map(id_map)

    return out



def compute_non_overlapping_peak_bounds(
    bin_centers: np.ndarray,
    y: np.ndarray,
    peak_indices: list[int],
    *,
    use_smooth_for_bounds: bool = True,
    smooth_window: int = 5,
) -> list[tuple[int, int]]:
    """
    For each peak index, compute a non-overlapping (left_idx, right_idx) range.

    Steps:
      1) Inflection bounds around each peak (concavity sign-change).
      2) Valley indices between adjacent peaks (argmin between peaks).
      3) Clip each peak's bounds to the valleys so ranges do not overlap.

    Returns list aligned with peak_indices: [(L0,R0), (L1,R1), ...]
    """
    if len(peak_indices) == 0:
        return []

    # Sort peaks (keep stable order)
    peaks = sorted(int(p) for p in peak_indices)

    # Optional light smoothing for derivative-based inflections only
    y_for_bounds = y.astype(float)
    if use_smooth_for_bounds and smooth_window > 1:
        w = int(smooth_window)
        if w % 2 == 0:
            w += 1
        kernel = np.ones(w, dtype=float) / float(w)
        y_for_bounds = np.convolve(y_for_bounds, kernel, mode="same")

    # 1) Inflection bounds
    inf_bounds: list[tuple[int, int]] = []
    for p in peaks:
        L, R = find_peak_bounds_by_inflection(bin_centers, y_for_bounds, p)
        inf_bounds.append((int(L), int(R)))

    # 2) Valleys between peaks (use raw y, not smoothed, to reflect true minima)
    valleys: list[int] = []
    for k in range(len(peaks) - 1):
        a = peaks[k]
        b = peaks[k + 1]
        if b <= a + 1:
            valleys.append(a)
            continue

        seg = y[a:b + 1]
        v_rel = int(np.argmin(seg))
        v = a + v_rel
        valleys.append(int(v))

    # 3) Clip bounds to enforce non-overlap
    final_bounds: list[tuple[int, int]] = []
    for k, p in enumerate(peaks):
        L_inf, R_inf = inf_bounds[k]

        # Right clip by valley to the next peak
        if k < len(valleys):
            R = min(R_inf, valleys[k])
        else:
            R = R_inf

        # Left clip by previous valley (strictly after it)
        if k > 0:
            L = max(L_inf, valleys[k - 1] + 1)
        else:
            L = L_inf

        # Ensure L <= p <= R
        L = int(np.clip(L, 0, p))
        R = int(np.clip(R, p, len(y) - 1))

        # If clipping makes it empty, force minimal range [p,p]
        if R < L:
            L = p
            R = p

        final_bounds.append((L, R))

    # Return aligned to the sorted peaks
    return final_bounds


def find_peak_bounds_by_inflection(bin_centers: np.ndarray, y: np.ndarray, peak_idx: int):
    """
    Given a peak index (local maximum) in y, return (left_idx, right_idx)
    where boundaries are the nearest NEGATIVE inflection points on each side.

    Negative inflection point criterion:
      - We look for where the second derivative crosses from negative to positive
        (concave down -> concave up) when moving away from the peak.

    Returns indices into bin_centers/y (inclusive bounds).
    """
    n = len(y)
    if n < 5:
        return 0, n - 1

    # Use central differences; assumes roughly uniform spacing
    # y'' at i approximated by y[i+1] - 2y[i] + y[i-1]
    y2 = np.zeros(n, dtype=float)
    y2[1:-1] = y[2:] - 2.0 * y[1:-1] + y[:-2]

    # ----- find left bound: nearest i<peak where y2 crosses (-) -> (+) -----
    left_bound = 0
    for i in range(peak_idx - 1, 1, -1):
        # crossing between i-1 and i
        if y2[i - 1] < 0.0 and y2[i] >= 0.0:
            left_bound = i
            break

    # ----- find right bound: nearest i>peak where y2 crosses (-) -> (+) -----
    right_bound = n - 1
    for i in range(peak_idx + 1, n - 2):
        # crossing between i and i+1
        if y2[i] < 0.0 and y2[i + 1] >= 0.0:
            right_bound = i + 1
            break

    # Sanity: ensure bounds wrap the peak
    left_bound = int(np.clip(left_bound, 0, peak_idx))
    right_bound = int(np.clip(right_bound, peak_idx, n - 1))

    return left_bound, right_bound


def find_peaks(bin_centers, percents, min_height=2.0):
    """
    Simple 1D peak finder: a peak is a local maximum where
    y[i] >= y[i-1] and y[i] >= y[i+1] and y[i] >= min_height.

    Returns a list of indices into bin_centers/percents.
    """
    peaks = []
    y = percents
    n = len(y)

    if n < 3:
        return peaks

    for i in range(1, n - 1):
        if y[i] >= y[i - 1] and y[i] >= y[i + 1] and y[i] >= min_height:
            peaks.append(i)

    return peaks


def write_peaks_txt(peaks_df: pd.DataFrame, txt_path: str) -> None:
    """
    Human-readable peak summary using inflection-bounded ranges.

    Structure:
        MODEL
          Peak N
            PeakBin
            PeakCenter
            Range [low, high]
            Range % of total
            Contributing pairs (sorted by PairCount)
    """
    with open(txt_path, "w", encoding="utf-8") as f:

        f.write("SURFACE CURVATURE PEAK SUMMARY (Inflection Defined)\n")
        f.write("=" * 80 + "\n\n")

        if peaks_df is None or peaks_df.empty:
            f.write("No peaks detected.\n")
            return

        # Group by model
        for model, model_df in peaks_df.groupby("Model", sort=False):

            f.write(f"\nMODEL: {model}\n")
            f.write("-" * 80 + "\n")

            # Group by peak
            for peak_id, peak_df in model_df.groupby("PeakIndex", sort=True):

                row0 = peak_df.iloc[0]

                peak_bin = int(row0["PeakBin"])
                peak_center = float(row0["PeakCenter"])

                range_left = int(row0["RangeBinLeft"])
                range_right = int(row0["RangeBinRight"])
                range_low = float(row0["RangeLow"])
                range_high = float(row0["RangeHigh"])
                range_percent = float(row0["RangePercentAllSurfs"])

                f.write(
                    f"\n  Peak {peak_id}\n"
                    f"    PeakBin        : {peak_bin}\n"
                    f"    PeakCenter     : {peak_center:.6f}\n"
                    f"    RangeBins      : {range_left} – {range_right}\n"
                    f"    RangeCurvature : [{range_low:.6f}, {range_high:.6f}]\n"
                    f"    Range%AllSurfs : {range_percent:.3f}%\n"
                )

                f.write("    Contributing Pairs:\n")

                # Sort by dominance in the peak-range
                peak_df = peak_df.sort_values("PairCount", ascending=False)

                for _, row in peak_df.iterrows():

                    pair = row["NormPair"]
                    count = int(row["PairCount"])
                    mu = float(row["CurvMean"])
                    sigma = float(row["CurvStd"])

                    peak_share = 100.0 * float(row["PairFracInRange"])
                    surface_type_pct = 100.0 * float(
                        row["PairFracOfAll"]) if "PairFracOfAll" in peak_df.columns else 0.0

                    f.write(
                        f"      {pair:10s} | "
                        f"Count={count:6d} | "
                        f"PeakShare={peak_share:6.2f}% | "
                        f"SurfaceType%={surface_type_pct:6.2f}% | "
                        f"μ={mu:.5f} | "
                        f"σ={sigma:.5f}\n"
                    )

        f.write("\nEND OF SUMMARY\n")


def normalize_atom_label(label: str) -> str:
    """
    Normalize individual atom labels so that things like
    H2'1, H2'2 -> H2'; NE/NE2 -> NE; NH1/NH2 -> NH; OG/OG1 -> OG.

    Rule: strip trailing digits at the end of the label.
    This does not touch labels like C2' (ends with a prime, not a digit).
    """
    label = label.strip()
    m = re.match(r"^(.*?)(\d+)$", label)
    if m:
        return m.group(1)  # drop trailing digits
    return label


def normalize_pair(pair: str) -> str:
    """
    Normalize a pair label and canonicalize order so that A-B == B-A.

    Example:
        'C2'-H2'2' -> 'C2'-H2'' (then sorted lexicographically)
    """
    left, right = pair.split("-")
    left_n = normalize_atom_label(left)
    right_n = normalize_atom_label(right)

    if left_n <= right_n:
        return f"{left_n}-{right_n}"
    else:
        return f"{right_n}-{left_n}"


def main(min_height=2.0, top_pairs=5, top_residues=5):
    """
    Analyze a *_curvature_bins.csv file, identify curvature peaks,
    and assign dominant surface types (normalized pairs) to each peak,
    including:
      * fraction of surfaces in the bin,
      * weighted mean/std of curvature,
      * residue breakdown (Residue-of-origin -> count).

    Parameters
    ----------
    min_height : float
        Minimum percent of total surfaces for a bin to be considered a peak.
    top_pairs : int
        Number of top normalized pairs to report per peak.
    top_residues : int
        Maximum number of residue types to list in the printed breakdown.
    """
    # File chooser for a single *_curvature_bins.csv
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)

    csv_path = filedialog.askopenfilename(
        title="Select *_curvature_bins.csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not csv_path:
        print("No file selected.")
        return

    df = pd.read_csv(csv_path)

    required_cols = {"Bin", "BinLow", "BinHigh", "Pair", "Count", "CurvMean", "CurvStd"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"Expected columns at least {required_cols}, got {set(df.columns)} in {csv_path}"
        )

    # Residue columns are needed for residue breakdown
    has_residues = {"Residue1", "Residue2"}.issubset(df.columns)

    # Add normalized pair column: CB-HB1, CB-HB2 -> CB-HB, etc.
    df["NormPair"] = df["Pair"].astype(str).apply(normalize_pair)

    # Total surfaces of each NormPair across the entire model/file
    pair_totals_all = (
        df.groupby("NormPair")["Count"]
        .sum()
        .astype(float)
        .to_dict()
    )

    # Aggregate per bin (total count per bin)
    bin_group = (
        df.groupby("Bin")
        .agg(
            BinLow=("BinLow", "first"),
            BinHigh=("BinHigh", "first"),
            TotalCount=("Count", "sum"),
        )
        .reset_index()
        .sort_values("Bin")
    )
    # --- Force a complete bin grid (fill missing bins with zeros) ---
    # Determine full bin id range
    bin_min = int(df["Bin"].min())
    bin_max = int(df["Bin"].max())
    full_bins = pd.DataFrame({"Bin": np.arange(bin_min, bin_max + 1, dtype=int)})

    # Merge to force missing bins to appear
    bin_group = full_bins.merge(bin_group, on="Bin", how="left")

    # Fill missing bin geometry and counts
    bin_group["TotalCount"] = bin_group["TotalCount"].fillna(0).astype(int)

    # BinLow/BinHigh may be missing for empty bins; infer using typical bin width
    # (assumes uniform binning, which your files appear to use)
    if bin_group["BinLow"].isna().any() or bin_group["BinHigh"].isna().any():
        # Use existing non-null bin edges to infer spacing
        edges = bin_group.dropna(subset=["BinLow", "BinHigh"]).sort_values("Bin")
        if len(edges) >= 2:
            step = float(edges["BinLow"].iloc[1] - edges["BinLow"].iloc[0])
        elif len(edges) == 1:
            step = float(edges["BinHigh"].iloc[0] - edges["BinLow"].iloc[0])
        else:
            step = 0.0

        # Fill BinLow by linear interpolation on Bin index
        bin_group["BinLow"] = bin_group["BinLow"].interpolate(method="linear")
        bin_group["BinHigh"] = bin_group["BinHigh"].interpolate(method="linear")

        # If still NaN at ends, extrapolate using step
        # Low end
        first_valid = bin_group["BinLow"].first_valid_index()
        if first_valid is not None:
            for k in range(first_valid - 1, -1, -1):
                bin_group.loc[k, "BinLow"] = float(bin_group.loc[k + 1, "BinLow"]) - step
                bin_group.loc[k, "BinHigh"] = float(bin_group.loc[k, "BinLow"]) + step

        # High end
        last_valid = bin_group["BinLow"].last_valid_index()
        if last_valid is not None:
            for k in range(last_valid + 1, len(bin_group)):
                bin_group.loc[k, "BinLow"] = float(bin_group.loc[k - 1, "BinLow"]) + step
                bin_group.loc[k, "BinHigh"] = float(bin_group.loc[k, "BinLow"]) + step

    # Recompute centers now that edges are complete
    bin_group["BinCenter"] = 0.5 * (bin_group["BinLow"].astype(float) + bin_group["BinHigh"].astype(float))

    # Recompute Percent (zeros stay zeros)
    total_surfs = bin_group["TotalCount"].sum()
    bin_group["Percent"] = 0.0
    if total_surfs > 0:
        bin_group["Percent"] = 100.0 * bin_group["TotalCount"] / float(total_surfs)

    total_surfs = bin_group["TotalCount"].sum()
    bin_group["Percent"] = 100.0 * bin_group["TotalCount"] / float(total_surfs)

    # Bin centers for plotting
    bin_group["BinCenter"] = 0.5 * (bin_group["BinLow"] + bin_group["BinHigh"])

    # Detect peaks
    bin_centers = bin_group["BinCenter"].to_numpy()
    percents = bin_group["Percent"].to_numpy()

    # Derive model name from file name: Ncp_curvature_bins.csv -> Ncp
    base_name = os.path.basename(csv_path)
    model_name = base_name.replace("_curvature_bins.csv", "").replace(".csv", "")

    print(f"\nModel: {model_name}")
    print(f"Total surfaces (from CSV counts): {int(total_surfs)}")
    
    # Find peaks above threshold (for saving/printing only)
    peak_indices = find_peaks(bin_centers, percents, min_height=min_height)
    if not peak_indices:
        print(f"No peaks found above {min_height:.2f}% in {os.path.basename(csv_path)}.")
        print("Plotting all data regardless of peak threshold...")
    else:
        print(f"Detected {len(peak_indices)} peaks (min height = {min_height:.2f}%)")

    # Compute NON-OVERLAPPING inflection+valley bounds for all peaks
    peak_bounds = compute_non_overlapping_peak_bounds(
        bin_centers,
        percents,
        peak_indices,
        use_smooth_for_bounds=True,
        smooth_window=5,
    )

    # Ensure peak_indices and bounds are aligned in sorted order
    peak_indices_sorted = sorted(int(p) for p in peak_indices)

    # Prepare peak summary rows for CSV (only for peaks above threshold)
    peak_rows = []

    for peak_idx_num, (i, (left_i, right_i)) in enumerate(zip(peak_indices_sorted, peak_bounds), start=1):

        # Use the precomputed NON-OVERLAPPING bounds
        left_i = int(left_i)
        right_i = int(right_i)

        # Peak bin (highest point)
        peak_bin_id = int(bin_group.iloc[i]["Bin"])

        peak_center = float(bin_group.iloc[i]["BinCenter"])

        # Range bounds (inclusive bin IDs)
        bin_id_left = int(bin_group.iloc[left_i]["Bin"])
        bin_id_right = int(bin_group.iloc[right_i]["Bin"])

        range_low = float(bin_group.iloc[left_i]["BinLow"])
        range_high = float(bin_group.iloc[right_i]["BinHigh"])

        # Total surfaces in the full dataset (already computed earlier; keep it)
        # total_surfs = ...

        # Range totals
        range_total = int(bin_group.iloc[left_i:right_i + 1]["TotalCount"].sum())
        range_percent_all = 100.0 * float(range_total) / float(total_surfs) if total_surfs > 0 else 0.0

        # All surfaces contributing to the PEAK RANGE (not just the peak bin)
        raw_range_df = df[(df["Bin"] >= bin_id_left) & (df["Bin"] <= bin_id_right)].copy()
        if raw_range_df.empty or range_total <= 0:
            continue

        # For each NormPair, we may have multiple underlying raw Pair/Residue combinations across the RANGE.
        groups = []
        for norm_pair, grp in raw_range_df.groupby("NormPair"):
            counts = grp["Count"].to_numpy(dtype=float)
            means = grp["CurvMean"].to_numpy(dtype=float)
            stds = grp["CurvStd"].to_numpy(dtype=float)

            N = float(counts.sum())
            if N <= 0:
                continue

            # Weighted mean curvature across the RANGE for this pair
            mu = float(np.sum(counts * means) / N)

            # Weighted variance: E[x^2] - (E[x])^2, where E[x^2] uses std^2 + mean^2
            ex2 = float(np.sum(counts * (stds ** 2 + means ** 2)) / N)
            var = max(ex2 - mu ** 2, 0.0)
            sigma = float(np.sqrt(var))

            frac_in_range = float(N) / float(range_total)

            total_pair_all = float(pair_totals_all.get(norm_pair, 0.0))
            frac_of_all_pair = (float(N) / total_pair_all) if total_pair_all > 0 else 0.0

            groups.append(
                {
                    "NormPair": norm_pair,
                    "PairCount": int(N),
                    "FracInRange": frac_in_range,
                    "FracOfAllPair": frac_of_all_pair,
                    "CurvMean": mu,
                    "CurvStd": sigma,
                }
            )

        grouped = pd.DataFrame(groups)
        if grouped.empty:
            continue

        grouped = grouped.sort_values("PairCount", ascending=False)

        print(
            f"\nPeak {peak_idx_num}: PeakBin {peak_bin_id}, "
            f"center = {peak_center:.5f}, "
            f"range = [{range_low:.5f}, {range_high:.5f}] "
            f"(bins {bin_id_left}–{bin_id_right}), "
            f"range % = {range_percent_all:.3f}%"
        )

        # Print top N normalized pairs with mean/std and residue breakdown (RANGE-based)
        for _, row in grouped.head(top_pairs).iterrows():
            pair = row["NormPair"]
            count = int(row["PairCount"])
            frac_range = 100.0 * float(row["FracInRange"])
            mu = float(row["CurvMean"])
            sigma = float(row["CurvStd"])

            residue_str = ""
            if has_residues:
                # All surfaces for this NormPair in this PEAK RANGE
                grp_np = raw_range_df[raw_range_df["NormPair"] == pair].copy()

                from collections import Counter
                counter = Counter()

                # Weight each residue by the number of surfaces (Count)
                for _, r in grp_np.iterrows():
                    cnt = int(r["Count"])

                    res1 = r.get("Residue1", None)
                    res2 = r.get("Residue2", None)

                    # Clean up NaNs
                    res1 = None if pd.isna(res1) else str(res1)
                    res2 = None if pd.isna(res2) else str(res2)

                    if res1 and res2:
                        if res1 == res2:
                            counter[res1] += cnt
                        else:
                            counter[res1] += cnt
                            counter[res2] += cnt
                    elif res1:
                        counter[res1] += cnt
                    elif res2:
                        counter[res2] += cnt

                if counter:
                    sorted_res = sorted(counter.items(), key=lambda x: (-x[1], x[0]))
                    parts = [f"{res} ({cnt})" for res, cnt in sorted_res[:top_residues]]
                    if len(sorted_res) > top_residues:
                        parts.append("...")
                    residue_str = " - " + ", ".join(parts)

            print(
                f"  {pair:8s} : {count:6d} "
                f"({frac_range:6.2f}% of this peak-range) "
                f"[μ = {mu:.5f}, σ = {sigma:.5f}]{residue_str}"
            )

        # Add all normalized pairs for this peak to output table (RANGE-based)
        for _, row in grouped.iterrows():
            peak_rows.append(
                {
                    "Model": model_name,
                    "PeakIndex": peak_idx_num,

                    # Peak anchor (highest bin)
                    "PeakBin": peak_bin_id,
                    "PeakCenter": peak_center,

                    # Inflection-bounded range
                    "RangeBinLeft": bin_id_left,
                    "RangeBinRight": bin_id_right,
                    "RangeLow": range_low,
                    "RangeHigh": range_high,
                    "RangePercentAllSurfs": range_percent_all,

                    # Pair contribution within the range
                    "NormPair": row["NormPair"],
                    "PairCount": int(row["PairCount"]),
                    "PairFracInRange": float(row["FracInRange"]),
                    "PairFracOfAll": float(row["FracOfAllPair"]),
                    "CurvMean": float(row["CurvMean"]),
                    "CurvStd": float(row["CurvStd"]),

                }
            )
        # -----------------------------------------------------------------

    # Export peak summary CSV next to the input file (only if peaks found above threshold)
    if peak_rows:
        peaks_df = pd.DataFrame(peak_rows)
        out_path = os.path.join(
            os.path.dirname(csv_path), f"{model_name}_curvature_peaks.csv"
        )
        peaks_df = merge_consecutive_almost_identical_peaks(
            peaks_df,
            bin_group=bin_group,
            topk=12,
            jaccard_thresh=0.85,
            count_rel_tol=0.05,
            min_pair_count=5,
        )

        peaks_df.to_csv(out_path, index=False)
        print(f"\nPeak summary exported to:\n  {out_path}")
        txt_path = os.path.splitext(out_path)[0] + ".txt"
        write_peaks_txt(peaks_df, txt_path)
        print(f"Peak summary text exported to:\n  {txt_path}")

    else:
        print(f"\nNo peaks above {min_height:.2f}% threshold to export.")

    # -------------------- PLOTTING (4A SCHEMA) --------------------
    # Optional: skip the first N bins (matches 4A’s "skip_bins" behavior)
    skip_bins = 10  # set to 10 if you want the same "ignore left tail" look as 4A
    start_idx = int(max(skip_bins, 0))

    x_segment = bin_centers[start_idx:]
    y_segment = percents[start_idx:]

    # 4A-style label + color
    model_upper = str(model_name).upper()
    letter = LETTER_CODES.get(model_upper, model_name)
    color = COLOR_MAP.get(model_upper, "k")

    xs = [x_segment]
    ys = [y_segment]
    labels = [f"{letter}"]
    colors = [color]
    # for i in range(len(xs[0])):
    #     print(xs[0][i], ys[0][i])
    line_plot(
        xs=xs,
        ys=ys,
        title="Surface Curvatures",
        x_label="Curvature",
        y_label="% of Surfs",
        legend_title="Model",
        labels=labels,
        title_size=35,
        x_label_size=30,
        y_label_size=30,
        colors=colors,
        tick_val_size=30,
        legend_orientation="Vertical",
        legend_label_size=10,
        tight_layout=True,
        Show=False,
        ylim=[0, 2.5],
        xlim=[0.025, 1.4]
    )

    # Only annotate peaks above threshold on the plot (keep your logic)
    # Only annotate peaks above the threshold on the plot
    ax = plt.gca()

    for peak_idx_num, i in enumerate(peak_indices, start=1):
        # If we skipped bins, shift indices for annotation
        if i < start_idx:
            continue

        x = float(bin_centers[i])
        y = float(percents[i])

        ax.text(
            x,
            y + 0.15,  # vertical offset in data units; tweak if needed
            str(peak_idx_num),
            ha="center",
            va="bottom",
            fontsize=10,
            rotation=0,
        )

    plt.show()
    # --------------------------------------------------------------


if __name__ == "__main__":
    # Example: min_height=2% of surfaces, top 5 normalized pairs per peak,
    # and top 5 residues in the breakdown.
    main(min_height=0.120, top_pairs=10, top_residues=20)
