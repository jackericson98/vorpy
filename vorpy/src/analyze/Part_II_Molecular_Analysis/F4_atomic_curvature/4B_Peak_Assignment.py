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
    Write a human-readable peak summary alongside the CSV output.

    Works with both schemas:
      - 4C combined peaks output (has: Model, Pair, PeakRank, Bin, BinCenter, BinLow, BinHigh, Height, Prominence, ...)
      - 4B peak assignment output (has: Model, NormPair, PeakIndex, Bin, BinCenter, BinLow, BinHigh, CurvMean, CurvStd, PairCount, ...)
    """
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("SURFACE CURVATURE PEAK SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        if peaks_df is None or peaks_df.empty:
            f.write("No peaks detected.\n")
            return

        cols = set(peaks_df.columns)

        # Figure out which columns exist
        model_col = "Model" if "Model" in cols else None
        pair_col = "Pair" if "Pair" in cols else ("NormPair" if "NormPair" in cols else None)
        peak_col = "PeakRank" if "PeakRank" in cols else ("PeakIndex" if "PeakIndex" in cols else None)

        # Helpers to safely fetch values
        def has(c: str) -> bool:
            return c in cols

        def fmt_float(v, nd=6) -> str:
            try:
                return f"{float(v):.{nd}f}"
            except Exception:
                return str(v)

        def fmt_int(v) -> str:
            try:
                return str(int(v))
            except Exception:
                return str(v)

        # Group by model if possible
        model_groups = [(None, peaks_df)]
        if model_col is not None:
            model_groups = list(peaks_df.groupby(model_col, sort=False))

        for model, model_df in model_groups:
            if model is not None:
                f.write(f"\nMODEL: {model}\n")
                f.write("-" * 70 + "\n")

            # Group by pair (Pair or NormPair) if possible
            pair_groups = [(None, model_df)]
            if pair_col is not None:
                pair_groups = list(model_df.groupby(pair_col, sort=False))

            for pair, pair_df in pair_groups:
                if pair is not None:
                    label = "Pair" if pair_col == "Pair" else "NormPair"
                    f.write(f"\n  {label}: {pair}\n")

                # Sort by peak index/rank if possible
                if peak_col is not None and peak_col in pair_df.columns:
                    pair_df = pair_df.sort_values(peak_col, kind="mergesort")

                for _, row in pair_df.iterrows():
                    # Common bin info
                    bin_id = fmt_int(row["Bin"]) if has("Bin") else "?"
                    center = fmt_float(row["BinCenter"], 6) if has("BinCenter") else "?"
                    low = fmt_float(row["BinLow"], 6) if has("BinLow") else "?"
                    high = fmt_float(row["BinHigh"], 6) if has("BinHigh") else "?"

                    peak_id = fmt_int(row[peak_col]) if (peak_col is not None and has(peak_col)) else "?"

                    f.write(
                        f"    Peak {peak_id} | Bin {bin_id} | "
                        f"Center = {center} | Range = [{low}, {high}]\n"
                    )

                    # 4C-style metrics
                    if has("Height") or has("Prominence"):
                        parts = []
                        if has("Height"):
                            parts.append(f"Height = {fmt_float(row['Height'], 3)}")
                        if has("Prominence"):
                            parts.append(f"Prominence = {fmt_float(row['Prominence'], 3)}")
                        f.write(f"      " + " | ".join(parts) + "\n")

                    if has("LeftBaseCenter") and has("RightBaseCenter"):
                        f.write(
                            f"      Bases: {fmt_float(row['LeftBaseCenter'], 6)} → "
                            f"{fmt_float(row['RightBaseCenter'], 6)}\n"
                        )

                    # 4B-style metrics
                    extra = []
                    if has("PairCount"):
                        extra.append(f"PairCount = {fmt_int(row['PairCount'])}")
                    if has("PairFracInBin"):
                        extra.append(f"FracInBin = {fmt_float(100.0 * float(row['PairFracInBin']), 2)}%")
                    if has("BinPercentAllSurfs"):
                        extra.append(f"Bin%AllSurfs = {fmt_float(row['BinPercentAllSurfs'], 3)}%")
                    if has("CurvMean"):
                        extra.append(f"μ = {fmt_float(row['CurvMean'], 6)}")
                    if has("CurvStd"):
                        extra.append(f"σ = {fmt_float(row['CurvStd'], 6)}")

                    if extra:
                        f.write("      " + " | ".join(extra) + "\n")

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

    # Prepare peak summary rows for CSV (only for peaks above threshold)
    peak_rows = []

    # Only process and save peaks above the threshold
    for peak_idx_num, i in enumerate(peak_indices, start=1):
        bin_id = int(bin_group.iloc[i]["Bin"])
        bin_low = float(bin_group.iloc[i]["BinLow"])
        bin_high = float(bin_group.iloc[i]["BinHigh"])
        bin_center = float(bin_group.iloc[i]["BinCenter"])
        bin_total = int(bin_group.iloc[i]["TotalCount"])
        bin_percent = float(bin_group.iloc[i]["Percent"])

        # Get all rows for this bin (includes Pair, Count, CurvMean, CurvStd, NormPair, Residue1/2)
        raw_bin_df = df[df["Bin"] == bin_id].copy()

        # For each NormPair, we may have multiple underlying raw Pair/Residue combinations.
        groups = []
        for norm_pair, grp in raw_bin_df.groupby("NormPair"):
            counts = grp["Count"].to_numpy(dtype=float)
            means = grp["CurvMean"].to_numpy(dtype=float)
            stds = grp["CurvStd"].to_numpy(dtype=float)

            N = counts.sum()
            if N <= 0:
                continue

            # Weighted mean curvature
            mu = float(np.sum(counts * means) / N)

            # Weighted variance: E[x^2] - (E[x])^2
            ex2 = np.sum(counts * (stds**2 + means**2)) / N
            var = max(ex2 - mu**2, 0.0)
            sigma = float(np.sqrt(var))

            frac_in_bin = N / float(bin_total)

            groups.append({
                "NormPair": norm_pair,
                "PairCount": int(N),
                "FracInBin": frac_in_bin,
                "CurvMean": mu,
                "CurvStd": sigma,
            })

        grouped = pd.DataFrame(groups)
        if grouped.empty:
            continue

        grouped = grouped.sort_values("PairCount", ascending=False)

        print(
            f"\nPeak {peak_idx_num}: Bin {bin_id}, "
            f"center = {bin_center:.5f}, range = [{bin_low:.5f}, {bin_high:.5f}], "
            f"bin % = {bin_percent:.3f}%"
        )

        # Print top N normalized pairs with mean/std and residue breakdown
        for _, row in grouped.head(top_pairs).iterrows():
            pair = row["NormPair"]
            count = int(row["PairCount"])
            frac_bin = 100.0 * float(row["FracInBin"])
            mu = row["CurvMean"]
            sigma = row["CurvStd"]

            residue_str = ""
            if has_residues:
                # All surfaces for this NormPair in this bin
                grp_np = raw_bin_df[raw_bin_df["NormPair"] == pair].copy()

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
                            # Intra-residue surface: count the residue once per surface
                            counter[res1] += cnt
                        else:
                            # Inter-residue surface: both residues participate
                            counter[res1] += cnt
                            counter[res2] += cnt
                    elif res1:
                        counter[res1] += cnt
                    elif res2:
                        counter[res2] += cnt

                if counter:
                    # Optionally, you can sanity-check:
                    # total_res_counts = sum(counter.values())
                    # print(f"  [DEBUG] {pair}: PairCount={count}, residue-participation={total_res_counts}")

                    # sort by count desc, then alphabetically
                    sorted_res = sorted(counter.items(), key=lambda x: (-x[1], x[0]))

                    # Build something like "M (14), I (12), T (10) ..."
                    parts = [
                        f"{res} ({cnt})" for res, cnt in sorted_res[:top_residues]
                    ]
                    if len(sorted_res) > top_residues:
                        parts.append("...")
                    residue_str = " - " + ", ".join(parts)

            print(
                f"  {pair:8s} : {count:6d} "
                f"({frac_bin:6.2f}% of this bin) "
                f"[μ = {mu:.5f}, σ = {sigma:.5f}]{residue_str}"
            )

        # Add all normalized pairs for this peak to output table
        for _, row in grouped.iterrows():
            peak_rows.append(
                {
                    "Model": model_name,
                    "PeakIndex": peak_idx_num,
                    "Bin": bin_id,
                    "BinCenter": bin_center,
                    "BinLow": bin_low,
                    "BinHigh": bin_high,
                    "BinPercentAllSurfs": bin_percent,
                    "NormPair": row["NormPair"],
                    "PairCount": int(row["PairCount"]),
                    "PairFracInBin": float(row["FracInBin"]),
                    "CurvMean": float(row["CurvMean"]),
                    "CurvStd": float(row["CurvStd"]),
                }
            )

    # Export peak summary CSV next to the input file (only if peaks found above threshold)
    if peak_rows:
        peaks_df = pd.DataFrame(peak_rows)
        out_path = os.path.join(
            os.path.dirname(csv_path), f"{model_name}_curvature_peaks.csv"
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
