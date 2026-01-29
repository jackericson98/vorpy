import os
import sys
from collections import Counter

import tkinter as tk
from tkinter import filedialog

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
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
    "BSA_20": "L"
}


COLOR_MAP = {
    "NA5": "#000000",           # Black
    "EDTA": "#0077FF",          # Bright Blue
    "DB1976": "#FF7F00",        # Vivid Orange
    "HAIRPIN": "#008F00",       # Strong Green
    "CAMBRIN": "#00C8FF",       # Cyan
    "BDNA": "#B300FF",          # Vivid Purple
    "HAMMERHEAD": "#FF0000",    # Bright Red
    "P53TET": "#FF00AA",        # Magenta
    "T4LP": "#00FF00",          # Neon Green
    "STREPTAVIDIN": "#A52A2A",  # Brown
    "NCP": "#8C8C8C",           # Neutral Gray
    "BSA": "#FFD700",           # Gold
    "BSA_20": "#FFD700"         # Gold
}


def _parse_ball_indices(balls_value):
    """
    Normalize the Balls entry to a pair of integer indices.

    Handles cases like:
      - '12,34'
      - '12 34'
      - [12, 34], (12, 34), np.array([...])
    """
    if isinstance(balls_value, str):
        # e.g. "12,34" or "12 34"
        if "," in balls_value:
            parts = balls_value.replace(" ", "").split(",")
        else:
            parts = balls_value.split()
        if len(parts) != 2:
            raise ValueError(f"Unexpected Balls string format: {balls_value!r}")
        return int(parts[0]), int(parts[1])

    if isinstance(balls_value, (list, tuple, np.ndarray)):
        if len(balls_value) != 2:
            raise ValueError(f"Unexpected Balls sequence length: {balls_value!r}")
        return int(balls_value[0]), int(balls_value[1])

    raise TypeError(f"Unsupported Balls type: {type(balls_value)} value={balls_value!r}")


def _load_pdb_index_map(pdb_path):
    """
    Parse a PDB file and build mappings:

      index 0 -> first ATOM/HETATM (#1 in the PDB)
      index 1 -> second ATOM/HETATM (#2), ...

    Returns:
      index_to_name:    dict[int -> atom name]
      index_to_coord:   dict[int -> (x, y, z)]
      index_to_resname: dict[int -> residue name (3-letter)]
    """
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    names = []
    coords = []
    resnames = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Standard PDB columns
                atom_name = line[12:16].strip()     # atom name
                res_name = line[17:20].strip()      # residue name
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                names.append(atom_name)
                coords.append((x, y, z))
                resnames.append(res_name)

    index_to_name = {i: names[i] for i in range(len(names))}
    index_to_coord = {i: coords[i] for i in range(len(coords))}
    index_to_resname = {i: resnames[i] for i in range(len(resnames))}
    return index_to_name, index_to_coord, index_to_resname



def main(excluded=None, bins=300, skip_bins=10):
    # Check the excluded list
    if excluded is None:
        excluded = []

    # Choose the directory containing molecular subfolders
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)
    # Get the folder
    folder = filedialog.askdirectory(title="Choose Molecular Data Directory")
    file_paths = []
    for subfolder in os.listdir(folder):
        if subfolder and subfolder[0] in excluded:
            continue

        model_folder = os.path.join(folder, subfolder)
        aw_logs = os.path.join(model_folder, "aw", "aw_logs.csv")

        # Optional safety check
        if not os.path.exists(aw_logs):
            print(f"Skipping {model_folder}: no aw_logs.csv found.")
            continue

        file_paths.append({
            "aw_logs": aw_logs,
            "model_folder": model_folder,
        })

    if not file_paths:
        print("No files selected.")
        return

    # model_data[model_name] = {
    #   "vals": np.ndarray of curvature values (finite only),
    #   "surfs": cleaned surfs DataFrame aligned with vals,
    #   "index_to_name": dict[int -> str],
    #   "index_to_coord": dict[int -> (x, y, z)],
    #   "curvature_col": str,
    #   "model_folder": str,
    # }
    model_data = {}

    for entry in file_paths:
        log_path = entry["aw_logs"]
        model_folder = entry["model_folder"]

        info = read_logs2(log_path, all_=False, surfs=True, balls=True)

        surfs = info["surfs"]
        model_name = info["data"]["name"].upper()
        print(model_name)

        print(f"\nProcessing file: {os.path.basename(log_path)}")
        print(f"Model name from info['data']['name']: {model_name}")
        print(f"Surface records: {len(surfs)}")
        print(f"DataFrame shape: {surfs.shape}")

        # Identify curvature column
        curvature_col = None
        for candidate in ["Mean Curvature", "Mean_Curvature", "mean_curvature"]:
            if candidate in surfs.columns:
                curvature_col = candidate
                break

        if curvature_col is None:
            print("  No mean curvature column found. Columns are:")
            print("  ", surfs.columns.tolist())
            continue

        print(f"  Using curvature column: {curvature_col}")

        vals_full = surfs[curvature_col].to_numpy()
        finite_mask = np.isfinite(vals_full)
        vals = vals_full[finite_mask]

        if vals.size == 0:
            print("  No finite curvature values found; skipping this file.")
            continue

        # Cleaned surfs aligned with finite vals
        surfs_clean = surfs.loc[finite_mask].reset_index(drop=True)

        # PDB path: folder/subfolder/model_name.pdb
        pdb_path = os.path.join(model_folder, f"{model_name}.pdb")

        try:
            index_to_name, index_to_coord, index_to_resname = _load_pdb_index_map(pdb_path)

        except FileNotFoundError as e:
            print(f"  {e}  Skipping this model.")
            continue

        # Store everything for this model
        model_data[model_name] = {
            "vals": vals,
            "surfs": surfs_clean,
            "index_to_name": index_to_name,
            "index_to_coord": index_to_coord,
            "index_to_resname": index_to_resname,
            "curvature_col": curvature_col,
            "model_folder": model_folder,
        }
    if not model_data:
        print("No valid curvature data found in selected files.")
        return

    # Build global bins across all models
    all_vals = np.concatenate([d["vals"] for d in model_data.values()])
    all_vals = all_vals[np.isfinite(all_vals)]

    curv_min = float(np.min(all_vals))
    curv_max = float(np.max(all_vals))
    increments = np.linspace(curv_min, curv_max, bins)

    # Bin centers
    _, bin_edges = np.histogram(all_vals, bins=increments)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    # Optionally skip the first N bins (for plotting only, not for CSV export)
    start_idx = skip_bins
    x_segment = bin_centers[start_idx:]

    xs = []
    ys = []
    labels = []
    colors = []

    # For each model, compute histogram and also print atom-pair groups per bin
    for model_name, data in model_data.items():
        vals = data["vals"]
        surfs_clean = data["surfs"]
        index_to_name = data["index_to_name"]
        model_folder = data["model_folder"]

        index_to_resname = data["index_to_resname"]

        # Surface area per surface
        if "Surface Area" not in surfs_clean.columns:
            raise KeyError(f"'Surface Area' column not found in surfs for {model_name}")
        sa_vals = surfs_clean["Surface Area"].to_numpy()

        # Per-surface residue labels
        res1_labels = np.empty(len(surfs_clean), dtype=object)
        res2_labels = np.empty(len(surfs_clean), dtype=object)

        counts, _ = np.histogram(vals, bins=increments)
        scaled = counts / float(len(vals))  # fraction
        scaled_pct = 100.0 * scaled

        xs.append(x_segment)
        ys.append(scaled_pct[start_idx:])

        letter = LETTER_CODES.get(model_name, "?")
        labels.append(f"{letter} - {model_name}")

        colors.append(COLOR_MAP.get(model_name, "k"))

        # Map each row to its bin index
        bin_indices = np.digitize(vals, increments) - 1
        bin_indices = np.clip(bin_indices, 0, len(increments) - 2)

        pair_labels = np.empty(len(surfs_clean), dtype=object)
        balls_array = surfs_clean["Balls"].to_numpy()

        for row_idx, balls_value in enumerate(balls_array):
            try:
                idx0, idx1 = _parse_ball_indices(balls_value)
            except Exception as e:
                print(f"  Warning: could not parse Balls={balls_value!r} in {model_name}: {e}")
                continue

            if idx0 not in index_to_name or idx1 not in index_to_name:
                print(
                    f"  Warning: atom index out of range for {model_name}: "
                    f"{idx0}, {idx1} (skipping this surface)"
                )
                continue

            name0 = index_to_name[idx0]
            name1 = index_to_name[idx1]
            res0 = index_to_resname.get(idx0)
            res1 = index_to_resname.get(idx1)

            # Sort names so CA-CB == CB-CA; keep residues aligned
            if name0 <= name1:
                pair_label = f"{name0}-{name1}"
                r1, r2 = res0, res1
            else:
                pair_label = f"{name1}-{name0}"
                r1, r2 = res1, res0

            pair_labels[row_idx] = pair_label
            res1_labels[row_idx] = r1
            res2_labels[row_idx] = r2

        print(f"\n=== {model_name} ({letter}) bin details ===")

        # Collect rows for CSV export, now including residue and surface area stats
        # Include ALL bins (skip_bins only affects plotting, not CSV export)
        rows = []

        for b in range(0, len(increments) - 1):
            # Include ALL surfaces in this bin (not just those with valid pair labels)
            mask_bin = (bin_indices == b)
            if not np.any(mask_bin):
                continue

            bin_low = increments[b]
            bin_high = increments[b + 1]

            # Build a small DataFrame for this bin: Pair, Residues, Curvature, Surface Area
            # Use placeholder for surfaces without valid pair labels
            pair_labels_bin = pair_labels[mask_bin].copy()
            pair_labels_bin[pair_labels_bin == None] = "UNKNOWN"
            
            bin_df = pd.DataFrame({
                "Pair": pair_labels_bin,
                "Residue1": res1_labels[mask_bin],
                "Residue2": res2_labels[mask_bin],
                "Curv": vals[mask_bin],
                "SurfArea": sa_vals[mask_bin],
            })

            # Group by atom pair + residue pair
            grouped = bin_df.groupby(["Pair", "Residue1", "Residue2"], dropna=False)

            for (pair, res1, res2), grp in grouped:
                curvs = grp["Curv"].to_numpy()
                areas = grp["SurfArea"].to_numpy()
                if curvs.size == 0:
                    continue

                count = int(curvs.size)
                mean_c = float(np.mean(curvs))
                std_c = float(np.std(curvs, ddof=0))
                mean_sa = float(np.mean(areas))
                std_sa = float(np.std(areas, ddof=0))

                rows.append({
                    "Bin": b + 1,
                    "BinLow": bin_low,
                    "BinHigh": bin_high,
                    "Pair": pair,
                    "Residue1": res1,
                    "Residue2": res2,
                    "Count": count,
                    "CurvMean": mean_c,
                    "CurvStd": std_c,
                    "SurfAreaMean": mean_sa,
                    "SurfAreaStd": std_sa,
                })

        # ---------- CSV EXPORT ----------
        import pandas as _pd

        df_out = _pd.DataFrame(rows)
        csv_path = os.path.join(model_folder, f"{model_name}_curvature_bins.csv")

        df_out.to_csv(csv_path, index=False)
        print(f"Exported CSV for {model_name}: {csv_path}")

    # Plot the distributions
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
        legend_orientation="horizontal",
        tight_layout=False,
    )

    plt.show()



if __name__ == "__main__":
    # Exclude models by first letter of subfolder, as before
    main(["A", "B", "C"], bins=150, skip_bins=10)
