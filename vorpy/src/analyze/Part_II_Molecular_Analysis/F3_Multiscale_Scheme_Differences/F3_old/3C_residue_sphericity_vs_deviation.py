import os
import sys
import warnings
import numpy as np
import tkinter as tk
from tkinter import filedialog
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"PIL|matplotlib\.backends\._backend_tk")

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.plot_templates.scatter import scatter
from vorpy.src.analyze.tools.batch.get_files import get_all_files
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.calculations.calcs import calc_sphericity


# ---------------------------- residue SA helper (yours) ----------------------------

PRINT_EVERY_RES = 100   # residue-level print frequency


def get_res_sa(logs):
    """
    Aggregate surface area by residue from logs['surfs'], counting only
    outer surfaces (where exactly one ball belongs to the residue).
    """
    res_data = {}
    res_key_dict = {}

    for _, atom in logs['atoms'].iterrows():
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
        if res_key not in res_data:
            res_data[res_key] = [atom['Index']]
        else:
            res_data[res_key].append(atom['Index'])
        res_key_dict[atom['Index']] = res_key

    res_sa_dict = {}

    for _, surf in logs['surfs'].iterrows():
        ball_1, ball_2 = surf['Balls']

        if ball_1 in res_key_dict:
            res_key = res_key_dict[ball_1]
        elif ball_2 in res_key_dict:
            res_key = res_key_dict[ball_2]
        else:
            continue

        if res_key not in res_sa_dict:
            res_sa_dict[res_key] = 0.0

        both_in_same_res = (ball_1 in res_data[res_key] and ball_2 in res_data[res_key])
        if both_in_same_res:
            continue

        res_sa_dict[res_key] += surf['Surface Area']

    return res_sa_dict


# ---------------------------- point builder + plotting wrapper ----------------------------

def build_residue_sphericity_vs_deviation_points(
    folder=None,
    exclude_keys=None,
    scheme="pow",                 # "pow" or "prm"
    deviation_kind="vol",         # "vol" or "sa"
    use_aw_sphericity=True,       # x-axis sphericity from AW (recommended)
    max_percent_diff=None,
):
    """
    Returns:
        x_vals, y_vals, records
    where each point is a residue across all systems.
    """

    outliers = []

    if exclude_keys is None:
        exclude_keys = []

    if folder is None:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()

    scheme = str(scheme).strip().lower()
    if scheme not in {"pow", "prm"}:
        raise ValueError("scheme must be 'pow' or 'prm'")

    deviation_kind = str(deviation_kind).strip().lower()
    if deviation_kind not in {"vol", "sa"}:
        raise ValueError("deviation_kind must be 'vol' or 'sa'")

    x_vals = []
    y_vals = []
    records = []

    subfolders = sorted(os.listdir(folder))
    n_sys = len(subfolders)

    for i_sys, subfolder in enumerate(subfolders, start=1):
        print(f"[{i_sys:02d}/{n_sys}] Processing system: {subfolder}")

        my_key = subfolder.split('_')[0]
        if my_key in exclude_keys:
            continue

        sub_path = os.path.join(folder, subfolder)

        try:
            aw_logs = read_logs2(os.path.join(sub_path, 'aw', 'aw_logs.csv'), all_=False, balls=True, surfs=True)
            pow_logs = read_logs2(os.path.join(sub_path, 'pow', 'pow_logs.csv'), all_=False, balls=True, surfs=True)
            prm_logs = read_logs2(os.path.join(sub_path, 'prm', 'prm_logs.csv'), all_=False, balls=True, surfs=True)
        except FileNotFoundError:
            continue

        aw_res_sa = get_res_sa(aw_logs)
        pow_res_sa = get_res_sa(pow_logs)
        prm_res_sa = get_res_sa(prm_logs)

        res_data = {}

        for _, atom in aw_logs['atoms'].iterrows():
            idx = atom['Index']

            pow_match = pow_logs['atoms'].loc[pow_logs['atoms']['Index'] == idx]
            prm_match = prm_logs['atoms'].loc[prm_logs['atoms']['Index'] == idx]

            if pow_match.empty or prm_match.empty:
                continue

            pow_atom = pow_match.iloc[0]
            prm_atom = prm_match.iloc[0]

            res_name = str(atom['Residue'])
            res_seq = int(atom['Residue Sequence'])
            res_key = f"{res_name}_{res_seq}"

            if res_key not in res_data:
                res_data[res_key] = {
                    "res_name": res_name,
                    "res_seq": res_seq,
                    "aw_vol": 0.0,
                    "pow_vol": 0.0,
                    "prm_vol": 0.0,
                }

            res_data[res_key]["aw_vol"] += float(atom["Volume"])
            res_data[res_key]["pow_vol"] += float(pow_atom["Volume"])
            res_data[res_key]["prm_vol"] += float(prm_atom["Volume"])

        n_res = len(res_data)
        n_added_before = len(x_vals)

        for i_res, (res_key, d) in enumerate(res_data.items(), start=1):

            if (i_res % PRINT_EVERY_RES) == 0 or i_res == n_res:
                print(f"    residues: {i_res} / {n_res}")

            aw_vol = float(d["aw_vol"])
            pow_vol = float(d["pow_vol"])
            prm_vol = float(d["prm_vol"])

            aw_sa = float(aw_res_sa.get(res_key, 0.0))
            pow_sa = float(pow_res_sa.get(res_key, 0.0))
            prm_sa = float(prm_res_sa.get(res_key, 0.0))

            # sphericity per scheme (guard against zeros)
            aw_sph = float(calc_sphericity(aw_vol, aw_sa)) if (aw_vol > 0.0 and aw_sa > 0.0) else np.nan
            pow_sph = float(calc_sphericity(pow_vol, pow_sa)) if (pow_vol > 0.0 and pow_sa > 0.0) else np.nan
            prm_sph = float(calc_sphericity(prm_vol, prm_sa)) if (prm_vol > 0.0 and prm_sa > 0.0) else np.nan


            # deviation vs AW (abs %)
            pow_vol_dev = (abs((pow_vol - aw_vol) / aw_vol) * 100.0) if aw_vol > 0.0 else np.nan
            prm_vol_dev = (abs((prm_vol - aw_vol) / aw_vol) * 100.0) if aw_vol > 0.0 else np.nan

            pow_sa_dev = (abs((pow_sa - aw_sa) / aw_sa) * 100.0) if aw_sa > 0.0 else np.nan
            prm_sa_dev = (abs((prm_sa - aw_sa) / aw_sa) * 100.0) if aw_sa > 0.0 else np.nan

            if scheme == "pow":
                dev = pow_vol_dev if deviation_kind == "vol" else pow_sa_dev
                sph_x = aw_sph if use_aw_sphericity else pow_sph
            else:
                dev = prm_vol_dev if deviation_kind == "vol" else prm_sa_dev
                sph_x = aw_sph if use_aw_sphericity else prm_sph

            if not np.isfinite(sph_x) or not np.isfinite(dev):
                continue

            # ---------------- sphericity sanity checks ----------------
            # Sphericity should be in (0, 1]. Anything above 1 is invalid -> store + skip.
            # Also skip NaNs/infs.
            if (not np.isfinite(sph_x)) or (sph_x <= 0.0) or (sph_x > 1.0):
                outliers.append({
                    "reason": "invalid_sphericity",
                    "system": my_key,
                    "res_key": res_key,
                    "scheme": scheme,
                    "deviation_kind": deviation_kind,
                    "use_aw_sphericity": bool(use_aw_sphericity),
                    "sph_x": float(sph_x) if np.isfinite(sph_x) else np.nan,
                    "aw_sph": float(aw_sph) if np.isfinite(aw_sph) else np.nan,
                    "pow_sph": float(pow_sph) if np.isfinite(pow_sph) else np.nan,
                    "prm_sph": float(prm_sph) if np.isfinite(prm_sph) else np.nan,
                    "aw_vol": aw_vol, "aw_sa": aw_sa,
                    "pow_vol": pow_vol, "pow_sa": pow_sa,
                    "prm_vol": prm_vol, "prm_sa": prm_sa,
                })
                continue

            if (max_percent_diff is not None) and (dev > float(max_percent_diff)):
                continue

            x_vals.append(float(sph_x))
            y_vals.append(float(dev))

            records.append({
                "system": my_key,
                "res_key": res_key,
                "res_name": d["res_name"],
                "res_seq": d["res_seq"],
                "x_sphericity": float(sph_x),
                "y_deviation": float(dev),
                "aw_sphericity": float(aw_sph) if np.isfinite(aw_sph) else np.nan,
                "pow_sphericity": float(pow_sph) if np.isfinite(pow_sph) else np.nan,
                "prm_sphericity": float(prm_sph) if np.isfinite(prm_sph) else np.nan,
                "pow_vol_dev": float(pow_vol_dev) if np.isfinite(pow_vol_dev) else np.nan,
                "pow_sa_dev": float(pow_sa_dev) if np.isfinite(pow_sa_dev) else np.nan,
                "prm_vol_dev": float(prm_vol_dev) if np.isfinite(prm_vol_dev) else np.nan,
                "prm_sa_dev": float(prm_sa_dev) if np.isfinite(prm_sa_dev) else np.nan,
                "aw_vol": aw_vol,
                "pow_vol": pow_vol,
                "prm_vol": prm_vol,
                "aw_sa": aw_sa,
                "pow_sa": pow_sa,
                "prm_sa": prm_sa,
                "scheme": scheme,
                "deviation_kind": deviation_kind,
                "use_aw_sphericity": bool(use_aw_sphericity),
            })
        n_added_after = len(x_vals)
        print(f"    done (points added: {n_added_after - n_added_before})\n")

    return np.asarray(x_vals, dtype=float), np.asarray(y_vals, dtype=float), records, outliers


def plot_residue_sphericity_vs_deviation(
    folder=None,
    exclude_keys=None,
    schemes=("pow", "prm"),        # plot 1 or 2 series
    deviation_kind="vol",
    use_aw_sphericity=True,
    max_percent_diff=None,
    title="Residue sphericity vs deviation",
    show=True,
    save=None,
):
    xs = []
    ys = []
    labels = []

    for sch in schemes:
        x, y, _, outliers = build_residue_sphericity_vs_deviation_points(
            folder=folder,
            exclude_keys=exclude_keys,
            scheme=sch,
            deviation_kind=deviation_kind,
            use_aw_sphericity=use_aw_sphericity,
            max_percent_diff=max_percent_diff,
        )

        # scatter() expects list-of-lists
        xs.append(list(x))
        ys.append(list(y))
        labels.append(sch.upper())

    x_label = "Residue sphericity (AW)" if use_aw_sphericity else "Residue sphericity (scheme)"
    y_label = f"Abs % {deviation_kind.upper()} deviation vs AW"

    scatter(
        xs=xs,
        ys=ys,
        labels=labels,
        title=title,
        x_axis_title=x_label,
        y_axis_title=y_label,
        Show=show,
        save=save,
        markers='o',
        marker_size=25,
        alpha=0.7,
        legend_title="Scheme",
        legend_loc="upper right",
        legend_bbox_to_anchor=(1.25, 0.97),
    )
    print(outliers)


# ---------------------------- example call ----------------------------

if __name__ == "__main__":
    plot_residue_sphericity_vs_deviation(
        folder=None,                 # pops dialog
        exclude_keys=["A", "B", "C"],
        schemes=("pow", "prm"),
        deviation_kind="vol",        # "vol" or "sa"
        use_aw_sphericity=True,
        max_percent_diff=200.0,
        title="Residue Sphericity vs Deviation (All Residues)",
        show=True,
        save=None,
    )
