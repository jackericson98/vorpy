import csv
import os
from os import path
import tkinter as tk
from tkinter import filedialog
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit

import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.inputs import read_pdb_line
from vorpy.src.calculations import calc_sphericity, calc_isoperimetric_quotient
from vorpy.src.calculations import sort_lists
from vorpy.src.analyze.tools.compare.read_logs import read_logs
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def select_folder():
    """Prompt user to select a higher-level folder using a GUI dialog."""
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder_path = filedialog.askdirectory(title='Choose the higher-level folder')
    return folder_path


def get_files(folder):
    """Traverse the folder and extract required PDB and log files.
    Logs are identified based on subfolder names containing 'aw' or 'pow'.
    """
    pdb_file, aw_logs, pow_logs = None, None, None

    for root_dir, sub_dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.pdb') and 'atoms' not in file and 'diff' not in file:
                pdb_file = path.join(root_dir, file)

        # Check sub-subfolders for logs
        for sub_dir in sub_dirs:
            if 'aw' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        aw_logs = path.join(root_dir, sub_dir, file)
            elif 'pow' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        pow_logs = path.join(root_dir, sub_dir, file)
    return pdb_file, aw_logs, pow_logs


def read_logs_data(pow_logs, aw_logs):
    """Read power logs and AW logs using the appropriate functions."""
    try:
        pow = read_logs2(pow_logs)
        aw = read_logs2(aw_logs)
        old_logs = False
    except Exception:
        pow = read_logs(pow_logs)
        aw = read_logs(aw_logs)
        old_logs = True
    return pow, aw, old_logs


def parse_pdb_file(pdb_file, pow, aw, old_logs):
    """Parse the PDB file and calculate differences and sphericities."""
    rads, diffs, abs_diffs = [], [], []
    aw_spcys, pow_spcys, spcy_diffs = [], [], []
    aw_isos, pow_isos, iso_diffs = [], [], []
    missing_rads = []

    with open(pdb_file, 'r') as pdb_reader:
        for i, line in enumerate(pdb_reader.readlines()):
            if i == 0:  # Skip the first line
                continue
            if i > 1001:
                continue

            # Extract data from PDB line
            npl = read_pdb_line(line)
            try:
                pow_atom = \
                    pow['atoms'].loc[pow['atoms']['Index' if not old_logs else 'num'] == i - 1].to_dict('records')[0]
                vor_atom = aw['atoms'].loc[aw['atoms']['Index' if not old_logs else 'num'] == i - 1].to_dict('records')[
                    0]
            except IndexError:
                missing_rads.append(npl['temperature_factor'])
                continue

            # Process radii and volume differences
            rads.append(npl['temperature_factor'])
            pow_vol, pow_sa = pow_atom['Volume' if not old_logs else 'volume'], pow_atom['Surface Area' if not old_logs else 'sa']
            aw_vol, aw_sa = vor_atom['Volume' if not old_logs else 'volume'], vor_atom['Surface Area' if not old_logs else 'sa']

            abs_vol_diff = (100 * abs(pow_vol - aw_vol) / aw_vol) if aw_vol else np.nan
            vol_diff = (100 * (pow_vol - aw_vol) / aw_vol) if aw_vol else np.nan
            if vol_diff >= 1000:  # Handle extreme cases
                vol_diff = np.nan

            diffs.append(vol_diff)
            abs_diffs.append(abs_vol_diff)

            # Calculate sphericity
            aw_sphericity = calc_sphericity(aw_vol, aw_sa)
            pow_sphericity = calc_sphericity(pow_vol, pow_sa)
            spcy_diff = (pow_sphericity - aw_sphericity) / aw_sphericity
            if 0 > aw_sphericity > 1 or 0 > pow_sphericity > 1 or abs(spcy_diff) > 1:
                aw_sphericity, pow_sphericity, spcy_diff = np.nan, np.nan, np.nan
            aw_spcys.append(aw_sphericity)
            pow_spcys.append(pow_sphericity)
            spcy_diffs.append(spcy_diff)

            # Calculate isoperimetric quotient
            aw_iso = calc_isoperimetric_quotient(aw_vol, aw_sa)
            pow_iso = calc_isoperimetric_quotient(pow_vol, pow_sa)
            iso_diff = (pow_iso - aw_iso) / aw_iso
            if aw_iso > 1 or pow_iso > 1 or abs(iso_diff) > 1:
                aw_iso, pow_iso, iso_diff = np.nan, np.nan, np.nan
            aw_isos.append(aw_iso)
            pow_isos.append(pow_iso)
            iso_diffs.append(iso_diff)

    return rads, diffs, abs_diffs, spcy_diffs, aw_spcys, pow_spcys, aw_isos, pow_isos, iso_diffs, missing_rads


def plot_spcy_results(rads, diffs, spcy_diffs, ax):
    """Plot the cumulative results using Matplotlib."""
    cmap = plt.cm.autumn

    norm = Normalize(vmin=min(rads), vmax=max(rads))
    sm = ScalarMappable(norm=norm, cmap=cmap)
    rads, spcy_diffs, diffs = sort_lists(rads, spcy_diffs, diffs, reverse=False)
    ax.plot([0, 0], [-200, 1000], c='k', alpha=0.4)
    ax.plot([-0.65, 0.65], [0, 0], c='k', alpha=0.4)
    ax.scatter(spcy_diffs, diffs, c=rads, cmap=cmap, alpha=0.15, s=1, marker='x')

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Ball Radius', fontsize=25)
    cbar.ax.tick_params(labelsize=20, width=2, length=12)
    # cbar.ax.set_yticklabels([f'{min(rads)}'])

    ax.set_xlabel('\u03A8 Diff', fontsize=30)
    ax.set_ylabel('Vol Diff', fontsize=30)
    ax.set_yticks([-100, 175, 450, 725, 1000])
    ax.set_xticks([-0.5, 0, 0.5])
    ax.set_xlim([-0.65, 0.65])
    ax.set_ylim([-200, 1000])
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=12)
    ax.set_title('\u03A8 vs Vol Diff', fontsize=30)


def plot_radii_vs_sphericity(rads, spcy_diffs, diffs, ax):
    """Plot radii vs sphericity with a color bar mapping to diffs, sorting points by radii."""
    # Sort data by radii (largest to smallest)
    sorted_data = sorted(zip(rads, spcy_diffs, diffs), key=lambda x: -x[0])
    rads_sorted, aw_spcys_sorted, diffs_sorted = zip(*sorted_data)

    cmap = plt.cm.cool_r
    norm = Normalize(vmin=-100, vmax=max(diffs_sorted))
    sm = ScalarMappable(norm=norm, cmap=cmap)

    ax.plot([-1, 6], [0, 0], c='k', alpha=0.25)
    ax.scatter(rads_sorted, aw_spcys_sorted, c=diffs_sorted, cmap=cmap, alpha=0.15, s=1, marker='x')

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Vol Diff', fontsize=25)
    cbar.ax.tick_params(labelsize=20, width=2, length=12)
    cbar.set_ticks(ticks=[-100, 175, 450, 725, 1000], labels=[-100, 175, 450, 725, 1000])

    ax.set_xlabel('Ball Radius', fontsize=30)
    ax.set_xlim([-1, 6])
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_ylabel('\u03A8 Diff', fontsize=30)
    ax.set_yticks([-0.5, 0, 0.5])
    ax.set_ylim(-0.65, 0.65)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=12)
    ax.set_title('Radii vs \u03A8 Diff', fontsize=30)


def display_metrics(diffs, spcy_diffs, rads, aw_spcys, pow_spcys, ax=None):
    """Display text-based metrics on a subplot."""
    num_points = len(spcy_diffs)
    spcy_less_than_zero = len([s for s in spcy_diffs if s < 0]) / num_points * 100
    mean_sphericity = np.nanmean(spcy_diffs)
    all_radii_mean = np.nanmean(rads)
    mean_neg_sphericity = np.nanmean([s for s in spcy_diffs if s < 0])
    mean_pos_sphericity = np.nanmean([s for s in spcy_diffs if s > 0])
    mean_radius_neg_sphericity = np.nanmean([r for r, s in zip(rads, spcy_diffs) if s < 0])

    mean_radius_pos_sphericity = np.nanmean([r for r, s in zip(rads, spcy_diffs) if s > 0])
    mean_aw_sphericity = np.nanmean(aw_spcys)
    mean_pow_sphericity = np.nanmean(pow_spcys)
    tot = len(spcy_diffs)
    q1 = len([0 for d, s in zip(diffs, spcy_diffs) if d > 0 > s])
    q2 = len([0 for d, s in zip(diffs, spcy_diffs) if d > 0 and s > 0])
    q3 = len([0 for d, s in zip(diffs, spcy_diffs) if d < 0 < s])
    q4 = len([0 for d, s in zip(diffs, spcy_diffs) if d < 0 and s < 0])

    metrics_text = (
        f"% \u03A8 < 0: {spcy_less_than_zero:.2f}%\n"
        f"Quadrant Count: Top Left: {q1} {100*q1/tot:.3f} %, Top Right: {q2} {100*q2/tot:.3f} %, Bottom Right {q3} {100*q3/tot:.3f} %, Bottom Left: {q4} {100*q4/tot:.3f} %\n"
        f"Mean Vol Diff: {np.nanmean(diffs):.3f}\n"
        f"Mean Abs Vol Diff: {np.nanmean([abs(_) for _ in diffs]):.3f}\n"
        f"Mean Abs Vol Diff (\u03A8 < 0): {np.nanmean([abs(_) for _, s in zip(diffs, spcy_diffs) if s < 0]):.3f}\n"
        f"Mean Abs Vol Diff (\u03A8 > 0): {np.nanmean([abs(_) for _, s in zip(diffs, spcy_diffs) if s > 0]):.3f}\n"
        f"Mean Abs Vol Diff (\u03A8 q1): {np.nanmean([s for s, d in zip(spcy_diffs, diffs) if d > 0 > s]):.3f}\n"
        f"Mean Abs Vol Diff (\u03A8 q2): {np.nanmean([s for s, d in zip(spcy_diffs, diffs) if d > 0 and s > 0]):.3f}\n"
        f"Mean Abs Vol Diff (\u03A8 q3): {np.nanmean([s for s, d in zip(spcy_diffs, diffs) if d < 0 < s]):.3f}\n"
        f"Mean Abs Vol Diff (\u03A8 q4): {np.nanmean([s for s, d in zip(spcy_diffs, diffs) if d < 0 and s < 0]):.3f}\n"
        f"Mean Vol Diff (r < 0): {np.nanmean([_ for _, s in zip(diffs, spcy_diffs) if s < 0]):.3f}\n"
        f"Mean Vol Diff (r > 0): {np.nanmean([_ for _, s in zip(diffs, spcy_diffs) if s > 0]):.3f}\n"
        
        
        f"Mean \u03A8: {mean_sphericity:.3f}\n"
        f"Mean \u03A8 (\u03A8 < 0): {mean_neg_sphericity:.3f}\n"
        f"Mean \u03A8 (\u03A8 > 0): {mean_pos_sphericity:.3f}\n"
        
        f"Mean Radius: {all_radii_mean:.3f}\n"
        f"Mean Radius (\u03A8 < 0): {mean_radius_neg_sphericity:.3f}\n"
        f"Mean Radius (\u03A8 > 0): {mean_radius_pos_sphericity:.3f}\n"
        f"Mean AW \u03A8: {mean_aw_sphericity:.3f}\n"
        f"Mean Power \u03A8: {mean_pow_sphericity:.3f}"
    )
    if ax is None:
        print(metrics_text)
    else:
        ax.axis('off')
        ax.text(0, 0.5, metrics_text, fontsize=20, ha='left', va='center', linespacing=1.5)


def process_all_subfolders(base_folder):
    """Process all subfolders within the base folder and combine results cumulatively."""
    cum_rads, cum_diffs, cum_abs_diffs, cum_iso_diffs, cum_spcy_diffs, cum_aw_spcys, cum_pow_spcys = [], [], [], [], [], [], []

    for folder_name in os.listdir(base_folder):
        folder_path = path.join(base_folder, folder_name)
        if path.isdir(folder_path):
            pdb_file, aw_logs, pow_logs = get_files(folder_path)
            if all([pdb_file, aw_logs, pow_logs]):
                pow, aw, old_logs = read_logs_data(pow_logs, aw_logs)
                results = parse_pdb_file(pdb_file, pow, aw, old_logs)
                rads, diffs, abs_diffs, spcy_diffs, aw_spcys, pow_spcys, *_ = results
                cum_rads.extend(rads)
                cum_diffs.extend(diffs)
                cum_abs_diffs.extend(abs_diffs)
                cum_spcy_diffs.extend(spcy_diffs)
                cum_aw_spcys.extend(aw_spcys)
                cum_pow_spcys.extend(pow_spcys)
                cum_iso_diffs.extend(results[7])
    return cum_rads, cum_diffs, cum_abs_diffs, cum_iso_diffs, cum_spcy_diffs, cum_aw_spcys, cum_pow_spcys


def main():
    """Main function to run the script."""
    base_folder = select_folder()
    print(base_folder)
    cum_rads, cum_diffs, cum_abs_diffs, cum_iso_diffs, cum_spcy_diffs, cum_aw_spcys, cum_pow_spcys = process_all_subfolders(
        base_folder)
    if not cum_rads:
        print("Error: No valid data found in the selected folder structure.")
        return
    ax3 = None
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    plot_spcy_results(cum_rads, cum_diffs, cum_spcy_diffs, ax1)
    plot_radii_vs_sphericity(cum_rads, cum_spcy_diffs, cum_diffs, ax2)
    display_metrics(cum_diffs, cum_spcy_diffs, cum_rads, cum_aw_spcys, cum_pow_spcys, ax3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()