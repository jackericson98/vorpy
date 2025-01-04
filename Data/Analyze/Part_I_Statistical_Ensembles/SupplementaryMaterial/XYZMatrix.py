import os
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
from System.system import System
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.sys_funcs.calcs.calcs import calc_dist
from scipy import stats


def get_syses(folder=None):
    """
    Gather normalized location data from simulation systems in a folder structure.

    Args:
        folder (str, optional): The root folder containing subfolders with system data. Defaults to prompting user.

    Returns:
        dict: A dictionary with keys as (cv, density) and values as lists of normalized locations.
    """
    if folder is None:
        # Prompt user to choose a folder if not provided
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory(title="Choose a data folder")

    loc_data = {}

    for subfolder in os.listdir(folder):
        split_subfolder = subfolder.split("_")
        try:
            sf_cv, sf_den = float(split_subfolder[1]), float(split_subfolder[3])
        except ValueError:
            continue

        # Get PDB, AW, and POW files
        pdb, aw, pow = get_files(os.path.join(folder, subfolder))

        try:
            my_sys = System(pdb, simple=True)
        except TypeError:
            print(f"Error loading system: {pdb}, {aw}, {pow}")
            continue

        foam_box = my_sys.data[0][2]
        norm_locs = []

        for _, ball in my_sys.balls.iterrows()[:1000]:
            my_loc = ball['loc']
            norm_loc = [my_loc[i] / foam_box for i in range(3)]
            norm_locs.append(norm_loc)

        if (sf_cv, sf_den) in loc_data:
            loc_data[(sf_cv, sf_den)] += norm_locs
        else:
            loc_data[(sf_cv, sf_den)] = norm_locs

    return loc_data


def distribution_of_overlaps(loc_data):
    """
    Plot the distribution of normalized locations for systems based on CV and density values.

    Args:
        loc_data (dict): A dictionary with keys as (cv, density) and values as lists of normalized locations.
    """
    cv_vals = sorted(set(key[0] for key in loc_data.keys()))
    density_vals = sorted(set(key[1] for key in loc_data.keys()), reverse=True)

    fig, axes = plt.subplots(len(density_vals), len(cv_vals), figsize=(20, 18), sharex=True, sharey=True)

    for i, density in enumerate(density_vals):
        for j, cv in enumerate(cv_vals):
            ax = axes[i, j]
            points = np.array(loc_data.get((cv, density), []))

            if len(points) > 0:
                scatter = ax.scatter(points[:, 0], points[:, 1], c=points[:, 2], cmap='gray', s=0.8, marker='x')
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_ticks([0, 1])
                cbar.ax.tick_params(labelsize=10)

            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticks([0, 1])
            ax.set_yticks([0, 1])

            if i == len(density_vals) - 1:
                ax.set_xlabel(f"CV={cv:.2f}", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"Density={density:.2f}", fontsize=10)

    fig.text(0.5, 0.02, 'CV Values', ha='center', fontsize=20)
    fig.text(0.02, 0.5, 'Density Values', va='center', rotation='vertical', fontsize=20)
    fig.suptitle("Distribution of Normalized Locations by CV and Density", fontsize=24)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    os.chdir('../../../..')
    loc_data = get_syses()
    distribution_of_overlaps(loc_data)
