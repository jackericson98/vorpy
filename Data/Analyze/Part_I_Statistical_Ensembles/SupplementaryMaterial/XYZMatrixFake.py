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
import random


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
    num_folders = len(os.listdir(folder))
    for k, subfolder in enumerate(os.listdir(folder)):

        print(f"\rFolder {k}/{num_folders} - {100 * (k/num_folders)}%", end="")
        split_subfolder = subfolder.split("_")
        try:
            sf_cv, sf_den = float(split_subfolder[1]), float(split_subfolder[3])
        except ValueError:
            continue
        except IndexError:
            continue

        # # Get PDB, AW, and POW files
        # pdb, aw, pow = get_files(os.path.join(folder, subfolder))

        # try:
        #     my_sys = System(pdb, simple=True)
        # except TypeError:
        #     print(f"Error loading system: {pdb}, {aw}, {pow}")
        #     continue
        # except IndexError:
        #     print(pdb)
        #     continue
        # except ValueError:
        #     print(pdb)
        #     continue

        # foam_box = float(my_sys.data[0][5][:-1])
        # norm_locs = []

        random_triples = [(random.random(), random.random(), random.random()) for _ in range(1000)]

        if (sf_cv, sf_den) in loc_data:
            loc_data[(sf_cv, sf_den)] += random_triples
        else:
            loc_data[(sf_cv, sf_den)] = random_triples
    return loc_data


# def distribution_of_overlaps(loc_data):
#     """
#     Plot the distribution of normalized locations for systems based on CV and density values.
#
#     Args:
#         loc_data (dict): A dictionary with keys as (cv, density) and values as lists of normalized locations.
#     """
#     cv_vals = sorted(set(key[0] for key in loc_data.keys()))
#     density_vals = sorted(set(key[1] for key in loc_data.keys()), reverse=True)
#
#     # Adjust the figure size for smaller subplots
#     fig, axes = plt.subplots(len(density_vals), len(cv_vals), figsize=(12, 9), sharex=True, sharey=True)
#
#     for i, density in enumerate(density_vals):
#         for j, cv in enumerate(cv_vals):
#             ax = axes[i, j]
#             points = np.array(loc_data.get((cv, density), []))
#
#             if len(points) > 0:
#                 scatter = ax.scatter(points[:, 0], points[:, 1], c=points[:, 2], cmap='gray', s=0.8, marker='x')
#
#             ax.set_xlim(0, 1)
#             ax.set_ylim(0, 1)
#             ax.set_xticks([0, 1])
#             ax.set_yticks([0, 1])
#
#             # Reduce font size for ticks to prevent overlap
#             ax.tick_params(axis='both', which='major', labelsize=8)
#
#     for ax, col in zip(axes[-1], cv_vals):
#         ax.set_xlabel("")  # Remove direct subplot labels, handled below
#
#     for ax, row in zip(axes[:, 0], density_vals):
#         ax.set_ylabel("")  # Remove direct subplot labels, handled below
#
#     # Add CV values to the bottom of the figure
#     for i, col in enumerate(cv_vals):
#         fig.text(0.13 + i * (0.87 / len(cv_vals)), 0.1, f"{col:.2f}", ha='center', fontsize=15)
#
#     for i, row in enumerate(density_vals[::-1]):
#         fig.text(0.05, 0.18 + i * (0.82 / len(density_vals)), f"{row:.2f}", va='center', rotation='horizontal', fontsize=15)
#
#     # Add overall axis labels
#     fig.text(0.525, 0.05, 'CV', ha='center', fontsize=20)
#     fig.text(0.02, 0.5, 'Density', va='center', rotation='vertical', fontsize=20)
#
#     # Adjust subplot spacing to prevent overlap
#     plt.tight_layout(pad=2.0, w_pad=1.0, h_pad=1.0)
#     fig.subplots_adjust(bottom=0.15, left=0.1, top=0.95, right=0.95)
#
#     plt.show()

def distribution_of_overlaps(loc_data):
    """
    Plot the distribution of normalized locations for systems based on CV and density values,
    including color bars for z-values on the right of the last column of subplots.

    Args:
        loc_data (dict): A dictionary with keys as (cv, density) and values as lists of normalized locations.
    """
    cv_vals = sorted(set(key[0] for key in loc_data.keys()))
    density_vals = sorted(set(key[1] for key in loc_data.keys()), reverse=True)

    # Adjust the figure size for smaller subplots
    fig, axes = plt.subplots(len(density_vals), len(cv_vals), figsize=(12, 9), sharex=True, sharey=True)

    for i, density in enumerate(density_vals):
        for j, cv in enumerate(cv_vals):
            ax = axes[i, j]
            points = np.array(loc_data.get((cv, density), []))

            if len(points) > 0:
                scatter = ax.scatter(points[:, 0], points[:, 1], c=points[:, 2], cmap='Greys', s=0.8, marker='x')

                # Add a color bar to the subplots in the last column
                if j == len(cv_vals) - 1:
                    cbar = fig.colorbar(scatter, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
                    cbar.set_label("z", fontsize=10, rotation='horizontal', va='center')
                    cbar.set_ticks([0, 1])
                # Add the y label
                if j == 0:
                    ax.set_ylabel('y', rotation='horizontal', va='center')
                # Add the x label
                if i == len(density_vals) - 1:
                    ax.set_xlabel('x', va='center')

            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticks([0, 1])
            ax.set_yticks([0, 1])

            # Reduce font size for ticks to prevent overlap
            ax.tick_params(axis='both', which='major', labelsize=8)

    # Add CV values to the bottom of the figure
    for i, col in enumerate(cv_vals):
        fig.text(0.16 + i * (0.83 / len(cv_vals)), 0.07, f"{col:.2f}", ha='center', fontsize=15)

    # Add density values to the left of the figure
    for i, row in enumerate(density_vals[::-1]):
        fig.text(0.05, 0.15 + i * (0.85 / len(density_vals)), f"{row:.2f}", va='center', rotation='horizontal', fontsize=15)

    # Add overall axis labels
    fig.text(0.525, 0.02, 'CV', ha='center', fontsize=20)
    fig.text(0.02, 0.5, 'Density', va='center', rotation='vertical', fontsize=20)

    # Adjust subplot spacing to prevent overlap
    plt.tight_layout(pad=2.0, w_pad=1.0, h_pad=1.0)
    fig.subplots_adjust(bottom=0.13, left=0.13, top=0.95, right=0.95)

    plt.show()



if __name__ == '__main__':
    os.chdir('../../../..')
    loc_data = get_syses()
    distribution_of_overlaps(loc_data)
