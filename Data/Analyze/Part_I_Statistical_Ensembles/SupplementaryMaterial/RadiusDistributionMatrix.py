import os
import time
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
from System.system import System
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.sys_funcs.calcs.calcs import calc_dist
from scipy import stats
from System.sys_funcs.calcs.calcs import get_time


def gamma(r, cv, mu=1):
    # Gamma parameters
    alpha = 1 / cv ** 2
    beta = alpha / mu  # To keep mean = 1
    gamma_dist = stats.gamma(a=alpha, scale=1 / beta)
    # Compute PDFs
    return gamma_dist.pdf(r)


def get_syses1(folder=None):
    # If the folder option isnt chosen prompt the user to choose a folder
    if folder is None:
        # Get the folder with all the logs
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes()
        folder = filedialog.askdirectory(title="Choose a data folder")

    # Create the densities data
    den_data = {}
    # Loop through the folders
    num_folders = len(os.listdir(folder))
    for k, subfolder in enumerate(os.listdir(folder)):
        print(f"Folder {k}/{num_folders} - {100 * (k / num_folders)}%")
        # Get the cv and density values
        split_subfolder = subfolder.split("_")
        try:
            sf_cv, sf_den = float(split_subfolder[1]), float(split_subfolder[3])
        except:
            continue
        # Get the pdb, aw, and pow
        pdb, aw, pow = get_files(folder + '/' + subfolder)
        # Get the logs and make a system
        try:
            my_sys = System(pdb, simple=True)
            if (sf_cv, sf_den) in den_data:
                den_data[(sf_cv, sf_den)].append([_ for _ in my_sys.balls['rad']])
            else:
                den_data[(sf_cv, sf_den)] = [[_ for _ in my_sys.balls['rad']]]
        except TypeError:
            print(pdb, aw, pow)
        except IndexError:
            print(pdb, aw, pow, subfolder)
    print(folder)
    return den_data


def get_syses(folder=None):
    """
    Efficiently process system data from a folder structure, organizing by cv and density.

    Args:
        folder (str, optional): The root folder containing subfolders with system data.
                               If not provided, prompts the user to select a folder.

    Returns:
        dict: A dictionary with keys as (cv, density) tuples and values as lists of radii data.
    """

    # Prompt user to choose a folder if not provided
    if folder is None:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory(title="Choose a data folder")

    den_data = {}

    # Get all subfolders once to avoid multiple calls to os.listdir
    subfolders = os.listdir(folder)
    num_folders = len(subfolders)
    start = time.perf_counter()
    for k, subfolder in enumerate(subfolders):
        print(f"\rProcessing Folder {k + 1}/{num_folders} - {100 * (k / num_folders):.2f}% - Time Elapsed = "
              f"{get_time(time.perf_counter() - start)}", end="")

        # Extract cv and density values
        split_subfolder = subfolder.split("_")
        try:
            sf_cv, sf_den = map(float, (split_subfolder[1], split_subfolder[3]))
        except (IndexError, ValueError):
            continue  # Skip invalid subfolders

        # Get file paths
        try:
            pdb, aw, pow = get_files(os.path.join(folder, subfolder))
        except Exception as e:
            print(f"Error getting files in {subfolder}: {e}")
            continue

        # Load system and extract radii
        try:
            my_sys = System(pdb, simple=True)
            radii = my_sys.balls['rad'].tolist()  # Efficiently convert column to list
            if (sf_cv, sf_den) in den_data:
                den_data[(sf_cv, sf_den)].append(radii)
            else:
                den_data[(sf_cv, sf_den)] = [radii]
        except (TypeError, IndexError) as e:
            print(f"Error processing system in {subfolder}: {e}")

    print(f"Processed folder: {folder}")
    return den_data


def distribution_of_overlaps(my_dict=None):

    cv_vals, density_vals = [], []
    for cv, den in my_dict:
        if cv not in cv_vals:
            cv_vals.append(cv)
        if den not in density_vals:
            density_vals.append(den)
    density_vals.sort(reverse=True)
    cv_vals.sort()

    fig, axes = plt.subplots(10, 11, figsize=(20, 18), sharex='all', sharey='all')

    for i, density in enumerate(density_vals):
        for j, cv in enumerate(cv_vals):
            # if cv not in {0.05, 0.1}:
            #     continue
            ax = axes[i, j]
            # Real Radii
            radii = []
            try:
                my_syses = my_dict[(cv, density)]

                for rads in my_syses:
                    radii += rads
            except KeyError:
                pass
            # Plot histogram
            data2 = ax.hist(radii, bins=20, alpha=0.5, color='blue', edgecolor='k', density=True)
            min_rad = min(radii)
            max_rad = max(radii)
            # Calculate bin width
            bin_edges = data2[1]
            bin_width = bin_edges[1] - bin_edges[0]
            total_area = sum(data2[0]) * bin_width

            # Generate x values for the PDF
            x_values = np.linspace(min_rad, max_rad, 100)

            # Scale gamma PDF to match the histogram area
            pdf_values = gamma(x_values, cv)
            scaled_pdf = pdf_values * total_area

            ax.set_ylim(bottom=0)

            # Create a secondary y-axis for the histogram
            ax2 = ax.twinx()

            # Plot the scaled gamma PDF
            ax.plot(x_values, scaled_pdf, label='Scaled PDF', color='red', linewidth=2)
            # ax2.set_ylabel('Probability', fontsize=25, color='red')
            # ax2.tick_params(axis='both', labelsize=20, colors='red')
            if cv == 1.0:
                ax2.set_yticks([0, 100], [0, 20000])
            else:
                ax2.set_yticks([])
            ax.set_yticks([])
            ax2.set_ylim(bottom=0, top=100)
            # ax2.tick_params(axis='y', colors='red')

    for ax, col in zip(axes[-1], cv_vals):
        ax.set_xlabel("")  # Remove direct subplot labels, handled below

    for ax, row in zip(axes[:, 0], density_vals):
        ax.set_ylabel("")  # Remove direct subplot labels, handled below

    # Add CV values to the bottom of the figure
    for i, col in enumerate(cv_vals):
        fig.text(0.13 + i * (0.815 / len(cv_vals)), 0.05, col, ha='center', fontsize=20)

    for i, row in enumerate(density_vals[::-1]):
        fig.text(0.05, 0.125 + i * (0.825 / len(density_vals)), row, va='center', rotation='horizontal', fontsize=20)

    fig.text(0.5, 0.02, 'CV Values', ha='center', fontsize=25)  # X-axis label for the figure
    fig.text(0.02, 0.5, 'Density Values', va='center', rotation='vertical', fontsize=25)  # Y-axis label for the figure

    # Add a main title for the entire figure
    fig.suptitle("Distribution of Ball Radii", fontsize=30)

    # Adjust layout to prevent overlapping labels
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.9)
    # plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    # Run the code
    os.chdir('../../../..')
    dicty = get_syses()
    distribution_of_overlaps(dicty)
