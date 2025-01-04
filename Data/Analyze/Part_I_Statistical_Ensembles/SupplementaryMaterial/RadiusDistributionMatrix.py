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


def gamma(r, cv, mu=1):
    # Gamma parameters
    alpha = 1 / cv ** 2
    beta = alpha / mu  # To keep mean = 1
    gamma_dist = stats.gamma(a=alpha, scale=1 / beta)
    # Compute PDFs
    return gamma_dist.pdf(r)

def get_syses(folder=None):
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
    for subfolder in os.listdir(folder):
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
            my_syses = my_dict[(cv, density)]
            radii = []
            for rads in my_syses:
                radii += rads
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

            # Update plot formatting
            # ax.set_xlabel('Radius', fontsize=25)
            # ax.set_ylabel('Count', fontsize=25, color='blue')
            # ax.set_xticks([round(_, 1) for _ in np.linspace(min_rad, max_rad, 4)])
            # ax.tick_params(axis='both', labelsize=20)
            # ax.tick_params(axis='y', colors='blue')
            ax.set_ylim(bottom=0)

            # Create a secondary y-axis for the histogram
            ax2 = ax.twinx()

            # Plot the scaled gamma PDF
            ax.plot(x_values, scaled_pdf, label='Scaled PDF', color='red', linewidth=2)
            # ax2.set_ylabel('Probability', fontsize=25, color='red')
            # ax2.tick_params(axis='both', labelsize=20, colors='red')
            if cv != 1.0:
                ax2.set_yticks([0, 100], [0, 20,000])
            ax.set_yticks([])
            ax2.set_ylim(bottom=0, top=100)
            # ax2.tick_params(axis='y', colors='red')


    for ax, col in zip(axes[-1], cv_vals):
        ax.set_xlabel("")  # Remove direct subplot labels, handled below

    for ax, row in zip(axes[:, 0], density_vals):
        ax.set_ylabel("")  # Remove direct subplot labels, handled below

    # Add CV values to the bottom of the figure
    for i, col in enumerate(cv_vals):
        fig.text(0.13 + i * (0.815 / len(cv_vals)), 0.05, col, ha='center', fontsize=15)

    for i, row in enumerate(density_vals[::-1]):
        fig.text(0.05, 0.125 + i * (0.825 / len(density_vals)), row, va='center', rotation='horizontal', fontsize=15)

    fig.text(0.5, 0.02, 'CV Values', ha='center', fontsize=20)  # X-axis label for the figure
    fig.text(0.02, 0.5, 'Density Values', va='center', rotation='vertical', fontsize=20)  # Y-axis label for the figure

    # Add a main title for the entire figure
    fig.suptitle("Distribution of radii by CV and Density", fontsize=20)

    # Adjust layout to prevent overlapping labels
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.9)
    # plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    # Run the code
    os.chdir('../../../..')
    print(os.getcwd())
    dicty = get_syses()
    print(dicty)
    distribution_of_overlaps(dicty)
