import os
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.system import System


def num_missing_balls(folder=None):
    # If the folder option isnt chosen prompt the user to choose a folder
    if folder is None:
        # Get the folder with all the logs
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes()
        folder = filedialog.askdirectory(title="Choose a data folder")

    # Create the densities data
    missing_balls = {}
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
            my_aw = read_logs2(aw)
        except TypeError:
            print(pdb, aw, pow)
            continue
        ther_balls_counter = 0
        # Loop through the balls
        for i, ball in my_aw['atoms'].iterrows():

            # Make sure the ball is complete
            if ball['Complete Cell?']:
                ther_balls_counter += 1
        # Create the missing balls counter
        if (sf_cv, sf_den) in missing_balls:
            missing_balls[(sf_cv, sf_den)][0] += ther_balls_counter
            missing_balls[(sf_cv, sf_den)][1] += 1000
        else:
            missing_balls[(sf_cv, sf_den)] = [ther_balls_counter, 1000]

    return missing_balls


def plot_heatmap(missing_balls):
    # Extract unique CV and density values from missing_balls
    cv_values = sorted({key[0] for key in missing_balls.keys()})
    density_values = sorted({key[1] for key in missing_balls.keys()})

    # Create a 2D array for percentages
    data = np.zeros((len(cv_values), len(density_values)))

    # Populate the data array with percentages from missing_balls
    for (cv, density), (ther_balls_counter, total_balls) in missing_balls.items():
        i = cv_values.index(cv)
        j = density_values.index(density)
        percentage = (ther_balls_counter / total_balls) * 100
        data[i, j] = percentage

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    heatmap = ax.imshow(data, cmap="coolwarm", aspect="auto")

    # Add color bar
    cbar = plt.colorbar(heatmap)
    cbar.set_label("Percentage (%)", fontsize=14)

    # Set axis labels and title
    ax.set_xlabel("Density Values", fontsize=14)
    ax.set_ylabel("CV Values", fontsize=14)
    ax.set_title("Heatmap of Percentages Across Density and CV Values", fontsize=16)

    # Set axis tick labels
    ax.set_xticks(np.arange(len(density_values)))
    ax.set_yticks(np.arange(len(cv_values)))
    ax.set_xticklabels([f"{v:.2f}" for v in density_values])
    ax.set_yticklabels([f"{v:.2f}" for v in cv_values])

    # Rotate x-axis tick labels for better readability
    plt.xticks(rotation=45)

    # Annotate each cell with the percentage value
    for i in range(len(cv_values)):
        for j in range(len(density_values)):
            ax.text(j, i, f"{data[i, j]:.1f}", ha="center", va="center", color="black")

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    m_blizzys = num_missing_balls()
    plot_heatmap(m_blizzys)
