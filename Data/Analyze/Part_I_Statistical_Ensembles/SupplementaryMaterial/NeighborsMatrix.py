import os
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
import matplotlib.pyplot as plt
import numpy as np


def get_neighbors_data(folder=None):
    # Get the folder
    if folder is None:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes()
        folder = filedialog.askdirectory()
    # Create the data dictionary
    data = {}
    # Go through the subfolders
    num_sub_folders = len(os.listdir(folder))
    for i, subfolder in enumerate(os.listdir(folder)):
        print(f'\rGathering Data from {subfolder} -> {100 * i / num_sub_folders:.3f} %', end="")

        my_pdb, my_aw, my_pow = get_files(folder + '/' + subfolder)
        try:
            aw_logs = read_logs2(my_aw, all_=False, balls=True)
            pow_logs = read_logs2(my_pow, all_=False, balls=True)
            rads, pow_naybs, aw_naybs = [], [], []
        except TypeError:
            continue
        for i, aw_ball in aw_logs['atoms'].iterrows():
            try:
                pow_ball = pow_logs['atoms'].iloc[i]
            except IndexError:
                continue
            rads.append(aw_ball['Radius'])
            aw_naybs.append(aw_ball['Number of Neighbors'])
            pow_naybs.append(pow_ball['Number of Neighbors'])
        data[subfolder] = {'rads': rads, 'aw': aw_naybs, 'pow': pow_naybs}

    return data


def plot_neighbor_data(data):
    perdic, cv_vals, density_vals = {}, [], []
    for file_name in data:
        curdic = data[file_name]
        split_file_name = file_name.split('_')
        try:
            cv, density, file_number = float(split_file_name[1]), float(split_file_name[3]), int(split_file_name[-1])
        except ValueError:
            continue

        my_my_data = [(curdic['rads'][i], curdic['pow'][i] - curdic['aw'][i]) for i in range(len(curdic['aw']))]
        if (cv, density) in perdic:
            perdic[(cv, density)] += my_my_data
        else:
            perdic[(cv, density)] = my_my_data

        if cv not in cv_vals:
            cv_vals.append(cv)
        if density not in density_vals:
            density_vals.append(density)

    density_vals.sort(reverse=True)
    cv_vals.sort()

    fig, axes = plt.subplots(len(density_vals), len(cv_vals), figsize=(20, 18), sharex='all', sharey='all')

    for i, density in enumerate(density_vals):
        for j, cv in enumerate(cv_vals):
            ax = axes[i, j]
            # Set labels on the left column and bottom row
            if j == 0:
                ax.set_ylabel("Nbor Diff")
            if i == len(density_vals) - 1:
                ax.set_xlabel("Radius", va='top')
            ax.set_ylim([-30, 30])
            try:
                my_data = perdic[(cv, density)]
                ax.plot([0, 5], [0, 0], c='k', linewidth=1)
                ax.scatter([_[0] for _ in my_data], [_[1] for _ in my_data], s=1, alpha=0.25, c='g')
            except KeyError:
                continue


    # Add CV values to the bottom of the figure
    for i, col in enumerate(cv_vals):
        fig.text(0.155 + i * (0.79 / len(cv_vals)), 0.07, col, ha='center', fontsize=20)

    for i, row in enumerate(density_vals[::-1]):
        fig.text(0.05, 0.18 + i * (0.76 / len(density_vals)), row, va='center', rotation='horizontal',
                 fontsize=20)
    # Additional figure-wide settings
    fig.text(0.5, 0.02, 'CV Values', ha='center', fontsize=25)
    fig.text(0.02, 0.5, 'Density Values', va='center', rotation='vertical', fontsize=25)
    fig.suptitle("Power - AW Neighbor Count Differences", fontsize=30)
    plt.subplots_adjust(left=0.135, bottom=0.15, top=0.9)

    plt.show()

if __name__ == '__main__':
    my_data = get_neighbors_data()
    plot_neighbor_data(my_data)
