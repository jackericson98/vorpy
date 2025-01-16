import os
import tkinter as tk
from tkinter import filedialog
import numpy as np
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from Data.Analyze.tools.compare.read_logs import read_logs
from System.system import System
import matplotlib.pyplot as plt


def get_net_neighbors(pdb, aw, pow):
    # Create the simple system
    my_sys = System(pdb, simple=True)
    # Get the logs
    old_logs = False
    try:
        aw_data = read_logs2(aw)
        pow_data = read_logs2(pow)
    except KeyError:
        old_logs = True
        aw_data = read_logs(aw)
        pow_data = read_logs(pow)

    # Create the lists for the data
    net_nborss, rads, nbor_per_diffs = [], [], []
    # Go through each ball
    for i, ball in my_sys.balls.iterrows():
        try:
            # Check for old logs
            if not old_logs:
                # Get the pow atom and the vor atom
                pow_atom, aw_atom = pow_data['atoms'].loc[pow_data['atoms']['Index'] == i - 1].to_dict('records')[0], \
                aw_data['atoms'].loc[aw_data['atoms']['Index'] == i - 1].to_dict('records')[0]

                # Calculate the difference in volume
                pow_nbors, aw_nbors = pow_atom['Number of Neighbors'], aw_atom['Number of Neighbors']

            else:
                # Get the pow atom and the vor atom
                pow_atom, aw_atom = pow_data['atoms'].loc[pow_data['atoms']['num'] == i - 1].to_dict('records')[0], \
                aw_data['atoms'].loc[aw_data['atoms']['num'] == i - 1].to_dict('records')[0]

                # Calculate the difference in volume
                pow_nbors, aw_nbors = len(pow_atom['neighbors']), len(aw_atom['neighbors'])
            # Get the data
            rads.append(ball['rad'])
            net_nborss.append(pow_nbors - aw_nbors)
            nbor_per_diffs.append(100 * (pow_nbors - aw_nbors) / aw_nbors)
        except IndexError:
            continue
    # Return the data
    return rads, net_nborss, nbor_per_diffs


def plot_data(rads, net_nbors, cv, den):
    # Get the average radius for each integer
    mean_rad_dict = {}
    for i, net_nbor in enumerate(net_nbors):
        if net_nbor in mean_rad_dict:
            mean_rad_dict[net_nbor].append(rads[i])
        else:
            mean_rad_dict[net_nbor] = [rads[i]]

    # Plot the net_nbors
    plt.scatter(rads, net_nbors, s=2, alpha=0.5)
    # Plot the average radius of each net neighbor difference
    plt.scatter([np.mean(mean_rad_dict[_]) for _ in mean_rad_dict], mean_rad_dict.keys(), s=10, c='r', marker='x')
    plt.title('# Power - # AW Neighbors', fontsize=30)
    plt.xlabel('Ball Radius', fontdict=dict(size=25))
    plt.ylabel(f'CV = {cv}, Density = {den}', fontdict=dict(size=25))

    plt.xticks(font=dict(size=20))
    plt.yticks(font=dict(size=20))
    plt.tick_params(axis='both', width=2, length=12)
    plt.tight_layout()
    # Show the plot
    plt.show()


if __name__ == '__main__':
    # Set the density and cv that we want
    my_cv = '0.05'
    my_den = '0.05'

    # Change the directoryh
    os.chdir('../../../..')
    # Get the folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    my_rads, my_nns, my_npds = [], [], []
    for subfolder in os.listdir(folder):
        # Get the density and the cv
        cv, den = subfolder.split('_')[1], subfolder.split('_')[3]
        if cv != my_cv and den != my_den:
            continue
        # Get the files
        my_pdb, my_aw, my_pow = get_files(subfolder)
        # Get the net neighbors and the corresponding radii
        _rads, _nns, _npds = get_net_neighbors(my_pdb, my_aw, my_pow)
        my_rads += _rads
        my_nns += _nns
        my_npds += _npds

    # Plot the data
    plot_data(my_rads, my_nns, my_cv, my_den)


