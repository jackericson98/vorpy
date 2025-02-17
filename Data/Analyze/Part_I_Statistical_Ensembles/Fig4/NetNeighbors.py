import os
import tkinter as tk
from tkinter import filedialog
import numpy as np
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from Data.Analyze.tools.compare.read_logs import read_logs
from System.system import System
import matplotlib.pyplot as plt
import random


def get_net_neighbors(pdb, aw, pow):
    # Create the simple system
    print(pdb)
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
    if old_logs:
        my_sys = System(pdb, simple=True)
        # Go through each ball
        for i, ball in my_sys.balls.iterrows():

            # Get the pow atom and the vor atom
            pow_atom, aw_atom = pow_data['atoms'].loc[pow_data['atoms']['num'] == i - 1].to_dict('records')[0], \
            aw_data['atoms'].loc[aw_data['atoms']['num'] == i - 1].to_dict('records')[0]

            # Calculate the difference in volume
            pow_nbors, aw_nbors = len(pow_atom['neighbors']), len(aw_atom['neighbors'])
            # Get the data
            rads.append(ball['rad'])
            net_nborss.append(pow_nbors - aw_nbors)
            nbor_per_diffs.append(100 * (pow_nbors - aw_nbors) / aw_nbors)
    else:
        for i, aw_atom in aw_data['atoms'].iterrows():

            # Check that the cell is complete
            if not aw_atom['Complete Cell?']:
                continue
            # Get the pow atom and the vor atom
            try:
                pow_atom = pow_data['atoms'].loc[pow_data['atoms']['Index'] == aw_atom['Index']].to_dict('records')[0]
            except IndexError:
                continue

            # Calculate the difference in volume
            pow_nbors, aw_nbors = pow_atom['Number of Neighbors'], aw_atom['Number of Neighbors']

            # Get the data
            rads.append(aw_atom['Radius'])
            net_nborss.append(pow_nbors - aw_nbors)
            nbor_per_diffs.append(100 * (pow_nbors - aw_nbors) / aw_nbors)

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

    # Filter out the net neighbors that don't have enough data points
    num_data, excludes = len(rads), []
    for _ in mean_rad_dict:
        if len(mean_rad_dict[_]) / num_data < 0.005:
            excludes.append(_)

    # Filter out the data_points that are bad
    filtered_rads, fltrd_nnbrs = zip(*[_ for _ in zip(rads, net_nbors) if _[1] not in excludes])

    # Calculate the number of elements to keep (10% of the data)
    num_to_keep = int(len(filtered_rads) * 0.1)

    # Generate random indices to keep
    indices_to_keep = random.sample(range(len(filtered_rads)), num_to_keep)

    # Create the reduced lists
    reduced_list1 = [filtered_rads[i] for i in indices_to_keep]
    reduced_list2 = [fltrd_nnbrs[i] for i in indices_to_keep]

    # Plot the net_nbors
    plt.scatter(reduced_list1, reduced_list2, s=2, alpha=0.8)
    # Plot the average radius of each net neighbor difference
    plt.scatter([np.mean(mean_rad_dict[_]) for _ in mean_rad_dict if _ not in excludes], [_ for _ in mean_rad_dict.keys() if _ not in excludes], s=50, c='r', marker='x')
    plt.title(f'Excess Power Neighbors\nCV = {cv}, Density = {den}', fontsize=30)
    plt.xlabel('Ball Radius', fontdict=dict(size=25))
    plt.ylabel('Pow - AW Nbors', fontdict=dict(size=25))
    x_ticks = [round(_, 1) for _ in np.linspace(min(filtered_rads), max(filtered_rads), 5)]
    # Store the full sorted y_ticks for later access
    all_y_ticks = sorted([_ for _ in mean_rad_dict.keys() if _ not in excludes])

    # # Determine the range for symmetric ticks around 0
    # if len(all_y_ticks) > 1:
    #     max_abs_tick = max(abs(min(all_y_ticks)), abs(max(all_y_ticks)))
    #     num_ticks = min(5, len(all_y_ticks))  # Limit to a maximum of 5 ticks
    #     step = (2 * max_abs_tick) // (num_ticks - 1) if num_ticks > 1 else 1
    #
    #     # Generate equal integer-spaced ticks symmetrically around 0
    #     half_range = (num_ticks - 1) // 2
    #     y_ticks = [i * step for i in range(-half_range, half_range + 1)]
    # else:
    #     y_ticks = all_y_ticks  # If only one tick, keep it as is

    plt.xticks(ticks=x_ticks, font=dict(size=20))
    plt.yticks(ticks=[-10, -5, 0, 5, 10], font=dict(size=20))
    plt.xlim([x_ticks[0] - 0.5 * (x_ticks[1] - x_ticks[0]), x_ticks[-1] + 0.5 * (x_ticks[1] - x_ticks[0])])
    plt.ylim([-12, 12])

    plt.tick_params(axis='both', width=2, length=12)
    plt.tight_layout()
    # Show the plot
    plt.show()


if __name__ == '__main__':
    # Set the density and cv that we want
    my_cv = '1.0'
    my_den = '0.5'

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
        try:
            cv, den = subfolder.split('_')[1], subfolder.split('_')[3]
        except IndexError:
            continue
        if cv != my_cv or den != my_den:
            continue
        # Get the files
        my_pdb, my_aw, my_pow = get_files(folder + '/' + subfolder)
        # Get the net neighbors and the corresponding radii
        try:
            _rads, _nns, _npds = get_net_neighbors(my_pdb, my_aw, my_pow)
        except TypeError:
            continue
        my_rads += _rads
        my_nns += _nns
        my_npds += _npds

    # Plot the data
    plot_data(my_rads, my_nns, my_cv, my_den)


