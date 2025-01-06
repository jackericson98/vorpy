import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from System.system import System
from System.sys_funcs.calcs.calcs import calc_dist
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2


def num_olaps(folder=None, cv=0.05):
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
        # Filter out the cv values we don't want
        if cv != sf_cv:
            continue
        # Get the pdb, aw, and pow
        pdb, aw, pow = get_files(folder + '/' + subfolder)
        # Get the logs and make a system
        try:
            my_sys = System(pdb, simple=True)
            my_aw = read_logs2(aw)
        except TypeError:
            print(pdb, aw, pow)
        # Create a list of overlap indices
        olap_ndxs = []
        # Loop through the balls
        for i, ball in my_aw['atoms'].iterrows():
            # Check if the ball is already accounted for
            if ball['Index'] in olap_ndxs:
                continue
            # Make sure the ball is complete
            if not ball['Complete Cell?']:
                continue
            # Get the ball location
            bloc = [ball['X'], ball['Y'], ball['Z']]
            brad = ball['Radius']
            # Go through the balls neighbors
            for neighbor in ball['Neighbors']:
                # No need to check if we know it overlaps
                if neighbor in olap_ndxs:
                    continue
                # Get the neighbors location
                nloc = my_sys.balls['loc'][neighbor]
                nrad = my_sys.balls['rad'][neighbor]
                # Calculate the distance
                if calc_dist(bloc, nloc) < nrad + brad:
                    # Add the ball and the neighbor and exit
                    olap_ndxs += [neighbor, ball['Index']]
                    break
        print(len(olap_ndxs), len(set(olap_ndxs)), olap_ndxs)
        # Add the overlap indexes to the density dictionary
        if sf_den in den_data:
            den_data[sf_den].append((len(set(olap_ndxs)) / 1000))
        else:
            den_data[sf_den] = [(len(set(olap_ndxs)) / 1000)]
    return den_data


def plot_olap_percent(olap_data, cv=0.05):
    # Get the data
    x, y = zip(*[(_, 100 * np.mean(olap_data[_])) for _ in olap_data])
    # Plot the bar
    x, y = zip(*sorted(zip(x, y)))
    print(x, y)
    plt.scatter(x, y)
    plt.ylim([0, 100])
    # Title
    plt.title(f"CV = {cv}")
    # Show
    plt.show()


if __name__ == '__main__':
    os.chdir('../../../..')
    # numy = num_olaps(cv=1.0)
    # print(numy)
    numy = {0.45: [0.732, 0.74, 0.728, 0.736, 0.754, 0.746, 0.744, 0.756, 0.738, 0.744, 0.756, 0.746, 0.738, 0.736, 0.742, 0.742, 0.73, 0.742, 0.734, 0.756],
            0.25: [0.614, 0.63, 0.59, 0.594, 0.618, 0.612, 0.598, 0.584, 0.598, 0.612, 0.618, 0.63, 0.636, 0.592, 0.594, 0.62, 0.594, 0.592, 0.582, 0.612],
            0.35: [0.72, 0.708, 0.704, 0.698, 0.708, 0.688, 0.706, 0.706, 0.696, 0.68, 0.684, 0.708, 0.696, 0.696, 0.708, 0.688, 0.686, 0.714, 0.708, 0.692],
            0.2: [0.528, 0.54, 0.536, 0.52, 0.546, 0.566, 0.532, 0.568, 0.524, 0.514, 0.528, 0.534, 0.568, 0.524, 0.562, 0.56, 0.522, 0.538, 0.552, 0.548],
            0.1: [0.352, 0.318, 0.342, 0.338, 0.34, 0.36, 0.332, 0.342, 0.376, 0.358, 0.354, 0.346, 0.348, 0.362, 0.344, 0.368, 0.374, 0.35, 0.35, 0.368],
            0.5: [0.768, 0.748, 0.776, 0.764, 0.772, 0.73, 0.744, 0.766, 0.762, 0.73, 0.77, 0.758, 0.752, 0.75, 0.768, 0.766, 0.766, 0.758, 0.764, 0.762],
            0.15: [0.45, 0.426, 0.422, 0.454, 0.446, 0.424, 0.45, 0.454, 0.428, 0.45, 0.434, 0.442, 0.44, 0.45, 0.456, 0.436, 0.502, 0.448, 0.43, 0.446],
            0.3: [0.666, 0.66, 0.678, 0.668, 0.656, 0.64, 0.666, 0.674, 0.652, 0.68, 0.636, 0.65, 0.642, 0.65, 0.66, 0.638, 0.676, 0.64, 0.666, 0.684],
            0.05: [0.178, 0.198, 0.182, 0.186, 0.196, 0.196, 0.178, 0.2, 0.214, 0.21, 0.192, 0.214, 0.188, 0.204, 0.212, 0.194, 0.22, 0.18, 0.214, 0.178],
            0.4: [0.736, 0.716, 0.716, 0.734, 0.72, 0.746, 0.746, 0.722, 0.736, 0.732, 0.718, 0.732, 0.734, 0.714, 0.736, 0.72, 0.736, 0.744, 0.714, 0.716]}
    # numy_1 = {0.25: [0.216, 0.2, 0.178, 0.204, 0.21, 0.214, 0.222, 0.21, 0.194, 0.198, 0.192, 0.232, 0.214, 0.18, 0.224, 0.2,
    #             0.22, 0.212, 0.2],
    #      0.4: [0.266, 0.274, 0.278, 0.25, 0.272, 0.036, 0.272, 0.276, 0.288, 0.276, 0.274, 0.25, 0.266, 0.282, 0.288, 0.244,
    #            0.27, 0.288, 0.284, 0.266],
    #      0.35: [0.264, 0.232, 0.254, 0.24, 0.26, 0.238, 0.244, 0.256, 0.294, 0.244, 0.252, 0.028, 0.254, 0.264, 0.252,
    #             0.274, 0.248, 0.248, 0.22, 0.218],
    #      0.15: [0.158, 0.158, 0.144, 0.164, 0.134, 0.152, 0.14, 0.17, 0.156, 0.138, 0.176, 0.17, 0.132, 0.166, 0.016, 0.14,
    #             0.178, 0.128, 0.016, 0.158],
    #      0.5: [0.308, 0.304, 0.312, 0.298, 0.306, 0.318, 0.316, 0.306, 0.302, 0.318, 0.316, 0.292, 0.302, 0.304, 0.288,
    #            0.326, 0.312, 0.278, 0.302, 0.306],
    #      0.3: [0.212, 0.242, 0.228, 0.246, 0.226, 0.232, 0.254, 0.24, 0.24, 0.22, 0.254, 0.212, 0.23, 0.23, 0.224, 0.22,
    #            0.226, 0.244, 0.228, 0.22],
    #      0.45: [0.256, 0.322, 0.27, 0.29, 0.31, 0.274, 0.28, 0.28, 0.306, 0.292, 0.268, 0.286, 0.29, 0.29, 0.294, 0.284,
    #             0.304, 0.286, 0.284, 0.278],
    #      0.2: [0.188, 0.164, 0.178, 0.188, 0.188, 0.174, 0.164, 0.2, 0.174, 0.188, 0.174, 0.182, 0.152, 0.174, 0.18, 0.2,
    #            0.174, 0.18, 0.184, 0.208],
    #      0.1: [0.102, 0.134, 0.106, 0.118, 0.116, 0.11, 0.12, 0.126, 0.128, 0.108, 0.108, 0.116, 0.13, 0.114, 0.126, 0.11,
    #            0.134, 0.13, 0.12, 0.112],
    #      0.05: [0.09, 0.062, 0.092, 0.084, 0.074, 0.076, 0.062, 0.104, 0.08, 0.092, 0.076, 0.068, 0.084, 0.066, 0.076, 0.07,
    #             0.004, 0.074, 0.082, 0.064]}
    plot_olap_percent(numy, 0.05)
