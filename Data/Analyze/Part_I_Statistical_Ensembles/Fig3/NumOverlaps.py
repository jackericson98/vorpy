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
        sf_cv, sf_den = float(split_subfolder[1]), float(split_subfolder[3])
        print(sf_cv)
        # Filter out the cv values we don't want
        if cv != sf_cv:
            continue
        # Get the pdb, aw, and pow
        pdb, aw, pow = get_files(subfolder)
        # Get the logs and make a system
        my_sys = System(pdb, simple=True)
        my_aw = read_logs2(aw)
        # Create a list of overlap indices
        olap_ndxs = []
        # Loop through the balls
        for i, ball in my_aw['atoms'].iterrows():
            # Check if the ball is already accounted for
            if ball['Index'] in olap_ndxs:
                continue
            # Make sure the ball is complete
            if not ball['Complete']:
                continue
            # Get the ball location
            bloc = [ball['X'], ball['Y'], ball['Z']]
            brad = ball['Radius']
            # Go through the balls neighbors
            for neighbor in ball['Neighbors']:
                # No need to check if we know it overlaps
                if neighbor in olap_ndxs:
                    continue
                # Get the neighbors loaction
                nloc = [my_aw['atoms']['X'][neighbor], my_aw['atoms']['Y'][neighbor], my_aw['atoms']['Z'][neighbor]]
                nrad = my_aw['atoms']['Radius'][neighbor]
                # Calculate the distance
                if calc_dist(bloc, nloc) > nrad + brad:
                    # Add the ball and the neighbor and exit
                    olap_ndxs += [neighbor, ball['Index']]
                    break
        # Add the overlap indexes to the density dictionary
        if sf_den in den_data:
            den_data[sf_den].append((len(set(olap_ndxs)) / len(my_aw['atoms']['Complete'])))
        else:
            den_data[sf_den] = [(len(set(olap_ndxs)) / len(my_aw['atoms']['Complete']))]
    return den_data


def plot_olap_percent(olap_data, cv=0.05):
    # Get the data
    x, y = zip(*[(_, 100 * np.mean(olap_data[_])) for _ in olap_data])
    # Plot the bar
    plt.bar(x, y)
    # Title
    plt.title(f"CV = {cv}")
    # Show
    plt.show()


if __name__ == '__main__':
    numy = num_olaps()
    print(numy)
    plot_olap_percent(numy)
