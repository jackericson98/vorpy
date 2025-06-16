


"""
This file plots the sphericity difference as a function of std neighbor distance and neighbor polydispersity
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools import read_logs2
from vorpy.src.inputs import read_pdb_simple
from vorpy.src.calculations import calc_aw_center, calc_pw_center, calc_dist


"""
Plotting script for average neighbor radius, as a function of ball radius. Would also like to color by % difference in pow and aw volume
"""

def select_file(title=None):
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title=title)
    return file_path


def get_data(pdb, aw, pw):
    """
    We need dictionaries with ball index as the key and radius, neighbors, and average radius as the values.
    """
    # Get the box size
    with open(pdb, 'r') as f:
        box_size = float(f.readline().split()[2])
    # Read the pdb file
    pdb_data = read_pdb_simple(pdb)

    # Get the logs from the file path
    aw_logs = read_logs2(aw, all_=False, balls=True)
    pw_logs = read_logs2(pw, all_=False, balls=True)
    
    # Initialize the dictionary
    ball_data = {k: {'radius': pdb_data[k]['radius'], 'loc': np.array([pdb_data[k]['x_coordinate'], pdb_data[k]['y_coordinate'], pdb_data[k]['z_coordinate']])} for k in pdb_data.keys()}

    # Go through the balls in the dictionary
    for i, ball in aw_logs['atoms'].iterrows():
        # Get the ball location
        loc = np.array([ball['X'], ball['Y'], ball['Z']])
        # Get the ball index
        ball_index = ball['Index']
        # Get the radius
        radius = ball['Radius']
        # Get the neighbors
        neighbors = ball['Neighbors']
        # Get the aw volume
        aw_volume = ball['Volume']
        # Get the aw sphericity
        aw_sphericity = ball['Sphericity']
        if aw_sphericity > 1:
            continue
        # Get the average neighbor radius
        aw_avg_neighbor_radius = np.mean([pdb_data[neighbor]['temperature_factor'] for neighbor in neighbors])
        # Add the ball data to the dictionary
        ball_data[ball_index] = {'radius': radius, 'loc': loc, 'aw_neighbors': neighbors, 'aw_volume': aw_volume, 'aw_sphericity': aw_sphericity, 'aw_avg_neighbor_radius': aw_avg_neighbor_radius}


    # Go through the balls in the pw dictionary
    for i, ball in pw_logs['atoms'].iterrows():
        # Get the ball index
        ball_index = ball['Index']
        # Skip if the ball index is not in the ball_data dictionary
        if ball_index not in ball_data or 'aw_volume' not in ball_data[ball_index]:
            continue
        # Get the ball_data
        pw_neighbors = ball['Neighbors']
        pw_volume = ball['Volume']
        pw_sphericity = ball['Sphericity']
        aw_volume = ball_data[ball_index]['aw_volume']
        # Get the % difference in volume
        vol_diff = 100 * (pw_volume - aw_volume) / aw_volume
        # Add the ball data to the dictionary
        ball_data[ball_index]['pw_neighbors'] = pw_neighbors
        ball_data[ball_index]['pw_volume'] = pw_volume
        ball_data[ball_index]['vol_diff'] = vol_diff
        ball_data[ball_index]['pw_sphericity'] = pw_sphericity
        ball_data[ball_index]['sphericity_diff'] = pw_sphericity - aw_sphericity
        pw_nbor_dists = [calc_pw_center(ball_data[ball_index]['radius'], ball_data[neighbor]['radius'], ball_data[ball_index]['loc'], ball_data[neighbor]['loc'])[0] for neighbor in pw_neighbors]
        aw_nbor_dists = [calc_aw_center(ball_data[ball_index]['radius'], ball_data[neighbor]['radius'], ball_data[ball_index]['loc'], ball_data[neighbor]['loc'])[0] for neighbor in ball_data[ball_index]['aw_neighbors']]
        ball_data[ball_index]['pw_nbor_dist_cv'] = np.std(pw_nbor_dists) / np.mean(pw_nbor_dists)
        ball_data[ball_index]['aw_nbor_dist_cv'] = np.std(aw_nbor_dists) / np.mean(aw_nbor_dists)
        ball_data[ball_index]['nbor_dist_cv_diff'] = ball_data[ball_index]['pw_nbor_dist_cv'] - ball_data[ball_index]['aw_nbor_dist_cv']
        
        # Get the average neighbor radius
        pw_avg_neighbor_radius = np.mean([pdb_data[neighbor]['temperature_factor'] for neighbor in pw_neighbors])
        # Add the average neighbor radius to the dictionary
        ball_data[ball_index]['pw_avg_neighbor_radius'] = pw_avg_neighbor_radius
        
    # Return the dictionary
    return {'box_size': box_size, 'balls': ball_data}


def plot_data(ball_data):
    """
    Plot the data with a line of best fit and points colored by volume difference.
    """

    plt.figure(figsize=(6, 5))

    # Get the data
    data = list(ball_data['balls'].values())
    
    # Extract x and y values
    y1 = [ball['sphericity_diff'] for ball in data if all(key in ball for key in ['aw_sphericity', 'pw_sphericity', 'pw_nbor_dist_cv', 'aw_nbor_dist_cv'])]
    x1 = [ball['nbor_dist_cv_diff'] for ball in data if all(key in ball for key in ['aw_sphericity', 'pw_sphericity', 'pw_nbor_dist_cv', 'aw_nbor_dist_cv'])]
    colors = [ball['radius'] for ball in data if all(key in ball for key in ['aw_sphericity', 'pw_sphericity', 'pw_nbor_dist_cv', 'aw_nbor_dist_cv'])]

    plt.figure(figsize=(6, 5))

    scatter1 = plt.scatter(x1, y1, c=colors, alpha=1.0, s=10)
    
    # Calculate and plot line of best fit
    z = np.polyfit(x1, y1, 1)
    p = np.poly1d(z)
    x_line = np.linspace(-0.15, 0.15, 100)
    plt.plot(x_line, p(x_line), 'r-', linewidth=2, alpha=0.8)
    
    # Calculate Spearman's rank correlation
    from scipy import stats
    spearman_corr, p_value = stats.spearmanr(x1, y1)
    
    # Add correlation coefficient to plot
    plt.text(0.05, 0.95, f'Spearman ρ = {spearman_corr:.3f}\np = {p_value:.3e}', 
             transform=plt.gca().transAxes, fontsize=16, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    scatter1.set_alpha(0.2)

    # Fit the data for aw_y with linear fit and plot it
    # z_aw = np.polyfit(x1, y1, 1)
    # aw_p = np.poly1d(z_aw)
    
    # window_size = 0.05
    # min_num_points = 3
    # x_fill = np.array([])
    # means = np.array([])
    # std_devs = np.array([])

    # # Loop through the data and get the standard deviations
    # for i in np.arange(-0.18 - window_size/2, 0.18 + window_size/2, window_size):
    #     # Get the window
    #     window_mask = (x1 >= i - window_size/2) & (x1 <= i + window_size/2)
    #     y_window_points = np.array(y1)[window_mask]

    #     # Only process if we have enough points
    #     if len(y_window_points) >= min_num_points:
    #         # Get the x value for filling between
    #         x_fill = np.append(x_fill, i + window_size/2)
    #         # Get the mean and standard deviation
    #         means = np.append(means, np.mean(y_window_points))
    #         std_devs = np.append(std_devs, np.std(y_window_points))
    # Plot the x and y axes
    plt.axhline(0, color='black', linewidth=1, zorder=0)
    plt.axvline(0, color='black', linewidth=1, zorder=0)

    # Plot the line and confidence interval
    # plt.plot(x_fill - window_size/2, means, 'r-', linewidth=4, alpha=0.5)
    # plt.fill_between(x_fill - window_size/2, means - 2*std_devs, means + 2*std_devs, color='red', alpha=0.1)


    # Add labels with different font sizes
    plt.ylabel('Sphericity Difference', fontsize=25)
    plt.xlabel('Neighbor Distance CV Difference', fontsize=25)
    plt.xlim(-0.3, 0.3)
    plt.ylim(-0.2, 0.2)
    plt.yticks([-0.15, 0.0, 0.15], fontsize=25)
    plt.xticks([-0.20, 0.0, 0.20], fontsize=25)
    plt.tick_params(axis='both', width=2, length=10)
    plt.title('Sphericity Difference vs.\nNeighbor Distance CV Difference', fontsize=30)
    # plt.legend(fontsize=20, markerscale=3)
    plt.tight_layout()
    # Show the plot
    plt.show()


if __name__ == '__main__':
    #get the folder they are in
    folder = filedialog.askdirectory(title='Select the folder')
    print(folder)
    # Get the data
    data = get_data(folder + '/balls.pdb', folder + '/aw/aw_logs.csv', folder + '/pow/pow_logs.csv')
    # Plot the data
    plot_data(data)