"""Calculates the volume of overlap between the power and the additively weighted surfaces"""
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
        # Get the average neighbor radius
        aw_avg_neighbor_radius = np.mean([pdb_data[neighbor]['temperature_factor'] for neighbor in neighbors])
        # Add the ball data to the dictionary
        ball_data[ball_index] = {'radius': radius, 'loc': loc, 'aw_neighbors': neighbors, 'aw_volume': aw_volume, 'aw_sphericity': aw_sphericity, 'aw_avg_neighbor_radius': aw_avg_neighbor_radius}
        # Get the closest neighbor
        ball_data[ball_index]['closest_neighbor'] = ball['Closest Neighbor']
        # Get the closest neighbor % overlap
        cn_rad, cn_loc = ball_data[ball['Closest Neighbor']]['radius'], ball_data[ball['Closest Neighbor']]['loc']
        cn_dist = calc_dist(loc, cn_loc)
        cn_overlap = max(100 * (-cn_dist + cn_rad + radius) / radius, 0)
        ball_data[ball_index]['cn_overlap'] = cn_overlap

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
        print(vol_diff)
        # Add the ball data to the dictionary
        ball_data[ball_index]['pw_neighbors'] = pw_neighbors
        ball_data[ball_index]['pw_volume'] = pw_volume
        ball_data[ball_index]['vol_diff'] = vol_diff
        
        
    # Return the dictionary
    return {'box_size': box_size, 'balls': ball_data}


def plot_data(ball_data):
    """
    Plot the data with a line of best fit and points colored by volume difference.
    """
    # Get the data
    data = list(ball_data['balls'].values())
    
    # Extract x and y values
    rads = [ball['radius'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    cn_overlap = [ball['cn_overlap'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    vol_diff = [np.log10(abs(ball['vol_diff']) + 1e-10) for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    num_neighbors = [len(ball['pw_neighbors']) for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    # print(vol_diff)
    # cn_overlap = [0 for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    # plt.figure(figsize=(5, 8))
    # Calculate average vol_diff for each integer number of neighbors
    # Count frequency of each neighbor number
    neighbor_counts = {}
    for n in num_neighbors:
        neighbor_counts[n] = neighbor_counts.get(n, 0) + 1
    
    
    # Plot histogram of neighbor counts on the second y-axis
    plt.bar(neighbor_counts.keys(), [_ / 100 for _ in neighbor_counts.values()], 
            alpha=0.5, color='gray', edgecolor='black', bottom=-1)
    # plt.set_ylabel('Count', fontsize=25)
    plt.tick_params(axis='y', labelsize=25, length=10, width=2)
    # plt.set_ylim(0, 100)

    
    # Calculate and plot averages using a 5-point window
    neighbor_avgs = {}
    for i in range(0, 75, 5):  # Step by 5 from 0 to 70
        # Get indices where num_neighbors is within the window
        indices = [j for j, n in enumerate(num_neighbors) if i <= n < i + 5]
        if indices:  # Only calculate if we have data for this window
            avg_vol_diff = np.mean([vol_diff[j] for j in indices])
            neighbor_avgs[i] = avg_vol_diff

    # Plot the averages
    x_vals = list(neighbor_avgs.keys())
    y_vals = list(neighbor_avgs.values())
    plt.plot(x_vals, y_vals, 'k-', linewidth=2, alpha=0.8)
    # plt.legend(fontsize=16)

    # Plot the scatter plot with the radii as colors
    scatter = plt.scatter(num_neighbors, vol_diff, c=cn_overlap, s=10, alpha=0.5, vmin=0, vmax=100)
    cbar = plt.colorbar(scatter, label='Closest Neighbor\nOverlap %', alpha=1.0)
    cbar.ax.tick_params(labelsize=25, length=10, width=2)
    cbar.ax.set_yticks([20, 50, 80])

    
    # Add labels with different font sizes
    plt.xlabel('Number of Neighbors', fontsize=25)
    plt.ylabel('Abs Vol Diff', fontsize=25)
    plt.ylim(-1, 4.5)
    # # plt.xlim(-0.25, 5)
    # plt.xlim(0.69, 1.35)
    plt.xticks([10, 30, 50], fontsize=25)
    # plt.xticks([0.8, 1.0, 1.2], fontsize=25)
    plt.yticks([0.0, 2.0, 4.0], [f'1', f'10\u00b3', f'10\u2075'], fontsize=25)
    plt.tick_params(axis='both', length=10, width=2)
    plt.grid(True)
    plt.title('Number of Neighbors &\nAverage Volume Difference', fontsize=30)
    # plt.legend(fontsize=16)
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