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
        ball_data[ball_index]['pw_neighbor_cv'] = np.std(pw_neighbors) / np.mean(pw_neighbors)
        ball_data[ball_index]['aw_neighbor_cv'] = np.std(ball_data[ball_index]['aw_neighbors']) / np.mean(ball_data[ball_index]['aw_neighbors'])
        ball_data[ball_index]['pw_volume'] = pw_volume
        ball_data[ball_index]['vol_diff'] = vol_diff
        ball_data[ball_index]['pw_sphericity'] = pw_sphericity
        ball_data[ball_index]['sphericity_diff'] = pw_sphericity - aw_sphericity
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
    # Get the data
    data = list(ball_data['balls'].values())
    
    # Extract x and y values
    x = [ball['radius'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])][::-1]
    y = [abs(ball['vol_diff']) for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])][::-1] 
    
    # Create scatter plot with colorbar
    plt.scatter(x, y, c='red', s=20, alpha=0.2)

    # Fit the data for aw_y with quadratic curve
    # aw_z = np.polyfit(x, aw_y, 2)
    # aw_p = np.poly1d(aw_z)
    # plt.plot(x, aw_p(x), 'r--', linewidth=2, label='AW fit')

    # # Fit the data for pw_y with quadratic curve
    # pw_z = np.polyfit(x, pw_y, 2)
    # pw_p = np.poly1d(pw_z)
    # plt.plot(x, pw_p(x), 'b--', linewidth=2, label='Pow fit')

    
    # Add labels with different font sizes
    plt.xlabel('Radius', fontsize=20)
    plt.ylabel('Volume Difference', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    # plt.legend(fontsize=16)
    plt.title('Volume Difference vs. Radius', fontsize=20)
    plt.tight_layout()
    # Show the plot
    plt.show()


if __name__ == '__main__':
    #get the folder they are in
    folder = filedialog.askdirectory(title='Select the folder')

    # Get the data
    data = get_data(folder + '/balls.pdb', folder + '/aw/aw_logs.csv', folder + '/pow/pow_logs.csv')
    # Plot the data
    plot_data(data)