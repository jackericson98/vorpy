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

def select_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    return file_path


def get_data(pdb, aw, pw):
    """
    We need dictionaries with ball index as the key and radius, neighbors, and average radius as the values.
    """
    # Read the pdb file
    pdb_data = read_pdb_simple(pdb)

    # Get the logs from the file path
    aw_logs = read_logs2(aw, all_=False, balls=True)
    pw_logs = read_logs2(pw, all_=False, balls=True)
    
    # Initialize the dictionary
    ball_data = {}


    # Go through the balls in the dictionary
    for i, ball in aw_logs['atoms'].iterrows():
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
        ball_data[ball_index] = {'radius': radius, 'aw_neighbors': neighbors, 'aw_volume': aw_volume, 'aw_sphericity': aw_sphericity, 'aw_avg_neighbor_radius': aw_avg_neighbor_radius}


    # Go through the balls in the pw dictionary
    for i, ball in pw_logs['atoms'].iterrows():
        # Get the ball index
        ball_index = ball['Index']
        # Skip if the ball index is not in the ball_data dictionary
        if ball_index not in ball_data:
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
        # Get the average neighbor radius
        pw_avg_neighbor_radius = np.mean([pdb_data[neighbor]['temperature_factor'] for neighbor in pw_neighbors])
        # Add the average neighbor radius to the dictionary
        ball_data[ball_index]['pw_avg_neighbor_radius'] = pw_avg_neighbor_radius
        
    # Return the dictionary
    return ball_data


def plot_data(ball_data):
    """
    Plot the data with a line of best fit and points colored by volume difference.
    """
    # Get the data
    data = list(ball_data.values())
    
    # Extract x and y values
    x = [ball['radius'] for ball in data]
    y = [ball['aw_avg_neighbor_radius'] for ball in data]
    colors = [np.log(abs(ball['vol_diff'])) for ball in data]

    for ball in data:
        print(ball['radius'], ball['aw_avg_neighbor_radius'], ball['vol_diff'], ball['aw_volume'], ball['pw_volume'])
    
    # Create scatter plot with colorbar
    scatter = plt.scatter(x, y, c=colors, cmap='coolwarm', alpha=0.5)
    plt.colorbar(scatter, label='Volume % Difference')
    
    # Add line of best fit
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    plt.plot(x, p(x), "r--", alpha=0.8)
    
    # Add labels
    plt.xlabel('Ball Radius')
    plt.ylabel('Average AW Neighbor Radius')
    plt.title('Ball Radius vs Average AW Neighbor Radius')
    
    # Show the plot
    plt.show()


if __name__ == '__main__':
    # Get the file path
    pdb_file = select_file()
    aw_file = select_file()
    pw_file = select_file()
    # Get the data
    ball_data = get_data(pdb_file, aw_file, pw_file)
    # Plot the data
    plot_data(ball_data)