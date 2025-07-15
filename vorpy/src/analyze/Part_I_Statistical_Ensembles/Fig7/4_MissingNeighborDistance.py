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
        exclusive_pw_neighbors = [neighbor for neighbor in pw_neighbors if neighbor not in ball_data[ball_index]['aw_neighbors']]
        exclusive_aw_neighbors = [neighbor for neighbor in ball_data[ball_index]['aw_neighbors'] if neighbor not in pw_neighbors]
        shared_neighbors = [neighbor for neighbor in pw_neighbors if neighbor in ball_data[ball_index]['aw_neighbors']]
        if len(exclusive_pw_neighbors) > 0:
            ball_data[ball_index]['pw_avg_nbor_dist'] = np.mean([calc_pw_center(ball_data[ball_index]['radius'], ball_data[neighbor]['radius'], ball_data[ball_index]['loc'], ball_data[neighbor]['loc'])[0] for neighbor in exclusive_pw_neighbors])
        else:
            ball_data[ball_index]['pw_avg_nbor_dist'] = 0
        if len(exclusive_aw_neighbors) > 0:
            ball_data[ball_index]['aw_avg_nbor_dist'] = np.mean([calc_aw_center(ball_data[ball_index]['radius'], ball_data[neighbor]['radius'], ball_data[ball_index]['loc'], ball_data[neighbor]['loc'])[0] for neighbor in exclusive_aw_neighbors])
        else:
            ball_data[ball_index]['aw_avg_nbor_dist'] = 0
        if len(shared_neighbors) > 0:
            ball_data[ball_index]['shared_avg_nbor_pw_dist'] = np.mean([calc_pw_center(ball_data[ball_index]['radius'], ball_data[neighbor]['radius'], ball_data[ball_index]['loc'], ball_data[neighbor]['loc'])[0] for neighbor in shared_neighbors])
            ball_data[ball_index]['shared_avg_nbor_aw_dist'] = np.mean([calc_aw_center(ball_data[ball_index]['radius'], ball_data[neighbor]['radius'], ball_data[ball_index]['loc'], ball_data[neighbor]['loc'])[0] for neighbor in shared_neighbors])
        else:
            ball_data[ball_index]['shared_avg_nbor_pw_dist'] = 0
            ball_data[ball_index]['shared_avg_nbor_aw_dist'] = 0
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
    x = [ball['radius'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    pw_avg_nbor_dist = [ball['pw_avg_nbor_dist'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])] 
    aw_avg_nbor_dist = [ball['aw_avg_nbor_dist'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    shared_avg_nbor_pw_dist = [ball['shared_avg_nbor_pw_dist'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]
    shared_avg_nbor_aw_dist = [ball['shared_avg_nbor_aw_dist'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])]

    # Calculate lines of best fit for all 4 lists
    z_pw = np.polyfit([_ for i, _ in enumerate(x) if pw_avg_nbor_dist[i] != 0], [_ for _ in pw_avg_nbor_dist if _ != 0], 1)
    z_aw = np.polyfit([_ for i, _ in enumerate(x) if aw_avg_nbor_dist[i] != 0], [_ for _ in aw_avg_nbor_dist if _ != 0], 1)
    z_shared_pw = np.polyfit([_ for i, _ in enumerate(x) if shared_avg_nbor_pw_dist[i] != 0], [_ for _ in shared_avg_nbor_pw_dist if _ != 0], 1)
    z_shared_aw = np.polyfit([_ for i, _ in enumerate(x) if shared_avg_nbor_aw_dist[i] != 0], [_ for _ in shared_avg_nbor_aw_dist if _ != 0], 1)

    # Create x values for the lines
    x_line = np.linspace(min(x), max(x), 100)

    # # Plot the lines of best fit    
    plt.plot(x_line, np.polyval(z_shared_pw, x_line), 'b-', label='Pow Shared', alpha=0.8)
    plt.plot(x_line, np.polyval(z_shared_aw, x_line), 'r-', label='AW Shared', alpha=0.8)
    plt.plot(x_line, np.polyval(z_pw, x_line), 'b--', label='Pow Only', alpha=0.8)
    plt.plot(x_line, np.polyval(z_aw, x_line), 'r--', label='AW Only', alpha=0.8)

    plt.scatter(x, shared_avg_nbor_pw_dist, c='blue', marker='o', label='Pow Shared', alpha=0.2, s=20)
    plt.scatter(x, shared_avg_nbor_aw_dist, c='red', marker='o', label='AW Shared', alpha=0.2, s=20)
    plt.scatter([_ for i, _ in enumerate(x) if pw_avg_nbor_dist[i] != 0], [_ for _ in pw_avg_nbor_dist if _ != 0], c='blue', marker='x', label='Pow Only', alpha=0.3, s=45)
    plt.scatter([_ for i, _ in enumerate(x) if aw_avg_nbor_dist[i] != 0], [_ for _ in aw_avg_nbor_dist if _ != 0], c='red', marker='x', label='AW Only', alpha=0.3, s=45)

    # plt.text(0.5, 0.5, f'Average Ball Radius: {np.mean(x):.2f}', fontsize=16, ha='center', va='center')
        
    # Add labels with different font sizes
    plt.xlabel('Radius', fontsize=25)
    plt.ylabel('Average Neighbor\nDistance', fontsize=25)
    plt.ylim(-0.5, 9)
    plt.xlim(-0.25, 5)
    # plt.xlim(0.69, 1.35)
    plt.xticks([1.0, 2.5, 4.0], fontsize=25)
    # plt.xticks([0.8, 1.0, 1.2], fontsize=25)
    plt.yticks([0.0, 4.0, 8.0], fontsize=25)

    plt.title('Neighbor Subsets &\nAverage Distance', fontsize=30)
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