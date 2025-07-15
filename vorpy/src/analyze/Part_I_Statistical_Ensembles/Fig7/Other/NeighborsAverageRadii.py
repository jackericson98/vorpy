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
    y = [ball['aw_avg_neighbor_radius'] for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])][::-1] 
    colors = [np.log(abs(ball['vol_diff'])) for ball in data if all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff'])][::-1]

        # Add curved line of best fit with smoothed local standard deviation interval
    x_new = np.array(x)  # Convert lists to numpy arrays
    y_new = np.array(y)
    
    # Sort data by x values to avoid interpolation issues
    sort_idx = np.argsort(x_new)
    x_new = x_new[sort_idx]
    y_new = y_new[sort_idx]
    
    # Remove duplicate x values by averaging corresponding y values
    unique_x_new, unique_indices_new = np.unique(x_new, return_index=True)
    unique_y_new = np.array([np.mean(y_new[x_new == val]) for val in unique_x_new])
    
    z = np.polyfit(unique_x_new, unique_y_new, 2)  # Use 2nd degree polynomial for curved fit
    p = np.poly1d(z)
    
    # Calculate local standard deviations using rolling window
    window_size = 20  # Adjust window size as needed
    y_pred_new = p(unique_x_new)
    std_smooth_new = []
    
    # Calculate standard deviation every 5 points for smoother interval
    step = 15
    x_smooth_new = unique_x_new[::step]
    
    for i in range(0, len(unique_x_new), step):
        # Get indices of points within window
        window_indices_new = np.where(np.abs(unique_x_new - unique_x_new[i]) <= (np.max(unique_x_new) - np.min(unique_x_new))/window_size)[0]
        # Calculate standard deviation of residuals in window
        window_residuals_new = unique_y_new[window_indices_new] - y_pred_new[window_indices_new]
        std_smooth_new.append(np.std(window_residuals_new))
    
    # Interpolate the standard deviations back to original x points
    from scipy.interpolate import interp1d
    f_new = interp1d(x_smooth_new, std_smooth_new, kind='cubic', bounds_error=False, fill_value=(std_smooth_new[0], std_smooth_new[-1]))
    local_std_new = f_new(unique_x_new)
    
    # Plot the curved fit and smoothed local standard deviation interval
    plt.plot(unique_x_new, p(unique_x_new), "r--", alpha=0.8, label='Curved Fit', linewidth=2)
    plt.fill_between(unique_x_new, p(unique_x_new) - 2*local_std_new, p(unique_x_new) + 2*local_std_new, color='red', alpha=0.35, label='±2σ Local Interval', linewidth=2)
    # plt.legend()

    for i in range(len(x)):
        print(x[i], y[i], colors[i])

    for ball in data:
        if not all(key in ball for key in ['radius', 'aw_avg_neighbor_radius', 'vol_diff']):
            continue
        print(ball['radius'], ball['aw_avg_neighbor_radius'], ball['vol_diff'], ball['aw_volume'], ball['pw_volume'])
    
    # Set font sizes
    plt.rcParams.update({'font.size': 20})
    min_val = np.min([np.log(abs(ball['vol_diff'])) for ball in data if all(key in ball for key in ['radius', 'vol_diff'])])
    max_val = np.max([np.log(abs(ball['vol_diff'])) for ball in data if all(key in ball for key in ['radius', 'vol_diff'])])
    
    # Create scatter plot with colorbar
    scatter = plt.scatter(x, y, c=colors, cmap='viridis', s=4, vmin=min_val, vmax=max_val)
    cbar = plt.colorbar(scatter, label='Volume % Diff', extend='both')
    
    # Set custom ticks to show log scale
    ticks = [-2, 0, 2, 4, 6]
    cbar.set_ticks(ticks)
    num_ticks = 5
    tick_locs = np.linspace(-2, 6, num_ticks)
    print(max_val, min_val)
    
    # Create labels with superscript for even indices
    labels = []
    uni_dict = {-3: '\u207B\u00B3', -2: '\u207B\u00B2', -1: '\u207B\u00B9', 0: '\u2070', 1: '\u00B9', 2: '\u00B2', 3: '\u00B3', 4: '\u2074', 5: '\u2075', 6: '\u2076', 7: '\u2077', 8: '\u2078', 9: '\u2079'}
    for i, tick in enumerate(ticks):
        labels.append(f'10{uni_dict[int(tick)]}')
        # labels.append(f'{tick}')
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(labels)
    

    
    # Add labels with different font sizes
    plt.xlabel('Radius', fontsize=20)
    plt.ylabel('Avg Neighbor Rad', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylim(0.5, 4.0)
    plt.title('Avg Neighbor Rad', fontsize=30)
    plt.tight_layout()
    # Show the plot
    plt.show()


if __name__ == '__main__':
    # Get the data
    data = get_data(select_file(title='Select the PDB file'), select_file(title='Select the AW log file'), select_file(title='Select the PW log file'))
    # Plot the data
    plot_data(data)