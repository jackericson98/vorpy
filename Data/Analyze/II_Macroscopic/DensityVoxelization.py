import csv
import math
import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
import numpy as np
from numpy import random, array
from System.sys_funcs.calcs.calcs import calc_dist
from System.sys_funcs.calcs.sorting import get_balls, box_search
from Data.Analyze.tools.plot_templates.histogram import histogram


# def record_density(bubbles, box, n_samples=100000):
#     """
#     Takes in bubbles and the box we are searching and finds a monte carlo density
#     """
#     # Start the count
#     count = 0
#     # Loop through the number of points to look for density/overlaps
#     for i in range(n_samples):
#         # Get the random point within the box
#         point = [random.uniform(min(box[0][i], box[1][i]), max(box[0][i], box[1][i])) for i in range(3)]
#         # Loop through the bubbles in the set and look for if any overlap
#         for j, bub in enumerate(bubbles):
#             # Check for overlap
#             if calc_dist(array(point), bub[0]) < bub[1]:
#                 count += 1
#                 break
#     return count / n_samples
#
#
# def get_boxes_density(sys, num_boxes, n_samples=10000):
#     # Get the number of splits per side
#     num_side_splits = int(np.cbrt(num_boxes)) + 1
#     # Get the large box
#     big_box = sys.groups[0].net.box['verts']
#     # Get the box side lengths
#     x_size, y_size, z_size = [(big_box[1][i] - big_box[0][i]) / num_side_splits for i in range(3)]
#     distance_reach = max(x_size, y_size, z_size)
#     max_ball = max(sys.groups[0].net.balls['rad'])
#     # Create the densities list
#     densities = []
#     count = 0
#     # Loop through the x
#     for i in range(num_side_splits):
#         for j in range(num_side_splits):
#             for k in range(num_side_splits):
#                 count += 1
#                 print("\r{} %".format(100 * round(count / num_side_splits ** 3, 5)), end="")
#                 # Create the box
#                 box = [[i * x_size + big_box[0][0], j * y_size + big_box[0][1], k * z_size + big_box[0][2]],
#                        [(i + 1) * x_size + big_box[0][0], (j + 1) * y_size + big_box[0][1], (k + 1) * z_size + big_box[0][2]]]
#                 # Get the initial box to stem from
#                 loc_sub_box = box_search([(box[0][l] + box[1][l]) / 2 for l in range(3)])
#                 # print(loc_sub_box)
#                 # Get the bubbles in the range of the box
#                 close_balls = get_balls([loc_sub_box], dist=distance_reach * max_ball)
#                 # Create the balls list
#                 balls = [(sys.balls['loc'][_], sys.balls['rad'][_]) for _ in close_balls]
#                 # Get the density
#                 densities.append(record_density(balls, box, n_samples))
#     return densities


import numpy as np
import math


def estimate_bounding_box_overlap(sphere_center, radius, sub_box_min, sub_box_max):
    """
    Estimate the fraction of the sphere that overlaps with the sub-box
    based on bounding box intersection.
    """
    # Initialize the overlap volume as 0
    overlap_volume = 0.0

    # Calculate the bounding box of the sphere
    sphere_box_min = [sphere_center[i] - radius for i in range(3)]
    sphere_box_max = [sphere_center[i] + radius for i in range(3)]

    # Find the intersection between the sphere's bounding box and the sub-box
    intersect_min = [max(sub_box_min[i], sphere_box_min[i]) for i in range(3)]
    intersect_max = [min(sub_box_max[i], sphere_box_max[i]) for i in range(3)]

    # Calculate the size of the intersection box
    intersect_sizes = [max(0, intersect_max[i] - intersect_min[i]) for i in range(3)]

    # If there is any overlap
    if all(size > 0 for size in intersect_sizes):
        # Estimate the volume of the intersection box
        intersect_volume = intersect_sizes[0] * intersect_sizes[1] * intersect_sizes[2]

        # Calculate the volume of the sphere
        sphere_volume = (4 / 3) * np.pi * radius ** 3

        # Estimate the fraction of the sphere's volume inside the sub-box
        overlap_volume = (intersect_volume / (2 * radius) ** 3) * sphere_volume

    return overlap_volume


def find_sub_box_density(spheres, box=None, num_divisions=27):
    # Find the subdivisions for the cube
    num_ax_divs = math.ceil(np.cbrt(num_divisions))

    # If no box is given, determine the min and max x, y, z
    if box is None:
        # Set up the mins and maxs variables
        mins, maxs = [np.inf, np.inf, np.inf], [-np.inf, -np.inf, -np.inf]
        # Loop through each of the spheres
        for loc, rad in spheres:
            # Loop through x, y, z
            for i in range(3):
                # Check if the min needs to be replaced
                if loc[i] - rad < mins[i]:
                    mins[i] = loc[i] - rad
                # Check if the max needs to be replaced
                if loc[i] + rad > maxs[i]:
                    maxs[i] = loc[i] + rad
        # set the box values
        box = mins, maxs

    # First get the sub box size
    sub_sizes = [(box[1][i] - box[0][i]) / num_ax_divs for i in range(3)]

    # Create the densities and indices list
    densities, ndxs = [], []

    # Volume of a sub-box
    sub_box_volume = sub_sizes[0] * sub_sizes[1] * sub_sizes[2]

    # Create the count variable
    count = 0

    # Now loop through each of the sub-boxes
    for i in range(num_ax_divs):
        for j in range(num_ax_divs):
            for k in range(num_ax_divs):

                # Print the progress
                count += 1
                # print("\r{} %".format(100 * round(count / num_ax_divs ** 3, 5)), end="")

                # Create the indices variable
                indices = (i, j, k)

                # Get the sub-box boundaries
                sub_box_min = [indices[m] * sub_sizes[m] + box[0][m] for m in range(3)]
                sub_box_max = [(indices[m] + 1) * sub_sizes[m] + box[0][m] for m in range(3)]
                # print(sub_box_min, sub_box_max)

                # Initialize the density for this sub-box
                density = 0.0

                # Loop through each sphere
                for (sphere_center, radius) in spheres:

                    # Calculate the overlap between the sphere and the sub-box
                    overlap_volume = estimate_bounding_box_overlap(sphere_center, radius, sub_box_min, sub_box_max)

                    # Add the proportional overlap volume to the density
                    density += overlap_volume / sub_box_volume

                    # Check for density larger than 1
                    if density >= 1:
                        density = 1.0
                        break

                # Add the density to the density list
                densities.append(density)
                # Add the indices to the ndxs list
                ndxs.append(indices)

    return ndxs, densities, box


if __name__ == '__main__':

    """
    We want to get the local density across a voxelized system, to know if there are pockets of extreme density and 
    pockets of low density
    """
    os.chdir('../../..')
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    logs_pdbs = {}
    folder = filedialog.askdirectory()
    for rroot, directories, files in os.walk(folder):
        for directory in directories:
            if 'aw' in directory or 'pow' in directory:
                continue
            logs_pdbs[directory] = {}
            for rrooot, dircs, filese in os.walk(rroot + '/' + directory):
                for file in filese:
                    if file[-3:] == 'pdb':
                        split_file = file.split('_')
                        logs_pdbs[directory]['pdb'] = rrooot + '/' + file
                        logs_pdbs[directory]['cv'] = split_file[3]
                        logs_pdbs[directory]['dens'] = split_file[1]

    def change_open_cell_density(dsty, n):
        if n < 10:
            a, b, c = 1.3379916075144123, 1.5342161903140346, -0.049539624149683714
        elif n < 100:
            a, b, c = 2.617546924616389, 0.5263832332311094, 0.023304967275296712
        elif n < 1000:
            a, b, c = 1.6102372828978275, 0.7263830756664165, 0.014383340152633836
        elif n < 10000:
            a, b, c = 1.333495416302907, 0.763842119925516, 0.014253487381526
        else:
            a, b, c = 1.2139394583920557, 0.7808992862460553, 0.014317251455209012
        return a * dsty ** 2 + b * dsty + c


    with open('densy_wensities.csv', 'w') as my_big_dawg_file:
        num_sub_boxes = 1000
        my_writer = csv.writer(my_big_dawg_file)
        my_writer.writerow(['File', 'Actual Density', 'Calculated Density (Unadjusted)', 'Calculated Density (Adjusted)', 'Radii CV', 'Number of Sub-Boxes', 'Box Max', 'Values'])
        for i, directory in enumerate(logs_pdbs):
            print("\r Calculating {}/{} --> {} %".format(i, len(logs_pdbs), 100 * round(i / len(logs_pdbs), 4)), end="")
            try:
                file = logs_pdbs[directory]['pdb']
            except KeyError:
                continue
            my_sys = System(file)
            locs = my_sys.balls['loc'].to_list()
            rads = my_sys.balls['rad'].to_list()
            ndxs, sub_densities, my_box = find_sub_box_density([*zip(locs, rads)], None, num_sub_boxes)

            vol = (my_box[1][0] - my_box[0][0]) * (my_box[1][1] - my_box[0][1]) * (my_box[1][2] - my_box[0][2])
            sub_box_vol = vol / len(sub_densities)
            tot_density = sum([sub_box_vol * _ for _ in sub_densities]) / vol
            adju_tot_dens = change_open_cell_density(float(logs_pdbs[directory]['dens']), 1000)
            my_writer.writerow([logs_pdbs[directory]['pdb'], logs_pdbs[directory]['dens'], tot_density, adju_tot_dens,
                                logs_pdbs[directory]['cv'], num_sub_boxes, max(my_box[1])] + sub_densities)



    #
    # my_sys = System(file)
    # my_sys.create_group(make_net=True)
    #

    # my_densities = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.003, 0.013, 0.002, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023, 0.012, 0.015, 0.006, 0.002, 0.005, 0.002, 0.0, 0.0, 0.0, 0.0, 0.005, 0.071, 0.034, 0.026, 0.039, 0.035, 0.006, 0.0, 0.0, 0.0, 0.0, 0.002, 0.018, 0.012, 0.007, 0.007, 0.004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.012, 0.025, 0.019, 0.005, 0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006, 0.015, 0.049, 0.055, 0.014, 0.059, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.001, 0.006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.016, 0.046, 0.063, 0.056, 0.004, 0.002, 0.0, 0.0, 0.0, 0.0, 0.071, 0.16, 0.389, 0.715, 0.345, 0.338, 0.015, 0.0, 0.0, 0.0, 0.0, 0.003, 0.403, 0.304, 0.561, 0.493, 0.355, 0.07, 0.0, 0.0, 0.0, 0.0, 0.058, 0.181, 0.23, 0.257, 0.211, 0.383, 0.0, 0.0, 0.0, 0.0, 0.0, 0.062, 0.43, 0.225, 0.503, 0.363, 0.107, 0.014, 0.0, 0.0, 0.0, 0.0, 0.008, 0.461, 0.401, 0.616, 0.263, 0.231, 0.055, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.063, 0.091, 0.028, 0.036, 0.008, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.077, 0.041, 0.045, 0.012, 0.001, 0.004, 0.0, 0.0, 0.0, 0.0, 0.03, 0.303, 0.473, 0.645, 0.133, 0.283, 0.004, 0.0, 0.0, 0.0, 0.0, 0.005, 0.456, 0.194, 0.397, 0.358, 0.419, 0.05, 0.0, 0.0, 0.0, 0.0, 0.008, 0.198, 0.396, 0.486, 0.494, 0.469, 0.014, 0.0, 0.0, 0.0, 0.0, 0.071, 0.403, 0.204, 0.421, 0.505, 0.483, 0.035, 0.0, 0.0, 0.0, 0.0, 0.077, 0.687, 0.627, 0.424, 0.575, 0.397, 0.015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.179, 0.083, 0.003, 0.018, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.005, 0.034, 0.014, 0.014, 0.042, 0.034, 0.0, 0.0, 0.0, 0.0, 0.0, 0.077, 0.238, 0.279, 0.579, 0.473, 0.355, 0.048, 0.0, 0.0, 0.0, 0.0, 0.036, 0.399, 0.623, 0.876, 0.574, 0.445, 0.054, 0.0, 0.0, 0.0, 0.0, 0.053, 0.298, 0.927, 0.78, 0.716, 0.448, 0.067, 0.0, 0.0, 0.0, 0.0, 0.084, 0.958, 0.461, 0.558, 0.416, 0.329, 0.069, 0.0, 0.0, 0.0, 0.0, 0.017, 0.304, 0.141, 0.328, 0.649, 0.241, 0.009, 0.0, 0.0, 0.0, 0.0, 0.003, 0.012, 0.003, 0.022, 0.018, 0.003, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.004, 0.064, 0.054, 0.042, 0.051, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.433, 0.464, 0.761, 0.724, 0.192, 0.023, 0.0, 0.0, 0.0, 0.0, 0.004, 0.279, 0.325, 0.799, 0.852, 0.417, 0.059, 0.0, 0.0, 0.0, 0.0, 0.038, 0.219, 0.399, 0.75, 0.65, 0.41, 0.029, 0.0, 0.0, 0.0, 0.0, 0.033, 0.606, 0.222, 0.461, 0.262, 0.566, 0.037, 0.0, 0.0, 0.0, 0.0, 0.031, 0.247, 0.221, 0.329, 0.333, 0.297, 0.006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007, 0.009, 0.049, 0.028, 0.034, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.017, 0.039, 0.019, 0.067, 0.015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.002, 0.2, 0.235, 0.266, 0.809, 0.417, 0.008, 0.0, 0.0, 0.0, 0.0, 0.061, 0.279, 0.374, 0.849, 0.4, 0.393, 0.002, 0.0, 0.0, 0.0, 0.0, 0.056, 0.348, 0.3, 0.609, 0.372, 0.191, 0.017, 0.0, 0.0, 0.0, 0.0, 0.13, 0.694, 0.372, 0.496, 0.286, 0.456, 0.039, 0.0, 0.0, 0.0, 0.0, 0.011, 0.428, 0.278, 0.632, 0.284, 0.315, 0.053, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.021, 0.03, 0.019, 0.031, 0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007, 0.0, 0.023, 0.049, 0.024, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.008, 0.106, 0.03, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.025, 0.055, 0.038, 0.029, 0.002, 0.0, 0.0, 0.0, 0.0, 0.014, 0.102, 0.022, 0.007, 0.058, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.009, 0.023, 0.036, 0.043, 0.007, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]



# def find_sub_box_density(spheres, box, num_divisions=27):
#     # Find the subdivisions for the cube
#     num_ax_divs = math.ceil(np.cbrt(num_divisions))
#     # First get the sub box size
#     sub_sizes = [(box[1][i] - box[0][i]) / num_ax_divs for i in range(3)]
#     # Create the densities and indices list
#     densities, ndxs = [], []
#     # Now loop through each of the boxes
#     for i in range(num_ax_divs):
#         for j in range(num_ax_divs):
#             for k in range(num_ax_divs):
#                 # Create the indices variable
#                 indices = i, j, k
#                 # Get the sub box locations
#                 sub_box = ([indices[m] * sub_sizes[m] + box[0][m] for m in range(3)],
#                            [(indices[m] + 1) * sub_sizes[m] + box[0][m] for m in range(3)])
#                 # Get the density of the sub box
#                 density = None  # ChatGPT ----> Replace None with the formula you come up with for sub box density from spheres
#                 # Add the density to the density list
#                 densities.append(density)
#                 # Add the indices to the ndxs list
#                 ndxs.append(indices)
#     return ndxs, densities










