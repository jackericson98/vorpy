import matplotlib.pyplot as plt
import numpy as np

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.network.fill import calc_surf_point
from vorpy.src.calculations.surf import calc_surf_func
from vorpy.src.calculations.calcs import calc_dist


def get_vor_points(l1, r1, l2, r2, res=0.1, height=5):
    """Must be opposite across the y axis"""
    l1, l2 = np.array(l1 + [0]), np.array(l2 + [0])
    func = calc_surf_func(l1, r1, l2, r2)
    x_pts, y_pts = [], []
    dist_between = calc_dist(l1, l2)
    x_center = l1[0] + 0.5 * (dist_between - r1 - r2) + r1

    if r1 == r2:
        x_pts = [x_center, x_center]
        y_pts = [height, - height]
        proj_xs, proj_ys = x_pts, y_pts
        return x_pts, y_pts, proj_xs, proj_ys
    direction = np.array([-r1, r1, 0])
    ray_hat = res * (direction / np.linalg.norm(direction))
    proj_xs, proj_ys = [], []
    my_pt = np.array([x_center, 0, 0])
    proj_pt = my_pt
    while my_pt[1] < height and len(x_pts) < 500:
        proj_pt = proj_pt + ray_hat
        my_pt = calc_surf_point([l1, l2], proj_pt, func)
        proj_xs.append(proj_pt[0])
        proj_ys.append(proj_pt[1])
        x_pts.append(my_pt[0])
        y_pts.append(my_pt[1])
    x_pts = x_pts[::-1] + x_pts
    y_pts = [-_ for _ in y_pts[::-1]] + y_pts
    return x_pts, y_pts, proj_xs, proj_ys


def plot_voronoi_partition(circle1, circle2, height=5):
    """
    Plots two circles and their Voronoi partition.

    Parameters:
    circle1 (tuple): (x, y, radius) of the first circle
    circle2 (tuple): (x, y, radius) of the second circle
    """
    fig, ax = plt.subplots()

    # Define the circles
    circle1_patch = plt.Circle((circle1[0], circle1[1]), circle1[2], facecolor='none', edgecolor='black')
    circle2_patch = plt.Circle((circle2[0], circle2[1]), circle2[2], facecolor='none', edgecolor='black')

    ax.add_patch(circle1_patch)
    ax.add_patch(circle2_patch)

    # Set the plot limits
    ax.set_xlim(min(circle1[0], circle2[0]) - max(circle1[2], circle2[2]) - 1,
                max(circle1[0], circle2[0]) + max(circle1[2], circle2[2]) + 1)
    ax.set_ylim(min([min([circle1[1], circle2[1]]) - max(circle1[2], circle2[2]) - 1, -height]),
                max([max(circle1[1], circle2[1]) + max(circle1[2], circle2[2]) + 1, height]))

    # Plot the surface points
    xs, ys, pxs, pys = get_vor_points(circle1[:2], circle1[2], circle2[:2], circle2[2], height=height)
    plt.plot(xs, ys, color='k')
    # plt.plot(pxs, pys)


    ax.set_axis_off()
    ax.set_aspect('equal')
    # plt.title('Voronoi Partition between Two Circles')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

# We want circle radius ratios 1:1, 1:3, 1:5, 1:10 and distances -r1, 0, +r1, +2r1


# Plot 1
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [1, 0, 1]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 2
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [2, 0, 1]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 3
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [3, 0, 1]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 4
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [4, 0, 1]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 5
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [6, 0, 0.33]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 6
circle1 = [0, 0, 2]  # (x, y, radius) of the first circle
circle2 = [8, 0, 6]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 7
circle1 = [0, 0, 2]  # (x, y, radius) of the first circle
circle2 = [10, 0, 6]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 8
circle1 = [0, 0, 2]  # (x, y, radius) of the first circle
circle2 = [12, 0, 6]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 9
circle1 = [0, 0, 1.5]  # (x, y, radius) of the first circle
circle2 = [7.5, 0, 7.5]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 10
circle1 = [0, 0, 1.5]  # (x, y, radius) of the first circle
circle2 = [9, 0, 7.5]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 11
circle1 = [0, 0, 1,5]  # (x, y, radius) of the first circle
circle2 = [10.5, 0, 7.5]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 12
circle1 = [0, 0, 1.5]  # (x, y, radius) of the first circle
circle2 = [12, 0, 7.5]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 13
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [10, 0, 10]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 14
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [11, 0, 10]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 15
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [12, 0, 10]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)

# Plot 16
circle1 = [0, 0, 1]  # (x, y, radius) of the first circle
circle2 = [13, 0, 10]  # (x, y, radius) of the second circle

plot_voronoi_partition(circle1, circle2, height=10)