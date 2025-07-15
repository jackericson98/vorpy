"""This script calculates the volume difference between the additively weighted and power surfaces as a function of radial distance"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.calculations import calc_dist, calc_aw_center, calc_pw_center, calc_tetra_vol, calc_surf_func
from vorpy.src.network import calc_surf_point, build_surf



def get_pow_vol(h, rad):
    return (1/3) * np.pi * rad**2 * h


def generate_circle_points(radius, resolution, x_center=0.0):
    # Compute the circumference
    circumference = 2 * np.pi * radius

    # Number of points based on resolution
    num_points = int(np.ceil(circumference / resolution))

    # Generate angles evenly spaced around the circle
    angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)

    # Circle lies in y-z plane, centered at (x_center, 0, 0)
    x = np.full(num_points, x_center)
    y = radius * np.cos(angles)
    z = radius * np.sin(angles)

    # Stack into Nx3 array of points
    points = np.column_stack((x, y, z))

    return points


def get_aw_vol(ball_locs, ball_rads, rad, func, pow_loc, aw_loc, res=0.05, ref_point=np.array([0, 0, 0])):
    # Generate the circle points
    circle_points = generate_circle_points(rad, res, pow_loc[0])

    # Project the points onto the surface
    proj_points = [calc_surf_point(ball_locs, circle_points[i], func) for i in range(len(circle_points))]
    
    # Build the surface
    surf = build_surf(ball_locs, ball_rads, None, res, 'aw', sfunc=func, perimeter=proj_points, surf_loc=aw_loc, surf_norm=np.array([1, 0, 0]))

    # Get the points and tris
    points, tris = surf[:2]

    # Loop through the tris to calculate the volume of the surface
    vol = 0
    for tri in tris:
        vol += calc_tetra_vol(ref_point, points[tri[0]], points[tri[1]], points[tri[2]])
    
    # Return the volume
    return vol


def compare_volumes(ball_locs, ball_rads, rad, func, pow_loc, aw_loc, res=0.05, ref_point=np.array([0, 0, 0])):
    # Get the volumes
    aw_vol = get_aw_vol(ball_locs, ball_rads, rad, func, pow_loc, aw_loc, res, ref_point)
    pw_vol = get_pow_vol(calc_dist(ref_point, pow_loc), rad)

    # Return the volume difference
    return 100 * (pw_vol - aw_vol) / aw_vol


def get_plotting_values(ball_locs, ball_rads, surf_res=0.05, res=0.05, max_rad=10, collect_intersect=True, ref_point=np.array([0, 0, 0])):
    # Get the center of the power surface
    pow_loc = calc_pw_center(*ball_rads, *ball_locs)[1]
    aw_loc = calc_aw_center(*ball_rads, *ball_locs)[1]

    # Calculate the surface function
    func = calc_surf_func(ball_locs[0], ball_rads[0], ball_locs[1], ball_rads[1])
    
    # Get the volume differences
    diffs = []
    rads = []

    # Get the volume difference for the first few points
    # for rad in np.arange(res/10, res, res/10):
    #     print(f"\rCalculating volume difference for {round(rad, 4)}/{max_rad}", end="")
    #     # Get the volume difference
    #     diffs.append(compare_volumes(ball_locs, ball_rads, rad, func, pow_loc, aw_loc, surf_res / 100))
    #     rads.append(rad)

    # Loop through the rads testing the volume difference as we go. 
    for i, rad in enumerate(np.arange(res, max_rad, res)):
        if i < 10:
            my_surf_res = surf_res / (10 - i)
        else:
            my_surf_res = surf_res
        print(f"\rCalculating volume difference for {round(rad, 4)}/{max_rad}", end="")
        # Get the volume difference
        diffs.append(compare_volumes(ball_locs, ball_rads, rad, func, pow_loc, aw_loc, my_surf_res, ref_point))
        rads.append(rad)

    # Get the intersection radius
    intersect_rad = np.sqrt(ball_rads[0]**2 - pow_loc[0]**2)

    if collect_intersect and intersect_rad > 0:
        # Get the volume difference at the intersection radius
        intersect_diff = compare_volumes(ball_locs, ball_rads, intersect_rad, func, pow_loc, aw_loc, surf_res, ref_point)
    else:
        intersect_diff = [0]

    # Return the volume difference
    return rads, diffs, (intersect_rad, intersect_diff)


def plot_volumes(balls_set, labels, colors, title, max_radial_dist=6, surf_res=0.05, radial_res=0.01, ref_point=np.array([0, 0, 0])):
    # Loop through the balls plotting each set
    for i, balls in enumerate(balls_set):
        print(f"\nPlotting {labels[i]}")
        rads, diffs, (intersect_rad, intersect_diff) = get_plotting_values([ball[0] for ball in balls], [ball[1] for ball in balls], surf_res, radial_res, max_radial_dist, ref_point=ref_point)
        
        # Extract values and create smoothed data for small radii
        values = [_ for _ in diffs]
        smoothed_values = values.copy()
        
        # Apply running window average for small radii (0 to 0.2)
        window_size = 5
        for j in range(len(rads)):
            if rads[j] <= 0.2:
                start_idx = max(0, j - window_size//2)
                end_idx = min(len(values), j + window_size//2 + 1)
                smoothed_values[j] = np.mean(values[start_idx:end_idx])
        my_x, my_y = [], []
        # Plot the smoothed data
        for k in range(len(rads)):
            # print(f"rads[k]: {rads[k]}, smoothed_values[k]: {smoothed_values[k]}")
            if smoothed_values[k] > 1000:
                continue
            else:
                my_x.append(rads[k])
                my_y.append(smoothed_values[k])
        plt.plot(my_x, my_y, color=colors[i], linewidth=3, label=labels[i])
        plt.scatter([intersect_rad], [intersect_diff], marker='o', color=colors[i], label=f'{labels[i]} Intersection', s=40)
    plt.legend(fontsize=10)
    plt.title(title, fontsize=30)
    plt.xlabel("+R", fontsize=25)
    plt.xlim(-0.5, 3.5)
    plt.xticks([0, 1.5, 3], fontsize=25)
    plt.ylim(-125, 625)
    plt.ylabel("Volume % Difference", fontsize=25)
    plt.yticks([0, 250, 500], fontsize=25)
    plt.tick_params(axis='both', length=10, width=2)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Define the ball pairs
    balls1 = [[(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([3.0, 0.0, 0.0]), 2.0)], [(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([2.5, 0.0, 0.0]), 2.0)], [(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([2.0, 0.0, 0.0]), 2.0)]]
    balls2 = [[(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([6.0, 0.0, 0.0]), 5.0)], [(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([5.5, 0.0, 0.0]), 5.0)], [(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([5.0, 0.0, 0.0]), 5.0)]]
    balls3 = [[(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([11.0, 0.0, 0.0]), 10.0)], [(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([10.5, 0.0, 0.0]), 10.0)], [(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([10.0, 0.0, 0.0]), 10.0)]]
    labels = ["0.0r Overlap", "0.5r Overlap", "1.0r Overlap"]
    colors = ['purple', 'green', 'orange']
    title1 = "Volume Difference"
    print("Plotting 1:2")
    plot_volumes(balls1, labels, colors, title1, ref_point=np.array([-2, 0, 0]), max_radial_dist=3)
