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
from vorpy.src.visualize import plot_surfs, plot_edges, plot_balls



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

    further_rad = 20
    # Generate the circle points
    circle_points = generate_circle_points(rad, res, pow_loc[0])

    # Project the points onto the surface
    aw_proj_points = [calc_surf_point(ball_locs, circle_points[i], func) for i in range(len(circle_points))]
    pow_proj_points = [[pow_loc[0], pnt[1], pnt[2]] for pnt in aw_proj_points]

    # Get the further circle points
    further_circle_points = generate_circle_points(further_rad, 4 *res, pow_loc[0])
    further_aw_proj_points = [calc_surf_point(ball_locs, further_circle_points[i], func) for i in range(len(further_circle_points))]
    further_pow_proj_points = [[pow_loc[0], pnt[1], pnt[2]] for pnt in further_aw_proj_points]
    
    # Build the surface
    aw_surf = build_surf(ball_locs, ball_rads, None, res, 'aw', sfunc=func, perimeter=aw_proj_points, surf_loc=aw_loc, surf_norm=np.array([1, 0, 0])) 
    pow_surf = build_surf(ball_locs, ball_rads, None, res, 'pow', sfunc=func, perimeter=pow_proj_points, surf_loc=pow_loc, surf_norm=np.array([1, 0, 0]))

    # Get the further projection surfaces
    further_aw_surf = build_surf(ball_locs, ball_rads, None, res, 'aw', sfunc=func, perimeter=further_aw_proj_points, surf_loc=aw_loc, surf_norm=np.array([1, 0, 0])) 
    further_pow_surf = build_surf(ball_locs, ball_rads, None, res, 'pow', sfunc=func, perimeter=further_pow_proj_points, surf_loc=pow_loc, surf_norm=np.array([1, 0, 0]))

    # Get the points and tris
    aw_points, aw_tris = aw_surf[:2]
    pow_points, pow_tris = pow_surf[:2]
    further_aw_points, further_aw_tris = further_aw_surf[:2]
    further_pow_points, further_pow_tris = further_pow_surf[:2]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plot_surfs(spnts=[aw_points], stris=[aw_tris], ax=ax, fig=fig, simp_color='red', alpha=0.2, colors=['red' for _ in range(len(aw_tris))])
    plot_surfs(spnts=[pow_points], stris=[pow_tris], ax=ax, fig=fig, simp_color='blue', alpha=0.2)
    plot_surfs(spnts=[further_aw_points], stris=[further_aw_tris], ax=ax, fig=fig, simp_color='red', alpha=0.01, colors=['red' for _ in range(len(further_aw_tris))])
    plot_surfs(spnts=[further_pow_points], stris=[further_pow_tris], ax=ax, fig=fig, simp_color='blue', alpha=0.01)

    plot_balls(alocs=ball_locs, arads=ball_rads, ax=ax, fig=fig, alpha=0.05, colors=['grey', 'grey'])

    # ax.scatter([pow_loc[0]], [pow_loc[1]], [pow_loc[2]], color='blue', s=100)
    # ax.scatter([aw_loc[0]], [aw_loc[1]], [aw_loc[2]], color='red', s=100)
    print(rad)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)

    ax.text(-3, 0.5, 0, f"Test Point", fontsize=20)
    ax.scatter([-2], [0], [0], color='black', s=30)
    ax.scatter([ball_locs[0][0]], [ball_locs[0][1]], [ball_locs[0][2]], color='black', s=30, marker='x')
    ax.text(ball_locs[0][0] - 0.4, ball_locs[0][1] + 0.2, ball_locs[0][2], f"S1", fontsize=20)
    ax.scatter([ball_locs[1][0]], [ball_locs[1][1]], [ball_locs[1][2]], color='black', s=30, marker='x')
    ax.text(ball_locs[1][0] + 0.2, ball_locs[1][1] + 0.2, ball_locs[1][2], f"S2", fontsize=20)


    for point in aw_proj_points[::5]:
        ax.plot([-2, point[0]], [0, point[1]], [0, point[2]], color='red', linewidth=0.5, linestyle='--', alpha=0.5)
    
    for point in pow_proj_points[3::5]:
        ax.plot([-2, point[0]], [0, point[1]], [0, point[2]], color='blue', linewidth=0.5, linestyle='--', alpha=0.5)

    ax.plot([1.4, 1.6], [aw_proj_points[int(len(aw_proj_points)/4)][1], aw_proj_points[int(len(aw_proj_points)/4)][1]], [aw_proj_points[int(len(aw_proj_points)/4)][2], aw_proj_points[int(len(aw_proj_points)/4)][2]], color='black', linewidth=1)
    ax.plot([1.5, 1.5], [0, aw_proj_points[int(len(aw_proj_points)/4)][1]], [0, aw_proj_points[int(len(aw_proj_points)/4)][2]], color='black', linewidth=1)
    ax.plot([1.4, 1.6], [0, 0], [0, 0], color='black', linewidth=1)
    ax.text(1.55, aw_proj_points[int(len(aw_proj_points)/4)][1] / 2, aw_proj_points[int(len(aw_proj_points)/4)][2] / 2, f"+D", fontsize=20)
    ax.plot([-10, 15], [0, 0], [0, 0], color='black', linewidth=0.5, alpha=0.5)
    plt.show()



    # Loop through the tris to calculate the volume of the surface
    vol = 0
    for tri in aw_tris:
        vol += calc_tetra_vol(ref_point, aw_points[tri[0]], aw_points[tri[1]], aw_points[tri[2]])
    
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
    plt.xticks([0, 2.5, 5], fontsize=25)
    plt.ylim(-125, 1100)
    plt.ylabel("Volume % Difference", fontsize=25)
    plt.yticks([0, 250, 500, 750, 1000], fontsize=25)
    plt.tick_params(axis='both', length=10, width=2)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Define the ball pairs
    balls1 = [[(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([2.5, 0.0, 0.0]), 2.0)]]
    balls2 = [[(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([5.0, 0.0, 0.0]), 5.0)]]
    balls3 = [[(np.array([0.0, 0.0, 0.0]), 1.0), (np.array([10.0, 0.0, 0.0]), 10.0)]]
    labels = ["0.0r Overlap", "0.5r Overlap", "1.0r Overlap"]
    colors = ['purple', 'green', 'orange']
    title1 = "Volume Difference"
    print("Plotting 1:10")
    plot_volumes(balls1, labels, colors, title1, ref_point=np.array([-2, 0, 0]), max_radial_dist=6, radial_res=1, surf_res=1)
