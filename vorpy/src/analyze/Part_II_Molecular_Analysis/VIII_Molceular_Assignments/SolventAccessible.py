import matplotlib.pyplot as plt
import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist, calc_com


# Function to draw circles
def draw_circle(ax, center, radius, color, alpha=1, fill=False, label=None, linewidth=1):
    circle = plt.Circle(center, radius, color=color, fill=fill, alpha=alpha, label=label, linewidth=linewidth)
    ax.add_artist(circle)


def get_sol_acc_points(mol_locs, mol_rads, sol_rad, res, d_theta):
    # Calculate the center of mass of the molecule
    mol_com = calc_com([_ for _ in mol_locs])

    # Set up the molecule double list
    mol_doubs = [(mol_locs[i], mol_rads[i]) for i in range(len(mol_locs))]

    # Set up the theta variable
    theta = 0

    # Define the array of sol ball locs
    sol_locs, touch_points, outer_points = [], [], []

    # Go through each of the angles
    while theta <= 2 * np.pi + d_theta:

        # Set the move distance to the sol_rad
        move_dist = sol_rad

        # Add the d_theta to theta
        theta += d_theta

        # Get the normalized ray
        norm_ray = np.array([np.cos(theta), np.sin(theta)])

        # Set up the circle distance variabble
        circ_dist = np.inf

        # Keep finding the location of the point until its close enough
        while abs(circ_dist) > res:
            # Add the scaled ray to the mol com and calcuklate if it overlaps
            scaled_ray = move_dist * norm_ray + mol_com
            # fig, ax = plt.subplots()
            # Create a list for all of the distances to see if we make it less than the threshold
            circ_dists = []

            # Go through each of the circles and chack for an overlap
            for circ in mol_doubs:
                # Calculate the distance
                circ_dists.append(calc_dist(circ[0], scaled_ray) - sol_rad - circ[1])

            # Find the closest circle
            circ_dist = min(circ_dists)

            move_dist = abs(move_dist - 0.5 * circ_dist)

        # Record the center of the sol circle
        sol_locs.append(scaled_ray)
        # get the closest circle
        close_circ = mol_doubs[circ_dists.index(circ_dist)]
        # Get the normalized ray between the two
        sol_ray = close_circ[0] - scaled_ray
        sol_ray_norm = sol_ray / np.linalg.norm(sol_ray)
        touch_points.append(scaled_ray + sol_rad * sol_ray_norm)
        outer_points.append(scaled_ray - sol_rad * sol_ray_norm)

    return sol_locs[1:], touch_points[1:], outer_points[1:]


if __name__ == '__main__':

    fig, ax = plt.subplots()

    # List of molecule radii
    mol_locs = [np.array([2, 0.5]), np.array([0, 1.5]), np.array([-1.4, -0.2]), np.array([-2.5, -1.5]),
                np.array([-2.4, -2.9]), np.array([-3.5, 0])]
    mol_rads = [2, 2.1, 2.3, 1.8, 1.9, 1.4]

    # Solvent accessible radii
    sol_rad = 1.4

    # Resolution
    res = 0.001

    # Create the d theta value
    d_theta = np.pi / 64

    sol_locs, touch_points, outer_points = get_sol_acc_points(mol_locs, mol_rads, sol_rad, res, d_theta)


    #
    # Convert lists to numpy arrays for easier handling
    solvent_accessible_points = np.array(sol_locs)
    vdw_points = np.array(touch_points)
    sol_ex_pts = np.array(outer_points)
    #
    ax.scatter(sol_ex_pts[:, 0], sol_ex_pts[:, 1], color='k', s=1)
    ax.fill_between(sol_ex_pts[:, 0], sol_ex_pts[:, 1], color='sandybrown', label='Probe Shell Volume')
    ax.fill_between(solvent_accessible_points[:, 0], solvent_accessible_points[:, 1], color='khaki', label='Solvent Accesible Volume')
    ax.fill_between(vdw_points[:, 0], vdw_points[:, 1], color='cornflowerblue', label='Molecular Volume')
    # Draw the molecule
    mol_label = 'van der Waals Volume'
    for center, radius in zip(mol_locs, mol_rads):
        draw_circle(ax, center, radius, 'mediumorchid', fill=True, alpha=1, label=mol_label)
        mol_label = None
    for center, radius in zip(mol_locs, mol_rads):
        draw_circle(ax, center, radius, 'k', fill=False, alpha=0.5)

    draw_circle(ax, solvent_accessible_points[-1], sol_rad, 'grey', alpha=0.8, fill=True, label='Probe')

    draw_circle(ax, solvent_accessible_points[-1], sol_rad, 'k', alpha=1)

    ax.scatter([solvent_accessible_points[-1][0]], [solvent_accessible_points[-1][1]], s=5, color='k')

    # Draw the solvent accessible and void volume points
    ax.scatter(solvent_accessible_points[:, 0], solvent_accessible_points[:, 1], s=0.5, color='k')

    ax.scatter(vdw_points[:, 0], vdw_points[:, 1], color='k', s=0.5)


    #
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 10])
    # Plot setup
    ax.set_aspect('equal', 'box')
    legend = ax.legend(loc='upper right', bbox_to_anchor=(1.25, 0.97), shadow=True, prop={'size': 20})
    ax.axis('off')
    plt.show()
