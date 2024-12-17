import time

import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist, calc_com, calc_angle_jit
from System.sys_funcs.calcs.edge import calc_circ, calc_edge_dir
from System.Network.edges.edge_project import edge_project
from System.sys_funcs.calcs.surf import calc_surf_func
from Visualize.mpl_visualize import plot_edges, plot_balls, plot_verts
import matplotlib.pyplot as plt


def edge_bad(e_points):
    """
    We need to find out if the edge points are out of order
    """
    for i in range(len(e_points) - 2):
        angle = calc_angle_jit(e_points[i + 1], e_points[i], e_points[i + 2])
        if angle < np.pi / 3 or angle > 5 * np.pi / 3:
            return False
    return True


def build_straight_edge(locs, rads, vlocs, res):
    # Get the location and radius of the circle inscribed between the edge atoms
    try:
        loc, rad = calc_circ(locs[0], locs[1], locs[2], rads[0], rads[1], rads[2])
    except TypeError:
        loc = calc_com([locs[0], locs[1], locs[2]])
        rad = calc_dist(loc, locs[0]) - rads[0]
    # Create the vals dictionary
    vals = {'loc': loc, 'rad': rad}
    # Determine the edge length
    edge_dist = calc_dist(vlocs[0], vlocs[1])
    # Divide the edge length by the resolution to find the number of points
    num_points = max(int(edge_dist / res) + 1, 3)
    # Create the new resolution to get even divisions of the edge
    new_res = edge_dist / num_points
    # Find the direction the edge heads
    edge_dir = vlocs[1] - vlocs[0]
    # Find the normalized vector between the vertices
    e_hat = edge_dir / np.linalg.norm(edge_dir)
    # Create the points
    e_points = [vlocs[0] + i * new_res * e_hat for i in range(num_points + 1)]
    # Return the edge
    return e_points, vals


def mid_edge_point(ep1, ep2, func, vmid, direction, new_direction=True):
    """
    Calculates the middle point between two edge points
    """
    # If the point is the first point we just need to move in the direction of the direction vector
    if new_direction:
        # Get the direction between the edges
        edir = ep2 - ep1
        # Get the distance between the points
        edist = calc_dist(ep1, ep2)
        # Get the ehat vector
        ehat = edir / edist
        # Get the middle point between the edge points
        proj_point = ep1 + 0.5 * edist * ehat
        # Get the normalized vector
        rn = proj_point - vmid
        # Normalize it
        direction = rn / np.linalg.norm(rn)
    # Project the point toward the projection point
    return edge_project(np.array(direction), np.array(vmid), np.array(func))


def build_edge(locs, rads, vlocs, res, blocs, brads, eballs, straight=False, vmid=None, dnorm=None, edub=False,
               edge_points1=None, edge_verts=None, redone_edge=False):
    """
    Build edge function. Takes in the locations and radii of the input balls and the vertices bounding the edge and
    outputs a fully resolved edge.
    """
    # If the edge is straight build the straight edge
    if straight or round(rads[0], 3) == round(rads[1], 3) == round(rads[2], 3):
        return build_straight_edge(locs, rads, vlocs, res)

    # Choose a curved surface to project onto. If the edge isn't straight at least 2 surfs are curved.
    if round(rads[0], 10) == round(rads[1], 10):
        func = calc_surf_func(locs[1], rads[1], locs[2], rads[2])
    else:
        func = calc_surf_func(locs[0], rads[0], locs[1], rads[1])

    # Get the edge direction
    edge_vals = None
    if vmid is None:
        edge_vals = calc_edge_dir(blocs, brads, eballs, vlocs, edub=edub)
        vmid, dnorm = edge_vals['vmid'], edge_vals['dnorm']

    if edge_points1 is not None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        plot_edges([edge_points1], fig, ax)
        plot_balls([blocs[_] for _ in eballs], [brads[_] for _ in eballs], fig=fig, ax=ax)
        ax.plot([edge_vals['vmid'][0], edge_vals['vmid'][0] + edge_vals['dnorm'][0]], [edge_vals['vmid'][1], edge_vals['vmid'][1] + edge_vals['dnorm'][1]], [edge_vals['vmid'][2], edge_vals['vmid'][2] + edge_vals['dnorm'][2]])
        plot_balls([edge_vals['loc']], [edge_vals['rad']], fig=fig, ax=ax, colors=['red'])
        plot_verts(vlocs, [1, 1], fig=fig, ax=ax, colors=['g', 'g'])
        print(edge_vals)
        print(vlocs)
        print(edge_verts)
        plt.show()

    # Check for the case 5
    if edge_vals is not None and edge_vals['case'] == 5:
        edge0 = build_edge(locs, rads, [vlocs[0], edge_vals['loc']], blocs, brads, eballs, res, straight,
                           edge_vals['vmid0'], edge_vals['dnorm0'])
        edge = build_edge(locs, rads, [edge_vals['loc'], edge_vals['loc2']], blocs, brads, eballs, res, straight,
                          edge_vals['vmid'], edge_vals['dnorm'])
        edge1 = build_edge(locs, rads, [edge_vals['loc2'], vlocs[1]], blocs, brads, eballs, res, straight,
                           edge_vals['vmid1'], edge_vals['dnorm1'])
        return edge0[0] + edge[0] + edge1[0], edge_vals
    # Create a catch for the time
    start = time.perf_counter()
    # Instantiate the edge points list with the vertices
    e_points = [*vlocs]
    # Main edge calculation loop
    while True:
        new_points_added = False  # Track whether new points are added in this iteration

        # Loop through the edge points
        i = 0
        while i < len(e_points) - 1:
            # Get the middle points
            ep1, ep2 = e_points[i], e_points[i + 1]

            # Check if the distance is greater than the resolution
            if calc_dist(ep1, ep2) > res:
                # Get the middle point
                mid_point = mid_edge_point(ep1, ep2, func, vmid, dnorm, new_direction=len(e_points) > 2)

                e_points.insert(i + 1, mid_point)  # Add the new midpoint
                new_points_added = True  # Mark that we added a new point

                i += 1  # Skip to the next segment
            i += 1

        # If no new points were added, the refinement is complete
        if not new_points_added:
            return e_points, edge_vals
        # Check the time
        if time.perf_counter() - start > 50:
            # If the edge is already being redone it is bad and there needs to be checked
            return build_edge(locs, rads, vlocs, res, blocs, brads, eballs, straight, vmid, -dnorm, edub,
                              redone_edge=True)


# Find projection values. Calculates the 181L end and projection points for the edge
def calc_edge_proj_pt(pv0, pv1, loc):
    # Get the projection point
    # Find the point in between the two vertex points
    r01 = pv1 - pv0  # Vector between vertices
    r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
    rn01 = r01 / r_mag  # Normal to the vector between the vertices
    pc01 = pv0 + 0.5 * rn01 * r_mag  # Center point

    # Determine if the theoretical center of the edge is inside the vertices or not
    dr = 1
    if np.sqrt(sum(np.square(loc - pv0))) < r_mag or np.sqrt(sum(np.square(loc - pv1))) < r_mag:
        dr = -1

    # Find the vector normal to the projection plane
    p_norm = dr * np.cross(loc - pc01, pv1 - pc01)
    # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
    r_pcr = - np.cross(p_norm, rn01)
    rn_pcr = r_pcr / np.linalg.norm(r_pcr)
    # Calculate the reference point
    return pc01 + 0.5 * r_mag * rn_pcr

# Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
def build_edge_old(locs, rads, vlocs, res, straight=None):
    # To ensure a better edge we cut the resolution in quarters
    res = res / 2
    # Check for straightness
    if straight is None:
        straight = False
        if rads[0] == rads[1] and rads[1] == rads[2]:
            straight = True
    # Get the location and radius of the circle inscribed between the edge atoms
    try:
        loc, rad = calc_circ(locs[0], locs[1], locs[2], rads[0], rads[1], rads[2])
    except TypeError:
        loc = calc_com([locs[0], locs[1], locs[2]])
        rad = calc_dist(loc, locs[0]) - rads[0]
    loc = np.array(loc)
    vals = {'loc': loc, 'rad': rad}
    # If the edge is straight return the bare minimum
    if straight:
        return vlocs, vals
    # Choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
    if round(rads[0], 10) == round(rads[1], 10):
        func = calc_surf_func(locs[1], rads[1], locs[2], rads[2])
    else:
        func = calc_surf_func(locs[0], rads[0], locs[1], rads[1])

    ################################################# Fill Edge ####################################################


    # Reset the edges points
    points = []
    # Typical case, no doublets
    pv0, pv1 = np.array(vlocs[0]), np.array(vlocs[1])
    # If the edge is completely straight add points in a line from pv0 to pv0 and return
    if straight or (rads[0] == rads[1] and rads[1] == rads[2]):
        # Get the vector between the two vectors and the number of point in the edge
        r = pv1 - pv0
        num_points = max(int(np.linalg.norm(r) / (4 * res)), 4)
        # Add the points
        for i in range(num_points + 1):
            points.append(pv0 + r * (i / num_points))
        return points, vals
    else:
        pa = calc_edge_proj_pt(pv0, pv1, loc)

    # Find the point in between the two vertex points
    r01 = pv1 - pv0  # Vector between vertices
    r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
    rn01 = r01 / r_mag  # Normal to the vector between the vertices
    # Find the number of points
    n = max(int(r_mag / res), 4)
    # Calculate the angle between the vertices and the reference point
    theta = calc_angle_jit(pa, pv0, pv1)
    # Add the first vertex to the list of points
    points = [pv0]
    # Find the edges points. Don't count the vertex
    for i in range(n + 1):
        if i == 0:
            A = 0.01 * theta / n
        elif i == 1:
            A = 0.99 * theta / n
        else:
            A = theta / n
        # Set pb to the previous point
        pb = points[-1]
        # Get the distance between pb and pa for c
        c = np.sqrt(sum(np.square(np.array(pb) - np.array(pa))))
        # Get the angle between pb, pa and pb + rno1
        B = calc_angle_jit(pb, pb + rn01, pa)
        # Get the last angle
        C = np.pi - B - A
        # Get the distance to our projection point or 'a' on our triangle
        a = np.sin(A) * c / np.sin(C)
        # Use that distance to project rn01 from pb to find our projection point or pc
        pc = pb + a * rn01
        # Get the vector from pa to pc
        r_ac = np.array(pc) - np.array(pa)
        r_nac = r_ac / np.linalg.norm(r_ac)
        # Project the vector onto the surface
        surf_point = edge_project(r_nac, pa, np.array(func))
        if surf_point is None:
            break
        points.append(surf_point)
    # Add the end point
    points.append(pv1)
    # Finally return the points
    return points, vals