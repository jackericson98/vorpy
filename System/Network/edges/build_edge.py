import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist, calc_com
from System.sys_funcs.calcs.edge import calc_circ, calc_edge_dir
from System.Network.edges.edge_project import edge_project
from System.sys_funcs.calcs.surf import calc_surf_func


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


def build_edge(locs, rads, vlocs, res, blocs, brads, eballs, straight=False, vmid=None, dnorm=None):
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
        edge_vals = calc_edge_dir(blocs, brads, eballs, vlocs)
        vmid, dnorm = edge_vals['vmid'], edge_vals['dnorm']

    # Check for the case 5
    if edge_vals is not None and edge_vals['case'] == 5:
        edge0 = build_edge(locs, rads, [vlocs[0], edge_vals['loc']], blocs, brads, eballs, res, straight,
                           edge_vals['vmid0'], edge_vals['dnorm0'])
        edge = build_edge(locs, rads, [edge_vals['loc'], edge_vals['loc2']], blocs, brads, eballs, res, straight,
                          edge_vals['vmid'], edge_vals['dnorm'])
        edge1 = build_edge(locs, rads, [edge_vals['loc2'], vlocs[1]], blocs, brads, eballs, res, straight,
                           edge_vals['vmid1'], edge_vals['dnorm1'])
        return edge0[0] + edge[0] + edge1[0], edge_vals

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
