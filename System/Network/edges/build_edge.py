import numpy as np
from System.Network.edges.edge_project import edge_project
from System.sys_funcs.calcs.calcs import calc_angle_jit, calc_com, calc_dist
from System.sys_funcs.calcs.surf import calc_surf_func
from System.sys_funcs.calcs.edge import calc_circ, calc_edge_dir


def check_edge_point(point, locs, rads):
    # Check 1: If a point is inside only one sphere.
    overlap_count = 0
    for i in range(3):
        if calc_dist(point, locs[i]) < rads[i]:
            overlap_count += 1
    if overlap_count in {1, 2}:
        return False
    return True


# Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
def build_edge(locs, rads, vlocs, res, blocs, brads, eballs, straight=None):
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
        edge_dist = calc_dist(vlocs[0], vlocs[1])
        num_points = max(int(edge_dist / res) + 1, 3)
        new_res = edge_dist / num_points
        edge_dir = vlocs[1] - vlocs[0]
        e_hat = edge_dir / np.linalg.norm(edge_dir)
        e_points = [vlocs[0]]
        last_point = vlocs[0]
        for i in range(num_points):
            last_point = last_point + e_hat * new_res
            e_points.append(last_point)
        return e_points, vals
    # Choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
    if round(rads[0], 10) == round(rads[1], 10):
        func = calc_surf_func(locs[1], rads[1], locs[2], rads[2])
    else:
        func = calc_surf_func(locs[0], rads[0], locs[1], rads[1])

    ################################################# Fill Edge ####################################################

    # Typical case, no doublets
    pv0, pv1 = np.array(vlocs[0]), np.array(vlocs[1])
    # If the edge is completely straight add points in a line from pv0 to pv0 and return
    edge_vals = calc_edge_dir(blocs, brads, eballs, vlocs)
    pa = edge_vals['vmid'] - 0.5 * edge_vals['vdist'] * edge_vals['dnorm']

    # Pull the vnorm value so we dont have to keep referencing the dictionary
    rn01 = edge_vals['vnorm']
    # Find the number of points
    n = max(int(edge_vals['vdist'] / res), 2)
    # Divide the vector between vertices
    increment = edge_vals['vdist'] / n

    # Add the first vertex to the list of points
    points = [pv0]
    # Find the edges points. Don't count the vertex
    for i in range(n):
        # Get the point along the vector between the vertices
        pc = pv0 + i * increment * rn01

        # Get the vector from pa to pc
        r_ac = np.array(pc) - np.array(pa)
        r_nac = r_ac / np.linalg.norm(r_ac)
        # Project the vector onto the surface
        surf_point = edge_project(r_nac, pa, np.array(func), locs, rads, points[-1], points[-2] if len(points) > 1 else None)
        if surf_point is None:
            break
        points.append(surf_point)
    # Add the end point
    points.append(pv1)
    # Finally return the points
    return points, vals


