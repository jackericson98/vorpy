import numpy as np
from System.Network.edges.edge_project import edge_project, calc_edge_proj_pt
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
def build_edge(locs, rads, vlocs, res, straight=None):
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
        surf_point = edge_project(r_nac, pa, np.array(func), points[-1], points[-2] if len(points) > 1 else None)
        if surf_point is None:
            break
        points.append(surf_point)
    # Add the end point
    points.append(pv1)
    # Finally return the points
    return points, vals


# Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
def build_edge_better(locs, rads, eballs, vlocs, res):
    # To ensure a better edge we cut the resolution in quarters
    res = res / 2

    # Calculate the edge direction
    edge_vals = calc_edge_dir(locs, rads, eballs, vlocs)

    # Choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
    if round(rads[0], 10) == round(rads[1], 10):
        func = calc_surf_func(locs[1], rads[1], locs[2], rads[2])
    else:
        func = calc_surf_func(locs[0], rads[0], locs[1], rads[1])

    ################################################# Fill Edge ####################################################

    # Typical case, no doublets
    pv0, pv1 = np.array(vlocs[0]), np.array(vlocs[1])
    # if the edge is just a normal edge or if the verts are outside the two locs half the distance between the verts
    if edge_vals['loc2'] is None:
        pa = edge_vals['vmid'] - 0.5 * edge_vals['vdist'] * edge_vals['dnorm']
    # Otherwise we want half the distance
    else:
        ldist = calc_dist(edge_vals['loc'], edge_vals['loc2'])
        pa = edge_vals['vmid'] - 0.5 * edge_vals['dnorm'] * (ldist - calc_dist(edge_vals['vmid'], edge_vals['loc']))

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
        surf_point = edge_project(r_nac, pa, np.array(func), points[-1], points[-2] if len(points) > 1 else None)
        if surf_point is None:
            break
        points.append(surf_point)
    # Add the end point
    points.append(pv1)
    # Finally return the points
    return points, edge_vals
