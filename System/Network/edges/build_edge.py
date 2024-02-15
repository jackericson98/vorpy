import numpy as np
from System.Network.edges.edge_project import edge_project, calc_edge_proj_pt
from System.sys_funcs.calcs.calcs import calc_angle_jit, calc_com, calc_dist
from System.sys_funcs.calcs.surf import calc_surf_func
from System.sys_funcs.calcs.circle import calc_circ


# Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
def build_edge(alocs, arads, vlocs, res, straight=None):
    # Check for straightness
    if straight is None:
        straight = False
        if arads[0] == arads[1] and arads[1] == arads[2]:
            straight = True
    # Choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
    if round(arads[0], 10) == round(arads[1], 10):
        func = calc_surf_func(alocs[1], arads[1], alocs[2], arads[2])
    else:
        func = calc_surf_func(alocs[0], arads[0], alocs[1], arads[1])

    ################################################# Fill Edge ####################################################

    # Get the location and radius of the circle inscribed between the edge atoms
    try:
        loc, rad = calc_circ(alocs[0], alocs[1], alocs[2], arads[0], arads[1], arads[2])
    except TypeError:
        loc = calc_com([alocs[0], alocs[1], alocs[2]])
        rad = calc_dist(loc, alocs[0]) - arads[0]
    loc = np.array(loc)
    vals = {'loc': loc, 'rad': rad}
    # Reset the edges points
    points = []
    # Typical case, no doublets
    pv0, pv1 = np.array(vlocs[0]), np.array(vlocs[1])
    # If the edge is completely straight add points in a line from pv0 to pv0 and return
    if straight or (arads[0] == arads[1] and arads[1] == arads[2]):
        # Get the vector between the two vectors and the number of point in the edge
        r = pv1 - pv0
        num_points = max(int(np.linalg.norm(r) / res), 4)
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
    return points, vals
