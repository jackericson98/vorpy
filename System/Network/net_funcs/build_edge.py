import numpy as np
from System.sys_funcs.calcs.calcs import calc_circ, calc_angle, calc_surf_func
from System.Network.net_objs.surface import Surface


# Find projection values. Calculates the correct end and projection points for the edge
def calc_edge_proj_pt(pv0, pv1, loc):
    # Get the projection point
    # Find the point in between the two vertex points
    r01 = pv1 - pv0  # Vector between vertices
    r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
    rn01 = r01 / r_mag  # Normal to the vector between the vertices
    pc01 = pv0 + 0.5 * rn01 * r_mag  # Center point

    # Determine if the theoretical center of the edge is inside the vertices or not
    dr = 1
    if np.sqrt(sum(np.square(np.array(loc) - np.array(pv0)))) < r_mag or \
            np.sqrt(sum(np.square(np.array(loc) - np.array(pv1)))) < r_mag:
        dr = -1

    # Find the vector normal to the projection plane
    p_norm = dr * np.cross(np.array(loc) - np.array(pc01), np.array(pv1) - np.array(pc01))
    # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
    r_pcr = - np.cross(p_norm, rn01)
    rn_pcr = r_pcr / np.linalg.norm(r_pcr)
    # Calculate the reference point
    return pc01 + 0.5 * r_mag * rn_pcr


# Project method. Projects a point onto the surface using a reference point
def edge_project(rn, pa, func, ep_1, ep_2=None):
    # Get the function values
    f = func
    # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
    a = f[0] * rn[0] ** 2 + f[1] * rn[1] ** 2 + f[2] * rn[2] ** 2 + f[3] * rn[0] * rn[1] + f[4] * rn[
        1] * rn[2] + f[5] * rn[2] * rn[0]
    b = 2 * f[0] * rn[0] * pa[0] + 2 * f[1] * rn[1] * pa[1] + 2 * f[2] * rn[2] * pa[2] + f[3] \
        * (rn[0] * pa[1] + rn[1] * pa[0]) + f[4] * (rn[1] * pa[2] + rn[2] * pa[1]) + f[5] \
        * (rn[2] * pa[0] + rn[0] * pa[2]) + f[6] * rn[0] + f[7] * rn[1] + f[8] * rn[2]
    c = f[0] * pa[0] ** 2 + f[1] * pa[1] ** 2 + f[2] * pa[2] ** 2 + f[3] * pa[0] * pa[1] + f[4] * pa[1] * pa[
        2] + f[5] * pa[2] * pa[0] + f[6] * pa[0] + f[7] * pa[1] + f[8] * pa[2] + f[9]
    # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
    # and add that point to our surface list of points
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        # Calculate the roots
        roots = np.roots([a, b, c])
        # If one root exists return it
        if len(roots) == 1:
            return pa + roots[0] * rn
        else:
            p1 = pa + min(roots) * rn
            p2 = pa + max(roots) * rn
        # If the point we are calculating is the first in the edge choose the one closest to the vertex
        if ep_2 is None:
            point = p1
            if np.sqrt(sum(np.square(np.array(p2) - np.array(ep_1)))) <= \
                    np.sqrt(sum(np.square(np.array(p1) - np.array(ep_1)))):
                point = p2
        # If we have 2 points to choose from, choose the one that makes the angle closer to 180
        else:
            point = p1
            if calc_angle(ep_1, ep_2, p2) >= calc_angle(ep_1, ep_2, p1):
                point = p2
        # Return the point we choose
        return point


# Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
def build_edge(edge, surf=None, res=None, straight=None):
    # Set the self defining straight value
    if straight is None:
        edge.straight = False
    # Get the location and radius of the circle inscribed between the edge atoms
    edge.loc, edge.rad = calc_circ(*[_.loc for _ in edge.atoms], *[_.rad for _ in edge.atoms])
    # Get the pvals
    # Typical case, no doublets
    edge.pv0, edge.pv1 = np.array(edge.verts[0].loc), np.array(edge.verts[1].loc)
    if not edge.straight:
        edge.pa = calc_edge_proj_pt(edge.pv0, edge.pv1, edge.loc)
    # Reset the edges points
    edge.points = []
    # Check to see if a minimum distance has been provided
    if res is None:
        # Get the network's minimum distance
        res = edge.net.surf_res
    # Check to see if a surface has been provided
    if surf is None:
        # Choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
        if round(edge.atoms[0].rad, 10) == round(edge.atoms[1].rad, 10):
            surf = Surface(edge.atoms[1:], edge.net)
        else:
            surf = Surface(edge.atoms[:2], edge.net)
    # Check to see if the surface has its function values
    if surf.func is None:
        surf.func = calc_surf_func(surf.atoms[0].loc, surf.atoms[0].rad, surf.atoms[1].loc, surf.atoms[1].rad)

    ################################################# Fill Edge ####################################################

    # If the edge is completely straight add points in a line from pv0 to pv0 and return
    if straight or (edge.atoms[0].rad == edge.atoms[1].rad and edge.atoms[1].rad == edge.atoms[2].rad):
        # Get the vector between the two vectors and the number of point in the edge
        r = edge.pv1 - edge.pv0
        num_points = max(int(np.linalg.norm(r) / edge.net.surf_res), 4)
        # Add the points
        for i in range(num_points + 1):
            edge.points.append(edge.pv0 + r * (i / num_points))
        return

    # Calculate the points

    # Find the point in between the two vertex points
    r01 = edge.pv1 - edge.pv0  # Vector between vertices
    r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
    rn01 = r01 / r_mag  # Normal to the vector between the vertices
    # Find the number of points
    n = max(int(r_mag / res), 4)
    # Calculate the angle between the vertices and the reference point
    theta = calc_angle(edge.pa, edge.pv0, edge.pv1)
    # Add the first vertex to the list of points
    edge.points = [edge.pv0.tolist()]
    # Find the edges points. Don't count the vertex
    for i in range(n + 1):
        if i == 0:
            A = 0.01 * theta / n
        elif i == 1:
            A = 0.99 * theta / n
        else:
            A = theta / n
        # Set pb to the previous point
        pb = edge.points[-1]
        # Get the distance between pb and pa for c
        c = np.sqrt(sum(np.square(np.array(pb) - np.array(edge.pa))))
        # Get the angle between pb, pa and pb + rno1
        B = calc_angle(pb, pb + rn01, edge.pa)
        # Get the last angle
        C = np.pi - B - A
        # Get the distance to our projection point or 'a' on our triangle
        a = np.sin(A) * c / np.sin(C)
        # Use that distance to project rn01 from pb to find our projection point or pc
        pc = pb + a * rn01
        # Get the vector from pa to pc
        r_ac = np.array(pc) - np.array(edge.pa)
        r_nac = r_ac / np.linalg.norm(r_ac)
        # Project the vector onto the surface
        surf_point = edge_project(r_nac, edge.pa, surf.func, edge.points[-1],
                                  edge.points[-2] if len(edge.points) > 1 else None)
        if surf_point is None:
            break
        edge.points.append(surf_point)
    # Add the end point
    edge.points.append(edge.pv1)
