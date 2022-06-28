import numpy as np


########################################################################################################################
"""Calculator functions"""


# Calculate distance function. Finds the distance between 2 points
def calc_dist(l1, l2):
    d = np.sqrt((l1[0]-l2[0])**2+(l1[1]-l2[2])**2+(l1[2]-l2[2])**2)
    return d


# Calculate center of mass function. Takes in a set of points and returns the coordinates of the com
def calc_com(points):
    # Se the running total for the x, y, z values to 0
    xtot, ytot, ztot = 0, 0, 0
    for point in points:
        xtot += point[0]
        ytot += point[1]
        ztot += point[2]
    return [xtot/len(points), ytot/len(points), ztot/len(points)]


# Calculate angle function. Finds the angle (in rads) between three points. The first being the common point
def calc_angle(p0, p1, p2=None):
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = np.array(p0), np.array(p1)
    else:
        v0, v1 = np.array(p1) - np.array(p0), np.array(p2) - np.array(p0)

    # Get the unit vectors
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)

    return np.arccos(np.clip(np.dot(n0, n1), -1.0, 1.0))


# Bisector function. Creates a bisector surface between 2 atoms
def calc_surf(a1, a2):

    # Set a1 to the smaller of the two atoms
    if a1.rad > a2.rad:
        a2, a1 = a1, a2

    # Grab the centers of the spheres
    x1, y1, z1 = a1.loc
    x2, y2, z2 = a2.loc

    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = a1.rad - a2.rad
    K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
    d = x1 - x2, y1 - y2, z1 - z2
    J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2

    # Instantiate/reset the hyperboloid coefficient vector lists
    ABC, DEF, GHI = [], [], []
    # Calculate hyperboloid coefficients
    for i in range(3):
        ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
        DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
        GHI.append(-8 * R ** 2 * a1.loc[i] - 4 * K * d[i])

    return ABC + DEF + GHI + [J] + [K] + [d]


# Calculate edge points function. Takes in an edge and a surface and updates the edge's points
def calc_edge_points(edge, surf):
    # Make sure that a0 is the atom with the smaller radius
    a0 = surf.atoms[0]
    if a0.rad > surf.atoms[1].rad:
        a0 = surf.atoms[1]
    # Find the angle made between the edges verts and the atom
    max_ang = calc_angle(a0.loc, edge.verts[0].loc, edge.verts[1].loc)
    num_points = int(np.degrees(max_ang))
    dtheta = max_ang/num_points
    # Make the first point in the edge its first vertex
    edge.points.append(edge.verts[0].loc)
    # Find the points along the edge, incrementing by angle
    for i in range(num_points):
        # Find the unit vector facing the COM from the previous point in the path
        r0 = np.array(edge.points[-1]) - np.array(edge.verts[1].loc)
        rn = r0 / np.linalg.norm(r0)
        # Get the angle between the path direction and center of the atom
        r_theta = calc_angle(edge.points[-1], rn, a0.loc)
        pn_1_theta = 180 - r_theta - dtheta
        # Using the law of sines, calculate the new sample point location
        b = np.sin(dtheta) * np.array(edge.points[-1]) / np.sin(pn_1_theta)
        new_samp = b[0] + edge.points[-1][0], b[1] + edge.points[-1][1], b[2] + edge.points[-1][2]
        # Get the new point location
        new_point = calc_spnt(surf, new_samp)
        edge.points.append(new_point)


# Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smaller of the surfaces 2 atoms through the point into the surface
def calc_spnt(surf, point):
    # Make sure that a0 is the atom with the smaller radius
    a0 = surf.atoms[0]
    if a0.rad > surf.atoms[1].rad:
        a0 = surf.atoms[1]

    # Grab the function
    if surf.func is None:
        surf.func = calc_surf(surf.atoms[0], surf.atoms[1])
    f = surf.func

    # Set up the unit vector
    vi = np.array(point) - np.array(a0.loc)
    vn = vi/np.linalg.norm(vi)
    # Find the location on the surface of the atom
    vi = np.array(a0.loc) + vn * a0.rad

    # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
    a = f[0] * vn[0] ** 2 + f[1] * vn[1] ** 2 + f[2] * vn[2] ** 2 + f[3] * vn[0] * vn[1] + f[4] * vn[1] * vn[2] + f[5] \
        * vn[2] * vn[0]
    b = 2 * f[0] * vn[0] * vi[0] + 2 * f[1] * vn[1] * vi[1] + 2 * f[2] * vn[2] * vi[2] + f[3] \
        * (vn[0] * vi[1] + vn[1] * vi[0]) + f[4] * (vn[1] * vi[2] + vn[2] * vi[1]) + f[5] \
        * (vn[2] * vi[0] + vn[0] * vi[2]) + f[6] * vn[0] + f[7] * vn[1] + f[8] * vn[2]
    c = f[0] * vi[0] ** 2 + f[1] * vi[1] ** 2 + f[2] * vi[2] ** 2 + f[3] * vi[0] * vi[1] + f[4] * vi[1] * vi[2] + \
        f[5] * vi[2] * vi[0] + f[6] * vi[0] + f[7] * vi[1] + f[8] * vi[2] + f[9]

    # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
    # and add that point to our surface list of points
    if b ** 2 - 4 * a * c > 0:
        mag = min(abs(np.roots([a, b, c])))
        return vi + mag*vn


# Edge trace function. Recursively goes around the edges of the surface and adds points for the vertices and the edges
def edge_trace(surf, en=None, vn=None):
    # First run check
    if en is None:
        en = surf.edges[0]
        vn = en.verts[0]
    # Check to see if en has been traced out yet
    if not en.points:
        calc_edge_points(en, surf)
    # Add the edge points and the vert points
    surf.vert_points.append(vn.loc)
    # Check to see if the points are ordered correctly
    if calc_dist(en.points[0], vn.loc) > calc_dist(en.points[-1], vn.loc):
        en.points.reverse()
    # Add the points to the edge points
    surf.edge_points += en.points
    # Find the next edge around the surface
    for em in surf.edges:
        # Check each of em's vertices against the lead vertex
        for i in range(2):
            if vn == em.verts[i] and em.verts[i].loc not in surf.vert_points:
                edge_trace(surf, em, em.verts[(i+1) % 2])
    return


# Make mesh function. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf, a0, density):
    # Check to see if the edges' points have been recorded yet
    if not surf.edge_points:
        edge_trace(surf)
    # Calculate the center of mass of the edge points
    com = calc_com(surf.edge_points)
    # Find where com maps on the surface
    com = calc_spnt(surf, com)
    # Make a list of paths
    paths = []
    max_angs = []
    # Create (and calculate the maximum angle for) each path
    for i in range(len(surf.edge_points)):
        paths.append([surf.edge_points[i]])
        max_angs.append(calc_angle(a0.loc, com, surf.edge_points[i]))
    # Get the maximum path information
    max_path_ndx = max_angs.index(max(max_angs))
    max_path = paths[max_path_ndx]
    # Decide how many rings based off of density
    num_rings = int(calc_dist(max_path[0], com) * density)
    # Create a list of dthetas
    dthetas = []
    for i in range(len(surf.edge_points)):
        dthetas.append(max_angs[i]/num_rings)

    # Build the surface ring by ring
    for i in range(num_rings):
        # Go through each path
        for j in range(len(paths)):
            # If the path has terminated continue to the next one
            if paths[j][-1] is None:
                continue
            # Find the unit vector facing the COM from the previous point in the path
            r0 = np.array(com) - np.array(paths[j][-1])
            rn = r0/np.linalg.norm(r0)
            # Get the angle between the path direction and center of the atom
            r_theta = calc_angle(paths[j][-1], rn, a0.loc)
            pn_1_theta = 180 - r_theta - dthetas[j]
            # Using the law of sines, calculate the new sample point location
            b = np.sin(dthetas[j]) * np.array(paths[j][-1]) / np.sin(pn_1_theta)
            new_samp = b[0] + paths[j][-1][0], b[1] + paths[j][-1][1], b[2] + paths[j][-1][2]
            # Get the new point location
            new_point = calc_spnt(surf, new_samp)

            # Find the most recent non-None point
            k = 1
            while paths[j-k][-1] is None:
                k += 1
            # Check to see if it is too close to the last point in the path before
            if calc_dist(new_point, paths[j-k][-1]) < 1/density:
                paths[j].append(None)
            else:
                paths[j].append(new_point)

    for path in paths:
        i = 0
        for i in range(len(path)):
            if path[i] is not None:
                surf.points.append(path[i])
