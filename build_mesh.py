import numpy as np


########################################################################################################################
"""Calculator functions"""


# Calculate distance function. Finds the distance between 2 points
def calc_dist(l1, l2):
    d = np.sqrt((l1[0]-l2[0])**2+(l1[1]-l2[1])**2+(l1[2]-l2[2])**2)
    return d


# Calculate center of mass function. Takes in a set of points and returns the coordinates of the com
def calc_com(points):
    # Set the running sum for the x, y, z values to 0
    xtot, ytot, ztot = 0, 0, 0
    for point in points:
        xtot = xtot + point[0]
        ytot = ytot + point[1]
        ztot = ztot + point[2]
    return xtot/len(points), ytot/len(points), ztot/len(points)


# Calculate angle function. Finds the angle (in rads) between three points. The first being the common point
def calc_angle(p0, p1, p2=None):
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = np.array(p0), np.array(p1)
    else:
        v0, v1 = np.array(p1) - np.array(p0), np.array(p2) - np.array(p0)
    # Get the unit vectors
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
    # Calculate the angle between the two vectors with catches for 180 and 0
    angle = np.arccos(np.clip(np.dot(n0, n1), -1.0, 1.0))
    return angle


# Bisector function. Creates a bisector surface between 2 atoms
def calc_surf(atoms):
    # Make sure that a0 is the atom with the smaller radius
    if atoms[0].rad > atoms[1].rad:
        atoms[0], atoms[1] = atoms[1], atoms[0]
    a0, a1 = atoms

    # Grab the centers of the spheres
    x1, y1, z1 = a0.loc
    x2, y2, z2 = a1.loc

    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = a0.rad - a1.rad
    K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
    d = x1 - x2, y1 - y2, z1 - z2
    J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2

    # Instantiate/reset the hyperboloid coefficient vector lists
    ABC, DEF, GHI = [], [], []
    # Calculate hyperboloid coefficients
    for i in range(3):
        ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
        DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
        GHI.append(-8 * R ** 2 * a0.loc[i] - 4 * K * d[i])

    return ABC + DEF + GHI + [J] + [K] + [d]


# Calculate edge points function. Takes in an edge and a surface and updates the edge's points.
def calc_edge_points(edge, surf, min_dist):
    # Get the location of the base atom
    pa = edge.atoms[0].loc
    # Get the locations of the vertices
    pv0 = np.array(edge.verts[0].loc)
    pv1 = np.array(edge.verts[1].loc)
    # Find the angle made between the edges vertices and the atom
    max_ang = calc_angle(pa, pv0, pv1)
    num_points = max(int(calc_dist(pv0, pv1) / min_dist), 100)
    # Set angle A to be the incremental angle decided by num points
    A = max_ang / num_points
    # Calculate each point along the way
    for i in range(1, num_points):
        # If the edge points are empty set pb to the start vertex. Else get the previous point in the path
        if not edge.points:
            pb = pv0
        else:
            pb = edge.points[-1]
        edge.points.append(find_next_point(pb, pv1, A, surf))


# Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smallest of the surfaces 2 atoms through the point into the surface
def calc_surf_point(surf, point):
    # Grab the function
    if surf.func is None:
        surf.func = calc_surf(surf)
    f = surf.func
    # Get the first atoms in the surfaces list of atoms
    a0 = surf.atoms[0]
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
    if round(b ** 2 - 4 * a * c, 2) >= 0:
        # Check to see of the point on the atom is inside the other atom
        if calc_dist(vi, surf.atoms[1].loc) - surf.atoms[1].rad < 0:
            mag = - abs(min(np.roots([a, b, c])))
        else:
            mag = min(abs(np.roots([a, b, c])))
        return vi + mag * vn


# Find next point function. Finds the
def find_next_point(pn_1, end, d_theta, surf):
    # Get the A angle
    A = d_theta
    # Get the smaller atom's location
    pa = surf.atoms[0].loc
    # Get the location of point b
    pb = pn_1
    # Get the distance between pb and pa
    c = calc_dist(pa, pb)
    # Get the angle between pa, pb and pv1
    B = calc_angle(pb, pa, end)
    # Get the last angle
    C = np.pi - A - B
    # Find a using the law of sines
    a = np.sin(A) * c / np.sin(C)
    # Find the intercept point by adding a to pb
    rn = end - pb
    rn_hat = rn / np.linalg.norm(rn)
    pc = pb + rn_hat * a
    # Calculate where the point intercepts the surface
    return calc_surf_point(surf, pc)


########################################################################################################################


def edge_trace1(surf, min_dist):
    for edge in surf.edges:
        # If the edge points exist already add them to the surfaces edge points and continue to the next edge
        if edge.points:
            surf.edge_points += edge.points
            continue
        calc_edge_points(edge, surf, min_dist)
        # Add the edge's points to the surface's edge points attribute
        surf.edge_points += edge.points


# Edge trace function
def edge_trace(surf, min_dist):
    # Go through each edge on the surface
    for edge in surf.edges:
        # If the edge points exist already continue to the next edge
        if not edge.points:
            calc_edge_points(edge, surf, min_dist)
    # Call the recursive edge tracing function
    edges = surf.edges
    # Get the first edge and vert
    v0, vert = edges[0].verts
    surf.edge_points = edges[0].points
    # Find the edge that is closest to vert and add the points accordingly
    v0_found = False
    while not v0_found:
        # Check the edges for similar vertices
        for edge in edges[1:]:
            # If the first vertex in the edge's vertex list is equal to vert, add the points
            if edge.verts[0] == vert:
                surf.edge_points += edge.points
                vert = edge.verts[1]
            elif edge.verts[1] == vert:
                surf.edge_points = edge.points[::-1]
                vert = edge.verts[0]
        if vert.loc == v0.loc:
            v0_found = True


# Make mesh function. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf, min_dist):
    # Set the atoms in the surface to make the smaller one listed first
    if surf.atoms[0].rad > surf.atoms[1].rad:
        surf.atoms = surf.atoms[1], surf.atoms[0]
    # Check to see if the edges' points have been recorded yet
    edge_trace1(surf, min_dist)
    # For each edge point set up a path list.
    paths = [[surf.edge_points[i]] for i in range(len(surf.edge_points))]
    # Grab the smallest of the 2 surface atoms' location
    pa = surf.atoms[0].loc
    # Calculate the center of mass point of the edge points and where it maps on the surface
    com = calc_surf_point(surf, calc_com(surf.edge_points))
    # Set up a list of end points
    ends = [com for i in range(len(paths))]
    # Get the angles between the edge points and the end points
    angs = []
    for i in range(len(paths)):
        # Calculate the angle for each path
        angs.append(calc_angle(pa, paths[i][0], ends[i]))
    # Get the maximum path
    max_path_ndx = angs.index(max(angs))
    max_path = paths[max_path_ndx][0]
    # Decide how many rings based off of the ellipticity and density
    num_rings = max(int(calc_dist(max_path, ends[max_path_ndx]) / min_dist), 10)
    # Get the incremental angle increases
    dthetas = [angs[i]/num_rings for i in range(len(angs))]
    # Set the pn_1 point to infinity
    pn_1 = [np.inf, np.inf, np.inf]
    num_paths = len(paths)
    # Go through ring by ring
    for j in range(num_rings):
        # Go through each of the remaining paths
        i = 0
        while i < num_paths:
            # Get the next point along the path
            pn = find_next_point(paths[i][-1], ends[i], dthetas[i], surf)
            # Check to see of the new point is too close to the previous point and the path has to end
            if calc_dist(pn, pn_1) < min_dist:
                # Add the path to the surfaces points and remove it from the paths list
                surf.points += paths.pop(i)[1:]
                ends.pop(i)
                dthetas.pop(i)
                num_paths -= 1
            else:
                # Set the pn_1 to pn and add it to the path
                pn_1 = pn
                paths[i].append(pn)
            # Increment i
            i += 1
    # Add the remaining paths to the surface
    for path in paths:
        surf.points += path[1:]



# Build meshes function. Runs make_mesh on all surfaces in the network
def build_meshes(sys, min_dist=None):
    # Set the minimum distance
    if min_dist is None:
        min_dist = calc_dist(sys.net.edges[0].verts[0].loc, sys.net.edges[0].verts[1].loc) / 30
    # Make each surface
    for surf in sys.net.surfs:
        make_mesh(surf, min_dist)
