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


# Calculate edge points function. Takes in an edge and a surface and updates the edge's points
def calc_edge_points(edge, surf):
    # Get the location of the base atom
    pa = edge.atoms[0].loc
    # Get the locations of the vertices
    pv0 = np.array(edge.verts[0].loc)
    pv1 = np.array(edge.verts[1].loc)
    # Add the first vertex to the edges points
    edge.points = [pv0]
    # Find the angle made between the edges vertices and the atom
    max_ang = calc_angle(pa, pv0, pv1)
    num_points = int(np.degrees(max_ang))  #############################################################################
    # Set angle A to be the incremental angle decided by num points
    A = max_ang / num_points
    # Go calculate each point along the way
    for i in range(num_points):
        # If the edge points are empty set pb to the start vertex
        if not edge.points:
            pb = pv0
        # Else get the location of the previous point
        else:
            pb = edge.points[-1]
        # Get the distance between pb and pa
        c = calc_dist(pa, pb)
        # Get the angle between pa, pb and pv1
        B = calc_angle(pb, pa, pv1)
        # Get the last angle
        C = np.pi - A - B
        # Find a using the law of sines
        a = np.sin(A) * c / np.sin(C)
        # Find the intercept point by adding a to pb
        rn = pv1 - pb
        rn_hat = rn/np.linalg.norm(rn)
        pc = pb + rn_hat * a
        # Calculate where the point intercepts the surface
        pn = calc_spnt(surf, pc)
        # Add the point to the edges list of points
        edge.points.append(pn)
    # Add the destination vertex point to the list of points
    edge.points.append(pv1)


# Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smaller of the surfaces 2 atoms through the point into the surface
def calc_spnt(surf, point):
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
    if b ** 2 - 4 * a * c > 0:
        mag = min(abs(np.roots([a, b, c])))
        return vi + mag*vn
    else:
        print(None)


# Calculate overlap points function. Used to calculate the points in and at the overlap of two intersecting spheres
def calc_olap_points(surf, com):
    # Not sure how to do this yet. Maybe calc_spoint works maybe it doesn't
    pass


########################################################################################################################

# Edge trace function
def edge_trace(surf):
    # Go through each edge on the surface
    for edge in surf.edges:
        # If the edge points exist already continue to the next edge
        if edge.points:
            continue
        calc_edge_points(edge, surf)
        # Add the edge's points to the surface's edge points attribute
        surf.edge_points += edge.points


# Make mesh function. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf, density=100):
    # Check to see if the edges' points have been recorded yet
    if not surf.edge_points:
        edge_trace(surf)
    # Grab the smaller of the 2 surface atoms
    a0 = surf.atoms[0]
    # Calculate the center of mass point of the edge points
    com = calc_com(surf.edge_points)
    # Find where com maps on the surface
    com = calc_spnt(surf, com)
    # Make a list of paths and angles
    paths = []
    ngl_lst = []
    # Create (and calculate the maximum angle for) each path
    for i in range(len(surf.edge_points)):
        paths.append(surf.edge_points[i])
        ngl_lst.append(calc_angle(a0.loc, com, surf.edge_points[i]))
    # Get the maximum path information
    max_path_ndx = ngl_lst.index(max(ngl_lst))
    max_path = paths[max_path_ndx]
    # Decide how many rings based off of density
    num_rings = int(calc_dist(max_path, com) * density/100)
    # Create a list of dthetas
    dthetas = []
    for i in range(len(surf.edge_points)):
        dthetas.append(ngl_lst[i]/num_rings)
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
            print(b)
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
        for i in range(len(path)):
            if path[i] is not None:
                surf.points.append(path[i])
    return surf


# Make mesh function. Calculates points on the surface
def make_mesh1(surf, density=0.0001):
    # Check to see if the edges' points have been recorded yet
    if not surf.edge_points:
        edge_trace(surf)
    # For each edge point set up a path list.
    paths = [[surf.edge_points[i]] for i in range(len(surf.edge_points))]
    print(paths)
    # Grab the smaller of the 2 surface atoms' location
    pa = surf.atoms[0].loc
    # Calculate the center of mass point of the edge points and where it maps on the surface
    com = calc_spnt(surf, calc_com(surf.edge_points))
    # Check to see if the atoms overlap
    if calc_dist(pa, surf.atoms[1].loc) - (surf.atoms[0].rad + surf.atoms[1].rad) < 0:
        ends = calc_olap_points(surf, com)
    else:
        ends = [com for i in range(len(paths))]
    # Get the angles between the edge points and the end points
    angs = []
    for i in range(len(paths)):
        # Calculate the angle for each path
        angs.append(calc_angle(pa, paths[i][0], ends[i]))
        # Get the maximum path information
    max_path_ndx = angs.index(max(angs))
    max_path = paths[max_path_ndx][0]
    # Decide how many rings based off of density
    num_rings = int(calc_dist(max_path, com) * 100)
    # Get the incremental angle increases
    dthetas = [angs[i]/num_rings for i in range(len(angs))]
    # Set the pn_1 point to infinity
    pn_1 = [np.inf, np.inf, np.inf]
    num_paths = len(paths)
    i = 0
    # Go through ring by ring
    for j in range(num_rings):
        # Go through each of the remaining paths
        while i < num_paths:
            # Get the A angle
            A = dthetas[i]
            # Get the location of point b
            pb = paths[i][-1]
            # Get the distance between pb and pa
            c = calc_dist(pa, pb)
            # Get the angle between pa, pb and pv1
            B = calc_angle(pb, pa, com)
            # Get the last angle
            C = np.pi - A - B
            # Find a using the law of sines
            a = np.sin(A) * c / np.sin(C)
            # Find the intercept point by adding a to pb
            rn = ends[i] - pb
            rn_hat = rn / np.linalg.norm(rn)
            pc = pb + rn_hat * a
            # Calculate where the point intercepts the surface
            pn = calc_spnt(surf, pc)
            # Check to see of the new point is too close to the previous point and the path has to end
            print(calc_dist(pn, pn_1))
            if calc_dist(pn, pn_1) < density:
                print(True)
                # Add the path to the surfaces points and remove it from the paths list
                surf.points += paths.pop(i)
                ends.pop(i)
                dthetas.pop(i)
                num_paths -= 1
            else:
                print(False)
                # Set the pn_1 to pn and add it to the path
                pn_1 = pn
                paths[i].append(pn)
            i += 1

    print(surf.points)
