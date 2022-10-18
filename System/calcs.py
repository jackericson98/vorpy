import numpy as np


# Calculate distance function. Takes in 2 points and returns the distance between them
def calc_dist(l1, l2):
    # Pythagorean theorem
    return np.sqrt((l1[0]-l2[0])**2+(l1[1]-l2[1])**2+(l1[2]-l2[2])**2)


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


# Calculate tetrahedron volume function. Calculated the volume of a tetrahedron defined by its vertices
def calc_tetra_vol(p0, p1, p2, p3):
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]
    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    vol = (1/6)*abs(np.dot(r03, np.cross(r01, r02)))
    return vol


# Calculate triangle are function. Takes in 3 points in 3 space and returns the area of the triangle created by them
def calc_tri(points):
    # Get the two triangles vectors
    AB = np.array(points[0]) - np.array(points[1])
    AC = np.array(points[0]) - np.array(points[2])
    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm((np.cross(AB, AC)))


# Calculate surface area function. Takes in a
def calc_sa(surf):
    sa = 0
    for tri in surf.tris:
        p0, p1, p2 = surf.points[tri[0]], surf.points[tri[1]], surf.points[tri[2]]
        sa += calc_tri([p0, p1, p2])
    # Return the total surface area
    return sa


# Calculate cell volume function. Grabs the points in a cell and calculates the volume made by the tetrahedrons
def calc_vol(atom):
    vol = 0
    # Go through each surface on the atom
    for surf in atom.surfs:
        for tri in surf.tris:
            p0, p1, p2, p3 = atom.loc, surf.points[tri[0]], surf.points[tri[1]], surf.points[tri[2]]
            vol += calc_tetra_vol(p0, p1, p2, p3)
    return vol


# Calculate center of mass function. Takes in a set of points and returns the coordinates of the com
def calc_com(atoms=None, points=None):
    if atoms:
        points = [atoms[i].loc for i in range(len(atoms))]
    # Set the running sum for the x, y, z values to 0
    tots = [0 for _ in range(len(points[0]))]
    for point in points:
        for i in range(len(points[0])):
            tots[i] += point[i]
    return [tots[i]/len(points) for i in range(len(points[0]))]


# Calculate edges center of mass function. Takes in a surface or edges and returns the center of mass of the points
def calc_edges_com(edges=None, surf=None, points=None):
    # Set up the points list
    myPoints = []
    # Check to see if edges were given
    if edges:
        for edge in edges:
            myPoints += edge.points
    # Check to see if a surface was given
    elif surf:
        for edge in surf.edges:
            myPoints += edge.points
    # Check to see if points were given
    elif points:
        myPoints = points
    else:
        print("Please give a surface or mesh!")\

    # Find the sum of the points
    x_tot, y_tot, z_tot = 0, 0, 0
    for point in myPoints:
        x_tot += point[0]
        y_tot += point[1]
        z_tot += point[2]
    # Return the center of mass of the edge points
    return [x_tot / len(myPoints), y_tot / len(myPoints), z_tot / len(myPoints)]


# Calculate circle function. Takes in 3 atoms, calculates the center and radius of inscribed circle and returns them
def calc_circ(atoms):
    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad
    # Get the relevant variables
    R2, R3 = atoms[1].rad, atoms[2].rad
    x2, y2, z2 = atoms[1].loc[0] - l1[0], atoms[1].loc[1] - l1[1], atoms[1].loc[2] - l1[2]
    x3, y3, z3 = atoms[2].loc[0] - l1[0], atoms[2].loc[1] - l1[1], atoms[2].loc[2] - l1[2]
    # Calculate coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (R1 - R2), R1 ** 2 - R2 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (R1 - R3), R1 ** 2 - R3 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3 = y2*z3 - z2*y3, z2*x3 - x2*z3, x2*y3 - y2*x3
    # More coefficients
    F = a3*b2*c1 - a2*b3*c1 - a3*b1*c2 + a1*b3*c2 + a2*b1*c3 - a1*b2*c3
    Fx0 = b3*c2*f1 - b2*c3*f1 - b3*c1*f2 + b1*c3*f2
    Fx1 = b3*c2*d1 - b2*c3*d1 - b3*c1*d2 + b1*c3*d2
    Fy0 = - a3*c2*f1 + a2*c3*f1 + a3*c1*f2 - a1*c3*f2
    Fy1 = - a3*c2*d1 + a2*c3*d1 + a3*c1*d2 - a1*c3*d2
    Fz0 = a3*b2*f1 - a2*b3*f1 - a3*b1*f2 + a1*b3*f2
    Fz1 = a3*b2*d1 - a2*b3*d1 - a3*b1*d2 + a1*b3*d2
    # Catch for F=0 (i.e. no circle exists)
    if F == 0:
        return
    # Find the radius of the tangential circle using the quadratic formula
    a = (Fx1 ** 2 + Fy1 ** 2 + Fz1 ** 2) / F ** 2 - 1
    b = 2 * (Fx0 * Fx1 + Fy0 * Fy1 + Fz0 * Fz1) / F ** 2 - 2 * R1
    c = (Fx0 ** 2 + Fy0 ** 2 + Fz0 ** 2) / F ** 2 - R1 ** 2
    # Calculate the discriminant.
    disc = b ** 2 - 4 * a * c
    # If the discriminant is negative then the tangential circle does not exist.
    if round(disc, 10) > 0:
        # Grab the two roots
        Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
        # If there is only one root return it
        if len(Rs) == 1:
            R = Rs[0]
        # If there are 2 roots choose between them
        else:
            # If the smaller of the two roots is negative return the other root
            if min(Rs) < 0:
                R = max(Rs)
            # If they're both positive, return the smaller of the two
            elif Rs[0] > 0 and Rs[1] > 0:
                R = min(Rs)
            # If they're both negative return
            else:
                return
        # Calculate the vertex based off of our coefficient values and the sphere's radius
        x = Fx0 / F + R * Fx1 / F + l1[0]
        y = Fy0 / F + R * Fy1 / F + l1[1]
        z = Fz0 / F + R * Fz1 / F + l1[2]
        return [[x, y, z], R]


# Calculate bisector value function. Takes in a bisector and point and returns the value of that point in regard to
# the bisector function. A point that is actually on the bisector surface should return a value of 0.
def calc_bisector_val(f, point):
    x, y, z = point
    val = f[0] * x ** 2 + f[1] * y ** 2 + f[2] * z ** 2 + \
          f[3] * x * y + f[4] * y * z + f[5] * z * x + \
          f[6] * x + f[7] * y + f[8] * z + f[9]
    return val


# Inverse Jacobian Function. Takes in 3 bisector functions and point in 3 space and returns the inverse Jacobian matrix.
# Only works with hyperboloid surfaces. Could be developed to more general case, but this is faster for now.
def inv_jac(funcs, point):
    # Get functions and point
    f1, f2, f3 = funcs
    x, y, z = point
    # Create the Jacobian Matrix
    jac_mat = np.array([[2 * f1[0] * x + f1[3] * y + f1[5] * z + f1[6], 2 * f1[1] * y + f1[3] * x + f1[4] * z + f1[7],
                         2 * f1[2] * z + f1[4] * y + f1[5] * x + f1[8]],
                        [2 * f2[0] * x + f2[3] * y + f2[5] * z + f2[6], 2 * f2[1] * y + f2[3] * x + f2[4] * z + f2[7],
                         2 * f2[2] * z + f2[4] * y + f2[5] * x + f2[8]],
                        [2 * f3[0] * x + f3[3] * y + f3[5] * z + f3[6], 2 * f3[1] * y + f3[3] * x + f3[4] * z + f3[7],
                         2 * f3[2] * z + f3[4] * y + f3[5] * x + f3[8]]])

    if jac_mat.shape[0] == jac_mat.shape[1] and np.linalg.matrix_rank(jac_mat) == jac_mat.shape[0]:
        # Calculate the inverse of this matrix and return it
        return np.linalg.inv(jac_mat)


"""Translator functions"""


# Rotate points function. Takes in a set of points and a vector and rotates the points and the vector so the v = [0,0,1]
def rotate_points(vec, points, reverse=False):
    # Get the vx, vy, vz vector components
    vx, vy, vz = vec
    # If vy or vz are zero we need a catch for divide by zero error.
    if round(vy, 2) == 0:
        phi = np.pi / 2
    else:
        phi = np.arctan(vx / vy)
    if round(vz, 2) == 0:
        theta = np.pi / 2
    else:
        theta = np.arctan(vy / vz)
    # If the points are to be sent back, provide the negative values for the angles
    if reverse:
        theta, phi = -theta, -phi
    # Get variables for sin(theta), cos(theta), sin(phi), cos(phi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    nps = []
    for p in points:
        px, py, pz = round(p[0], 7), round(p[1], 7), round(p[2], 7)
        # Multiplying the x, y rotation matrices gives the following:
        npx = px * cp - py * sp
        npy = px * ct * sp + py * ct * cp - pz * st
        npz = px * st * sp + py * st * cp + pz * ct
        # Add the new points to the list
        nps.append([npx, npy, npz])
    return nps


"""System checks"""


# Check surf function. Takes in a set of atoms and a list of surfs and returns the corresponding surf or None if no surf
def check_surf(s_atoms, surf_list):
    # Go through each surf in the surf list
    for surf in surf_list:
        # Check if the given atoms correspond to the atoms in the surf
        if s_atoms.issubset(surf.atoms):
            # Return the surf
            return surf
    return


# Check edge function. Takes in a set of atoms and a list of edges and returns the corresponding edge or None if no edge
def check_edge(e_atoms, edge_list):
    # Go through each edge in the edge list
    for edge in edge_list:
        # Skip for doublets
        if edge.doublet:
            continue
        # Check if the given atoms correspond to the atoms in the edge
        if e_atoms.issubset(edge.atoms):
            # Return the edge
            return edge
    return


# Check vert function. Takes in a set of atoms and a list of verts and returns the corresponding edge or None if no vert
def check_vert(v_atoms, vert_list):
    # Go through each edge in the edge list
    for vert in vert_list:
        # Check if the given atoms correspond to the atoms in the edge
        if set(v_atoms).issubset(vert.atoms):
            # Return the edge
            return vert
    return


# Search vertices function. Searches a list of indices of atoms sorted by smallest atom and where the vertex would be
def search_verts(my_vert_ndxs, vert_ndx):
    # If the length of the test list is equal to 0 return the next index
    if len(my_vert_ndxs) <= 1:
        # If there exists one vertex already and the new vertex is less than the old vertex return 1
        if len(my_vert_ndxs) > 0 and vert_ndx > my_vert_ndxs[0]:
            return 1
        # Otherwise, return 0
        return 0
    # Get the middle of the list of vertices
    mid_list_ndx = len(my_vert_ndxs) // 2
    # If the search element (my_list) is greater than the test element (test_lol) search the lower half of test_lol
    if vert_ndx > my_vert_ndxs[mid_list_ndx]:
        my_vert_ndx = search_verts(my_vert_ndxs[mid_list_ndx:], vert_ndx)
        return my_vert_ndx + len(my_vert_ndxs[:mid_list_ndx])
    # If the search element (my_list) is less than the test element (test_lol) search the upper half of test_lol
    elif vert_ndx < my_vert_ndxs[mid_list_ndx]:
        my_vert_ndx = search_verts(my_vert_ndxs[:mid_list_ndx], vert_ndx)
        return my_vert_ndx
    # If the search element (my_list) is greater than the test element (test_lol) search the lower half of test_lol
    elif vert_ndx == my_vert_ndxs[mid_list_ndx]:
        return mid_list_ndx
