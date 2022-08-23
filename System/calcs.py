import numpy as np
import os


# Calculate distance function. Takes in 2 points and returns the distance between them
def calc_dist(l1, l2):
    d = np.sqrt((l1[0]-l2[0])**2+(l1[1]-l2[1])**2+(l1[2]-l2[2])**2)
    return d


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
def calc_atoms_com(atoms):
    # Set the running sum for the x, y, z values to 0
    xtot, ytot, ztot = 0, 0, 0
    for atom in atoms:
        xtot = xtot + atom.loc[0]
        ytot = ytot + atom.loc[1]
        ztot = ztot + atom.loc[2]
    return xtot/len(atoms), ytot/len(atoms), ztot/len(atoms)


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


# Sort by distance function. Sorts all atoms in the System by distance from COM of given atoms
def sortbyDist(atoms, net):
    # Find the point closest to each of the atoms
    loc = [0, 0, 0]
    for i in range(len(atoms)):
        loc = loc[0] + atoms[i].loc[0], loc[1] + atoms[i].loc[1], loc[2] + atoms[i].loc[2]
    loc = [loc[0]/len(atoms), loc[1]/len(atoms), loc[2]/len(atoms)]
    # Initialize the lists
    dist_list, atom_list = [], []
    # Go through all the atoms in the molecules
    for atom in net.atoms:
        # Don't include the atoms in our list of atom
        if atom in atoms:
            continue
        # Get the distance between the atoms and subtract their radii
        dist = calc_dist(loc, atom.loc) - atom.rad
        dist_list.append(dist)
        atom_list.append(atom)
    # Selection sort the atom list based off their distances from the point
    for i in range(len(dist_list)):
        low_in = i
        for j in range(i+1, len(dist_list)):
            if dist_list[low_in] > dist_list[j]:
                low_in = j
                dist_list[i], dist_list[low_in] = dist_list[low_in], dist_list[i]
                atom_list[i], atom_list[low_in] = atom_list[low_in], atom_list[i]
    # Return a list with the length specified
    return atom_list


"""Translator functions"""


# Rotate points function. Takes in a set of points and a vector and rotates the points and the vector so the v = [0,0,1]
def rotate_points(vec, points):
    # Get the vx, vy, vz vector components
    vx, vy, vz = vec
    # If x and y are 0 no transform is needed
    if vy == 0 == vx:
        return points
    elif vz == 0 == vy:
        theta = 0
        phi = np.pi/2
    elif vz == 0 == vx:
        phi = 0
        theta = np.pi/2
    elif vz == 0:
        theta = np.pi / 2
        phi = np.arctan(vy/vx)
    else:
        theta = np.arctan(vx / vz)
        phi = np.arctan(vy / vz)
    # Get variables for sin(theta), cos(theta), sin(phi), cos(phi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    nps = []
    for p in points:
        px, py, pz = round(p[0], 7), round(p[1], 7), round(p[2], 7)
        # Multiplying the x, y rotation matrices gives the following:
        npx = px * cp + pz * sp
        npy = px * st * sp + py * ct - pz * st * cp
        npz = - px * ct * sp + py * st + pz * ct * cp
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
        if v_atoms.issubset(vert.atoms):
            # Return the edge
            return vert
    return
