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


# Sort by distance function. Sorts all atoms in the System by distance from COM of given atoms
def sortbyDist(atoms, net, length=None):
    # If the length of the returned list is not specified return the whole list
    if length is None:
        length = len(net.atoms)
    # Find the point closest to each of the atoms
    loc = [0, 0, 0]
    for i in range(len(atoms)):
        f = i + 1
        loc = loc[0] + atoms[i].loc[0] / f, loc[1] + atoms[i].loc[1] / f, loc[2] + atoms[i].loc[2] / f
    # Initialize the lists
    dist_list = []
    atom_list = []
    # Go through all the atoms in the molecules
    for atom2 in net.atoms:
        # Don't include the atoms in our list of atom
        if atom2 in atoms:
            continue
        # Get the distance between the atoms and subtract their radii
        dist = calc_dist(loc, atom2.loc) - atom2.rad
        dist_list.append(dist)
        atom_list.append(atom2)
    # Selection sort the atom list based off their distances from the point
    for i in range(len(dist_list)):
        low_in = i
        for j in range(i+1, len(dist_list)):
            if dist_list[low_in] > dist_list[j]:
                low_in = j
                dist_list[i], dist_list[low_in] = dist_list[low_in], dist_list[i]
                atom_list[i], atom_list[low_in] = atom_list[low_in], atom_list[i]

    # Return a list with the length specified
    return atom_list[:length]


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
    # If the discriminant is negative then the tangential sphere does not exist.
    if disc > 0:
        circs = []
        Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
        Rs.sort()
        # Go through each circle and gather its points
        for R in Rs:
            # Calculate the vertex based off of our coefficient values and the sphere's radius
            x = Fx0 / F + R * Fx1 / F + l1[0]
            y = Fy0 / F + R * Fy1 / F + l1[1]
            z = Fz0 / F + R * Fz1 / F + l1[2]
            # Add the circle to the circle array
            circs.append([[x, y, z], R])
        return circs[0]
    # Catch for negative discriminant
    else:
        return


"""Translator functions"""


# Rotate points function. Takes in a set of points and a vector and returns a set of rotated points about the origin
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
