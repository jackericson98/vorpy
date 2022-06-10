"""Imports"""
import numpy as np
from objects import Atom, Edge, Vertex, Surface

"""Operations:
    1. build_network: System -> Network
    2. first_site: Network -> Vertex
    3. get_circle(
    3. make_edges: Vertex, Network -> Network
    4. calc_circle: atoms -> loc, rad
    5. calc_vertex: atoms -> Vertex
    """


########################################################################################################################
"""Calculator functions"""


# Sort by distance function. Sorts all atoms in the system by distance from COM of given atoms
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
        dist = (np.sqrt((loc[0]-atom2.loc[0])**2 + (loc[1]-atom2.loc[1])**2 + (loc[2]-atom2.loc[2])**2)) - atom2.rad
        dist_list.append(dist)
        atom_list.append(atom2)

    # Sort the list of atoms by the distances from the given atom
    atom_list = [atom for _, atom in sorted(zip(dist_list, atom_list))]
    # Return a list with the length specified
    return atom_list[:length]


# Move sphere function. Takes in a location, an Atom object and a direction and updates the Atom's location
def move(loc, atom, to_home=False):
    # Change whether we are adding or subtracting the location to the sphere's location.
    d = 1
    if not to_home:
        d = -1
    # Update the atom's location
    atom.loc[0] = atom.loc[0] + d * loc[0]
    atom.loc[1] = atom.loc[1] + d * loc[1]
    atom.loc[2] = atom.loc[2] + d * loc[2]


# Calculate circle function. Takes in 3 atoms, calculates the center and radius of inscribed circle and returns them
def calc_circ(atoms):
    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad
    # Move each sphere to surround the base sphere now located at the origin
    for sphere in atoms:
        move(l1, sphere)
    # Get the relevant variables
    R2, R3 = atoms[1].rad, atoms[2].rad
    x2, y2, z2 = atoms[1].loc
    x3, y3, z3 = atoms[2].loc
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
    # Find the radius of the tangential circle using the quadratic formula
    a = (Fx1 ** 2 + Fy1 ** 2 + Fz1 ** 2) / F ** 2 - 1
    b = 2 * (Fx0 * Fx1 + Fy0 * Fy1 + Fz0 * Fz1) / F ** 2 - 2 * R1
    c = (Fx0 ** 2 + Fy0 ** 2 + Fz0 ** 2) / F ** 2 - R1 ** 2
    # Calculate the discriminant.
    disc = b ** 2 - 4 * a * c
    # If the discriminant is negative then the tangential sphere does not exist.
    if disc > 0:
        circs = []
        Rs = [R for R in np.roots([a, b, c]) if np.isreal(R) and R > 0]
        Rs.sort()
        # Go through each circle and gather its points
        for R in Rs:
            # Calculate the vertex based off of our coefficient values and the sphere's radius
            x = Fx0 / F + R * Fx1 / F + l1[0]
            y = Fy0 / F + R * Fy1 / F + l1[1]
            z = Fz0 / F + R * Fz1 / F + l1[2]
            # Add the circle to the circle array
            circs.append([[x, y, z], R])

        return circs
    else:
        # Catch for negative discriminant
        print("Negative discriminant!")
        return


# Calculate vertex function. Takes in 4 atoms, calculates the center and radius of the inscribed sphere and returns them
def calc_vertex(atoms):

    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad

    # Set the radii and x, y, z values for the 3 spheres
    R2, R3, R4 = atoms[1].rad, atoms[2].rad, atoms[3].rad
    x2, y2, z2 = atoms[1].loc[0] - l1[0], atoms[1].loc[1] - l1[1], atoms[1].loc[2] - l1[1]
    x3, y3, z3 = atoms[2].loc[0] - l1[0], atoms[2].loc[1] - l1[1], atoms[2].loc[2] - l1[1]
    x4, y4, z4 = atoms[3].loc[0] - l1[0], atoms[3].loc[1] - l1[1], atoms[3].loc[2] - l1[1]

    # Calculate our system of linear equations coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (R1 - R2), R1 ** 2 - R2 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (R1 - R3), R1 ** 2 - R3 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3, d3, f3 = 2 * x4, 2 * y4, 2 * z4, 2 * (R1 - R4), R1 ** 2 - R4 ** 2 + x4 ** 2 + y4 ** 2 + z4 ** 2

    A, B, C, d, f = [a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f3]

    # Calculate the ranks of the matrices
    ABC_rank = np.linalg.matrix_rank([A, B, C])
    m_rank = np.linalg.matrix_rank([A, B, C, d])
    f_rank = np.linalg.matrix_rank([A, B, C, d, f])

    # Calculate the F values
    F, F10, F11, F20, F21, F30, F31 = np.linalg.det([A, B, C]), np.linalg.det([f, B, C]), \
                                      np.linalg.det([[-d[0], -d[1], -d[2]], B, C]), np.linalg.det([A, f, C]), \
                                      np.linalg.det([A, [-d[0], -d[1], -d[2]], C]), np.linalg.det([A, B, f]), \
                                      np.linalg.det([A, B, [-d[0], -d[1], -d[2]]])
    # Instantiate our root arrays
    xs, ys, zs, Rs = [], [], [], []
    verts = []
    # Case 1:
    if ABC_rank == 3 and m_rank == 3 and f_rank == 3:
        # Calculate the radius polynomial coefficients
        a = (F11 ** 2 + F21 ** 2 + F31 ** 2) / F ** 2 - 1
        b = 2 * (F10 * F11 + F20 * F21 + F30 * F31) / F ** 2 - 2 * R1
        c = (F10 ** 2 + F20 ** 2 + F30 ** 2) / F ** 2 - R1 ** 2
        # If the discriminant is positive, find the real positive roots of the quadratic
        if -4*a*c + b**2 > 0:
            Rs = [R for R in np.roots([a, b, c]) if np.isreal(R) and R > 0]
        # Instantiate the verts array
        verts = []
        # Go through each radius and calculate the vertex
        for R in Rs:
            x, y, z = F10/F + R*F11/F, F20/F + R*F21/F, F30/F + R*F31/F
            # Move the vertex back to the actual location of the atoms
            verts.append([[x + l1[0], y + l1[1], z + l1[2]], R])

    # Case 2:
    elif ABC_rank == 2 and m_rank == 3 and f_rank == 3:
        # Case 2 subcases:
        # Case 2.1
        if np.linalg.matrix_rank([A, B, d]) == 3:
            # Calculate the z value polynomial coefficients
            a = F**2 + F11**2 + F21**2 - F31**2
            b = 2*(F10*F11 + F20*F21 - F30*F31 - F*F31*R1)
            c = F10**2 + F20**2 - (F30 + F*R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                zs = [z for z in np.roots([a, b, c]) if np.isreal(z) and z > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for z in zs:
                x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l1[0], y + l1[1], z + l1[2]], R])

        # Case 2.2
        elif np.linalg.matrix_rank([A, d, C]) == 3:
            # Calculate the z value polynomial coefficients
            a = F ** 2 + F11 ** 2 - F21 ** 2 + F31 ** 2
            b = 2 * (F10 * F11 - F20 * F21 + F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F30 ** 2 - (F20 + F * R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                ys = [y for y in np.roots([a, b, c]) if np.isreal(y) and y > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for y in ys:
                x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l1[0], y + l1[1], z + l1[2]], R])

        # Case 2.3
        elif np.linalg.matrix_rank([d, B, C]):
            # Calculate the z value polynomial coefficients
            a = F ** 2 + F11 ** 2 + F21 ** 2 - F31 ** 2
            b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F20 ** 2 - (F30 + F * R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                xs = [x for x in np.roots([a, b, c]) if np.isreal(x) and x > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for x in xs:
                R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l1[0], y + l1[1], z + l1[2]], R])

    return verts


########################################################################################################################
"""Gather information functions"""

# Doublet making function. Takes in a doublet and adds the two vertices tot herir respective places.
def doublet(verts, net):
    return v0


# Find vertex function. Takes in an edge, network and finds the other site along that edge.
def find_vertex(edge, net, a0=None):
    # Set up our vertex
    vert, doublet = None, None
    # Grab the 50 closest atoms
    prox_list = sortbyDist(edge.atoms, net, length=50)
    # Go through each of the closest atoms
    for a3 in prox_list:
        # Calculate the value of the vertex made from the edge atoms and our test atom
        verts = calc_vertex([edge.atoms + a3])
        # Check the different vertex cases: None, singlet, doublet
        if len(verts) == 0:
            continue
        elif len(verts) == 1:
            vert = verts[0]
        elif len(verts) == 2:
            # Calculate the distance between
            d1 = np.sqrt((verts[0][0] - a0.loc[0])**2 + (verts[0][1] - a0.loc[1])**2 + (verts[0][2] - a0.loc[2])**2)
            d2 = np.sqrt((verts[1][0] - a0.loc[0])**2 + (verts[1][1] - a0.loc[1])**2 + (verts[1][2] - a0.loc[2])**2)
            vert = verts[0]
            if d2 < d1:
                vert = verts[1]
            doublet = vert


    if vert == doublet:
        vert = doublet(vert, net)
    return vert


# Get initial vertex function. Finds an optimal starting vertex for the network.
def get_v1(net):
    # Grab the first atom in the network. This will be replaced later with an optimized
    a0 = net.atoms[0]
    # Find it's closest neighbor
    neighbors0 = sortbyDist([a0], net)
    a1 = neighbors0[1]
    # Find the smallest circle you can make with a0, a1 and a third atom
    a2 = find_circle(a0, a1, net)[0]
    # Find the first site
    e0 = Edge([a0, a1, a2], None)
    v1 = find_vertex(e0, net)
    for i in range(4):
        print(net.atoms.index(v1.atoms[i]))

    return v1


########################################################################################################################
"""Recursive network finding function"""


# Find edges function. Recursively traces out the network and records vertex locations,
def find_edges(vertex, net):

    # Find the edges
    for edge in vertex.edges:
        # If the edge already has two vertices skip it
        if len(edge.verts) == 2:
            continue

        # Find the smallest vertex that can be made with the atoms in the edge that is not the vertex
        vn = find_vertex(edge, net, vertex)

        net.verts.append(vn)
        edge.verts.append(vn)
        find_edges(vn, net)
        print(vn.loc)



########################################################################################################################


# Build Network function. Takes in a system, runs as a shell for the recursive next_site function and returns a Network
def build_network(mol):
    # Grab the network object from the molecule.
    net = mol.net
    # Find the first vertex
    v1 = get_v1(net)
    # Add the vertex to the edge of the vertex

    find_edges(v1, net)

    return myNet