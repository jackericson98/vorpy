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
    print(atoms[0].loc, atoms[1].loc, atoms[2].loc)
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
    # Catch for F=0 (i.e. no circle exists)
    if F == 0:
        return [[0, 0, 0], np.inf]
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
            if R < 0:
                print('negative radius')
            else:
                print('positive radius')
            x, y, z = F10/F + R*F11/F, F20/F + R*F21/F, F30/F + R*F31/F
            # Move the vertex back to the actual location of the atoms
            verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R))

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
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R))

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
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R))

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
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R))

    if verts == []:
        return None
    else:
        return verts[0]


# Calculate edge function. Chases an edge toward the next vertex
def calc_edge(edge, net, dt=None):
    # Find the closest neighbors
    neighbors = sortbyDist(edge.atoms, net, length=50)
    # Get the old vertex
    v0 = edge.verts[0]
    # Estimate a working dt
    if dt is None:
        dt = edge.atoms[0].rad / 20

    an = None
    # Find the atom from the v0 that is not in the edge's atom list
    for atom in v0.atoms:
        # Set our changing atom to the outlier atom found above and change the radius by 5%
        if atom not in edge.atoms:
            an = atom

    # Calculate the bottleneck of the edge
    bn = calc_circ(edge.atoms)[0][1]
    # Adjust the size of the atom to fit through the bottleneck
    if bn < 1.05*an.rad:
        an.rad = 0.95*bn

    # Find the vertex between the edge atoms and the adjusted atom
    vn = calc_vertex(edge.atoms+[an])

    # If we get a None vertex the shrink went too far. Keep increasing radius until a vertex is found.
    while vn is None:
        an.rad = an.rad * 1.01
        vn = calc_vertex(edge.atoms+[an])
    # If the radius is larger than the bottleneck, continue and hope that the other side of the edge will be able to
    if vn.rad > bn:
        return

    # Find the initial direction by getting the vector between the new vertex formed after the atom got smaller
    dr = np.array([v0.loc[0] - vn.loc[0], v0.loc[1] - vn.loc[1], v0.loc[2] - vn.loc[2]])

    elen = 0
    vfound = False
    vert = None
    # Keep adding points to the edge until the next vertex is found or the edge left the network
    while not vfound:

        # Normalize the direction vector
        dr_mag = np.sqrt(dr.dot(dr))
        dr = dr / dr_mag

        # Add up the length of the edge
        elen += dr_mag
        if elen > net.rad:
            edge.verts.append(None)
            return None

        # Record vns location before changing it
        vn_1 = vn
        # Move the atom along the direction of the edge by dt increments
        print(dt)
        an.loc = an.loc[0] + dt*dr[0], an.loc[1] + dt*dr[1], an.loc[2] + dt*dr[2]
        # Calculate the new vertex
        vn = calc_vertex(edge.atoms + [an])
        # Add the vertex location to the edges points
        edge.points.append(vn.loc)
        # Find the new move direction by finding the direction from vn-1 to vn
        dr = np.array([vn.loc[0] - vn_1.loc[0], vn.loc[1] - vn_1.loc[1], vn.loc[2] - vn_1.loc[2]])

        # Check to see if we have passed a vertex
        for vert in neighbors:
            # Calculate the vectors between the vertex and the new and old edge points
            d1 = np.array([vn_1.loc[0] - vert.loc[0], vn_1.loc[0] - vert.loc[0], vn_1.loc[0] - vert.loc[0]])
            d2 = np.array([vn.loc[0] - vert.loc[0], vn.loc[0] - vert.loc[0], vn.loc[0] - vert.loc[0]])
            # Check to see if the vertex is in between the new and old edge points
            if np.sqrt(d1.dot(d1)) <= dr_mag and np.sqrt(d2.dot(d2)) <= dr_mag:
                # If so, we have found our vert and exit
                vfound = True

    edge.verts.append(vert)
    net.verts.append(vert)

    return vert


# Check vertex function. Used to see if the new vertex is either better or worse than the old or not allowed.
def check_vertex(vert, myVert):
    return myVert


########################################################################################################################
"""Gather information functions"""


# Doublet making function. Takes in a doublet and adds the two vertices to their respective places.
def doublet(verts, net):
    return verts[0]


# Find circle function. Finds the smallest circle between the two given atoms and every other atom and return that atom
def find_circle(a0, a1, net, num_checks=12):
    # Instantiate variables
    a2 = None
    neighbors = sortbyDist([a0, a1], net)
    rad = np.inf
    # Go through the num_checks closest atoms and find the smallest circle
    for atom in neighbors[1:num_checks]:
        # Check if any are the same as the new atom
        if a0 == atom or a1 == atom:
            continue
        # Calculate the radius of the circle made by the three atoms
        new_rad = calc_circ([a0, a1, atom])[0][1]
        # Check the new radius against the smallest found and make it the smallest if it is
        if new_rad < rad:
            rad = new_rad
            a2 = atom
    # Return the atom found to have the smallest circle
    return a2


# Get initial vertex function. Finds an optimal starting vertex for the network.
def find_v0(net):
    # Grab the first atom in the network. This will be replaced later with an optimized
    a0 = net.atoms[0]
    # Find it's closest neighbor
    neighbors0 = sortbyDist([a0], net)
    a1 = neighbors0[1]
    # Find the smallest circle you can make with a0, a1 and a third atom
    a2 = find_circle(a0, a1, net)
    # Find the first site by choosing the smallest interstitial sphere that can be made with the 3 atoms and the 50
    # closest atoms
    neighbors2 = sortbyDist([a0, a1, a2], net)
    r = np.inf
    myVert, my_an = None, None
    # Go through the closest atoms to the triplet and find the smallest vertex that can be made with a neighbor
    for an in neighbors2:
        # Calculate the vertex and check if None
        vert = calc_vertex([a0, a1, a2, an])
        if vert is None:
            continue
        # Check if the vn has a smaller radius than the current smallest radius
        if vert.rad < r:
            r = vert.rad
            myVert = vert
            my_an = an
    # Add connections to the network
    myVert.atoms = a0, a1, a2, my_an
    net.verts.append(myVert)
    return myVert


########################################################################################################################
"""Recursive network finding function"""


# Find edges function. Recursively traces out the network and records vertex locations,
def find_edges(vertex, net, nedges=4, recurse=True):
    # Create the edge objects or grab them from the network and connect them
    for i in range(nedges):
        # Go through each iteration of atom combinations to make the edges.
        myAtoms = [vertex.atoms[i], vertex.atoms[(i+1) % 4], vertex.atoms[(i+2) % 4]]

        # If there aren't edges, it is the first pass and the edges don't need to be checked.
        if not net.edges:
            vertex.edges.append(Edge(myAtoms, vertex))
            continue

        # Check the atoms against each edge in the networks' atoms
        for edge in net.edges:
            # Reset the counter
            like_atoms = 0
            # Check each edge atom against my atoms
            for atom in edge.atoms:
                # If an atom is a match, increment the counter
                if atom in myAtoms:
                    like_atoms += 1
            # If the number of like atoms is less than 3 (i.e. not a match) create an edge object
            if like_atoms < 3:
                vertex.edges.append(Edge([myAtoms], vertex))
            # Else add the edge to the vertex's edge list and the vert to the edges list
            else:
                vertex.edges.append(edge)
                edge.verts.append(vertex)
    # Find the ends to the edges
    for edge in vertex.edges:
        # If the edge already has two vertices skip it
        if len(edge.verts) >= 2:
            continue
        # Calculate the vertex
        vn = calc_edge(edge, net)

        print(vn, vertex)
        if vn is None:
            edge.verts.append(None)
            continue
        # Spread to other sites or just find this one edge
        if recurse:
            find_edges(vn, net)
        else:
            return edge


########################################################################################################################


# Build Network function. Takes in a system, runs as a shell for the recursive next_site function and returns a Network
def build_network(mySys):
    # Find the first vertex
    v0 = find_v0(mySys.net)
    # Initiate the recursive network finding algorithm on the network and the first vertex
    find_edges(v0, mySys.net)

    # Return the completed network
    return mySys.net
