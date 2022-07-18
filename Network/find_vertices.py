from System.system import Edge, Atom, Vertex
from Network.calculators import *
from Presentation.Visualize.visualize import plot_atoms, plot_verts
import matplotlib.pyplot as plt

########################################################################################################################
"""Finding functions"""


# Doublet making function. Takes in a doublet and adds the two vertices to their respective places.
def doublet(verts, net):
    return verts[0]


# Find v0 function. Finds the first vertex in the network
def find_v0(sys):
    # Find the center of mass of the atoms
    com = calc_com(sys.atoms)
    # First choose an appropriate initial atom based of com proximity
    min_dist = np.inf
    a0 = None
    # Go through each atom determining if it is closer to the com
    for atom in sys.atoms:
        # Set the new com distance
        com_dist = calc_dist(atom.loc, com)
        # If is less than the current closest atom's distance to the center of mass update the variables
        if com_dist < min_dist:
            min_dist = com_dist
            a0 = atom
    # Find the set of atoms with the minimum distance between surfaces
    min_dist = np.inf
    a1 = None
    # Go through each atom determining the atom with the minimum distance between it and a0's surfaces
    for atom in sys.atoms:
        # Skip a0
        if atom == a0:
            continue
        # Set the new atom distances
        a_dist = calc_dist(a0.loc, atom.loc) - (a0.rad + atom.rad)
        # If the new atom distance is less than the previous minimum distance update the variables
        if a_dist < min_dist:
            min_dist = a_dist
            a1 = atom
    # Find the set of atoms with the minimum inscribed circle
    min_rad = np.inf
    a2 = None
    # Go through each other atom to determine the smallest circle that can be made with our 2 atoms and a third
    for atom in sys.atoms:
        # Skip a0, a1
        if atom == a0 or atom == a1:
            continue
        # Calculate the circle made with the 3 atoms
        circ = calc_circ([a0, a1, atom])
        # If the radius of the circle that is made with the 3 atoms is smaller than the previous smallest radius replace
        if circ and abs(circ[0][1]) < min_rad:
            min_rad = abs(circ[0][1])
            a2 = atom
    # Find the set of atoms with the minimum inscribed sphere
    min_rad = np.inf
    myVert = None
    # Go through each other atom to determine the smallest possible inscribed sphere
    for atom in sys.atoms:
        # Skip a0, a1, a2
        if atom == a0 or atom == a1 or atom == a2:
            continue
        # Get the vertex made from the atoms
        vert = calc_vert([a0, a1, a2] + [atom])
        # If the radius of the inscribed
        if vert and vert.rad < min_rad:
            min_rad = vert.rad
            myVert = vert
    # Return the vertex
    return myVert


# Find site function. When given an edge, this function returns the closest next vertex in the System
def find_site(edge, net):
    # Get the edges location and radius
    circs = calc_circ(edge.atoms)
    if circs is None:
        return
    edge.loc, edge.rad = circs[0][0], abs(circs[0][1])
    # Get the edge's direction
    edge.dir = calc_dir(edge)
    # Grab the edge's previous vertex
    vn_1 = edge.verts[0]
    min_val = np.inf
    vn = None
    # Go through each atom in the System
    for atom in net.atoms:
        # Check for edge atoms and backwards atom by making sure atom is not in vn_1's list of atoms
        if {atom}.issubset(vn_1.atoms):
            continue
        # Set up an atoms list
        atoms = [atom] + edge.atoms
        # Calculate the vertex
        vert = calc_vert(atoms)
        if vert is None:
            continue
        # Find the relative distance between the verts along the edge
        dist = calc_rel_dist(vn_1, vert, edge)
        if dist < min_val:
            min_val = dist
            vn = vert
    return vn


# Find site function (#2)> Takes in an edge and finds the smallest vertex attached to it
def find_site1(edge, net):
    myVert = None
    # Loop through the atoms to see if they over
    for atom in net.atoms:
        # This filters out any of the atoms in the edge or the remaining atom from the previous vertex
        if {atom}.issubset(edge.verts[0].atoms):
            continue
        # Calculate the vertex with atom
        vert = calc_vert(edge.atoms + [atom])
        if vert is None:
            continue
        # Check if the vertex overlaps with any of the networks atoms
        overlap = False
        for a_test in net.atoms:
            if {a_test}.issubset(edge.atoms + [atom]):
                continue
            if round(calc_dist(a_test.loc, vert.loc) - (a_test.rad + vert.rad), 7) < 0:
                overlap = True
                break
        if not overlap:
            print("vert found")
            myVert = vert
            break
    return myVert


# Find network function. Keeps searching the network until all verts are found
def find_network(sys):
    # Find the first vertex in the System
    v0 = find_v0(sys)
    # Add v0 to the System
    sys.net.verts.append(v0)
    # Set up the vertex stack
    vert_stack = [v0]
    # While the verts stack is not empty
    while vert_stack:
        # Get the vertex from the top of the stack
        vert = vert_stack.pop()
        # Set up the edge stack
        e_stack = [Edge([vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]], [vert]) for i in range(4)]
        # While the edge stack is not empty
        while e_stack:
            # Get the edge from the top of the stack
            edge = e_stack.pop()
            # Find the next site in the network
            myVert = find_site(edge, sys.net)
            # If the vertex is none continue
            if myVert is None:
                print("None vert")
                continue
            # If the vertex exists in the network add the vertex to the edge and move on to the next edge in the stack
            found_vert = check_vert(set(myVert.atoms), sys.net.verts)
            if not found_vert:
                vert_stack.append(myVert)
                sys.net.verts.append(myVert)
