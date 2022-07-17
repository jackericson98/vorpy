"""Imports"""
import numpy as np
from System.system import Edge, Vertex
from calculators import *


########################################################################################################################
"""Gather information functions"""


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
        new_rad = calc_circ([a0, a1, atom])
        if new_rad:
            new_rad = new_rad[0][1]
        else:
            continue
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
    a1 = neighbors0[0]
    # Find the smallest circle you can make with a0, a1 and a third atom
    a2 = find_circle(a0, a1, net)
    # Find the first site by choosing the smallest interstitial sphere that can be made with the 3 atoms and the 50
    # closest atoms
    neighbors2 = sortbyDist([a0, a1, a2], net)
    r = np.inf
    myVert, my_an = None, None
    # Go through the closest atoms to the triplet and find the smallest vertex that can be made with a neighbor
    for an in neighbors2[:15]:
        if an.loc == a0.loc or an.loc == a1.loc or an.loc == a2.loc:
            continue
        # Calculate the vertex and check if None
        vert = calc_vert([a0, a1, a2, an])
        # Don't worry about None vert types
        if vert is None:
            continue
        # Check if the vn has a smaller radius than the current smallest radius
        if vert.rad < r:
            r = vert.rad
            myVert = vert
            my_an = an
    # Add connections to the network
    net.verts.append(myVert)
    return myVert


########################################################################################################################
"""Recursive network finding function"""


# Find edges function. Recursively traces out the network and records vertex locations,
def find_edges(vertex, net):
    # Create 4 edges with the 4 combinations of atoms that can be created
    for i in range(4):
        myEdge = Edge([vertex.atoms[i], vertex.atoms[(i + 1) % 4], vertex.atoms[(i + 2) % 4]], [vertex])
        etest = check_edge(myEdge, net)
        if etest:
            vertex.edges.append(etest)
        else:
            vertex.edges.append(myEdge)
    # Create the edge objects or grab them from the network and connect them
    for edge in vertex.edges:
        # Check to see if the edge exists. If it does move on the next edge in the vertex
        net_edge = check_edge(edge, net)
        if net_edge:
            continue
        # Create a vertex
        vn = calc_edge(edge, net)
        # If the vertex is None give the edge a None vertex and continue to the next edge
        if vn is None:
            edge.verts.append(Vertex([np.inf, np.inf, np.inf], np.inf))
            net.edges.append(edge)
            continue
        # Check the vertex to see if it exists in the network
        net_vert = check_vert(vn, net)
        # If it does, add the vertex to the edge and the edge to the network
        if net_vert:
            edge.verts.append(net_vert)
            net.edges.append(edge)
        # If both the edge and the vertex do not exist in the network, we have a true new site
        else:
            edge.verts.append(vn)
            net.edges.append(edge)
            net.verts.append(vn)
            find_edges(vn, net)
    return


########################################################################################################################


# Build Network function. Takes in a System, runs as a shell for the recursive next_site function and returns a Network
def build_network(mySys):
    # Find the first vertex
    v0 = find_v0(mySys.net)
    # Initiate the recursive network finding algorithm on the network and the first vertex
    find_edges(v0, mySys.net)
    return mySys
