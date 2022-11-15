from System.calcs import *
# from System.atom import Atom
# from System.Network.vertex import Vertex
from System.Network.edge import Edge
from System.Network.surface import Surface
# from Visualize.visualize import *


############################################## Doublets ################################################################


def find_doublet_edges(net, vert):
    """
    Doublet function used to create the edges surrounding and inside of vertices
    :return:
    """
    # Grab the two vertices
    v0, v1, atoms = vert, vert.doublet, vert.atoms
    # Set up a variable for doublet edges
    edges = []
    # Find what type of doublet it is (i.e. the # of edges) by counting the number of "free" inscribed circles
    for i in range(4):
        # Get the current combination of atoms to test doubletness of the edge
        edge_atoms = [net.atoms.index(_) for _ in [atoms[i], atoms[(i + 1) % 4], atoms[(i + 2) % 4]]]
        edges.append(edge_atoms)
    # Add the connecting edges
    for edge_atoms in edges:
        # Create a vertices list
        verts = []
        # Go through the other vertices looking for a matching set of atoms
        for myNdx in net.vert_ndxs:
            # Check to see if the current check vertex is the doublet itself
            if myNdx == v0.ndx:
                continue
            # Count the number of vertices that the edge atoms share 3 atoms with
            num_shared_atoms = len([0 for _ in edge_atoms if _ in myNdx])
            # If the edge connects with 2 other vertices it is a doublet connecting edge
            if num_shared_atoms == 3:
                verts.append(net.verts[net.vert_ndxs.index(myNdx)])
        # If there are 2 vertices it is a connecting edge
        if len(verts) == 2:
            for vert in verts:
                # Pair the edges to their correct vertices
                d00, d01 = calc_dist(v0.loc, vert.loc), calc_dist(v1.loc, vert.loc)
                if d00 < d01:
                    net.edges.append(Edge(atoms=[net.atoms[_] for _ in edge_atoms], net=net, verts=[v0, vert]))
                else:
                    net.edges.append(Edge(atoms=[net.atoms[_] for _ in edge_atoms], net=net, verts=[v1, vert]))
        # If there are no vertices that contain all three edge atoms (other than the doublet) it is a pure doublet edge
        elif len(verts) == 0:
            myEdge = Edge(atoms=[net.atoms[_] for _ in edge_atoms], verts=[v0, v1], net=net, doublet=True)
            net.edges.append(myEdge)


def make_objects(net):
    """
    Checks the atoms of the vertices for patterns and creates edges and surfaces
    :param net: Network object to pull information from
    """
    # Reset the network's list of edges and surfaces for a clean slate
    net.edges = []
    net.surfs = []

    ################################################# Create the edges #################################################

    vert_dubs = []
    # Go through the vertices in the network searching for potential edges
    for vert1 in net.verts:
        if vert1.doublet in vert_dubs:
            continue

        if vert1.doublet is not None:
            find_doublet_edges(net, vert1)
            vert_dubs.append(vert1)
        else:
            # Set up the edges
            my_edges = []
            # Check every combination of vert atoms as a potential edge
            for i in range(4):
                # Grab the atoms
                atoms = [vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]]
                # If the edge has been found before, continue
                if check_edge(atoms, net.edges):
                    continue
                verts = []
                # Find the possible verts (the original vert and the new vert)
                for vert2 in net.verts:
                    if len([0 for _ in atoms if _ in vert2.atoms]) == 3:
                        verts.append(vert2)
                # If the number of valid vertices for the edge is 1
                if len(verts) == 1:
                    continue
                # Create the edge
                my_edge = Edge(list(atoms), net, verts)
                # Add the edge to the System
                my_edges.append(my_edge)
            # Add the edges to the network's list of edges
            net.edges += my_edges

    ################################################### Create the surfaces ############################################

    # Go through the edges in the network
    for edge1 in net.edges:
        # Go through the edge's atoms combinations
        for i in range(3):
            atoms = {edge1.atoms[i], edge1.atoms[(i + 1) % 3]}
            # If the surface has been found before continue
            if check_surf(atoms, net.surfs):
                continue
            # Put together a list of edges that have our atoms
            edges = []
            # Go through the edges in the system
            for edge2 in net.edges:
                # If the surface's atoms are in the edge add it
                if atoms.issubset(edge2.atoms):
                    edges.append(edge2)
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in net.verts:
                # If the surface's atoms are shared with the vertex, add it to the list
                if atoms.issubset(vert2.atoms):
                    verts.append(vert2)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):
                my_surf = Surface(list(atoms), verts=verts, net=net, edges=edges)
                net.surfs.append(my_surf)


# Connect network function.
def connect(net):
    """
    Takes in a disconnected network of atoms, vertices, surfaces and edges and connects it
    :param net: network to connect
    """

    ################################################# Connect the atoms ################################################

    # Go through the atoms in the network adding vertices, edges and surfaces
    for atom in net.atoms:

        # Reset the atom's vert list
        atom.verts = []
        # Go through the verts in the network
        for vert in net.verts:
            # If the atom is in the vertices atoms add the vertex to the atom's list of vertices
            if {atom}.issubset(vert.atoms):
                atom.verts.append(vert)

        # Reset the atom's edge list
        atom.edges = []
        # Go through the edges in the network
        for edge in net.edges:
            # If the atom is in the edge's list of atoms add the edge to the atoms list of edges
            if {atom}.issubset(edge.atoms):
                atom.edges.append(edge)

        # Reset the atom's surf list
        atom.surfs = []
        # Go through the surfs in the network
        for surf in net.surfs:
            # If the atom is in the surfs list of atoms add the surf to the atoms list of surfs
            if {atom}.issubset(surf.atoms):
                atom.surfs.append(surf)

    ############################################# Connect the vertices #################################################

    # Go through the vertices in the network adding the edges and surfaces that share atoms
    for vert in net.verts:

        # Reset the vertexes edge list
        vert.edges = []
        # Go through the edges in the network
        for edge in net.edges:
            # If the edges atoms are in the vertices atoms add it to the vertex
            if set(edge.atoms).issubset(vert.atoms):
                vert.edges.append(edge)

        # Reset the vertexes surf list
        vert.surfs = []
        # Go through the surfaces in the network
        for surf in net.surfs:
            # If the surfaces atoms are in the vertexes atoms add it to the vertex
            if set(surf.atoms).issubset(vert.atoms):
                vert.surfs.append(surf)

    ########################################## Connect the Edges #######################################################

    # Add the surfs to the edges
    for edge in net.edges:
        # Reset the edges surf list
        edge.surfs = []
        # Go through the surfaces in the network
        for surf in net.surfs:
            # If the surfaces atoms are in the edges atoms add it to the edge
            if set(surf.atoms).issubset(edge.atoms):
                edge.surfs.append(surf)


def build(net):
    """
    Takes in a broken network of vertices and spits out a connected network

    :param net: network to build
    """
    # Make the surface and edge objects
    make_objects(net)
    # Connect the network
    connect(net)
