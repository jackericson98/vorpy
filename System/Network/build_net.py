from System.calcs import *
# from System.atom import Atom
# from System.Network.vertex import Vertex
from System.Network.edge import Edge
from System.Network.surface import Surface
# from Visualize.visualize import *


############################################## Doublets ################################################################


def doublify(net):
    """
    Finds all doublet edges throughout the network and adds them
    :return:
    """
    # Find the doublets, separate them and connect their edges
    net.doublets = []
    for vert in net.verts:
        # If we come across a doublet add it to the list
        if vert.doublet is not None:
            net.doublets.append(vert)

    # Go through the doublets
    for dub in net.doublets:

        # Add the doublet to the network
        net.verts.insert(net.verts.index(dub) + 1, dub.doublet)

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's atoms
        con_verts = []
        for vert in net.verts:
            if len([0 for _ in dub.ndx if _ in vert.ndx]) == 3:
                con_verts.append(vert)
        # Divide the connecting outer vertices between the two doublet vertices
        dub_verts, dub_dub_verts = [], []
        for vert in con_verts:
            if calc_dist(vert.loc, dub.loc) < calc_dist(vert.loc, dub.doublet.loc):
                dub_verts.append(vert)
            else:
                dub_dub_verts.append(vert)

        # Create the edge objects for each of the vertices connected to the primary doublet vertex
        dub.edges = []
        for vert in dub_verts:
            # Create the edge from the atoms in both dub and vert and add it to the network and each vertex
            my_edge = Edge([net.atoms[_] for _ in [_ for _ in vert.ndx if _ in dub.ndx]], net=net, verts=[dub, vert])
            net.edges.append(my_edge)
            dub.edges.append(my_edge)
            vert.edges.append(my_edge)

        # Create the edge objects for each of the vertices connected to the secondary doublet vertex
        dub.doublet.edges = []
        for vert in dub_dub_verts:
            # Create the edge from the atoms in both dub.doublet and vert and add it to the network and each vertex
            my_edge = Edge(atoms=[net.atoms[_] for _ in [_ for _ in vert.ndx if _ in dub.doublet.ndx]], net=net,
                           verts=[dub.doublet, vert])
            net.edges.append(my_edge)
            dub.doublet.edges.append(my_edge)
            vert.edges.append(my_edge)

        ########################################## Create the inner edges ##########################################

        potential_edges = [[dub.ndx[i], dub.ndx[(i + 1) % 4], dub.ndx[(i + 2) % 4]] for i in range(4)]
        for ndx in potential_edges:
            ndx.sort()
        known_edges_ndxs = [edge.ndx for edge in dub.edges + dub.doublet.edges]

        # Gather the other combinations of atoms and create the remaining inner atoms
        edges = [Edge([net.atoms[_] for _ in ndx], net=net, verts=[dub, dub.doublet], doublet=True)
                 for ndx in potential_edges if ndx not in known_edges_ndxs]

        # Add the edges to the network and the doublet vertices
        net.edges += edges
        dub.edges += edges
        dub.doublet.edges += edges


def make_objects(net):
    """
    Checks the atoms of the vertices for patterns and creates edges and surfaces
    :param net: Network object to pull information from
    """
    # Reset the network's list of edges and surfaces for a clean slate
    net.surfs = []
    net.edges = []

    # Fill in the doublets and set their outer edges
    doublify(net)

    ################################################# Create the edges #################################################

    # Go through the vertices in the network searching for potential edges
    for vert1 in net.verts:
        if vert1.doublet is not None:
            continue
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
                if vert2.doublet is None and len([0 for _ in atoms if _ in vert2.atoms]) == 3:
                    verts.append(vert2)
            # If the number of valid vertices for the edge is 1
            if len(verts) == 1:
                continue
            # Create the edge
            my_edge = Edge(atoms=atoms, net=net, verts=verts)
            # Add the edge to the System
            my_edges.append(my_edge)
        # Add the edges to the network's list of edges
        net.edges += my_edges

    ################################################### Create the surfaces ############################################

    # Go through the edges in the network
    for edge1 in net.edges:
        # Go through the edge's atoms combinations
        for i in range(3):
            atoms = [edge1.atoms[i], edge1.atoms[(i + 1) % 3]]
            atom_ndxs = [net.atoms.index(atom) for atom in atoms]
            atom_ndxs.sort()
            # If the surface has been found before continue
            if check_surf(atoms, net.surfs):
                continue
            # Put together a list of edges that have our atoms
            edges = []
            # Go through the edges in the system
            for edge2 in net.edges:
                # If the surface's atoms are in the edge add it
                if set(atoms).issubset(edge2.atoms):
                    edges.append(edge2)
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in net.verts:
                # If the surface's atoms are shared with the vertex, add it to the list
                if set(atoms).issubset(vert2.atoms):
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
