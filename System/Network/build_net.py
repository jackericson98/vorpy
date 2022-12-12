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
    :param: net - Network object
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
        dub_ndx = net.verts.index(dub)
        net.verts.insert(dub_ndx + 1, dub.doublet)
        net.vert_ndxs.insert(dub_ndx + 1, dub.ndx)

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's atoms
        con_verts = []
        for vert in net.verts:
            try:
                if len([0 for _ in dub.ndx if _ in vert.ndx]) == 3:
                    con_verts.append(vert)
            except AttributeError:
                pass
        # Divide the connecting outer vertices between the two doublet vertices
        dub_verts, dub_dub_verts = [], []
        for vert in con_verts:
            try:
                if calc_dist(vert.loc, dub.loc) < calc_dist(vert.loc, dub.doublet.loc):
                    dub_verts.append(vert)
                else:
                    dub_dub_verts.append(vert)
            except AttributeError:
                pass
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
                           verts=[dub.doublet, vert], surfs=[])
            net.edges.append(my_edge)
            dub.doublet.edges.append(my_edge)
            vert.edges.append(my_edge)

        ########################################## Create the inner edges ##########################################


        potential_edges = [[dub.ndx[i], dub.ndx[(i + 1) % 4], dub.ndx[(i + 2) % 4]] for i in range(4)]
        for ndx in potential_edges:
            ndx.sort()
        known_edges_ndxs = [edge.ndx for edge in dub.edges + dub.doublet.edges]

        # Gather the other combinations of atoms and create the remaining inner atoms
        edges = [Edge(atoms=[net.atoms[_] for _ in ndx], net=net, verts=[dub, dub.doublet], doublet=True, surfs=[])
                 for ndx in potential_edges if ndx not in known_edges_ndxs]

        # Add the edges to the network and the doublet vertices
        net.edges += edges
        dub.edges += edges
        dub.doublet.edges += edges


def build(net):
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
        # Check every combination of vert atoms as a potential edge
        for i in range(4):
            # Grab the atoms
            atoms = [vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]]
            # Get the atoms indices
            atom_ndxs = [net.atoms.index(_) for _ in atoms]
            atom_ndxs.sort()
            # Get the edge
            edge_ndx = ndx_search(net.edge_ndxs, atom_ndxs)
            # If the edge has been found before, continue
            if len(net.edge_ndxs) > edge_ndx and atom_ndxs == net.edge_ndxs[edge_ndx]:
                continue
            # Create the vertices list
            verts = []
            # Find the possible verts (the original vert and the new vert)
            for vert2 in net.verts:
                if vert2.doublet is None and len([0 for _ in atom_ndxs if _ in vert2.ndx]) == 3:
                    verts.append(vert2)
            # If the number of valid vertices for the edge is 1
            if len(verts) == 1:
                continue
            # Create the edge
            my_edge = Edge(atoms=atoms, net=net, verts=verts, surfs=[])
            # Add the edge to the System
            net.edges.insert(edge_ndx, my_edge)
            net.edge_ndxs.insert(edge_ndx, my_edge.ndx)
            # Edd the edge to it's objects
            for atom in atoms:
                atom.edges.append(my_edge)
            for vert in verts:
                vert.edges.append(my_edge)


    ################################################### Create the surfaces ############################################

    # Go through the edges in the network
    for edge1 in net.edges:
        # Go through the edge's atoms combinations
        for i in range(3):
            # Get the atoms and their sorted list of ndxs
            atoms = [edge1.atoms[i], edge1.atoms[(i + 1) % 3]]
            atom_ndxs = [net.atoms.index(atom) for atom in atoms]
            atom_ndxs.sort()
            # If the surface has been found before continue
            surf_ndx = ndx_search(net.surf_ndxs, atom_ndxs)
            # If the edge has been found before, continue
            if len(net.surf_ndxs) > surf_ndx and atom_ndxs == net.surf_ndxs[surf_ndx]:
                continue
            # Put together a list of edges that have our atoms
            edges = []
            # Go through the edges in the system
            for edge2 in net.edges:
                # If the surface's atoms are in the edge add it
                if len([0 for ndx in atom_ndxs if ndx in edge2.ndx]) == 2:
                    edges.append(edge2)
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in net.verts:
                # If the surface's atoms are shared with the vertex, add it to the list
                if len([0 for ndx in atom_ndxs if ndx in vert2.ndx]) == 2:
                    verts.append(vert2)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):
                # Create the surface
                my_surf = Surface(atoms=atoms, verts=verts, net=net, edges=edges)
                net.surfs.insert(surf_ndx, my_surf)
                net.surf_ndxs.insert(surf_ndx, atom_ndxs)
                # Add the surface to its objects
                for atom in atoms:
                    atom.surfs.append(my_surf)
                for vert in verts:
                    vert.surfs.append(my_surf)
                for edge in edges:
                    if edge.surfs is None:
                        edge.surfs = []
                    edge.surfs.append(my_surf)
