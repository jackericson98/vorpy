from System.sys_funcs.calcs import *
# from System.atom import Atom
# from System.Network.vertex import Vertex
from System.Network.net_objs.edge import Edge
from System.Network.net_objs.surface import Surface
# from Visualize.visualize import *


############################################## Doublets ################################################################


def doublify(net, get_edges=True):
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

        if not get_edges:
            return

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's atoms
        con_verts = []
        for vert in net.verts:
            if len([0 for _ in dub.ndx if _ in vert.ndx]) == 3:
                con_verts.append(vert)
        # Divide the connecting outer vertices between the two doublet vertices
        dub_verts, dub_dub_verts = [], []
        for vert in con_verts:
            # Decide between the two sides of the doublet for the outer vertex
            if calc_dist(vert.loc, dub.loc) < calc_dist(vert.loc, dub.doublet.loc):
                dub_verts.append(vert)
            else:
                dub_dub_verts.append(vert)
        # Create the edge objects for each of the vertices connected to the primary doublet vertex
        dub.edges = []
        for vert in dub_verts:
            # Create the edge from the atoms in both dub and vert and add it to the network and each vertex
            my_edge = Edge([net.atoms[_] for _ in [_ for _ in vert.ndx if _ in dub.ndx]], net=net, verts=[dub, vert],
                           surfs=[])
            edge_ndx = ndx_search(net.edge_ndxs, my_edge.ndx)
            net.edges.insert(edge_ndx, my_edge)
            net.edge_ndxs.insert(edge_ndx, my_edge.ndx)
            dub.edges.append(my_edge)
            vert.edges.append(my_edge)

        # Create the edge objects for each of the vertices connected to the secondary doublet vertex
        dub.doublet.edges = []
        for vert in dub_dub_verts:
            # Create the edge from the atoms in both dub.doublet and vert and add it to the network and each vertex
            my_edge = Edge(atoms=[net.atoms[_] for _ in [_ for _ in vert.ndx if _ in dub.doublet.ndx]], net=net,
                           verts=[dub.doublet, vert], surfs=[])
            edge_ndx = ndx_search(net.edge_ndxs, my_edge.ndx)
            net.edges.insert(edge_ndx, my_edge)
            net.edge_ndxs.insert(edge_ndx, my_edge.ndx)
            dub.doublet.edges.append(my_edge)
            vert.edges.append(my_edge)

        ########################################## Create the inner edges ##########################################

        # Create a list of every edge possibility
        potential_edges = [[dub.ndx[i], dub.ndx[(i + 1) % 4], dub.ndx[(i + 2) % 4]] for i in range(4)]
        for ndx in potential_edges:
            ndx.sort()
        # Get the connecting edge
        known_edges_ndxs = [edge.ndx for edge in dub.edges + dub.doublet.edges]
        # Gather the other combinations of atoms and create the remaining inner atoms
        edges = [Edge(atoms=[net.atoms[_] for _ in ndx], net=net, verts=[dub, dub.doublet], doublet=True, surfs=[])
                 for ndx in potential_edges if ndx not in known_edges_ndxs]

        # Add the edges to the network and the doublet vertices
        for edge in edges:
            edge_ndx = ndx_search(net.edge_ndxs, edge.ndx)
            net.edges.insert(edge_ndx, edge)
            net.edge_ndxs.insert(edge_ndx, edge.ndx)

        dub.edges += edges
        dub.doublet.edges += edges


def build(net, get_edges=True, get_surfs=True):
    """
    Checks the atoms of the vertices for patterns and creates edges and surfaces
    :param get_surfs:
    :param get_edges:
    :param net: Network object to pull information from
    """
    # Reset the network's list of edges and surfaces for a clean slate
    if get_edges:
        net.edges = []
    if get_surfs:
        net.surfs = []

    # Fill in the doublets and set their outer edges
    doublify(net, get_edges=get_edges)
    # Add the vertices to the atoms
    for vert in net.verts:
        for atom in vert.atoms:
            atom.verts.append(vert)

    ################################################# Create the edges #################################################

    # Escape if only finding the vertices
    if not get_edges:
        return
    # Go through the vertices in the network searching for potential edges
    for i in range(len(net.verts)):
        print("\rConnecting Network: {:.2f} %".format(min(100.0, 100 * (len(net.edges)) / (2 * len(net.verts)))), end="")
        vert1 = net.verts[i]
        if vert1.doublet is not None:
            continue
        # Check every combination of vert atoms as a potential edge
        for j in range(4):
            # Grab the atoms
            atoms = [vert1.atoms[j], vert1.atoms[(j + 1) % 4], vert1.atoms[(j + 2) % 4]]
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
                if vert2.doublet is None and len([0 for _ in atom_ndxs if _ in vert2.ndx]) == 3 and vert2.ndx not in [_.ndx for _ in verts]:
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

    # Escape if only connecting the edges
    if not get_surfs:
        return

    # Go through the edges in the network
    for i in range(len(net.edges)):
        print("\rConnecting Network: {:.2f} %".format(min(100.0, 100 * (i + len(net.edges) - 1) / (2 * len(net.verts)))), end="")
        edge1 = net.edges[i]
        # Go through the edge's atoms combinations
        for j in range(3):
            # Get the atoms and their sorted list of ndxs
            atoms = [edge1.atoms[j], edge1.atoms[(j + 1) % 3]]
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
