import time

from System.sys_funcs.calcs.calcs import *
from System.Network.net_objs.edge import Edge
from System.Network.net_objs.surface import Surface


############################################## Doublets ################################################################


def connect_vert(net, vert):
    con_verts = []
    for vert1 in net.verts:
        if len([0 for _ in vert.ndx if _ in vert1.ndx]) == 3:
            con_verts.append(vert1)
    return con_verts


def doublify(net, get_edges=True):
    """
    Finds all doublet edges throughout the network and adds them
    :param: net - Network object
    :return:
    """
    # Find the doublets, separate them and connect their edges
    doublets = []
    net.dub_ndxs = []
    for vert in net.verts:
        # If we come across a doublet add it to the list
        if vert.doublet is not None:
            doublets.append(vert)
    # Go through the doublets
    for dub in doublets:
        # Add the doublet to the network
        dub_ndx = ndx_search(net.vert_ndxs, dub.ndx)
        net.verts.insert(dub_ndx + 1, dub.doublet)
        net.vert_ndxs.insert(dub_ndx + 1, dub.ndx)
        net.dub_ndxs += [dub_ndx, dub_ndx + 1]

        if not get_edges:
            return

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's atoms
        con_verts = connect_vert(net, dub)
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
        edges = [Edge(atoms=[net.atoms[_] for _ in ndx], net=net, verts=[dub, dub.doublet], surfs=[])
                 for ndx in potential_edges if ndx not in known_edges_ndxs]

        # Add the edges to the network and the doublet vertices
        for edge in edges:
            edge_ndx = ndx_search(net.edge_ndxs, edge.ndx)
            net.edges.insert(edge_ndx, edge)
            net.edge_ndxs.insert(edge_ndx, edge.ndx)

        dub.edges += edges
        dub.doublet.edges += edges


def build(net, my_group, from_scratch=True):
    """
    Checks the atoms of the vertices for patterns and creates edges and surfaces
    :param my_group:
    :param from_scratch:
    :param net: Network object to pull information from
    """
    # Reset the network's list of edges and surfaces for a clean slate
    if from_scratch or (net.edges is None or len(net.edges) == 0):
        net.edges = []
    if from_scratch or (net.surfs is None or len(net.surfs) == 0):
        net.surfs = []

    # Fill in the doublets and set their outer edges
    doublify(net)
    net.sort_verts(my_group)
    # Add the vertices to the atoms
    for vert in net.verts:
        for atom in vert.atoms:
            atom.verts.append(vert)

    ################################################# Create the edges #################################################

    # Go through the vertices in the network searching for potential edges
    for i, vert1 in enumerate(net.verts):
        my_time = time.perf_counter() - net.my_time
        h, m, s = get_time(my_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * (0.5 * len(net.edges)) / (3/2 * len(net.verts)))), end="")
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
                check_ndx = set(vert2.ndx)
                if vert2.doublet is None and len([0 for _ in atom_ndxs if _ in check_ndx]) == 3 and \
                        vert2.ndx not in [_.ndx for _ in verts]:
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
    for i, edge1 in enumerate(net.edges):
        my_time = time.perf_counter() - net.my_time
        h, m, s = get_time(my_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * (len(net.surfs) + 0.5 * (len(net.edges))) / ((3/2) * len(net.verts)))), end="")
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
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in net.verts:
                # If the surface's atoms are shared with the vertex, add it to the list
                if len([0 for ndx in atom_ndxs if ndx in vert2.ndx]) == 2:
                    verts.append(vert2)
            # Put together a list of edges that have our atoms
            edges = []
            # Go through the edges in the system
            for vert in verts:
                for edge2 in vert.edges:
                    if edge2 in edges:
                        continue
                    # If the surface's atom s are in the edge add it
                    if len([0 for ndx in atom_ndxs if ndx in edge2.ndx]) == 2:
                        edges.append(edge2)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):
                no_surf = False
                # Check to see if the surface is worth adding
                for vert in verts:
                    if len(vert.edges) <= 2:
                        no_surf = True
                if no_surf:
                    continue
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


def filter_edges(net):
    # Set up the edge list
    bad_edges = []
    # Go through the edges in the network
    for i, edge in enumerate(net.edges):
        # Check the number of surfaces and vertices
        if len(edge.surfs) <= 1 or len(edge.verts) <= 1:
            bad_edges.append(i)

    # Carefully remove the edges from the network
    for i in reversed(bad_edges):
        # Get the edge from the network
        edge = net.edges[i]
        # Remove the edge from the edges vertices
        for vert in edge.verts:
            vert.edges.remove(edge)
        # Remove the edge from the surfaces
        for surf in edge.surfs:
            surf.edges.remove(edge)
        # Remove the edge from the network
        net.edges.pop(i)
        net.edge_ndxs.pop(i)
