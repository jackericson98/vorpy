import time
from System.sys_funcs.calcs.calcs import *


############################################## Doublets ################################################################


def doublify(vatoms, vlocs, vdubs):
    """
    Finds all doublet edges throughout the network and adds them
    :param: net - Network object
    :return:
    """
    eatoms, everts, esurfs = [], [], []
    # Go through the doublets
    for i in range(len(vatoms)):

        # Skip the verts that aren't doublets
        if i >= len(vatoms) - 1 or vdubs[i + 1] != 1:
            continue

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's atoms
        con_verts = []
        for j in range(len(vatoms)):
            if len([0 for _ in vatoms[i] if _ in vatoms[j]]) == 3:
                con_verts.append(j)

        # Divide the connecting outer vertices between the two doublet vertices
        dub_verts, dub_dub_verts = [], []
        for j in con_verts:
            # Decide between the two sides of the doublet for the outer vertex
            if calc_dist(np.array(vlocs[j]), np.array(vlocs[i])) < calc_dist(np.array(vlocs[j]), np.array(vlocs[i + 1])):
                dub_verts.append(j)
            else:
                dub_dub_verts.append(j)

        known_edges = []
        # Create the edge objects for each of the vertices connected to the primary doublet vertex
        for j in dub_verts:
            # Create the edge from the atoms in both dub and vert and add it to the network and each vertex
            edge_atoms = [_ for _ in vatoms[i] if _ in vatoms[j]]
            edge_ndx = ndx_search(eatoms, edge_atoms)
            eatoms.insert(edge_ndx, edge_atoms)
            everts.insert(edge_ndx, [i, j])
            esurfs.insert(edge_ndx, [])
            known_edges.append(edge_atoms)

        # Create the edge objects for each of the vertices connected to the secondary doublet vertex
        for j in dub_dub_verts:
            # Create the edge from the atoms in both dub.doublet and vert and add it to the network and each vertex
            edge_atoms = [_ for _ in vatoms[i] if _ in vatoms[j]]
            edge_ndx = ndx_search(eatoms, edge_atoms)
            eatoms.insert(edge_ndx, edge_atoms)
            everts.insert(edge_ndx, [i + 1, j])
            esurfs.insert(edge_ndx, [])
            known_edges.append(edge_atoms)

        ########################################## Create the inner edges ##########################################

        # Create a list of every edge possibility
        potential_edges = [[vatoms[i][k], vatoms[i][(k + 1) % 4], vatoms[i][(k + 2) % 4]] for k in range(4)]
        for ndx in potential_edges:
            ndx.sort()

        # Gather the other combinations of atoms and create the remaining inner atoms
        inner_edges = [ndx for ndx in potential_edges if ndx not in known_edges]

        # Add the edges to the network and the doublet vertices
        for edge in inner_edges:
            edge_ndx = ndx_search(eatoms, edge)
            eatoms.insert(edge_ndx, edge)
            everts.insert(edge_ndx, [i, i + 1])
    # Return the partial lists
    return eatoms, everts


def get_build_edges(averts, vatoms, vlocs, vdubs, start_time):
    # Get the doublet edges
    eatoms, everts = doublify(vatoms, vlocs, vdubs)

    # Go through the vertices in the network searching for potential edges
    for i, vert1 in enumerate(vatoms):
        # Print the time and process
        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * (0.5 * len(eatoms)) / (3 / 2 * len(vatoms)))),
              end="")

        # If the vertex is a doublet it has its edges already, so skip
        if vdubs[i] == 1 or (i + 1 < len(vdubs) and vdubs[i + 1] == 1):
            continue

        # Go through the atoms in the vertex looking for shared atoms
        for atom in vert1:
            # Go through the vertices in each atom
            for j in averts[atom]:
                # Get the atoms for vert2
                vert2 = vatoms[j]
                # Check the number of shared atoms between vert1 and vert2
                shared_atoms = [_ for _ in vert1 if _ in vert2]
                # Check if this edge is real
                if len(shared_atoms) == 3:
                    # Get the index of the edge in the edge list
                    edge_ndx = ndx_search(eatoms, shared_atoms)
                    # Check if we have found this edge before
                    if edge_ndx >= len(eatoms) or eatoms[edge_ndx] != shared_atoms:
                        # Add the edges atoms and the edges verts to their respective lists
                        eatoms.insert(edge_ndx, shared_atoms)
                        everts.insert(edge_ndx, [i, j])
    # Return the edge's atoms and verts
    return eatoms, everts


def add_build_edges(num_atoms, eatoms, num_verts, everts):
    # Create the empty atom list of edges
    aedges = [[i for i in range(0)] for _ in range(num_atoms)]
    # Create the empty vertex list of edges
    vedges = [[i for i in range(0)] for _ in range(num_verts)]
    # Go through the edges in the network
    for i, edge_atoms in enumerate(eatoms):
        # Go through the atoms in the edge
        for j in edge_atoms:
            # Add the edge to each atom
            aedges[j].append(i)
        # Get the edges vertices
        edge_verts = everts[i]
        # Go through the verts in the edge
        for j in edge_verts:
            vedges[j].append(i)
    # Return the newly filled in lists
    return aedges, vedges


def get_build_surfs1(vatoms, vedges, eatoms, start_time):
    # Set up the surface lists
    satoms, sverts, sedges = [], [], []

    # Go through the edges in the network
    for i, edge1 in enumerate(eatoms):

        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2),
                      min(100.0, 100 * (len(satoms) + 0.5 * (len(eatoms))) / ((3 / 2) * len(vatoms)))), end="")

        # Go through the edge's atoms combinations
        for j in range(3):
            # Get the atoms and their sorted list of ndxs
            atoms = [edge1[j], edge1[(j + 1) % 3]]
            atom_ndxs = atoms[:]
            atom_ndxs.sort()
            # If the surface has been found before continue
            surf_ndx = ndx_search(satoms, atom_ndxs)
            # If the edge has been found before, continue
            if len(satoms) > surf_ndx and atom_ndxs == satoms[surf_ndx]:
                continue

            # Limit the list of verts to possible vertices
            # max_vert_ndx =

            # Put together a list of verts that have our atoms
            verts = []
            for k, vert2 in enumerate(vatoms):
                # If the surface's atoms are shared with the vertex, add it to the list
                if len([0 for ndx in atom_ndxs if ndx in vert2]) == 2:
                    verts.append(k)

            # Put together a list of edges that have our atoms
            edges = []
            # Go through the edges in the system
            for k, edge2 in enumerate(eatoms):
                # If the surface's atom s are in the edge add it
                if len([0 for ndx in atom_ndxs if ndx in edge2]) == 2:
                    edges.append(k)

            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):

                no_surf = False
                # Check to see if the surface is worth adding
                for vert_ndx in verts:
                    if len(vedges[vert_ndx]) <= 2:
                        no_surf = True
                if no_surf:
                    continue
                incomplete = False
                for vert in verts:
                    if len(vedges[vert]) > 3:
                        incomplete = True
                if incomplete:
                    continue

                satoms.insert(surf_ndx, atom_ndxs)
                sedges.insert(surf_ndx, edges)
                sverts.insert(surf_ndx, verts)
    return satoms, sverts, sedges


def get_build_surfs(averts, aedges, vatoms, vedges, eatoms, start_time):
    # Set up the surface lists
    satoms, sverts, sedges = [], [], []

    # Go through the edges in the network
    for i, edge1 in enumerate(eatoms):

        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2),
                      min(100.0, 100 * (len(satoms) + 0.5 * (len(eatoms))) / ((3 / 2) * len(vatoms)))), end="")
        # Get the possible surfs from the edge's atoms
        test_surfs = [edge1[:2], edge1[1:], edge1[::2]]
        # Go through each possible surface for the edge
        for test_surf in test_surfs:
            # If the surface has been found before continue
            surf_ndx = ndx_search(satoms, test_surf)
            # If the surface has been found before, continue
            if len(satoms) > surf_ndx and test_surf == satoms[surf_ndx]:
                continue
            # Set up the surf edges and surf verts lists
            surf_edges, surf_verts = [], []
            # Go through the atoms in the surface looking for edge candidates
            for atom in test_surf:
                # Get the atom's edges
                for edge in aedges[atom]:
                    # Get the edges atoms
                    edge2 = eatoms[edge]
                    # If the number of shared atoms is 2 add the edge to the test surf's list of edges
                    if len([_ for _ in edge2 if _ in test_surf]) == 2 and edge not in surf_edges:
                        # Add the edge
                        surf_edges.append(edge)
                # Get the atom's vertices
                for vert in averts[atom]:
                    # Get the vertices atoms
                    vert2 = vatoms[vert]
                    # If the number of shared atoms is 2 add the edge to the test surf's list of edges
                    if len([_ for _ in vert2 if _ in test_surf]) == 2 and vert not in surf_verts:
                        # Add the edge
                        surf_verts.append(vert)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(surf_verts) == len(surf_edges):

                no_surf = False
                # Check to see if the surface is worth adding
                for vert_ndx in surf_verts:
                    if len(vedges[vert_ndx]) <= 2:
                        no_surf = True
                if no_surf:
                    continue

                satoms.insert(surf_ndx, test_surf)
                sedges.insert(surf_ndx, surf_edges)
                sverts.insert(surf_ndx, surf_verts)
    return satoms, sverts, sedges


def add_build_surfs(num_atoms, satoms, num_verts, sverts, num_edges, sedges):
    # Atoms
    asurfs = [[] for _ in range(num_atoms)]
    for i, surf_atoms in enumerate(satoms):
        for j in surf_atoms:
            asurfs[j] += [i]

    # Verts
    vsurfs = [[] for _ in range(num_verts)]
    for i, surf_verts in enumerate(sverts):
        for j in surf_verts:
            vsurfs[j] += [i]

    # Edges
    esurfs = [[] for _ in range(num_edges)]
    for i, surf_edges in enumerate(sedges):
        for j in surf_edges:
            esurfs[j] += [i]

    return asurfs, vsurfs, esurfs


def build(vatoms, vlocs, vdubs, num_atoms, my_time):
    """
    Checks the atoms of the vertices for patterns and creates edges and surfaces
    :param my_group:
    :param from_scratch:
    :param net: Network object to pull information from
    """

    # Create the lists
    averts, aedges = [[] for _ in range(num_atoms)], [[] for _ in range(num_atoms)]

    # Add the vertices to the atoms
    for i, ndx in enumerate(vatoms):
        for j in ndx:
            averts[j].append(i)

    ################################################# Create the edges #################################################

    # Fill in the doublets and set their outer edges
    eatoms, everts = get_build_edges(averts, vatoms, vlocs, vdubs, my_time)

    # Add the edges to their atoms and vertices
    aedges, vedges = add_build_edges(num_atoms, eatoms, len(vatoms), everts)

    ################################################### Create the surfaces ############################################

    # Get the surfaces
    satoms, sverts, sedges = get_build_surfs(averts, aedges, vatoms, vedges, eatoms, my_time)

    # Add the surface objects to their correct indices
    asurfs, vsurfs, esurfs = add_build_surfs(num_atoms, satoms, len(vatoms), sverts, len(eatoms), sedges)

    # Package the lists neatly for easier parsing
    atom_lists = {'averts': averts, 'aedges': aedges, 'asurfs': asurfs}
    vert_lists = {'vedges': vedges, 'vsurfs': vsurfs}
    edge_lists = {'eatoms': eatoms, 'everts': everts, 'esurfs': esurfs}
    surf_lists = {'satoms': satoms, 'sverts': sverts, 'sedges': sedges}

    # Return the dictionary lists
    return atom_lists, vert_lists, edge_lists, surf_lists
