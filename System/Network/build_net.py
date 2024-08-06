import time
import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist, get_time
from System.sys_funcs.calcs.sorting import ndx_search


############################################## Doublets ################################################################


def doublify(vballs, vlocs, vdubs):
    """
    Finds all doublet edges throughout the network and adds them
    :param: net - Network object
    :return:
    """
    eballs, everts, esurfs = [], [], []
    # Go through the doublets
    for i in range(len(vballs)):

        # Skip the verts that aren't doublets
        if i >= len(vballs) - 1 or vdubs[i + 1] != 1:
            continue

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's atoms
        con_verts = []
        for j in range(len(vballs)):
            if len([0 for _ in vballs[i] if _ in vballs[j]]) == 3:
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
            edge_balls = [_ for _ in vballs[i] if _ in vballs[j]]
            edge_ndx = ndx_search(eballs, edge_balls)
            eballs.insert(edge_ndx, edge_balls)
            everts.insert(edge_ndx, [i, j])
            esurfs.insert(edge_ndx, [])
            known_edges.append(edge_balls)

        # Create the edge objects for each of the vertices connected to the secondary doublet vertex
        for j in dub_dub_verts:
            # Create the edge from the atoms in both dub.doublet and vert and add it to the network and each vertex
            edge_balls = [_ for _ in vballs[i] if _ in vballs[j]]
            edge_ndx = ndx_search(eballs, edge_balls)
            eballs.insert(edge_ndx, edge_balls)
            everts.insert(edge_ndx, [i + 1, j])
            esurfs.insert(edge_ndx, [])
            known_edges.append(edge_balls)

        ########################################## Create the inner edges ##########################################

        # Create a list of every edge possibility
        potential_edges = [[vballs[i][k], vballs[i][(k + 1) % 4], vballs[i][(k + 2) % 4]] for k in range(4)]
        for ndx in potential_edges:
            ndx.sort()

        # Gather the other combinations of atoms and create the remaining inner atoms
        inner_edges = [ndx for ndx in potential_edges if ndx not in known_edges]

        # Add the edges to the network and the doublet vertices
        for edge in inner_edges:
            edge_ndx = ndx_search(eballs, edge)
            eballs.insert(edge_ndx, edge)
            everts.insert(edge_ndx, [i, i + 1])
    # Return the partial lists
    return eballs, everts


def get_build_edges(b_verts, v_balls, vlocs, vdubs, start_time):
    # Get the doublet edges
    e_balls, e_verts = doublify(v_balls, vlocs, vdubs)

    # Go through the vertices in the network searching for potential edges
    for i, vert1 in enumerate(v_balls):
        # Print the time and process
        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * (0.5 * len(e_balls)) / (3 / 2 * len(v_balls)))),
              end="")

        # If the vertex is a doublet it has its edges already, so skip
        if vdubs[i] == 1 or (i + 1 < len(vdubs) and vdubs[i + 1] == 1):
            continue

        # Go through the atoms in the vertex looking for shared atoms
        for ball in vert1:
            # Go through the vertices in each atom
            for j in b_verts[ball]:
                # Get the atoms for vert2
                vert2 = v_balls[j]
                # Check the number of shared atoms between vert1 and vert2
                shared_balls = [_ for _ in vert1 if _ in vert2]
                # Check if this edge is real
                if len(shared_balls) == 3:
                    # Get the index of the edge in the edge list
                    edge_ndx = ndx_search(e_balls, shared_balls)
                    # Check if we have found this edge before
                    if edge_ndx >= len(e_balls) or e_balls[edge_ndx] != shared_balls:
                        # Add the edges atoms and the edges verts to their respective lists
                        e_balls.insert(edge_ndx, shared_balls)
                        e_verts.insert(edge_ndx, [i, j])
    # Return the edge's atoms and verts
    return e_balls, e_verts


def add_build_edges(num_balls, e_balls, num_verts, e_verts):
    # Create the empty atom list of edges
    b_edges = [[i for i in range(0)] for _ in range(num_balls)]
    # Create the empty vertex list of edges
    v_edges = [[i for i in range(0)] for _ in range(num_verts)]
    # Go through the edges in the network
    for i, edge_balls in enumerate(e_balls):
        # Go through the atoms in the edge
        for j in edge_balls:
            # Add the edge to each atom
            b_edges[j].append(i)
        # Get the edges vertices
        edge_verts = e_verts[i]
        # Go through the verts in the edge
        for j in edge_verts:
            v_edges[j].append(i)
    # Return the newly filled in lists
    return b_edges, v_edges


def get_build_surfs1(v_balls, v_edges, e_balls, start_time):
    # Set up the surface lists
    s_balls, s_verts, s_edges = [], [], []

    # Go through the edges in the network
    for i, edge1 in enumerate(e_balls):

        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2),
                      min(100.0, 100 * (len(s_balls) + 0.5 * (len(e_balls))) / ((3 / 2) * len(v_balls)))), end="")

        # Go through the edge's atoms combinations
        for j in range(3):
            # Get the atoms and their sorted list of ndxs
            balls = [edge1[j], edge1[(j + 1) % 3]]
            ball_ndxs = balls[:]
            ball_ndxs.sort()
            # If the surface has been found before continue
            surf_ndx = ndx_search(s_balls, ball_ndxs)
            # If the edge has been found before, continue
            if len(s_balls) > surf_ndx and ball_ndxs == s_balls[surf_ndx]:
                continue

            # Limit the list of verts to possible vertices
            # max_vert_ndx =

            # Put together a list of verts that have our atoms
            verts = []
            for k, vert2 in enumerate(v_balls):
                # If the surface's atoms are shared with the vertex, add it to the list
                if len([0 for ndx in ball_ndxs if ndx in vert2]) == 2:
                    verts.append(k)

            # Put together a list of edges that have our atoms
            edges = []
            # Go through the edges in the system
            for k, edge2 in enumerate(e_balls):
                # If the surface's atom s are in the edge add it
                if len([0 for ndx in ball_ndxs if ndx in edge2]) == 2:
                    edges.append(k)

            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):

                no_surf = False
                # Check to see if the surface is worth adding
                for vert_ndx in verts:
                    if len(v_edges[vert_ndx]) <= 2:
                        no_surf = True
                if no_surf:
                    continue
                incomplete = False
                for vert in verts:
                    if len(v_edges[vert]) > 3:
                        incomplete = True
                if incomplete:
                    continue

                s_balls.insert(surf_ndx, ball_ndxs)
                s_edges.insert(surf_ndx, edges)
                s_verts.insert(surf_ndx, verts)
    return s_balls, s_verts, s_edges


def get_build_surfs(b_verts, b_edges, v_balls, v_edges, e_balls, start_time):
    # Set up the surface lists
    s_balls, s_verts, s_edges = [], [], []

    # Go through the edges in the network
    for i, edge1 in enumerate(e_balls):

        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2),
                      min(100.0, 100 * (len(s_balls) + 0.5 * (len(e_balls))) / ((3 / 2) * len(v_balls)))), end="")
        # Get the possible surfs from the edge's atoms
        test_surfs = [edge1[:2], edge1[1:], edge1[::2]]
        # Go through each possible surface for the edge
        for test_surf in test_surfs:
            # If the surface has been found before continue
            surf_ndx = ndx_search(s_balls, test_surf)
            # If the surface has been found before, continue
            if len(s_balls) > surf_ndx and test_surf == s_balls[surf_ndx]:
                continue
            # Set up the surf edges and surf verts lists
            surf_edges, surf_verts = [], []
            # Go through the atoms in the surface looking for edge candidates
            for ball in test_surf:
                # Get the atom's edges
                for edge in b_edges[ball]:
                    # Get the edges atoms
                    edge2 = e_balls[edge]
                    # If the number of shared atoms is 2 add the edge to the test surf's list of edges
                    if len([_ for _ in edge2 if _ in test_surf]) == 2 and edge not in surf_edges:
                        # Add the edge
                        surf_edges.append(edge)
                # Get the atom's vertices
                for vert in b_verts[ball]:
                    # Get the vertices atoms
                    vert2 = v_balls[vert]
                    # If the number of shared atoms is 2 add the edge to the test surf's list of edges
                    if len([_ for _ in vert2 if _ in test_surf]) == 2 and vert not in surf_verts:
                        # Add the edge
                        surf_verts.append(vert)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(surf_verts) == len(surf_edges):

                no_surf = False
                # Check to see if the surface is worth adding
                for vert_ndx in surf_verts:
                    if len(v_edges[vert_ndx]) <= 2:
                        no_surf = True
                if no_surf:
                    continue

                s_balls.insert(surf_ndx, test_surf)
                s_edges.insert(surf_ndx, surf_edges)
                s_verts.insert(surf_ndx, surf_verts)
    return s_balls, s_verts, s_edges


def add_build_surfs(num_balls, s_balls, num_verts, s_verts, num_edges, s_edges):
    # Atoms
    b_surfs = [[] for _ in range(num_balls)]
    for i, surf_balls in enumerate(s_balls):
        for j in surf_balls:
            b_surfs[j] += [i]

    # Verts
    v_surfs = [[] for _ in range(num_verts)]
    for i, surf_verts in enumerate(s_verts):
        for j in surf_verts:
            v_surfs[j] += [i]

    # Edges
    e_surfs = [[] for _ in range(num_edges)]
    for i, surf_edges in enumerate(s_edges):
        for j in surf_edges:
            e_surfs[j] += [i]

    return b_surfs, v_surfs, e_surfs


def build(v_balls, v_locs, v_dubs, num_balls, my_time):
    """
    Checks the atoms of the vertices for patterns and creates edges and surfaces
    :param my_group:
    :param from_scratch:
    :param net: Network object to pull information from
    """

    # Create the lists
    b_verts, b_edges = [[] for _ in range(num_balls)], [[] for _ in range(num_balls)]

    # Add the vertices to the atoms
    for i, ndx in enumerate(v_balls):
        for j in ndx:
            b_verts[j].append(i)

    ################################################# Create the edges #################################################

    # Fill in the doublets and set their outer edges
    e_balls, e_verts = get_build_edges(b_verts, v_balls, v_locs, v_dubs, my_time)

    # Add the edges to their atoms and vertices
    b_edges, v_edges = add_build_edges(num_balls, e_balls, len(v_balls), e_verts)

    ################################################### Create the surfaces ############################################

    # Get the surfaces
    s_balls, s_verts, s_edges = get_build_surfs(b_verts, b_edges, v_balls, v_edges, e_balls, my_time)

    # Add the surface objects to their 181L indices
    b_surfs, v_surfs, e_surfs = add_build_surfs(num_balls, s_balls, len(v_balls), s_verts, len(e_balls), s_edges)

    # Package the lists neatly for easier parsing
    ball_lists = {'averts': b_verts, 'aedges': b_edges, 'asurfs': b_surfs}
    vert_lists = {'vedges': v_edges, 'vsurfs': v_surfs}
    edge_lists = {'eatoms': e_balls, 'everts': e_verts, 'esurfs': e_surfs}
    surf_lists = {'satoms': s_balls, 'sverts': s_verts, 'sedges': s_edges}

    # Return the dictionary lists
    return ball_lists, vert_lists, edge_lists, surf_lists
