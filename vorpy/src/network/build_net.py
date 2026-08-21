import time
import numpy as np
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import get_time
from vorpy.src.calculations import ndx_search


def _record_timing(timings, name, start):
    """Add elapsed wall time to a timing bucket."""
    if timings is not None:
        timings[name] = timings.get(name, 0.0) + (time.perf_counter() - start)


def _print_build_timings(timings, counts):
    """Print a compact timing summary for network construction."""
    total = timings.get('total', 0.0)
    phases = [
        ('Setup', 'setup'),
        ('Ball -> vertex index', 'ball_vertex_index'),
        ('Doublet edges', 'doublets'),
        ('Regular edge search', 'regular_edges'),
        ('Edge filtering', 'edge_filter'),
        ('Edge adjacency', 'edge_adjacency'),
        ('Surface construction', 'surfaces'),
        ('Interface validation', 'validation'),
        ('Surface adjacency', 'surface_adjacency'),
        ('Packaging', 'packaging'),
    ]
    print('\n\n' + '=' * 70)
    print('NETWORK BUILD TIMING')
    print('=' * 70)
    for label, key in phases:
        seconds = timings.get(key, 0.0)
        pct = 100.0 * seconds / total if total else 0.0
        print(f'{label:<24} {seconds:>10.4f} s  {pct:>6.2f} %')
    print('-' * 70)
    print(f'{"TOTAL":<24} {total:>10.4f} s  {100.0 if total else 0.0:>6.2f} %')
    print('\nNETWORK BUILD SIZE')
    print(f'Balls:     {counts.get("balls", 0):,}')
    print(f'Vertices:  {counts.get("verts", 0):,}')
    print(f'Doublets:  {counts.get("doublets", 0):,}')
    print(f'Edges:     {counts.get("edges", 0):,}')
    print(f'Surfaces:  {counts.get("surfs", 0):,}')
    if counts.get('edge_candidates') is not None:
        print(f'Edge candidate visits: {counts["edge_candidates"]:,}')
    if counts.get('surface_candidates') is not None:
        print(f'Surface candidates:    {counts["surface_candidates"]:,}')
    print('=' * 70)


def spans_interface(ball_indices, iface_grps):
    """
    Return True when ball_indices contains at least one defining ball from
    each side of an interface.

    iface_grps must contain exactly two collections of network-local ball
    indices.
    """
    if iface_grps is None:
        return True

    if len(iface_grps) != 2:
        raise ValueError(
            "iface_grps must contain exactly two ball-index collections."
        )

    ball_set = set(ball_indices)
    group1 = set(iface_grps[0])
    group2 = set(iface_grps[1])

    return bool(ball_set & group1) and bool(ball_set & group2)


def belongs_to_group(ball_indices, group):
    """Return True when at least one defining ball belongs to the requested group."""
    if group is None:
        return True
    return any(ball in group for ball in ball_indices)


############################################## Doublets ################################################################


def doublify(v_balls, v_locs, v_dubs, timings=None):
    """
    Construct edges associated with doublet vertices.

    Doublet vertices represent two geometric locations defined by the same
    four balls. Connecting vertices are assigned to the nearer doublet
    location, and any remaining three-ball combinations form the internal
    edges between the two doublet locations.

    Parameters
    ----------
    v_balls : list
        Four defining ball indices for each vertex.
    v_locs : list
        Vertex locations.
    v_dubs : list
        Doublet flags for each vertex.

    Returns
    -------
    tuple
        ``(e_balls, e_verts)`` containing edge-defining ball indices and
        the vertex indices connected by each edge.
    """
    timing_start = time.perf_counter()
    e_balls, e_verts, e_surfs = [], [], []
    # Go through the doublets
    for i in range(len(v_balls)):

        # Skip the verts that aren't doublets
        if i >= len(v_balls) - 1 or v_dubs[i + 1] != 1:
            continue

        ################################################ Create the outer edges ########################################

        # Find all vertices that match edges with the doublet's balls
        con_verts = []
        for j in range(len(v_balls)):
            if len([0 for _ in v_balls[i] if _ in v_balls[j]]) == 3:
                con_verts.append(j)

        # Divide the connecting outer vertices between the two doublet vertices
        dub_verts, dub_dub_verts = [], []
        for j in con_verts:
            # Decide between the two sides of the doublet for the outer vertex
            if calc_dist(np.array(v_locs[j]), np.array(v_locs[i])) < calc_dist(np.array(v_locs[j]), np.array(v_locs[i + 1])):
                dub_verts.append(j)
            else:
                dub_dub_verts.append(j)

        known_edges = []
        # Create the edge objects for each of the vertices connected to the primary doublet vertex
        for j in dub_verts:
            # Create the edge from the balls in both dub and vert and add it to the network and each vertex
            edge_balls = [_ for _ in v_balls[i] if _ in v_balls[j]]
            edge_ndx = ndx_search(e_balls, edge_balls)
            e_balls.insert(edge_ndx, edge_balls)
            e_verts.insert(edge_ndx, [i, j])
            e_surfs.insert(edge_ndx, [])
            known_edges.append(edge_balls)

        # Create the edge objects for each of the vertices connected to the secondary doublet vertex
        for j in dub_dub_verts:
            # Create the edge from the balls in both dub.doublet and vert and add it to the network and each vertex
            edge_balls = [_ for _ in v_balls[i] if _ in v_balls[j]]
            edge_ndx = ndx_search(e_balls, edge_balls)
            e_balls.insert(edge_ndx, edge_balls)
            e_verts.insert(edge_ndx, [i + 1, j])
            e_surfs.insert(edge_ndx, [])
            known_edges.append(edge_balls)

        ########################################## Create the inner edges ##########################################

        # Create a list of every edge possibility
        potential_edges = [[v_balls[i][k], v_balls[i][(k + 1) % 4], v_balls[i][(k + 2) % 4]] for k in range(4)]
        for ndx in potential_edges:
            ndx.sort()

        # Gather the other combinations of balls and create the remaining inner balls
        inner_edges = [ndx for ndx in potential_edges if ndx not in known_edges]

        # Add the edges to the network and the doublet vertices
        for edge in inner_edges:
            edge_ndx = ndx_search(e_balls, edge)
            e_balls.insert(edge_ndx, edge)
            e_verts.insert(edge_ndx, [i, i + 1])
    # Return the partial lists
    _record_timing(timings, 'doublets', timing_start)
    return e_balls, e_verts


def get_build_edges(b_verts, v_balls, v_locs, v_dubs, start_time, timings=None, counts=None):
    """
    Construct network edges from shared three-ball vertex definitions.

    Doublet-specific edges are created first. Remaining edges are discovered
    by finding pairs of vertices that share exactly three defining balls.

    Parameters
    ----------
    b_verts : list
        Vertex indices associated with each ball.
    v_balls : list
        Defining ball indices for each vertex.
    v_locs : list
        Vertex locations.
    v_dubs : list
        Doublet flags for each vertex.
    start_time : float
        Build start time used for progress reporting.

    Returns
    -------
    tuple
        ``(e_balls, e_verts)`` containing edge definitions and their vertices.
    """
    # Get the doublet edges
    e_balls, e_verts = doublify(v_balls, v_locs, v_dubs, timings=timings)

    regular_edge_start = time.perf_counter()
    edge_candidate_visits = 0

    # Go through the vertices in the network searching for potential edges
    for i, vert1 in enumerate(v_balls):
        # Print the time and process
        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * (0.5 * len(e_balls)) / (3 / 2 * len(v_balls)))),
              end="")

        # If the vertex is a doublet it has its edges already, so skip
        if v_dubs[i] == 1 or (i + 1 < len(v_dubs) and v_dubs[i + 1] == 1):
            continue

        # Go through the balls in the vertex looking for shared balls
        for ball in vert1:
            # Go through the vertices in each ball
            for j in b_verts[ball]:
                edge_candidate_visits += 1
                # Get the balls for vert2
                vert2 = v_balls[j]
                # Check the number of shared balls between vert1 and vert2
                shared_balls = [_ for _ in vert1 if _ in vert2]
                # Check if this edge is real
                if len(shared_balls) == 3:
                    # Get the index of the edge in the edge list
                    edge_ndx = ndx_search(e_balls, shared_balls)
                    # Check if we have found this edge before
                    if edge_ndx >= len(e_balls) or e_balls[edge_ndx] != shared_balls:
                        # Add the edges balls and the edges verts to their respective lists
                        e_balls.insert(edge_ndx, shared_balls)
                        e_verts.insert(edge_ndx, [i, j])
    # Return the edge's balls and verts
    _record_timing(timings, 'regular_edges', regular_edge_start)
    if counts is not None:
        counts['edge_candidates'] = edge_candidate_visits
    return e_balls, e_verts


def add_build_edges(num_balls, e_balls, num_verts, e_verts):
    """Create ball-to-edge and vertex-to-edge adjacency lists."""
    # Create the empty ball list of edges
    b_edges = [[] for _ in range(num_balls)]
    # Create the empty vertex list of edges
    v_edges = [[] for _ in range(num_verts)]
    # Go through the edges in the network
    for i, edge_balls in enumerate(e_balls):
        # Go through the balls in the edge
        for j in edge_balls:
            # Add the edge to each ball
            b_edges[j].append(i)
        # Get the edges vertices
        edge_verts = e_verts[i]
        # Go through the verts in the edge
        for j in edge_verts:
            v_edges[j].append(i)
    # Return the newly filled in lists
    return b_edges, v_edges


def get_build_surfs(b_verts, b_edges, v_balls, v_edges, e_balls, start_time, group=None,
                    interface=False, iface_grps=None, timings=None, counts=None):
    """
    Construct valid two-ball surfaces from the connected edge network.

    Each three-ball edge defines three possible two-ball surfaces. Candidate
    surfaces are retained only when their edge and vertex topology is complete.
    Group networks retain surfaces containing at least one group ball, while
    interface networks retain only surfaces spanning both interface groups.

    Parameters
    ----------
    b_verts : list
        Vertex indices associated with each ball.
    b_edges : list
        Edge indices associated with each ball.
    v_balls : list
        Defining ball indices for each vertex.
    v_edges : list
        Edge indices associated with each vertex.
    e_balls : list
        Three defining ball indices for each edge.
    start_time : float
        Build start time used for progress reporting.
    group : collection, optional
        Ball indices belonging to a normal group network.
    interface : bool, optional
        Whether the network represents an interface.
    iface_grps : tuple, optional
        Two ball-index collections defining the interface sides.

    Returns
    -------
    tuple
        ``(s_balls, s_verts, s_edges)`` describing surface definitions and
        their associated vertices and edges.
    """
    surface_start = time.perf_counter()
    surface_candidates = 0

    # Set up the surface lists
    s_balls, s_verts, s_edges = [], [], []

    # Go through the edges in the network
    for i, edge1 in enumerate(e_balls):

        the_time = time.perf_counter() - start_time
        h, m, s = get_time(the_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: connecting network: {:.2f} %"
              .format(int(h), int(m), round(s, 2),
                      min(100.0, 100 * (len(s_balls) + 0.5 * (len(e_balls))) / ((3 / 2) * len(v_balls)))), end="")
        # Get the possible surfs from the edge's balls
        test_surfs = [edge1[:2], edge1[1:], edge1[::2]]
        # Go through each possible surface for the edge
        for test_surf in test_surfs:
            surface_candidates += 1
            # Check if the surface is in the interface or not
            if interface and not spans_interface(test_surf, iface_grps):
                continue
            if not interface and group is not None and not belongs_to_group(test_surf, group):
                continue

            # If the surface has been found before continue
            surf_ndx = ndx_search(s_balls, test_surf)
            # If the surface has been found before, continue
            if len(s_balls) > surf_ndx and test_surf == s_balls[surf_ndx]:
                continue
            # Set up the surf edges and surf verts lists
            surf_edges, surf_verts = [], []
            # Go through the balls in the surface looking for edge candidates
            for ball in test_surf:
                # Get the ball's edges
                for edge in b_edges[ball]:
                    # Get the edges balls
                    edge2 = e_balls[edge]
                    # If the number of shared balls is 2 add the edge to the test surf's list of edges
                    if len([_ for _ in edge2 if _ in test_surf]) == 2 and edge not in surf_edges:
                        # Add the edge
                        surf_edges.append(edge)
                # Get the ball's vertices
                for vert in b_verts[ball]:
                    # Get the vertices balls
                    vert2 = v_balls[vert]
                    # If the number of shared balls is 2 add the edge to the test surf's list of edges
                    if len([_ for _ in vert2 if _ in test_surf]) == 2 and vert not in surf_verts:
                        # Add the edge
                        surf_verts.append(vert)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(surf_verts) == len(surf_edges):

                if interface:
                    # In an interface-only polygon, each surface vertex should
                    # have exactly two edges belonging to this surface.
                    surf_edge_set = set(surf_edges)

                    invalid_surface = False

                    for vert_ndx in surf_verts:
                        surface_degree = sum(
                            edge_ndx in surf_edge_set
                            for edge_ndx in v_edges[vert_ndx]
                        )

                        if surface_degree != 2:
                            invalid_surface = True
                            break

                    if invalid_surface:
                        continue

                else:
                    # Existing complete-network requirement.
                    no_surf = False

                    for vert_ndx in surf_verts:
                        if len(v_edges[vert_ndx]) <= 2:
                            no_surf = True
                            break

                    if no_surf:
                        continue

                s_balls.insert(surf_ndx, test_surf)
                s_edges.insert(surf_ndx, surf_edges)
                s_verts.insert(surf_ndx, surf_verts)
    _record_timing(timings, 'surfaces', surface_start)
    if counts is not None:
        counts['surface_candidates'] = surface_candidates
    return s_balls, s_verts, s_edges


def add_build_surfs(num_balls, s_balls, num_verts, s_verts, num_edges, s_edges):
    """
    Create reverse adjacency lists linking balls, vertices, and edges to surfaces.

    Returns
    -------
    tuple
        ``(b_surfs, v_surfs, e_surfs)``.
    """
    # balls
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


def build(v_balls, v_locs, v_dubs, num_balls, my_time,
          group=None, interface=False, iface_grps=None):
    """
    Connect solved vertices into the network's edges and surfaces.

    Vertices sharing three defining balls form edges, and connected edges
    sharing two defining balls form surfaces. Normal group networks retain
    topology containing at least one requested group ball. Interface networks
    retain only topology spanning both interface groups.

    Parameters
    ----------
    v_balls : list
        Four defining ball indices for each vertex.
    v_locs : list
        Vertex locations.
    v_dubs : list
        Doublet flags for each vertex.
    num_balls : int
        Total number of balls in the parent system.
    my_time : float
        Build start time used for progress reporting.
    group : collection, optional
        Ball indices defining a normal group network.
    interface : bool, optional
        Whether interface-only topology should be constructed.
    iface_grps : tuple, optional
        Two ball-index collections defining the interface.

    Returns
    -------
    tuple
        Dictionaries containing ball, vertex, edge, and surface adjacency data.
    """

    build_start = time.perf_counter()
    timings = {}
    counts = {'balls': num_balls, 'verts': len(v_balls), 'doublets': sum(1 for dub in v_dubs if dub == 1)}

    setup_start = time.perf_counter()
    if group is not None:
        group = set(group)

    if interface and iface_grps is None:
        raise ValueError("Interface network construction requires iface_grps.")

    if iface_grps is not None:
        iface_grps = tuple(set(group_indices) for group_indices in iface_grps)
        overlap = iface_grps[0] & iface_grps[1]
        if overlap:
            raise ValueError(f"Interface groups overlap by {len(overlap)} balls.")
    _record_timing(timings, 'setup', setup_start)

    # Create the ball -> vertex index.
    ball_vertex_start = time.perf_counter()
    b_verts, b_edges = [[] for _ in range(num_balls)], [[] for _ in range(num_balls)]
    for i, ndx in enumerate(v_balls):
        for j in ndx:
            b_verts[j].append(i)
    _record_timing(timings, 'ball_vertex_index', ball_vertex_start)

    ################################################# Create the edges #################################################

    # Fill in the doublets and regular edges.
    e_balls, e_verts = get_build_edges(b_verts, v_balls, v_locs, v_dubs, my_time, timings=timings, counts=counts)

    edge_filter_start = time.perf_counter()
    if interface:
        retained_edges = [(edge_balls, edge_verts) for edge_balls, edge_verts in zip(e_balls, e_verts)
                          if spans_interface(edge_balls, iface_grps)]
        e_balls = [edge_balls for edge_balls, _ in retained_edges]
        e_verts = [edge_verts for _, edge_verts in retained_edges]
    elif group is not None:
        retained_edges = [(edge_balls, edge_verts) for edge_balls, edge_verts in zip(e_balls, e_verts)
                          if belongs_to_group(edge_balls, group)]
        e_balls = [edge_balls for edge_balls, _ in retained_edges]
        e_verts = [edge_verts for _, edge_verts in retained_edges]
    _record_timing(timings, 'edge_filter', edge_filter_start)
    counts['edges'] = len(e_balls)

    edge_adjacency_start = time.perf_counter()
    b_edges, v_edges = add_build_edges(num_balls, e_balls, len(v_balls), e_verts)
    _record_timing(timings, 'edge_adjacency', edge_adjacency_start)

    ################################################### Create the surfaces ############################################

    s_balls, s_verts, s_edges = get_build_surfs(b_verts, b_edges, v_balls, v_edges, e_balls, my_time, group=group,
                                                interface=interface, iface_grps=iface_grps, timings=timings, counts=counts)
    counts['surfs'] = len(s_balls)

    validation_start = time.perf_counter()
    if interface:
        invalid_edges = [edge_balls for edge_balls in e_balls if not spans_interface(edge_balls, iface_grps)]
        invalid_surfs = [surf_balls for surf_balls in s_balls if not spans_interface(surf_balls, iface_grps)]
        if invalid_edges or invalid_surfs:
            raise RuntimeError("Interface topology validation failed.\n"
                               f"Invalid edges: {len(invalid_edges)}\n"
                               f"Invalid surfaces: {len(invalid_surfs)}")
    _record_timing(timings, 'validation', validation_start)

    surface_adjacency_start = time.perf_counter()
    b_surfs, v_surfs, e_surfs = add_build_surfs(num_balls, s_balls, len(v_balls), s_verts, len(e_balls), s_edges)
    _record_timing(timings, 'surface_adjacency', surface_adjacency_start)

    packaging_start = time.perf_counter()
    ball_lists = {'verts': b_verts, 'edges': b_edges, 'surfs': b_surfs}
    vert_lists = {'edges': v_edges, 'surfs': v_surfs}
    edge_lists = {'balls': e_balls, 'verts': e_verts, 'surfs': e_surfs}
    surf_lists = {'balls': s_balls, 'verts': s_verts, 'edges': s_edges}
    _record_timing(timings, 'packaging', packaging_start)

    timings['total'] = time.perf_counter() - build_start
    _print_build_timings(timings, counts)
    return ball_lists, vert_lists, edge_lists, surf_lists