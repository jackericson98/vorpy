import time
import numpy as np
from itertools import combinations
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
    if counts.get('unique_surface_keys') is not None:
        print(f'Unique surface keys:   {counts["unique_surface_keys"]:,}')
    if counts.get('vertex_surface_candidates') is not None:
        print(f'Vertex-surface keys:   {counts["vertex_surface_candidates"]:,}')
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


def doublify(b_verts, v_balls, v_locs, v_dubs, timings=None, counts=None):
    """Construct doublet edges while preserving duplicate 3-ball edge definitions."""
    timing_start = time.perf_counter()
    e_balls, e_verts, e_surfs = [], [], []
    doublet_candidate_visits = 0

    for i in range(len(v_balls)):
        if i >= len(v_balls) - 1 or v_dubs[i + 1] != 1:
            continue

        # Any vertex sharing exactly three balls with this doublet must belong to
        # at least one of its four ball->vertex lists. Sorting reproduces the
        # original full-network j iteration order without scanning every vertex.
        candidates = set()
        for ball in v_balls[i]:
            candidates.update(b_verts[ball])

        con_verts = []
        for j in sorted(candidates):
            doublet_candidate_visits += 1
            if len([0 for ball in v_balls[i] if ball in v_balls[j]]) == 3:
                con_verts.append(j)

        dub_verts, dub_dub_verts = [], []
        for j in con_verts:
            if calc_dist(np.array(v_locs[j]), np.array(v_locs[i])) < calc_dist(np.array(v_locs[j]), np.array(v_locs[i + 1])):
                dub_verts.append(j)
            else:
                dub_dub_verts.append(j)

        # IMPORTANT: keep list insertion semantics here. Doublet topology can
        # legitimately contain multiple edges with identical three-ball keys but
        # different endpoint vertices, so these must NOT be collapsed in a dict.
        known_edges = []
        for j in dub_verts:
            edge_balls = [ball for ball in v_balls[i] if ball in v_balls[j]]
            edge_ndx = ndx_search(e_balls, edge_balls)
            e_balls.insert(edge_ndx, edge_balls)
            e_verts.insert(edge_ndx, [i, j])
            e_surfs.insert(edge_ndx, [])
            known_edges.append(edge_balls)

        for j in dub_dub_verts:
            edge_balls = [ball for ball in v_balls[i] if ball in v_balls[j]]
            edge_ndx = ndx_search(e_balls, edge_balls)
            e_balls.insert(edge_ndx, edge_balls)
            e_verts.insert(edge_ndx, [i + 1, j])
            e_surfs.insert(edge_ndx, [])
            known_edges.append(edge_balls)

        potential_edges = [[v_balls[i][k], v_balls[i][(k + 1) % 4], v_balls[i][(k + 2) % 4]] for k in range(4)]
        for ndx in potential_edges:
            ndx.sort()

        inner_edges = [ndx for ndx in potential_edges if ndx not in known_edges]
        for edge in inner_edges:
            edge_ndx = ndx_search(e_balls, edge)
            e_balls.insert(edge_ndx, edge)
            e_verts.insert(edge_ndx, [i, i + 1])

    _record_timing(timings, 'doublets', timing_start)
    if counts is not None:
        counts['doublet_candidates'] = doublet_candidate_visits
        counts['doublet_edges'] = len(e_balls)
    return e_balls, e_verts


def get_build_edges(b_verts, v_balls, v_locs, v_dubs, start_time, net=None, timings=None, counts=None):
    """Construct regular edges from indexed 3-ball combinations while retaining exact doublet multiplicity."""
    e_balls, e_verts = doublify(b_verts, v_balls, v_locs, v_dubs, timings=timings, counts=counts)
    regular_edge_start = time.perf_counter()

    # Index every vertex by each of its four possible 3-ball definitions. Store
    # the fourth ball too: two vertices define a regular edge only when they share
    # exactly three balls, not all four.
    triple_index = {}
    for vi, vert in enumerate(v_balls):
        vert_set = set(vert)
        for triple in combinations(sorted(vert), 3):
            fourth = next(ball for ball in vert_set if ball not in triple)
            triple_index.setdefault(triple, []).append((vi, fourth))

    # Any key already produced by doublify must be considered existing exactly as
    # in the original ndx_search check. We preserve all duplicate doublet entries
    # in e_balls/e_verts; this set is only for regular-edge existence testing.
    existing_keys = {tuple(edge) for edge in e_balls}
    regular_edges = {}
    regular_edge_keys = 0

    # i order matches the original loop. For a given triple, the first candidate j
    # with a different fourth ball is the first vertex that shares exactly 3 balls.
    for i, vert in enumerate(v_balls):
        if net is not None and i % 1000 == 0:
            percentage = 10.0 + 40.0 * (i + 1) / len(v_balls)
            net.update_progress("Connecting network", percentage)
        if v_dubs[i] == 1 or (i + 1 < len(v_dubs) and v_dubs[i + 1] == 1):
            continue

        vert_set = set(vert)
        for triple in combinations(sorted(vert), 3):
            regular_edge_keys += 1
            if triple in existing_keys or triple in regular_edges:
                continue

            fourth_i = next(ball for ball in vert_set if ball not in triple)
            for j, fourth_j in triple_index[triple]:
                if fourth_j != fourth_i:
                    regular_edges[triple] = [i, j]
                    break

    # Merge regular edges with the duplicate-preserving doublet lists, then sort
    # by the 3-ball definition. Python's stable sort keeps duplicate doublet edges
    # in their existing relative order.
    combined = [(tuple(edge), verts) for edge, verts in zip(e_balls, e_verts)]
    combined.extend((edge, verts) for edge, verts in regular_edges.items())
    combined.sort(key=lambda item: item[0])
    e_balls = [list(edge) for edge, _ in combined]
    e_verts = [verts for _, verts in combined]

    _record_timing(timings, 'regular_edges', regular_edge_start)
    if counts is not None:
        counts['edge_candidates'] = regular_edge_keys
        counts['unique_edge_keys'] = len(triple_index)
        counts['regular_edges'] = len(regular_edges)
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


def get_build_surfs(b_verts, b_edges, v_balls, v_edges, e_balls, start_time, net=None, group=None,
                    interface=False, iface_grps=None, timings=None, counts=None):
    """Construct surfaces from direct 2-ball topology indices while preserving original validity rules."""
    surface_start = time.perf_counter()

    # Build each candidate surface's edge membership once. Enumerating e_balls
    # in edge-index order preserves the same surf_edges order as the original
    # b_edges scan. Duplicate doublet edges remain distinct because we store
    # edge indices, not just three-ball keys.
    surf_edge_map = {}
    surface_candidates = 0
    for edge_ndx, edge_balls in enumerate(e_balls):
        for surf in combinations(edge_balls, 2):
            surface_candidates += 1
            key = tuple(sorted(surf))
            surf_edge_map.setdefault(key, []).append(edge_ndx)

    # Build vertex membership for only surface keys that actually occur in the
    # edge network. Enumerating vertices in index order preserves surf_verts
    # ordering from the original b_verts scan.
    surf_vert_map = {key: [] for key in surf_edge_map}
    vertex_surface_candidates = 0
    for vert_ndx, vert_balls in enumerate(v_balls):
        for surf in combinations(vert_balls, 2):
            vertex_surface_candidates += 1
            key = tuple(sorted(surf))
            if key in surf_vert_map:
                surf_vert_map[key].append(vert_ndx)

    s_balls, s_verts, s_edges = [], [], []
    keys = sorted(surf_edge_map)
    for n, key in enumerate(keys):
        if net is not None and n % 5000 == 0 and keys:
            percentage = 60.0 + 30.0 * (n + 1) / len(keys)
            net.update_progress("Building Topology", percentage)

        test_surf = list(key)
        if interface and not spans_interface(test_surf, iface_grps):
            continue
        if not interface and group is not None and not belongs_to_group(test_surf, group):
            continue

        surf_edges = surf_edge_map[key]
        surf_verts = surf_vert_map[key]
        if len(surf_verts) != len(surf_edges):
            continue

        if interface:
            surf_edge_set = set(surf_edges)
            invalid_surface = False
            for vert_ndx in surf_verts:
                surface_degree = sum(edge_ndx in surf_edge_set for edge_ndx in v_edges[vert_ndx])
                if surface_degree != 2:
                    invalid_surface = True
                    break
            if invalid_surface:
                continue
        else:
            no_surf = False
            for vert_ndx in surf_verts:
                if len(v_edges[vert_ndx]) <= 2:
                    no_surf = True
                    break
            if no_surf:
                continue

        s_balls.append(test_surf)
        s_edges.append(surf_edges)
        s_verts.append(surf_verts)

    _record_timing(timings, 'surfaces', surface_start)
    if counts is not None:
        counts['surface_candidates'] = surface_candidates
        counts['unique_surface_keys'] = len(surf_edge_map)
        counts['vertex_surface_candidates'] = vertex_surface_candidates
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
          group=None, interface=False, iface_grps=None, net=None):
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

    if net is not None:
        net.update_progress("Building Topology", 0.0)

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
    if net is not None:
        net.update_progress("Building Topology", 10.0)

    # Fill in the doublets and regular edges.
    e_balls, e_verts = get_build_edges(b_verts, v_balls, v_locs, v_dubs, my_time, net=net, timings=timings, counts=counts)

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
    if net is not None:
        net.update_progress("Building Topology", 60.0)
    _record_timing(timings, 'edge_adjacency', edge_adjacency_start)

    ################################################### Create the surfaces ############################################

    if net is not None:
        net.update_progress("Building Topology", 60.0)

    s_balls, s_verts, s_edges = get_build_surfs(b_verts, b_edges, v_balls, v_edges, e_balls, my_time, group=group,
                                                interface=interface, iface_grps=iface_grps, timings=timings,
                                                net=net, counts=counts)
    counts['surfs'] = len(s_balls)
    if net is not None:
        net.update_progress("Building Topology", 90.0)

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
    if net is not None:
        net.update_progress("Building Topology", 97.0)

    _record_timing(timings, 'surface_adjacency', surface_adjacency_start)

    packaging_start = time.perf_counter()
    ball_lists = {'verts': b_verts, 'edges': b_edges, 'surfs': b_surfs}
    vert_lists = {'edges': v_edges, 'surfs': v_surfs}
    edge_lists = {'balls': e_balls, 'verts': e_verts, 'surfs': e_surfs}
    surf_lists = {'balls': s_balls, 'verts': s_verts, 'edges': s_edges}
    _record_timing(timings, 'packaging', packaging_start)

    timings['total'] = time.perf_counter() - build_start

    if net is not None:
        net.update_progress("Building Topology", 100.0)

    # Timing is always collected; -v / net.verbose controls only printing.
    if net is not None and net.settings.get('verbose', False):
        _print_build_timings(timings, counts)

    return ball_lists, vert_lists, edge_lists, surf_lists