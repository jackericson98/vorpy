import time
import numpy as np

from vorpy.src.network.build_net import (
    _minimum_distance_edge_pairing,
    get_build_edges,
    add_build_edges,
    get_build_surfs,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_ball_vertex_index(v_balls, num_balls=None):
    """Return the ball -> vertex adjacency used by build_net."""
    if num_balls is None:
        num_balls = max(max(v) for v in v_balls) + 1

    b_verts = [[] for _ in range(num_balls)]
    for vert_ndx, balls in enumerate(v_balls):
        for ball in balls:
            b_verts[ball].append(vert_ndx)
    return b_verts


def _normalized_pairs(pairs):
    """Normalize undirected vertex pairs for order-independent assertions."""
    return {tuple(sorted(pair)) for pair in pairs}


def _ala78_ca_c_regression_case():
    """
    Synthetic geometric realization of the B/ALA78 CA-C topology that exposed
    the AW multi-edge bug.

    The important topology is preserved from the real failure:

        shared surface balls: 2590 (CA), 2596 (C)

        repeated 3-ball key:
            [2590, 2596, 2598]

        four vertices containing that key:
            old v25 -> [2590, 2596, 2598, 2599]
            old v26 -> [2590, 2591, 2596, 2598]
            old v28 -> [2587, 2590, 2596, 2598]
            old v31 -> [2590, 2596, 2598, 41358]

    The vertices are stored contiguously here so the unit test does not depend
    on unrelated vertices from the full Streptavidin solve.

    Their synthetic coordinates are arranged around the expected CA-C surface
    cycle.  That makes the two geometrically local segments for the repeated
    triple:

        old v28 <-> old v26
        old v25 <-> old v31
    """

    # old diagnostic vertex index -> local unit-test vertex index
    old_to_local = {
        0: 0,
        1: 1,
        2: 2,
        11: 3,
        12: 4,
        13: 5,
        15: 6,
        17: 7,
        24: 8,
        25: 9,
        26: 10,
        28: 11,
        31: 12,
    }

    v_balls = [
        [2586, 2588, 2590, 2596],   # old v0
        [2586, 2587, 2590, 2596],   # old v1
        [2588, 2589, 2590, 2596],   # old v2
        [2590, 2592, 2595, 2596],   # old v11
        [2590, 2592, 2593, 2596],   # old v12
        [2590, 2595, 2596, 2597],   # old v13
        [2589, 2590, 2596, 2597],   # old v15
        [2590, 2593, 2596, 41358],  # old v17
        [2590, 2591, 2596, 2599],   # old v24
        [2590, 2596, 2598, 2599],   # old v25
        [2590, 2591, 2596, 2598],   # old v26
        [2587, 2590, 2596, 2598],   # old v28
        [2590, 2596, 2598, 41358],  # old v31
    ]

    # Expected boundary cycle around the CA-C surface:
    # old 0 -> 1 -> 28 -> 26 -> 24 -> 25 -> 31 -> 17
    #      -> 12 -> 11 -> 13 -> 15 -> 2 -> 0
    old_cycle = [0, 1, 28, 26, 24, 25, 31, 17, 12, 11, 13, 15, 2]

    # Place the cycle on a circle.  The repeated-triple pairings that should
    # be selected are therefore spatially adjacent.
    v_locs = [None] * len(v_balls)
    n = len(old_cycle)
    for cycle_ndx, old_ndx in enumerate(old_cycle):
        theta = 2.0 * np.pi * cycle_ndx / n
        local_ndx = old_to_local[old_ndx]
        v_locs[local_ndx] = np.array(
            [np.cos(theta), np.sin(theta), 0.0],
            dtype=float,
        )

    v_dubs = [0] * len(v_balls)

    expected_repeated_pairs = {
        tuple(sorted((old_to_local[28], old_to_local[26]))),
        tuple(sorted((old_to_local[25], old_to_local[31]))),
    }

    return {
        "v_balls": v_balls,
        "v_locs": v_locs,
        "v_dubs": v_dubs,
        "old_to_local": old_to_local,
        "expected_repeated_pairs": expected_repeated_pairs,
        "surface_pair": [2590, 2596],
        "repeated_triple": [2590, 2596, 2598],
    }


# ---------------------------------------------------------------------------
# _minimum_distance_edge_pairing
# ---------------------------------------------------------------------------

def test_minimum_distance_edge_pairing_normal_two_vertex_case():
    """A normal 3-ball edge definition with two endpoints produces one edge."""
    entries = [(0, 10), (1, 11)]
    v_locs = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
    ]

    pairs = _minimum_distance_edge_pairing(entries, v_locs)

    assert _normalized_pairs(pairs) == {(0, 1)}


def test_minimum_distance_edge_pairing_rejects_same_fourth_ball():
    """Two duplicate four-ball definitions must not be connected as an edge."""
    entries = [(0, 10), (1, 10)]
    v_locs = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
    ]

    pairs = _minimum_distance_edge_pairing(entries, v_locs)

    assert pairs == []


def test_minimum_distance_edge_pairing_rejects_odd_endpoint_count():
    """Repeated regular-edge definitions require an even number of endpoints."""
    entries = [(0, 10), (1, 11), (2, 12)]
    v_locs = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
    ]

    pairs = _minimum_distance_edge_pairing(entries, v_locs)

    assert pairs == []


def test_aw_multiedge_pairing_creates_two_segments():
    """
    Four AW vertices sharing one 3-ball key must produce two edge segments,
    not collapse to a single edge.
    """
    case = _ala78_ca_c_regression_case()

    triple = set(case["repeated_triple"])
    entries = []

    for vert_ndx, balls in enumerate(case["v_balls"]):
        if triple.issubset(balls):
            fourth = next(ball for ball in balls if ball not in triple)
            entries.append((vert_ndx, fourth))

    assert len(entries) == 4

    pairs = _minimum_distance_edge_pairing(entries, case["v_locs"])

    assert len(pairs) == 2
    assert _normalized_pairs(pairs) == case["expected_repeated_pairs"]


def test_aw_multiedge_pairing_is_entry_order_invariant():
    """
    Reordering the endpoint entries must not change the geometric edge pairing.
    """
    case = _ala78_ca_c_regression_case()

    triple = set(case["repeated_triple"])
    entries = []
    for vert_ndx, balls in enumerate(case["v_balls"]):
        if triple.issubset(balls):
            fourth = next(ball for ball in balls if ball not in triple)
            entries.append((vert_ndx, fourth))

    expected = case["expected_repeated_pairs"]

    orders = [
        entries,
        list(reversed(entries)),
        [entries[2], entries[0], entries[3], entries[1]],
        [entries[1], entries[3], entries[0], entries[2]],
    ]

    for ordered_entries in orders:
        pairs = _minimum_distance_edge_pairing(ordered_entries, case["v_locs"])
        assert _normalized_pairs(pairs) == expected


# ---------------------------------------------------------------------------
# get_build_edges
# ---------------------------------------------------------------------------

def test_get_build_edges_preserves_multiple_segments_for_same_three_ball_key():
    """
    Regression for the original AW failure: one 3-ball key may legitimately
    correspond to multiple disconnected edge segments.
    """
    case = _ala78_ca_c_regression_case()

    v_balls = case["v_balls"]
    v_locs = case["v_locs"]
    v_dubs = case["v_dubs"]

    num_balls = max(max(v) for v in v_balls) + 1
    b_verts = _build_ball_vertex_index(v_balls, num_balls=num_balls)

    e_balls, e_verts = get_build_edges(
        b_verts=b_verts,
        v_balls=v_balls,
        v_locs=v_locs,
        v_dubs=v_dubs,
        start_time=time.time(),
        net=None,
    )

    repeated_key = tuple(case["repeated_triple"])
    repeated_segments = [
        e_verts[i]
        for i, balls in enumerate(e_balls)
        if tuple(balls) == repeated_key
    ]

    assert len(repeated_segments) == 2
    assert _normalized_pairs(repeated_segments) == case["expected_repeated_pairs"]


def test_get_build_edges_normal_three_ball_key_still_has_one_segment():
    """The AW fix must not duplicate ordinary two-endpoint edge definitions."""
    v_balls = [
        [0, 1, 2, 3],
        [0, 1, 2, 4],
    ]
    v_locs = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
    ]
    v_dubs = [0, 0]

    b_verts = _build_ball_vertex_index(v_balls, num_balls=5)

    e_balls, e_verts = get_build_edges(
        b_verts=b_verts,
        v_balls=v_balls,
        v_locs=v_locs,
        v_dubs=v_dubs,
        start_time=time.time(),
        net=None,
    )

    target_key = (0, 1, 2)
    segments = [
        e_verts[i]
        for i, balls in enumerate(e_balls)
        if tuple(balls) == target_key
    ]

    assert len(segments) == 1
    assert _normalized_pairs(segments) == {(0, 1)}


# ---------------------------------------------------------------------------
# Surface-closure regression
# ---------------------------------------------------------------------------

def test_ala78_ca_c_multiedge_surface_closes():
    """
    Full topology regression for the surface that originally disappeared.

    Before the fix, B/ALA78 CA(2590)-C(2596) contained 13 vertices but only
    12 edges because [2590, 2596, 2598] was collapsed to one edge segment.

    After the fix the surface must contain:
        13 vertices
        13 edges
        degree 2 at every boundary vertex

    This fixture intentionally contains only the local CA-C boundary topology,
    not the complete surrounding Streptavidin network.  We therefore construct
    the candidate as a two-group interface surface.  That exercises the same
    edge/surface closure machinery without triggering the normal group-network
    guard that rejects synthetic vertices having no additional surrounding
    edges.
    """
    case = _ala78_ca_c_regression_case()

    v_balls = case["v_balls"]
    v_locs = case["v_locs"]
    v_dubs = case["v_dubs"]

    num_balls = max(max(v) for v in v_balls) + 1
    b_verts = _build_ball_vertex_index(v_balls, num_balls=num_balls)

    e_balls, e_verts = get_build_edges(
        b_verts=b_verts,
        v_balls=v_balls,
        v_locs=v_locs,
        v_dubs=v_dubs,
        start_time=time.time(),
        net=None,
    )

    b_edges, v_edges = add_build_edges(
        num_balls=num_balls,
        e_balls=e_balls,
        num_verts=len(v_balls),
        e_verts=e_verts,
    )

    s_balls, s_verts, s_edges = get_build_surfs(
        b_verts=b_verts,
        b_edges=b_edges,
        v_balls=v_balls,
        v_edges=v_edges,
        e_balls=e_balls,
        start_time=time.time(),
        net=None,
        interface=True,
        iface_grps=({2590}, {2596}),
    )

    target_pair = case["surface_pair"]
    target_indices = [
        i for i, balls in enumerate(s_balls)
        if balls == target_pair
    ]

    assert len(target_indices) == 1

    surf_ndx = target_indices[0]
    surf_verts = s_verts[surf_ndx]
    surf_edges = s_edges[surf_ndx]

    assert len(surf_verts) == 13
    assert len(surf_edges) == 13

    surf_edge_set = set(surf_edges)
    degrees = {
        vert_ndx: sum(
            edge_ndx in surf_edge_set
            for edge_ndx in v_edges[vert_ndx]
        )
        for vert_ndx in surf_verts
    }

    assert set(degrees.values()) == {2}


def test_ala78_ca_c_repeated_triple_is_part_of_closed_surface():
    """
    Verify that both segments generated from the repeated key are actually
    retained in the reconstructed CA-C boundary.
    """
    case = _ala78_ca_c_regression_case()

    v_balls = case["v_balls"]
    v_locs = case["v_locs"]
    v_dubs = case["v_dubs"]

    num_balls = max(max(v) for v in v_balls) + 1
    b_verts = _build_ball_vertex_index(v_balls, num_balls=num_balls)

    e_balls, e_verts = get_build_edges(
        b_verts=b_verts,
        v_balls=v_balls,
        v_locs=v_locs,
        v_dubs=v_dubs,
        start_time=time.time(),
        net=None,
    )

    b_edges, v_edges = add_build_edges(
        num_balls=num_balls,
        e_balls=e_balls,
        num_verts=len(v_balls),
        e_verts=e_verts,
    )

    s_balls, s_verts, s_edges = get_build_surfs(
        b_verts=b_verts,
        b_edges=b_edges,
        v_balls=v_balls,
        v_edges=v_edges,
        e_balls=e_balls,
        start_time=time.time(),
        net=None,
        interface=True,
        iface_grps=({2590}, {2596}),
    )

    surf_ndx = s_balls.index(case["surface_pair"])
    surf_edge_indices = set(s_edges[surf_ndx])

    repeated_key = tuple(case["repeated_triple"])
    repeated_edge_indices = {
        i
        for i, balls in enumerate(e_balls)
        if tuple(balls) == repeated_key
    }

    assert len(repeated_edge_indices) == 2
    assert repeated_edge_indices.issubset(surf_edge_indices)
