from types import SimpleNamespace

import json
import numpy as np
import pandas as pd
import pytest

from vorpy.src.network.network import Network
from vorpy.src.network.perimeter import build_perimeter, canonicalize_edge_orientation
from vorpy.src.network.aw_geometry_reuse import (
    AWGeometryReuseError,
    AWReuseDecision,
    assess_aw_reuse,
    build_aw_interface_view,
    build_solve_universe_key,
    certify_aw_interface_completeness,
    collect_aw_solver_provenance,
    compare_aw_interface_view,
    compare_solved_aw_networks,
    diagnose_aw_interface_difference,
    extract_aw_geometry,
)


def _settings():
    return {
        "net_type": "aw",
        "surf_res": 0.25,
        "max_vert": 5.0,
        "box_size": 1.5,
        "num_splits": 2,
        "build_type": "all",
        "foam_box": None,
    }


def _network(vertex_balls, *, interface=None, radii=None):
    locations = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
    ]
    if radii is None:
        radii = [1.0] * len(locations)
    balls = pd.DataFrame({
        "loc": locations,
        "rad": radii,
        "complete": [True] * len(locations),
        "stable_id": [f"atom-{index}" for index in range(len(locations))],
    })
    verts = pd.DataFrame({
        "balls": [list(value) for value in vertex_balls],
        "loc": [(0.1 + index, 0.2, 0.3) for index in range(len(vertex_balls))],
        "rad": [2.0] * len(vertex_balls),
        "dub": [0] * len(vertex_balls),
    })
    system = SimpleNamespace(
        frame_count=0,
        update_progress=lambda **kwargs: None,
    )
    if interface is None:
        group = list(range(len(locations)))
        iface_grps = None
    else:
        group = sorted(set(interface[0]) | set(interface[1]))
        iface_grps = interface
    network = Network(
        locs=locations,
        rads=radii,
        group=group,
        iface_grps=iface_grps,
        group_name="synthetic",
        settings=_settings(),
        balls=balls,
        verts=verts,
        system=system,
    )
    network.metrics["vert"] = 0.0
    network.connect()
    return network


def _local_surface_network(*, incomplete=(), interface=False):
    """Two curved, four-sided interface patches sharing a junction edge."""
    locations = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0), (0.0, 1.0, 1.0), (0.0, 0.0, 1.0),
        (0.5, 0.5, 0.2), (100.0, 100.0, 100.0),
    ]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [1.0] * len(locations),
        "complete": [index not in set(incomplete) for index in range(len(locations))],
        "stable_id": [f"local-atom-{index}" for index in range(len(locations))],
    })
    vertex_locations = locations[:6]
    verts = pd.DataFrame({
        "balls": [[0, 1, 2, 3]] * 6,
        "loc": vertex_locations,
        "rad": [0.2] * 6,
        "dub": [0] * 6,
        "edges": [
            [0, 3, 6], [0, 1], [1, 2],
            [2, 3, 4], [4, 5], [5, 6],
        ],
        "surfs": [[0, 1], [0], [0], [0, 1], [1], [1]],
    })
    edge_vertices = [[0, 1], [1, 2], [2, 3], [3, 0], [3, 4], [4, 5], [5, 0]]
    edge_balls = [
        [0, 1, 2], [1, 2, 3], [0, 2, 3], [0, 1, 3],
        [0, 2, 3], [0, 1, 2], [0, 1, 3],
    ]
    edge_points = [
        [(0.0, 0.0, 0.0), (0.5, -0.05, 0.0), (1.0, 0.0, 0.0)],
        [(1.0, 0.0, 0.0), (1.05, 0.5, 0.0), (1.0, 1.0, 0.0)],
        [(1.0, 1.0, 0.0), (0.5, 1.05, 0.0), (0.0, 1.0, 0.0)],
        [(0.0, 1.0, 0.0), (-0.05, 0.5, 0.0), (0.0, 0.0, 0.0)],
        [(0.0, 1.0, 0.0), (0.0, 1.0, 0.5), (0.0, 1.0, 1.0)],
        [(0.0, 1.0, 1.0), (0.0, 0.5, 1.05), (0.0, 0.0, 1.0)],
        [(0.0, 0.0, 1.0), (0.0, -0.05, 0.5), (0.0, 0.0, 0.0)],
    ]
    edges = pd.DataFrame({
        "balls": edge_balls,
        "verts": edge_vertices,
        "surfs": [[0], [0], [0], [0, 1], [1], [1], [1]],
        "points": edge_points,
        "length": [1.0] * 7,
    })
    surfaces = pd.DataFrame({
        "balls": [[0, 2], [1, 3]],
        "verts": [[0, 1, 2, 3], [3, 4, 5, 0]],
        "edges": [[0, 1, 2, 3], [3, 4, 5, 6]],
        "points": [
            [locations[index] for index in (0, 1, 2, 3)],
            [locations[index] for index in (3, 4, 5, 0)],
        ],
        "tris": [[(0, 1, 2), (0, 2, 3)], [(0, 1, 2), (0, 2, 3)]],
        "sa": [1.0, 1.0],
        "orientation": [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
    })
    system = SimpleNamespace(frame_count=0, update_progress=lambda **kwargs: None)
    return Network(
        locs=locations,
        rads=[1.0] * len(locations),
        group=[0, 1, 2, 3] if interface else list(range(len(locations))),
        iface_grps=({0, 1}, {2, 3}) if interface else None,
        group_name="local_surface",
        settings=_settings(),
        balls=balls,
        verts=verts,
        edges=edges,
        surfs=surfaces,
        system=system,
    )


def _production_box_network(
        *, interface=None, max_vert=40.0, radii=None, calculate_curvature=False
):
    """Small AW system solved through Network.build, not assembled tables."""
    locations = np.array([
        [0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [-20.0, 0.0, 0.0],
        [0.0, 20.0, 0.0], [0.0, -20.0, 0.0], [0.0, 0.0, 20.0],
        [0.0, 0.0, -20.0],
    ])
    radii = [1.0] * len(locations) if radii is None else list(radii)
    settings = _settings()
    settings.update({
        "max_vert": max_vert,
        "surf_scheme": "int_mean_curv",
        "aw_provenance": True,
    })
    if interface is None:
        group = list(range(len(locations)))
        iface_grps = None
    else:
        iface_grps = tuple(set(side) for side in interface)
        group = sorted(set().union(*iface_grps))
    network = Network(
        locs=locations.tolist(),
        rads=radii,
        group=group,
        iface_grps=iface_grps,
        group_name="production_box",
        settings=settings,
        system=None,
    )
    network.build(calculate_curvature=calculate_curvature)
    return network


def test_aw_perimeter_is_invariant_to_edge_order_and_direction():
    corners = {
        "a": np.array([-1.0, -1.0, 0.0]),
        "b": np.array([1.0, -1.0, 0.0]),
        "c": np.array([1.0, 1.0, 0.0]),
        "d": np.array([-1.0, 1.0, 0.0]),
    }
    edges = [
        ([corners["a"], corners["b"]], ("a", "b")),
        ([corners["b"], corners["c"]], ("b", "c")),
        ([corners["c"], corners["d"]], ("c", "d")),
        ([corners["d"], corners["a"]], ("d", "a")),
    ]
    shuffled = [
        (list(reversed(points)), tuple(reversed(endpoints)))
        for points, endpoints in reversed(edges)
    ]

    first = build_perimeter(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [1.0, 1.0],
        [points for points, _ in edges], edge_endpoints=[endpoints for _, endpoints in edges],
    )
    second = build_perimeter(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [1.0, 1.0],
        [points for points, _ in shuffled], edge_endpoints=[endpoints for _, endpoints in shuffled],
    )

    assert np.array_equal(np.asarray(first[0]), np.asarray(second[0]))
    assert np.array_equal(first[1], second[1])
    assert np.array_equal(first[2], second[2])


def test_aw_edge_orientation_reverses_directional_metadata():
    points = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
    metadata = {
        "loc": np.array([2.0, 0.0, 0.0]),
        "loc2": np.array([3.0, 0.0, 0.0]),
        "rad": 1.0,
        "rad2": 2.0,
        "vnorm": np.array([1.0, 0.0, 0.0]),
        "dnorm": np.array([0.0, 1.0, 0.0]),
    }

    ordered, reversed_metadata, was_reversed = canonicalize_edge_orientation(
        points, ("z", "a"), metadata
    )

    assert was_reversed is True
    assert np.array_equal(ordered, list(reversed(points)))
    assert np.array_equal(reversed_metadata["vnorm"], [-1.0, 0.0, 0.0])
    assert np.array_equal(reversed_metadata["dnorm"], [0.0, -1.0, 0.0])
    assert np.array_equal(reversed_metadata["loc"], metadata["loc2"])
    assert reversed_metadata["rad"] == metadata["rad2"]


def test_aw_perimeter_rejects_multiple_boundary_cycles_explicitly():
    square_a = [
        ([[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0]], ("a", "b")),
        ([[1.0, -1.0, 0.0], [1.0, 1.0, 0.0]], ("b", "c")),
        ([[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]], ("c", "d")),
        ([[-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0]], ("d", "a")),
    ]
    square_b = [
        ([[3.0, 3.0, 0.0], [4.0, 3.0, 0.0]], ("e", "f")),
        ([[4.0, 3.0, 0.0], [4.0, 4.0, 0.0]], ("f", "g")),
        ([[4.0, 4.0, 0.0], [3.0, 4.0, 0.0]], ("g", "h")),
        ([[3.0, 4.0, 0.0], [3.0, 3.0, 0.0]], ("h", "e")),
    ]
    edges = square_a + square_b
    with pytest.raises(ValueError, match="Multiple boundary cycles"):
        build_perimeter(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [1.0, 1.0],
            [[np.asarray(point) for point in points] for points, _ in edges],
            edge_endpoints=[endpoints for _, endpoints in edges],
        )


def test_full_aw_extraction_materializes_interface_without_mutating_source():
    full = _network(((0, 1, 2, 4), (0, 1, 2, 5)))
    before = {
        name: frame.copy(deep=True)
        for name, frame in (
            ("balls", full.balls),
            ("verts", full.verts),
            ("edges", full.edges),
            ("surfs", full.surfs),
        )
    }

    store = extract_aw_geometry(full, measure_memory=True)
    view = build_aw_interface_view(store, {0, 1}, {2, 3})
    materialized = view.materialize(measure_memory=True)

    assert len(materialized.verts) == 2
    assert len(materialized.edges) == len(full.edges)
    assert len(materialized.surfs) == len(full.surfs)
    assert materialized.group_a == (0, 1)
    assert materialized.group_b == (2, 3)
    assert store.metrics.seconds >= 0.0
    assert materialized.metrics.seconds >= 0.0
    assert store.metrics.peak_tracemalloc_bytes is not None
    assert materialized.metrics.peak_tracemalloc_bytes is not None

    for name, frame in before.items():
        assert getattr(full, name).equals(frame)


def test_globally_incomplete_but_locally_complete_interface_can_be_certified():
    full = _local_surface_network(incomplete=(7,))
    direct = _local_surface_network(interface=True)
    store = extract_aw_geometry(full, allow_incomplete=True)
    view = build_aw_interface_view(store, {0, 1}, {2, 3})

    certificate = certify_aw_interface_completeness(view)
    assessment = assess_aw_reuse(view, direct_network=direct)

    assert not certificate.globally_complete
    assert certificate.proven
    assert 7 not in certificate.relevant_generators
    assert certificate.boundary_edges == 0
    assert assessment.decision is AWReuseDecision.EXACT_REUSE
    assert assessment.comparison.equivalent
    assert assessment.comparison.counts == {
        "vertices": {"extracted": 6, "direct": 6},
        "edges": {"extracted": 7, "direct": 7},
        "surfaces": {"extracted": 2, "direct": 2},
    }


def test_missing_required_interface_vertex_rejects_topology():
    full = _local_surface_network()
    direct = _local_surface_network(interface=True)
    direct.verts = direct.verts.iloc[:-1].reset_index(drop=True)
    direct.edges = direct.edges.iloc[0:0].copy()
    direct.surfs = direct.surfs.iloc[0:0].copy()
    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})

    assessment = assess_aw_reuse(view, direct_network=direct)

    assert assessment.decision is AWReuseDecision.REJECT_TOPOLOGY
    assert "vertices" in assessment.reasons


def test_missing_required_interface_edge_rejects_topology():
    full = _local_surface_network()
    direct = _local_surface_network(interface=True)
    direct.edges = direct.edges.iloc[0:0].copy()
    direct.surfs = direct.surfs.iloc[0:0].copy()
    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})

    assessment = assess_aw_reuse(view, direct_network=direct)

    assert assessment.decision is AWReuseDecision.REJECT_TOPOLOGY
    assert "edges" in assessment.reasons


def test_open_surface_boundary_rejects_incomplete_reuse():
    full = _local_surface_network()
    full.surfs.at[0, "edges"] = [0, 1, 2]
    store = extract_aw_geometry(full, allow_incomplete=True)
    view = build_aw_interface_view(store, {0, 1}, {2, 3})

    assessment = assess_aw_reuse(view)

    assert assessment.decision is AWReuseDecision.REJECT_INCOMPLETE
    assert "open_interface_surface_boundary" in assessment.reasons
    assert assessment.certificate.boundary_edges > 0


def test_incomplete_relevant_competitor_rejects_incomplete_reuse():
    full = _local_surface_network(incomplete=(6,))
    store = extract_aw_geometry(full, allow_incomplete=True)
    view = build_aw_interface_view(store, {0, 1}, {2, 3})

    assessment = assess_aw_reuse(view)

    assert assessment.decision is AWReuseDecision.REJECT_INCOMPLETE
    assert "incomplete_relevant_competitor" in assessment.reasons
    assert 6 in assessment.certificate.relevant_generators


def test_changed_environment_rejects_reuse_explicitly():
    full = _local_surface_network()
    direct = _local_surface_network(interface=True)
    direct.balls.loc[6, "rad"] = 1.25
    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})

    assessment = assess_aw_reuse(view, direct_network=direct)

    assert assessment.decision is AWReuseDecision.REJECT_ENVIRONMENT
    assert "competitor_universe" in assessment.reasons


def test_equivalence_without_direct_builder_remains_unresolved():
    full = _local_surface_network()
    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})

    assessment = assess_aw_reuse(view)

    assert assessment.decision is AWReuseDecision.UNRESOLVED
    assert assessment.reasons == ("direct_interface_equivalence_not_checked",)


def test_materialization_preserves_surface_mesh_orientation_and_incidence():
    locations = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [1.0] * 4,
        "complete": [True] * 4,
        "stable_id": [f"surface-atom-{index}" for index in range(4)],
    })
    verts = pd.DataFrame({
        "balls": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "loc": [(0.25, 0.25, 0.25), (0.25, 0.25, -0.25)],
        "rad": [2.0, 2.0],
        "dub": [0, 0],
        "edges": [[0, 1], [0, 1]],
        "surfs": [[0], [0]],
    })
    edges = pd.DataFrame({
        "balls": [[0, 1, 2], [0, 2, 3]],
        "verts": [[0, 1], [0, 1]],
        "surfs": [[0], [0]],
        "points": [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]],
        "length": [1.0, 1.0],
    })
    surfaces = pd.DataFrame({
        "balls": [[0, 2]],
        "verts": [[0, 1]],
        "edges": [[0, 1]],
        "points": [[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]],
        "tris": [[(0, 1, 0)]],
        "sa": [3.5],
        "orientation": [(1.0, 0.0, 0.0)],
    })
    system = SimpleNamespace(frame_count=0, update_progress=lambda **kwargs: None)
    full = Network(
        locs=locations,
        rads=[1.0] * 4,
        group=list(range(4)),
        settings=_settings(),
        balls=balls,
        verts=verts,
        edges=edges,
        surfs=surfaces,
        system=system,
    )
    direct = Network(
        locs=locations,
        rads=[1.0] * 4,
        group=[0, 1, 2, 3],
        iface_grps=({0, 1}, {2, 3}),
        settings=_settings(),
        balls=balls.copy(deep=True),
        verts=verts.copy(deep=True),
        edges=edges.copy(deep=True),
        surfs=surfaces.copy(deep=True),
        system=system,
    )

    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})
    materialized = view.materialize()
    comparison = compare_aw_interface_view(view, direct)

    assert comparison.equivalent
    assert len(materialized.surfs) == 1
    assert materialized.surfs.iloc[0]["sa"] == 3.5
    assert tuple(materialized.surfs.iloc[0]["orientation"]) == (1.0, 0.0, 0.0)
    assert materialized.surfs.iloc[0]["verts"] == [0, 1]
    assert materialized.surfs.iloc[0]["edges"] == [0, 1]


def test_materialization_avoids_a_second_network_solve(monkeypatch):
    full = _network(((0, 1, 2, 4), (0, 1, 2, 5)))
    view = build_aw_interface_view(
        extract_aw_geometry(full),
        {0, 1},
        {2, 3},
    )

    def forbidden_build(*args, **kwargs):
        raise AssertionError("interface materialization must not call Network.build")

    monkeypatch.setattr(Network, "build", forbidden_build)
    materialized = view.materialize()

    assert len(materialized.verts) == 2


def test_extracted_interface_matches_direct_builder_when_topology_context_matches():
    vertices = ((0, 1, 2, 4), (0, 1, 2, 5))
    full = _network(vertices)
    direct = _network(vertices, interface=({0, 1}, {2, 3}))

    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})
    comparison = compare_aw_interface_view(view, direct)

    assert comparison.equivalent
    assert comparison.same_competitor_universe
    assert comparison.compatible_geometry_context
    assert comparison.source_scope == "full"
    assert comparison.direct_scope == "interface"
    assert comparison.mismatches == {}
    assert comparison.counts["vertices"] == {"extracted": 2, "direct": 2}


def test_active_non_interface_competitor_is_reported_as_topology_mismatch():
    # This competitor is part of the full solve universe.  Changing its
    # radius in the direct network must invalidate reuse even when the
    # selected interface tables happen to have the same topology.
    vertices = (
        (0, 1, 2, 4),
        (0, 1, 2, 5),
        (2, 4, 5, 6),
    )
    full = _network(vertices)
    changed_radii = [1.0] * 7
    changed_radii[6] = 1.25
    direct = _network(
        vertices[:2],
        interface=({0, 1}, {2, 3}),
        radii=changed_radii,
    )

    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})
    comparison = compare_aw_interface_view(view, direct)

    assert not comparison.same_competitor_universe
    assert not comparison.compatible_geometry_context
    assert not comparison.equivalent
    assert "competitor_universe" in comparison.mismatches


def test_interface_comparison_reports_topology_mismatch_without_reuse():
    vertices = ((0, 1, 2, 4), (0, 1, 2, 5))
    full = _network(vertices)
    direct = _network(vertices, interface=({0, 1}, {2, 3}))
    direct.edges = direct.edges.iloc[0:0].copy()

    view = build_aw_interface_view(extract_aw_geometry(full), {0, 1}, {2, 3})
    comparison = compare_aw_interface_view(view, direct)

    assert comparison.same_competitor_universe
    assert comparison.compatible_geometry_context
    assert not comparison.equivalent
    assert "edges" in comparison.mismatches


def test_solve_keys_are_deterministic_and_scope_is_part_of_identity():
    full = _network(((0, 2, 4, 5), (1, 2, 4, 5)))
    direct = _network(((0, 2, 4, 5), (1, 2, 4, 5)), interface=({0, 1}, {2, 3}))

    full_key_a = build_solve_universe_key(full, frame_id=4)
    full_key_b = build_solve_universe_key(full, frame_id=4)
    direct_key = build_solve_universe_key(direct, frame_id=4)
    changed_frame_key = build_solve_universe_key(full, frame_id=5)
    full.settings["surf_res"] = 0.5
    changed_resolution_key = build_solve_universe_key(full, frame_id=4)

    assert full_key_a == full_key_b
    assert full_key_a.digest != direct_key.digest
    assert full_key_a.competitor_digest == direct_key.competitor_digest
    assert full_key_a.geometry_digest == direct_key.geometry_digest
    assert full_key_a.frame_id == "4"
    assert changed_frame_key.digest != full_key_a.digest
    assert changed_resolution_key.digest != full_key_a.digest

    full.settings["surf_res"] = 0.25
    full_store = extract_aw_geometry(full)
    direct_store = extract_aw_geometry(direct)
    assert [record.primitive_id for record in full_store.vertices] == [
        record.primitive_id for record in direct_store.vertices
    ]


def test_incomplete_source_fails_closed():
    network = _network(((0, 2, 4, 5), (1, 2, 4, 5)))
    network.balls.loc[6, "complete"] = False

    with pytest.raises(AWGeometryReuseError, match="incomplete"):
        extract_aw_geometry(network)


def test_non_aw_source_fails_closed():
    network = _network(((0, 2, 4, 5), (1, 2, 4, 5)))
    network.settings["net_type"] = "pow"

    with pytest.raises(AWGeometryReuseError, match="non-AW"):
        extract_aw_geometry(network)


def test_production_builder_provenance_is_opt_in_and_fail_closed(monkeypatch):
    monkeypatch.setattr("vorpy.src.network.find_net_verts.write_verts", lambda net: None)
    network = _production_box_network()

    payload = network.aw_solver_provenance
    assert payload["enabled"] is True
    assert len(payload["canonical_generators"]) == 7
    assert payload["vertex_search"]["status"] in {
        "traversal_exhausted", "truncated", "no_seed", "unobserved",
    }
    assert payload["completeness"]["completeness_proven"] is False
    assert payload["edge_construction"]["completeness_proven"] is False
    assert payload["surface_construction"]["completeness_proven"] is False
    assert json.dumps(payload, allow_nan=False)


def test_independently_solved_production_networks_are_compared_observationally(monkeypatch):
    monkeypatch.setattr("vorpy.src.network.find_net_verts.write_verts", lambda net: None)
    full = _production_box_network()
    direct = _production_box_network(interface=({1}, {0, 2, 3, 4, 5, 6}))
    before = {
        name: frame.copy(deep=True)
        for name, frame in (("verts", full.verts), ("edges", full.edges), ("surfs", full.surfs))
    }

    result = compare_solved_aw_networks(full, direct)

    assert result.assessment.decision in set(AWReuseDecision)
    assert result.assessment.comparison is not None
    assert result.assessment.comparison.counts["vertices"]["extracted"] == len(result.view.vertex_ids)
    assert result.assessment.comparison.counts == {
        "vertices": {"extracted": 4, "direct": 4},
        "edges": {"extracted": 4, "direct": 4},
        "surfaces": {"extracted": 1, "direct": 1},
    }
    assert set(result.assessment.comparison.mismatches) == {
        "edges", "surfaces",
    }
    assert result.assessment.comparison.equivalent is False
    assert result.assessment.decision is AWReuseDecision.REJECT_TOPOLOGY
    assert result.full_provenance.completeness["completeness_proven"] is False
    for name, frame in before.items():
        assert getattr(full, name).equals(frame)


def test_production_first_difference_report_preserves_exact_evidence(monkeypatch):
    monkeypatch.setattr("vorpy.src.network.find_net_verts.write_verts", lambda net: None)
    full = _production_box_network()
    direct = _production_box_network(interface=({1}, {0, 2, 3, 4, 5, 6}))
    view = build_aw_interface_view(
        extract_aw_geometry(full, allow_incomplete=True),
        {1},
        {0, 2, 3, 4, 5, 6},
    )

    report = diagnose_aw_interface_difference(view, direct)
    extracted_surface = report["stages"]["surfaces"]["extracted_records"][0]
    direct_surface = report["stages"]["surfaces"]["direct_records"][0]
    assert np.max(np.abs(
        np.asarray(extracted_surface["points"], dtype=float)
        - np.asarray(direct_surface["points"], dtype=float)
    )) <= 1.0e-10
    assert np.array_equal(extracted_surface["tris"], direct_surface["tris"])
    assert extracted_surface["sa"] == direct_surface["sa"]
    assert report["first_divergence"] == "edge_sample_roundoff"
    assert report["first_physical_geometry_difference"] == "none"
    assert report["same_competitor_universe"] is True
    assert report["compatible_geometry_context"] is True
    assert report["stages"]["vertices"]["canonical_equivalent"] is True
    assert report["stages"]["edges"]["canonical_equivalent"] is False
    assert report["stages"]["edges"]["topology_equivalent"] is True
    assert report["stages"]["edges"]["geometry_equivalent"] is False
    assert report["stages"]["edges"]["incidence_equivalent"] is True
    assert report["stages"]["edges"]["physical_geometry_equivalent"] is False
    assert report["stages"]["edges"]["physical_geometry_within_tolerance"] is True
    assert report["stages"]["edges"]["physical_incidence_equivalent"] is True
    assert report["stages"]["edges"]["representation_only"] is True
    assert report["stages"]["surfaces"]["canonical_equivalent"] is False
    assert report["stages"]["surfaces"]["topology_equivalent"] is True
    assert report["stages"]["surfaces"]["geometry_equivalent"] is True
    assert report["stages"]["surfaces"]["incidence_equivalent"] is False
    assert report["stages"]["surfaces"]["physical_geometry_equivalent"] is True
    assert report["stages"]["surfaces"]["physical_geometry_within_tolerance"] is True
    assert report["stages"]["surfaces"]["physical_incidence_equivalent"] is True
    assert report["stages"]["surfaces"]["representation_only"] is True
    assert report["stages"]["edges"]["extracted_count"] == 4
    assert report["stages"]["edges"]["direct_count"] == 4
    assert report["stages"]["surfaces"]["extracted_count"] == 1
    assert report["stages"]["surfaces"]["direct_count"] == 1
    assert len(report["stages"]["vertices"]["extracted_records"]) == 4
    assert {"identity", "balls", "loc", "rad"} <= set(
        report["stages"]["vertices"]["extracted_records"][0]
    )
    assert {"identity", "balls", "verts", "points", "vals"} <= set(
        report["stages"]["edges"]["extracted_records"][0]
    )
    assert {"identity", "balls", "verts", "edges", "points", "tris", "orientation"} <= set(
        report["stages"]["surfaces"]["direct_records"][0]
    )
    assert len(report["generator_identities"]["extracted"]) == 7
    assert report["generator_identities"]["extracted"] == report["generator_identities"]["direct"]
    assert report["discarded_candidates"]["direct_solver_provenance"]
    assert json.dumps(report, sort_keys=True, allow_nan=False)


def test_production_aw_mesh_and_curvature_match_when_enabled(monkeypatch):
    monkeypatch.setattr("vorpy.src.network.find_net_verts.write_verts", lambda net: None)
    full = _production_box_network(calculate_curvature=True)
    direct = _production_box_network(
        interface=({1}, {0, 2, 3, 4, 5, 6}), calculate_curvature=True
    )
    full_surface = next(
        row for _, row in full.surfs.iterrows()
        if set(row["balls"]) == {0, 1}
    )
    direct_surface = direct.surfs.iloc[0]

    assert np.array_equal(full_surface["tris"], direct_surface["tris"])
    assert np.array_equal(
        np.asarray(full_surface["points"], dtype=float),
        np.asarray(direct_surface["points"], dtype=float),
    )
    for column in (
        "sa", "mean_curv", "avg_mean_curv", "gauss_curv", "avg_gauss_curv",
        "int_mean_curv", "int_mean_curv_sq", "int_gauss_curv", "surf_energy",
    ):
        assert np.allclose(
            np.asarray(full_surface[column], dtype=float),
            np.asarray(direct_surface[column], dtype=float),
        )


def test_production_active_competitor_and_incomplete_search_never_certify_reuse(monkeypatch):
    monkeypatch.setattr("vorpy.src.network.find_net_verts.write_verts", lambda net: None)
    full = _production_box_network()
    changed = _production_box_network(
        interface=({1}, {0, 2, 3, 4, 5, 6}),
        radii=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0],
    )
    changed_result = compare_solved_aw_networks(full, changed)
    assert changed_result.assessment.comparison is not None
    assert not changed_result.assessment.comparison.same_competitor_universe
    assert changed_result.assessment.decision is AWReuseDecision.REJECT_ENVIRONMENT

    incomplete = _production_box_network(
        interface=({1}, {0, 2, 3, 4, 5, 6}),
        max_vert=0.1,
    )
    incomplete_result = compare_solved_aw_networks(full, incomplete)
    assert incomplete_result.interface_provenance.vertex_search["completeness_proven"] is False
    assert incomplete_result.interface_provenance.vertex_search["status"] == "no_seed"
    assert incomplete_result.assessment.decision is AWReuseDecision.REJECT_ENVIRONMENT
    assert "vertices" in incomplete_result.assessment.reasons
