from itertools import combinations
from types import SimpleNamespace

import pandas as pd

from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.interfaces import (
    aw_curvature_summary,
    build_full_interface,
    power_curvature_summary,
)


def _tetra_network(scheme="prm"):
    pairs = list(combinations(range(4), 2))
    triples = list(combinations(range(4), 3))
    balls = pd.DataFrame({
        "loc": [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)],
        "rad": [1., 1., 1., 1.],
        "complete": [True] * 4,
    })
    pair_edges = {
        pair: [i for i, tri in enumerate(triples) if set(pair).issubset(tri)]
        for pair in pairs
    }
    tri_surfaces = {
        tri: [i for i, pair in enumerate(pairs) if set(pair).issubset(tri)]
        for tri in triples
    }
    surfaces = pd.DataFrame({
        "balls": [list(pair) for pair in pairs],
        "edges": [pair_edges[pair] for pair in pairs],
        "verts": [[0, 1] for _ in pairs],
        "sa": [1.0] * len(pairs),
    })
    edges = pd.DataFrame({
        "balls": [list(tri) for tri in triples],
        "verts": [[0, 1] for _ in triples],
        "surfs": [tri_surfaces[tri] for tri in triples],
    })
    verts = pd.DataFrame({
        "balls": [list(range(4)), list(range(4))],
        "loc": [(0.5, 0.5, 0.5), (0.5, 0.5, -0.5)],
        "edges": [list(range(4)), list(range(4))],
        "surfs": [list(range(6)), list(range(6))],
    })
    return SimpleNamespace(settings={"net_type": scheme}, balls=balls,
                           surfs=surfaces, edges=edges, verts=verts,
                           group_name="tetra")


def test_full_bicolor_interface_maps_surface_geometry_and_topology_for_each_scheme(tmp_path):
    for scheme in ("prm", "pow", "aw"):
        network = _tetra_network(scheme)
        dual = build_dual(network)
        interface = build_full_interface(network, dual, {0, 1}, {2, 3})

        assert interface.scheme == scheme
        assert interface.selection_mode == "full"
        assert interface.alpha is None
        assert len(interface.candidate_surfaces) == 4
        assert len(interface.selected_surfaces) == 4
        assert interface.interface_atoms == frozenset(range(4))
        assert interface.area_A2 == 4.0
        assert interface.component_count == 1
        assert {edge.incidence_status for edge in interface.edges} == {"shared"}
        assert len(interface.interface_edges) == 4
        assert len(interface.interface_vertices) == 2
        assert interface.euler_characteristic == 2
        assert interface.metadata["alpha_or_M_filter_applied"] is False

        interface.export(tmp_path / scheme)
        assert (tmp_path / scheme / "interface_surfaces.csv").is_file()
        assert (tmp_path / scheme / "interface_summary.json").is_file()


def test_existing_curvature_results_map_to_common_schema_without_recalculation():
    aw_result = SimpleNamespace(
        selection_mode="raw", interface_area=100.0, surface_curvature=2.0,
        edge_curvature={"signed": -4.0, "unsigned": 10.0,
                        "positive": 3.0, "negative": -7.0},
    )
    aw = aw_curvature_summary(aw_result)
    assert aw.selection_mode == "full"
    assert aw.smooth_H_A == 2.0
    assert aw.raw_signed_edge_A_rad == -4.0
    assert aw.conventional_edge_H_A == -2.0
    assert aw.combined_H_A == 0.0
    assert aw.H_per_area_inv_A == 0.0
    assert aw.orientation.startswith("interface normal A-to-B")

    power_result = SimpleNamespace(
        edge_totals={"signed": -6.0, "unsigned": 12.0,
                     "positive": 3.0, "negative": -9.0},
        smooth_surface_curvature=0.0,
    )
    power = power_curvature_summary(power_result, area_A2=20.0)
    assert power.selection_mode == "cazals_alpha0_power"
    assert power.smooth_H_A == 0.0
    assert power.raw_signed_edge_A_rad == -6.0
    assert power.conventional_edge_H_A == -3.0
    assert power.combined_H_A == -3.0
    assert power.cancellation_fraction == 0.5
