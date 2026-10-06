from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex
from vorpy.src.geometry.interfaces import build_alpha_interface


def _network(scheme="pow", reverse_rows=False):
    ids = [10, 13, 16]
    rows = [
        {"loc": (0.0, 0.0, 0.0), "rad": 1.0, "complete": True,
         "chain_name": "A", "res_seq": 1, "res_name": "GLY", "name": "CA", "num": 101},
        {"loc": (2.0, 0.0, 0.0), "rad": 1.0, "complete": True,
         "chain_name": "B", "res_seq": 2, "res_name": "ALA", "name": "CB", "num": 102},
        {"loc": (0.0, 2.0, 0.0), "rad": 1.0, "complete": True,
         "chain_name": "B", "res_seq": 3, "res_name": "SER", "name": "OG", "num": 103},
    ]
    order = [2, 1, 0] if reverse_rows else [0, 1, 2]
    balls = pd.DataFrame([rows[i] for i in order], index=[ids[i] for i in order])
    return SimpleNamespace(
        settings={"net_type": scheme}, balls=balls,
        surfs=pd.DataFrame([
            {"balls": [10, 13], "edges": [7], "verts": [4, 5], "sa": 3.0},
            {"balls": [10, 16], "edges": [7], "verts": [4, 5], "sa": 4.0},
            {"balls": [13, 16], "edges": [7], "verts": [4, 5], "sa": 5.0},
        ], index=[20, 21, 22]),
        edges=pd.DataFrame([
            {"balls": [10, 13, 16], "verts": [4, 5], "surfs": [20, 21, 22], "length": 2.0},
        ], index=[7]),
        verts=pd.DataFrame([
            {"balls": [10, 13, 16, 19], "loc": (0.0, 0.0, 1.0), "surfs": [20, 21], "edges": [7]},
            {"balls": [10, 13, 16, 19], "loc": (0.0, 0.0, -1.0), "surfs": [20, 21], "edges": [7]},
        ], index=[4, 5]),
    )


def _filtration(network, scheme="pow", second_pair_birth=1.0, triangle=None):
    records = {dimension: {} for dimension in range(4)}
    for generator in (10, 13, 16):
        record = AlphaSimplex(0, (generator,), 0.0, 0.0, "test", True)
        records[0][(generator,)] = record
    for pair, birth in (((10, 13), 0.0), ((10, 16), second_pair_birth)):
        record = AlphaSimplex(1, pair, birth, birth, "test", True)
        records[1][pair] = record
    if triangle is not None:
        records[2][(10, 13, 16)] = AlphaSimplex(
            2, (10, 13, 16), triangle, triangle, "test", triangle is not None,
        )
    units = "A" if scheme == "aw" else "A^2"
    filtration = AlphaFiltration(
        network, scheme, records,
        {"native_alpha_units": units, "filtration_kind": "test filtration"},
    )
    return filtration


def test_alpha_mapping_selects_same_network_surface_and_restricts_incidence():
    network = _network()
    interface = build_alpha_interface(
        network, build_dual(network), _filtration(network), 0.0, {10}, {13, 16},
    )
    assert [(row.generator_ids, row.feature_id) for row in interface.selected_surfaces] == [((10, 13), 20)]
    assert interface.area_A2 == pytest.approx(3.0)
    excluded = next(row for row in interface.surfaces if row.feature_id == 21)
    assert excluded.selected is False
    assert excluded.status == "excluded_alpha"
    assert interface.selection_mappings[0]["mapped"] is True
    assert interface.selection_mappings[0]["surface_id"] == 20
    assert interface.selection_mappings[0]["stable_atom_i"] == "A|1|GLY|CA"
    assert interface.selected_surfaces[0].stable_atom_i == "A|1|GLY|CA"
    assert interface.edges[0].classification == "boundary_selected"
    assert interface.edges[0].reason.startswith("physical edge bounds one alpha-selected surface")
    assert interface.alpha_incidence_audit[0]["dual_triangle_status"] == "unresolved"
    assert interface.alpha_incidence_audit[0]["selected_surface_ids"] == (20,)
    assert {row.feature_id for row in interface.interface_vertices} == {4, 5}
    assert interface.curvature.smooth_H_A == 0.0


def test_alpha_selected_surface_subset_and_stable_ids_survive_row_reordering():
    forward = _network()
    reverse = _network(reverse_rows=True)
    results = []
    for network in (forward, reverse):
        result = build_alpha_interface(
            network, build_dual(network), _filtration(network), 0.0, {10}, {13, 16},
        )
        results.append(result)
        full_pairs = {row.generator_ids for row in result.candidate_surfaces}
        assert {row.generator_ids for row in result.selected_surfaces} <= full_pairs
        assert result.area_A2 == pytest.approx(sum(row.area_A2 for row in result.selected_surfaces))
    assert results[0].selection_mappings[0]["stable_atom_i"] == results[1].selection_mappings[0]["stable_atom_i"]
    assert results[0].selection_mappings[0]["surface_id"] == results[1].selection_mappings[0]["surface_id"]


def test_primitive_alpha_maps_to_planar_voronoi_surface():
    network = _network("prm")
    result = build_alpha_interface(
        network, build_dual(network), _filtration(network, "prm"),
        0.0, {10}, {13, 16},
    )
    assert result.scheme == "prm"
    assert len(result.selected_surfaces) == 1
    assert result.selected_surfaces[0].smooth_H_A == 0.0
    assert result.curvature.smooth_H_A == 0.0


def test_two_selected_surfaces_make_a_physical_interior_edge():
    network = _network()
    result = build_alpha_interface(
        network, build_dual(network),
        _filtration(network, second_pair_birth=-1.0, triangle=None),
        0.0, {10}, {13, 16},
    )
    assert len(result.selected_surfaces) == 2
    assert result.edges[0].classification == "interior_selected"
    assert result.edges[0].dual_triangle_present_at_alpha is False
    assert result.edges[0].dual_triangle_status == "unresolved"


def test_aw_alpha_smooth_term_comes_from_existing_aw_analyzer(monkeypatch):
    network = _network("aw")
    network.surfs["int_mean_curv_by_ball"] = [None, None, None]
    selected = SimpleNamespace(surface_id=20, integrated_mean_curvature=0.75)
    fake_result = SimpleNamespace(selected_surfaces=[selected], selected_edges=[])
    monkeypatch.setattr(
        "vorpy.src.analyze.aw_interface_curvature.analyze_aw_interface_curvature",
        lambda *args, **kwargs: fake_result,
    )
    interface = build_alpha_interface(
        network, build_dual(network), _filtration(network, "aw"), 0.0,
        {10}, {13, 16},
    )
    assert interface.curvature.smooth_H_A == pytest.approx(0.75)
    assert interface.curvature.raw_signed_edge_A_rad == pytest.approx(0.0)
    assert interface.curvature.combined_H_A == pytest.approx(0.75)


def test_aw_group_reversal_preserves_area_unsigned_and_flips_signed(monkeypatch):
    network = _network("aw")

    def fake_analyzer(_network, group_a, group_b, **_kwargs):
        sign = 1.0 if 10 in group_a else -1.0
        return SimpleNamespace(
            selected_surfaces=[SimpleNamespace(
                surface_id=20, integrated_mean_curvature=sign * 0.75,
            )],
            selected_edges=[],
        )

    monkeypatch.setattr(
        "vorpy.src.analyze.aw_interface_curvature.analyze_aw_interface_curvature",
        fake_analyzer,
    )
    dual = build_dual(network)
    forward = build_alpha_interface(
        network, dual, _filtration(network, "aw"), 0.0, {10}, {13, 16},
    )
    reverse = build_alpha_interface(
        network, dual, _filtration(network, "aw"), 0.0, {13, 16}, {10},
    )
    assert forward.area_A2 == reverse.area_A2 == pytest.approx(3.0)
    assert forward.curvature.raw_unsigned_edge_A_rad == reverse.curvature.raw_unsigned_edge_A_rad
    assert forward.curvature.combined_H_A == pytest.approx(-reverse.curvature.combined_H_A)


def test_alpha_mapping_export_includes_required_audit_files(tmp_path):
    network = _network()
    interface = build_alpha_interface(
        network, build_dual(network), _filtration(network), 0.0, {10}, {13, 16},
    )
    output = interface.export(tmp_path)
    for name in (
        "interface_summary.json", "interface_surfaces.csv", "interface_edges.csv",
        "curvature_summary.csv", "alpha_incidence_audit.csv",
        "alpha_surface_mapping.csv",
    ):
        assert (output / name).is_file()
