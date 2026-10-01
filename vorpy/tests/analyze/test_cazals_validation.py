from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest

from vorpy.src.analyze.cazals_validation import (
    classify_same_group_pair,
    compnd_chain_mapping,
    facet_area_statistics,
    interface_atom_asymmetry,
    interface_atom_ids,
    interface_components,
    power_facet_polygon_area,
    main as cazals_validation_main,
)
from vorpy.src.analyze.cazals_edge_forensics import (
    audit_final_interface_polygons,
    classify_interface_edge,
    condition_beta_audit_rows,
    select_interior_finite_edges,
    summarize_edge_records,
)
from vorpy.src.analyze.power_interface_curvature import (
    PowerVertex,
    cazals_signed_beta,
)


def test_cazals_validation_entrypoint_has_explicit_system_choices():
    with pytest.raises(SystemExit) as error:
        cazals_validation_main(["--system", "2kai-typo"])
    assert error.value.code == 2


def test_selected_facet_atom_sets_are_unique_and_ratio_is_symmetric():
    facets = [SimpleNamespace(generator_ids=(0, 4)),
              SimpleNamespace(generator_ids=(1, 4)),
              SimpleNamespace(generator_ids=(0, 5))]
    a, b = interface_atom_ids(facets, {0, 1, 2})
    assert a == {0, 1}
    assert b == {4, 5}
    assert interface_atom_asymmetry(a, b) == pytest.approx(1.0)
    assert interface_atom_asymmetry({0, 1, 2, 3}, {4, 5}) == pytest.approx(2.0)


def test_bounded_power_facet_polygon_area_and_area_summary():
    points = np.array([[0., 0., -1.], [0., 0., 1.],
                       [1., 0., 0.], [0., 1., 0.],
                       [-1., 0., 0.], [0., -1., 0.]])
    tets = {(0, 1, 2, 3), (0, 1, 3, 4), (0, 1, 4, 5), (0, 1, 2, 5)}
    positions = {(0, 1, 2, 3): (1., 1., 0.),
                 (0, 1, 3, 4): (-1., 1., 0.),
                 (0, 1, 4, 5): (-1., -1., 0.),
                 (0, 1, 2, 5): (1., -1., 0.)}
    vertices = {tet: PowerVertex(tet, positions[tet], 0., 0.) for tet in tets}
    result = SimpleNamespace(points=points, full_simplices={3: tets}, power_vertices=vertices)
    facet = SimpleNamespace(generator_ids=(0, 1), incident_regular_triangles=(
        (0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 1, 5)))
    area = power_facet_polygon_area(result, facet)
    assert area["status"] == "finite_bounded"
    assert area["polygon_vertex_count"] == 4
    assert area["area_A2"] == pytest.approx(4.0)
    summary = facet_area_statistics([area["area_A2"], 0.5, 12.0])
    assert summary["VIA_A2"] == pytest.approx(16.5)
    assert summary["mean_A2"] == pytest.approx(5.5)
    assert summary["median_A2"] == pytest.approx(4.0)
    assert summary["fraction_below_1_A2"] == pytest.approx(1 / 3)
    assert summary["fraction_above_10_A2"] == pytest.approx(1 / 3)


def test_facet_components_use_shared_edges_and_scc_area_threshold():
    f1, f2, f3, f4 = ((0, 4), (1, 4), (2, 5), (3, 5))
    result = SimpleNamespace(
        final_facets=[SimpleNamespace(generator_ids=f) for f in (f1, f2, f3, f4)],
        included_edges=[SimpleNamespace(bicolor_facets=(f1, f2))],
    )
    components = interface_components(result, {f1: 4., f2: 4., f3: 0.5, f4: 0.5})
    assert len(components) == 3
    assert components[0]["facets"] == (f1, f2)
    assert components[0]["via_fraction"] == pytest.approx(8 / 9)
    assert components[0]["significant"]
    assert sum(component["significant"] for component in components) == 1


def test_scc_threshold_includes_exactly_seven_point_five_percent():
    large, boundary, small = (0, 5), (1, 5), (2, 6)
    result = SimpleNamespace(
        final_facets=[SimpleNamespace(generator_ids=f) for f in (large, boundary, small)],
        included_edges=[],
    )
    components = interface_components(result, {large: 92.5, boundary: 7.5, small: 0.0})
    significant = [component for component in components if component["significant"]]
    assert len(significant) == 2
    assert sorted(component["via_fraction"] for component in significant) == pytest.approx([0.075, 0.925])


def test_supplied_pdb_component_identity_is_parsed_from_compnd():
    data = Path(__file__).parents[2] / "data"
    two_kai = compnd_chain_mapping(data / "2KAI.pdb")
    one_udi = compnd_chain_mapping(data / "1UDI.pdb")
    assert any(x["molecule"] == "KALLIKREIN A" and "A" in x["chains"] for x in two_kai)
    assert any(x["molecule"] == "KALLIKREIN A" and "B" in x["chains"] for x in two_kai)
    assert any("TRYPSIN INHIBITOR" in x["molecule"] and "I" in x["chains"] for x in two_kai)
    assert any(x["molecule"] == "URACIL-DNA GLYCOSYLASE" and "E" in x["chains"] for x in one_udi)
    assert any("INHIBITOR PROTEIN" in x["molecule"] and "I" in x["chains"] for x in one_udi)


def test_known_thirty_degree_beta_lbeta_reversal_and_ordering():
    angle = np.deg2rad(30.0)
    points = np.array([[0., 0., 0.],
                       [np.cos(angle / 2), np.sin(angle / 2), 0.],
                       [np.cos(angle / 2), -np.sin(angle / 2), 0.]])
    _, _, beta = cazals_signed_beta((0, 1, 2), points, {0}, {1, 2})
    assert beta == pytest.approx(np.pi / 6)
    assert 5.0 * beta == pytest.approx(5.0 * np.pi / 6)
    assert 0.5 * 5.0 * beta == pytest.approx(5.0 * np.pi / 12)
    assert cazals_signed_beta((2, 0, 1), points, {0}, {1, 2})[2] == pytest.approx(beta)

    _, _, reversed_beta = cazals_signed_beta((2, 1, 0), points, {1, 2}, {0})
    assert reversed_beta == pytest.approx(-beta)
    assert reversed_beta * 5.0 == pytest.approx(-5.0 * beta)
    assert abs(reversed_beta) == pytest.approx(abs(beta))
    assert abs(reversed_beta * 5.0) == pytest.approx(abs(beta * 5.0))
    assert 5.0 == pytest.approx(5.0)


def test_covalent_pair_classification_is_topology_membership_not_distance():
    bonds = {(3, 9)}
    assert classify_same_group_pair((9, 3), bonds) == "covalent"
    assert classify_same_group_pair((3, 10), bonds) == "noncovalent"


def test_final_facet_polygon_edges_are_classified_as_boundary_or_interior():
    points = np.array([[0., 0., -1.], [0., 0., 1.],
                       [1., 0., 0.], [0., 1., 0.],
                       [-1., 0., 0.], [0., -1., 0.]])
    tets = {(0, 1, 2, 3), (0, 1, 3, 4), (0, 1, 4, 5), (0, 1, 2, 5)}
    positions = {(0, 1, 2, 3): (1., 1., 0.), (0, 1, 3, 4): (-1., 1., 0.),
                 (0, 1, 4, 5): (-1., -1., 0.), (0, 1, 2, 5): (1., -1., 0.)}
    vertices = {tet: PowerVertex(tet, positions[tet], 0., 0.) for tet in tets}
    facet = SimpleNamespace(generator_ids=(0, 1), incident_regular_triangles=(
        (0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 1, 5)))
    result = SimpleNamespace(points=points, full_simplices={3: tets},
                             power_vertices=vertices, final_facets=[facet])
    rows, issues = audit_final_interface_polygons(result)
    assert not issues
    assert len(rows) == 4
    assert {row["classification"] for row in rows} == {"BOUNDARY"}
    assert all(row["incident_selected_facet_count"] == 1 for row in rows)
    assert classify_interface_edge(2) == "INTERIOR"
    assert classify_interface_edge(3) == "NONMANIFOLD"
    assert classify_interface_edge(0, geometry_status="unbounded_boundary_power_edge") == "UNBOUNDED"


def test_polygon_edge_adjacency_confirms_two_selected_facet_polygons_share_a_segment():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                       [0., 0., -1.], [0., 0., 2.], [0., 2., 0.],
                       [2., 0., 0.]])
    t0, t1 = (0, 1, 2, 3), (0, 1, 2, 4)
    t2, t3 = (0, 1, 4, 5), (0, 1, 3, 5)
    t4, t5 = (0, 2, 4, 6), (0, 2, 3, 6)
    positions = {t0: (0., 0., 0.), t1: (0., 0., 1.),
                 t2: (0., 1., 1.), t3: (0., 1., 0.),
                 t4: (1., 0., 1.), t5: (1., 0., 0.)}
    vertices = {tet: PowerVertex(tet, pos, 0., 0.) for tet, pos in positions.items()}
    f01 = SimpleNamespace(generator_ids=(0, 1), incident_regular_triangles=(
        (0, 1, 2), (0, 1, 4), (0, 1, 5), (0, 1, 3)))
    f02 = SimpleNamespace(generator_ids=(0, 2), incident_regular_triangles=(
        (0, 1, 2), (0, 2, 4), (0, 2, 6), (0, 2, 3)))
    result = SimpleNamespace(points=points, full_simplices={3: set(positions)},
                             power_vertices=vertices, final_facets=[f01, f02])
    rows, issues = audit_final_interface_polygons(result)
    assert not issues
    shared = [row for row in rows if row["triangle"] == (0, 1, 2)]
    assert len(shared) == 1
    assert shared[0]["incident_facets"] == ((0, 1), (0, 2))
    assert shared[0]["classification"] == "INTERIOR"
    assert shared[0]["length"] == pytest.approx(1.0)
    assert sum(row["classification"] == "BOUNDARY" for row in rows) == 6


def test_edge_set_filtering_exposes_the_thirteen_triangle_alpha_discrepancy():
    rows = [{"classification": "INTERIOR", "geometry_status": "finite",
             "triangle_alpha_zero": i < 851, "length": 1.0, "beta": 0.1,
             "length_beta": 0.1} for i in range(864)]
    assert len(select_interior_finite_edges(rows)) == 864
    triangle_selected = select_interior_finite_edges(rows, require_triangle_alpha=True)
    assert len(triangle_selected) == 851
    assert len(rows) - len(triangle_selected) == 13
    summary = summarize_edge_records(triangle_selected)
    assert summary["C_signed"] == pytest.approx(85.1)
    assert summary["s_H_signed_deg"] == pytest.approx(np.degrees(0.1))


def test_condition_b_audit_retains_m_over_r_values_and_acceptance():
    facet = SimpleNamespace(alpha_zero_selected=True, generator_ids=(3, 8),
                            smaller_expanded_radius=2.5,
                            largest_orthogonal_ball_radius=8.0,
                            m_over_r=3.2, condition_beta_accepted=True,
                            condition_beta_status="accepted",
                            incident_regular_triangles=())
    excluded = SimpleNamespace(alpha_zero_selected=False, generator_ids=(1, 8),
                               smaller_expanded_radius=2.0,
                               largest_orthogonal_ball_radius=None, m_over_r=None,
                               condition_beta_accepted=None,
                               condition_beta_status="not_evaluated_not_alpha0")
    result = SimpleNamespace(facets=[facet, excluded], condition_beta_m=5.0,
                             full_simplices={3: set()}, power_vertices={})
    rows = condition_beta_audit_rows(result, {3: "A:foo", 8: "B:bar"})
    assert len(rows) == 1
    assert rows[0]["m_over_r"] == pytest.approx(3.2)
    assert rows[0]["threshold_M"] == 5.0
    assert rows[0]["accepted"] is True
