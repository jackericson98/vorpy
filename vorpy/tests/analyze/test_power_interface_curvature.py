import math

import numpy as np
import pytest

from vorpy.src.analyze.power_interface_curvature import (
    PowerVertex,
    _closure_issues,
    _condition_beta_for_facet,
    analyze_power_interface_curvature,
    analyze_power_regular_complex,
    calculate_power_vertex,
    cazals_signed_beta,
    clip_finite_power_edge_to_balls,
    oriented_power_facet_normal,
    power_edge_geometry,
)


def test_equal_radius_power_vertex_has_expected_equal_power_residuals():
    points = np.array([
        [0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.],
    ])
    weights = np.ones(4)
    vertex = calculate_power_vertex((3, 1, 0, 2), points, weights)
    assert vertex.tetrahedron == (0, 1, 2, 3)
    assert vertex.position == pytest.approx((1., 1., 1.))
    assert vertex.power_value == pytest.approx(2.)
    assert vertex.residual < 1e-12


def test_unequal_weight_power_vertex_satisfies_all_four_equations():
    points = np.array([
        [0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.],
    ])
    weights = np.array([1., 4., 2.25, 6.25])
    vertex = calculate_power_vertex((0, 1, 2, 3), points, weights)
    powers = [np.dot(vertex.position - points[i], vertex.position - points[i]) - weights[i]
              for i in range(4)]
    assert max(powers) - min(powers) < 1e-12


def test_alpha_clipping_solves_exact_quadratic_sublevel_segment():
    points = np.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    length, status = clip_finite_power_edge_to_balls(
        (0, 1, 2), (-2., .5, .5), (2., .5, .5), points, np.ones(3)
    )
    assert status == "clipped_to_alpha_balls"
    assert length == pytest.approx(2.0 * np.sqrt(0.5))


def test_cazals_beta_uses_singleton_angle_and_source_sign_convention():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    singleton, group, beta = cazals_signed_beta((2, 0, 1), points, {0}, {1, 2})
    assert singleton == 0
    assert group == "A"
    assert beta == pytest.approx(math.pi / 2)
    reverse_singleton, reverse_group, reverse_beta = cazals_signed_beta(
        (1, 2, 0), points, {1, 2}, {0}
    )
    assert reverse_singleton == 0
    assert reverse_group == "B"
    assert reverse_beta == pytest.approx(-beta)
    # Reordering simplex IDs cannot affect the geometric result.
    assert cazals_signed_beta((1, 0, 2), points, {0}, {1, 2})[2] == pytest.approx(beta)
    normal = oriented_power_facet_normal((0, 1), points, {0}, {1, 2})
    reversed_normal = oriented_power_facet_normal((1, 0), points, {1, 2}, {0})
    assert reversed_normal == pytest.approx(-normal)


def test_condition_beta_uses_maximum_power_vertex_radius_and_m_over_r():
    edge = (0, 1)
    triangles = ((0, 1, 2), (0, 1, 3))
    triangle_tets = {
        triangles[0]: ((0, 1, 2, 4), (0, 1, 2, 5)),
        triangles[1]: ((0, 1, 3, 4), (0, 1, 3, 5)),
    }
    vertices = {}
    # Fill all four endpoint records; max power=9 -> m=3, r=2, ratio=1.5.
    vertices[(0, 1, 2, 4)] = PowerVertex((0, 1, 2, 4), (0., 0., 0.), 4., 0.)
    vertices[(0, 1, 2, 5)] = PowerVertex((0, 1, 2, 5), (0., 0., 1.), 9., 0.)
    vertices[(0, 1, 3, 4)] = PowerVertex((0, 1, 3, 4), (0., 1., 0.), 4., 0.)
    vertices[(0, 1, 3, 5)] = PowerVertex((0, 1, 3, 5), (0., 1., 1.), 1., 0.)
    status, radius, ratio, accepted, _ = _condition_beta_for_facet(
        edge, triangles, triangle_tets, vertices, np.array([2., 3., 1., 1.]), 5.
    )
    assert status == "accepted"
    assert radius == pytest.approx(3.)
    assert ratio == pytest.approx(1.5)
    assert accepted

    vertices[(0, 1, 2, 5)] = PowerVertex((0, 1, 2, 5), (0., 0., 1.), 144., 0.)
    status, radius, ratio, accepted, _ = _condition_beta_for_facet(
        edge, triangles, triangle_tets, vertices, np.array([2., 3., 1., 1.]), 5.
    )
    assert status == "rejected_m_over_r"
    assert ratio == pytest.approx(6.)
    assert not accepted


def test_unbounded_facet_and_dual_edge_remain_explicitly_classified():
    status, radius, ratio, accepted, reason = _condition_beta_for_facet(
        (0, 1), ((0, 1, 2),), {(0, 1, 2): ((0, 1, 2, 3),)},
        {(0, 1, 2, 3): PowerVertex((0, 1, 2, 3), (0., 0., 0.), 1., 0.)},
        np.ones(4), 5.,
    )
    assert status == "rejected_unbounded"
    assert radius is None and math.isinf(ratio) and not accepted
    assert "ray" in reason

    edge_status, coordinates, length, edge_reason = power_edge_geometry(
        (0, 1, 2), ((0, 1, 2, 3),),
        {(0, 1, 2, 3): PowerVertex((0, 1, 2, 3), (1., 2., 3.), 0., 0.)},
    )
    assert edge_status == "unbounded_boundary_power_edge"
    assert coordinates == ((1., 2., 3.),)
    assert length is None and "ray" in edge_reason


def test_full_and_alpha_complex_closure_validator_reports_missing_faces():
    bad = {
        0: {(0,), (1,), (2,)},
        1: {(0, 1), (0, 2)},
        2: {(0, 1, 2)},
        3: set(),
    }
    assert _closure_issues(bad) == [((0, 1, 2), (1, 2))]


def test_regular_complex_entry_point_preserves_alpha_selection_separately():
    points = np.array([
        [0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.], [0., 0., -2.],
    ])
    tetrahedra = {(0, 1, 2, 3), (0, 1, 2, 4)}
    full = {d: set() for d in range(4)}
    for tet in tetrahedra:
        for d in range(4):
            from itertools import combinations
            full[d].update(tuple(sorted(face)) for face in combinations(tet, d + 1))
    filtration = {simplex: -1.0 for simplex in full[0]}
    filtration.update({simplex: 1.0 for simplex in full[1]})
    filtration[(0, 1)] = -0.1
    filtration.update({simplex: 2.0 for simplex in full[2] | full[3]})
    result = analyze_power_regular_complex(
        points, np.ones(5), {0}, {1, 2, 3, 4}, full, filtration
    )
    assert result.full_counts[1] == len(full[1])
    assert result.alpha_counts[1] == 1
    assert result.alpha_counts[0] == 5
    assert len(result.all_bicolor_facets) > len(result.alpha_zero_facets)
    assert result.smooth_surface_curvature == 0.0


def test_gudhi_pipeline_keeps_full_and_alpha_zero_complex_separate(tmp_path):
    pytest.importorskip("gudhi")
    points = np.array([
        [0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.],
    ])
    result = analyze_power_interface_curvature(
        points, np.ones(4), {0}, {1, 2, 3}, alpha=0.0
    )
    assert all(result.alpha_counts[d] <= result.full_counts[d] for d in range(4))
    assert not _closure_issues(result.full_simplices)
    assert not _closure_issues(result.alpha_simplices)
    assert result.smooth_surface_curvature == 0.0
    result.export(tmp_path)
    assert (tmp_path / "power_interface_surfaces.csv").exists()
    assert (tmp_path / "power_interface_edges.csv").exists()
    assert (tmp_path / "power_simplex_filtration.csv").exists()
