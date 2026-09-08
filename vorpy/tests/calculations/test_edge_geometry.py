import numpy as np
import pytest

from vorpy.src.calculations.edge_geometry import (
    LineEdgeGeometry,
    ParametricEdgeGeometry,
    canonical_edge_endpoints,
    line_geometry_from_endpoints,
    normalize_edge_construction,
    AdditivelyWeightedTrisectorBranch,
    match_aw_trisector_to_endpoints,
    validate_edge_samples,
)


def test_line_geometry_is_exact():
    edge = LineEdgeGeometry([1, 2, 3], [4, 6, 3])
    assert edge.point(0.0) == pytest.approx([1, 2, 3])
    assert edge.point(1.0) == pytest.approx([4, 6, 3])
    assert edge.arc_length() == pytest.approx(5.0)
    assert edge.curvature(0.5) == pytest.approx(0.0)
    assert edge.integrated_curvature() == pytest.approx(0.0)
    assert edge.integrated_curvature_squared() == pytest.approx(0.0)


def test_circle_parameterization_has_analytic_integrals():
    radius = 3.0
    edge = ParametricEdgeGeometry(
        point=lambda t: [radius*np.cos(t), radius*np.sin(t), 0],
        tangent=lambda t: [-radius*np.sin(t), radius*np.cos(t), 0],
        second_derivative=lambda t: [-radius*np.cos(t), -radius*np.sin(t), 0],
        t_min=0,
        t_max=np.pi/2,
        conic_type="circle",
    )
    assert edge.arc_length() == pytest.approx(radius*np.pi/2)
    assert edge.integrated_curvature() == pytest.approx(np.pi/2)
    assert edge.integrated_curvature_squared() == pytest.approx(np.pi/(2*radius))


def test_line_is_rigid_transform_invariant():
    start = np.array([0.2, -1.0, 3.0])
    end = np.array([2.4, 4.0, -0.5])
    theta = 0.7
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ])
    translation = np.array([8.0, -2.0, 1.5])
    original = LineEdgeGeometry(start, end)
    transformed = LineEdgeGeometry(
        rotation@start + translation, rotation@end + translation
    )
    assert transformed.arc_length() == pytest.approx(original.arc_length())


def test_canonical_endpoint_orientation_ignores_discovery_order():
    locations = [np.array([4.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
    forward = canonical_edge_endpoints([9, 3], locations)
    reverse = canonical_edge_endpoints([3, 9], locations[::-1])
    assert forward[0][0] == 3
    assert forward[1][0] == 9
    assert forward[0][1] == pytest.approx(reverse[0][1])
    assert forward[1][1] == pytest.approx(reverse[1][1])

    line_forward = line_geometry_from_endpoints([9, 3], locations)
    line_reverse = line_geometry_from_endpoints([3, 9], locations[::-1])
    assert line_forward.point(0.0) == pytest.approx(line_reverse.point(0.0))
    assert line_forward.unit_tangent(0.5) == pytest.approx(
        line_reverse.unit_tangent(0.5)
    )


def test_current_curved_metadata_normalizes_without_mutating_input():
    raw = {
        "case": "3a",
        "loc": np.array([1.0, 2.0, 3.0]),
        "rad": 2.5,
        "loc2": None,
        "rad2": None,
        "vmid": np.array([0.5, 0.0, 0.0]),
        "pnorm": np.array([0.0, 0.0, 1.0]),
        "dnorm": np.array([0.0, 1.0, 0.0]),
        "pa": np.array([0.5, -1.0, 0.0]),
        "outside": True,
    }
    normalized = normalize_edge_construction(raw)
    assert normalized.case == "3a"
    assert normalized.center == pytest.approx([1.0, 2.0, 3.0])
    assert normalized.radius == pytest.approx(2.5)
    assert normalized.secondary_center is None
    assert normalized.outside is True
    assert set(raw) == {
        "case", "loc", "rad", "loc2", "rad2", "vmid", "pnorm",
        "dnorm", "pa", "outside",
    }


def test_straight_edge_metadata_with_only_center_and_radius_is_supported():
    normalized = normalize_edge_construction({"loc": [1, 2, 3], "rad": 4})
    assert normalized.case is None
    assert normalized.center == pytest.approx([1, 2, 3])
    assert normalized.radius == pytest.approx(4)
    assert normalized.midpoint is None


def test_aw_trisector_points_have_equal_clearance_to_all_generators():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    edge = AdditivelyWeightedTrisectorBranch(
        locations, radii, rho_min=2.0, rho_max=4.0, branch=1
    )
    for rho in np.linspace(2.0, 4.0, 9):
        assert edge.clearance_residuals(rho) == pytest.approx(
            np.zeros(3), abs=1e-12
        )
        assert edge.parameter_for_point(edge.point(rho)) == pytest.approx(rho)


def test_aw_trisector_derivatives_match_centered_finite_differences():
    edge = AdditivelyWeightedTrisectorBranch(
        [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
        [1., 1.5, 2.], 2.0, 4.0, branch=-1
    )
    rho = 3.0
    step = 1e-5
    numerical_first = (edge.point(rho + step) - edge.point(rho - step))/(2*step)
    numerical_second = (
        edge.point(rho + step) - 2*edge.point(rho) + edge.point(rho - step)
    )/step**2
    assert edge.tangent(rho) == pytest.approx(numerical_first, rel=1e-8, abs=1e-8)
    assert edge.second_derivative(rho) == pytest.approx(
        numerical_second, rel=2e-5, abs=2e-5
    )


def test_aw_trisector_is_rigid_transform_invariant():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    theta = 0.41
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.],
        [np.sin(theta), np.cos(theta), 0.],
        [0., 0., 1.],
    ])
    translation = np.array([7., -3., 2.])
    original = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=1
    )
    transformed = AdditivelyWeightedTrisectorBranch(
        locations@rotation.T + translation, radii, 2.0, 4.0, branch=1
    )
    assert transformed.arc_length() == pytest.approx(original.arc_length(), rel=1e-12)
    assert transformed.integrated_curvature() == pytest.approx(
        original.integrated_curvature(), rel=1e-12
    )
    assert transformed.integrated_curvature_squared() == pytest.approx(
        original.integrated_curvature_squared(), rel=1e-12
    )


def test_aw_rho_parameterization_rejects_turning_point_interval():
    with pytest.raises(ValueError, match="turning point"):
        AdditivelyWeightedTrisectorBranch(
            [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
            [1., 1.5, 2.], 1.0, 2.0, branch=1
        )


def test_aw_endpoint_match_reproduces_canonical_vertices_and_samples():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    reference = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=-1
    )
    start, end = reference.point(4.0), reference.point(2.0)
    match = match_aw_trisector_to_endpoints(
        locations, radii, [3, 9], [start, end]
    )
    assert match.status == "matched"
    assert match.geometry.point(0.0) == pytest.approx(start, abs=1e-12)
    assert match.geometry.point(1.0) == pytest.approx(end, abs=1e-12)
    samples = [reference.point(rho) for rho in np.linspace(2.0, 4.0, 11)]
    validation = validate_edge_samples(match.geometry, samples)
    assert validation.count == 11
    assert validation.maximum_error < 1e-12


def test_aw_endpoint_match_is_independent_of_discovery_order():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    reference = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=1
    )
    points = [reference.point(2.2), reference.point(3.8)]
    forward = match_aw_trisector_to_endpoints(
        locations, radii, [5, 8], points
    )
    reverse = match_aw_trisector_to_endpoints(
        locations, radii, [8, 5], points[::-1]
    )
    assert forward.status == reverse.status == "matched"
    for t in np.linspace(0.0, 1.0, 7):
        assert forward.geometry.point(t) == pytest.approx(
            reverse.geometry.point(t), abs=1e-12
        )


def test_aw_endpoint_match_reports_cross_branch_without_approximation():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    positive = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=1
    )
    negative = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=-1
    )
    match = match_aw_trisector_to_endpoints(
        locations, radii, [1, 2], [positive.point(2.0), negative.point(4.0)]
    )
    assert match.status == "cross_branch"
    assert match.geometry is None


def test_aw_endpoint_match_reports_equal_rho_as_turning_point_case():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    positive = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=1
    )
    negative = AdditivelyWeightedTrisectorBranch(
        locations, radii, 2.0, 4.0, branch=-1
    )
    match = match_aw_trisector_to_endpoints(
        locations, radii, [1, 2], [positive.point(3.0), negative.point(3.0)]
    )
    assert match.status == "turning_point"
    assert match.geometry is None
