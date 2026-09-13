import numpy as np
import pytest

from vorpy.src.calculations import (
    LineEdgeGeometry,
    canonical_edge_endpoints,
    line_geometry_from_endpoints,
    normalize_edge_construction,
    AdditivelyWeightedTrisectorBranch,
    AdditivelyWeightedTrisectorConic,
    match_aw_trisector_to_endpoints,
    match_aw_trisector_conic_to_samples,
    validate_edge_samples,
    aw_edge_point_geometry,
    boundary_turning_angle,
    aw_cell_turning_angle,
    aw_cell_edge_mean_curvature,
    aw_edge_mean_curvatures
)

from vorpy.src.calculations.edge_geometry import (
    ParametricEdgeGeometry,
    signed_geodesic_curvature,
    integrated_geodesic_curvature,
    aw_face_edge_geodesic_curvature,
    integrated_aw_face_edge_geodesic_curvature,
    network_surface_boundary_cycle,
    aw_face_edge_boundary_orientation,
oriented_aw_face_edge_geodesic_curvature
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


def test_nonsingular_aw_hyperbola_crosses_turning_point_smoothly():
    edge = AdditivelyWeightedTrisectorConic(
        [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
        [1., 1.5, 2.], -1.0, 1.0, component=1
    )
    assert edge.conic_type == "hyperbola"
    assert edge.hyperbola_orientation == "rho"
    assert edge.clearance_residuals(0.0) == pytest.approx(
        np.zeros(3), abs=1e-12
    )
    assert np.all(np.isfinite(edge.tangent(0.0)))
    assert np.linalg.norm(edge.tangent(0.0)) > 0.0


def test_nonsingular_aw_match_uses_samples_to_cross_rho_turning_point():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    reference = AdditivelyWeightedTrisectorConic(
        locations, radii, -0.7, 1.1, component=1
    )
    parameters = np.linspace(-0.7, 1.1, 19)
    samples = [reference.point(value) for value in parameters]
    match = match_aw_trisector_conic_to_samples(
        locations, radii, [8, 3],
        [samples[-1], samples[0]], samples[::-1],
    )
    assert match.status == "matched_nonsingular"
    assert match.geometry.point(0.0) == pytest.approx(samples[0], abs=1e-12)
    assert match.geometry.point(1.0) == pytest.approx(samples[-1], abs=1e-12)
    validation = validate_edge_samples(match.geometry, samples)
    assert validation.maximum_error < 1e-12


def test_nonsingular_aw_match_handles_equal_rho_endpoints():
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.])
    reference = AdditivelyWeightedTrisectorConic(
        locations, radii, -1.0, 1.0, component=1
    )
    samples = [
        reference.point(value) for value in np.linspace(-1.0, 1.0, 21)
    ]
    match = match_aw_trisector_conic_to_samples(
        locations, radii, [2, 7],
        [samples[0], samples[-1]], samples,
    )
    assert match.status == "matched_nonsingular"
    assert match.geometry.arc_length() == pytest.approx(
        reference.arc_length(), rel=1e-12
    )


def test_nonsingular_aw_parabola_is_regular_at_normal_coordinate_zero():
    theta = 0.5
    locations = [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]]
    radii = [
        5.0,
        5.0 + 4.0*np.cos(theta),
        5.0 + 5.0*np.sin(theta),
    ]
    edge = AdditivelyWeightedTrisectorConic(
        locations, radii, -1.0, 1.0
    )
    assert edge.conic_type == "parabola"
    assert edge.clearance_residuals(0.0) == pytest.approx(
        np.zeros(3), abs=1e-12
    )
    assert np.linalg.norm(edge.tangent(0.0)) > 0.0


def test_aw_edge_point_geometry_returns_exact_local_vectors():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    indices = [10, 20, 30]
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    result = aw_edge_point_geometry(
        edge_geometry=edge,
        t=0.5,
        generator_indices=indices,
        generator_locations=locations,
        cell_index=10,
    )

    assert result.point == pytest.approx([1.0, 1.0, 0.0])
    assert result.tangent == pytest.approx([0.0, 0.0, 1.0])

    assert result.first_normal == pytest.approx(
        [1.0, 0.0, 0.0]
    )

    assert result.second_normal == pytest.approx(
        [0.0, 1.0, 0.0]
    )

    assert result.cell_index == 10
    assert result.neighbor_indices == (20, 30)


def test_aw_edge_point_geometry_changes_with_cell_perspective():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    indices = [10, 20, 30]
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    cell_10 = aw_edge_point_geometry(
        edge, 0.5, indices, locations, 10
    )

    cell_20 = aw_edge_point_geometry(
        edge, 0.5, indices, locations, 20
    )

    assert cell_10.first_normal == pytest.approx(
        -cell_20.first_normal
    )


def test_aw_edge_tangent_is_orthogonal_to_both_incident_normals():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    result = aw_edge_point_geometry(
        edge,
        0.5,
        [10, 20, 30],
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
        ],
        10,
    )

    assert np.dot(
        result.tangent,
        result.first_normal
    ) == pytest.approx(0.0, abs=1e-12)

    assert np.dot(
        result.tangent,
        result.second_normal
    ) == pytest.approx(0.0, abs=1e-12)


def test_aw_edge_point_geometry_is_rigid_transform_covariant():
    start = np.array([1.0, 1.0, -2.0])
    end = np.array([1.0, 1.0, 2.0])

    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    theta = 0.43
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta),  np.cos(theta), 0.0],
        [0.0,            0.0,           1.0],
    ])

    translation = np.array([4.0, -3.0, 7.0])

    original = aw_edge_point_geometry(
        LineEdgeGeometry(start, end),
        0.37,
        [10, 20, 30],
        locations,
        10,
    )

    transformed = aw_edge_point_geometry(
        LineEdgeGeometry(
            rotation @ start + translation,
            rotation @ end + translation,
        ),
        0.37,
        [10, 20, 30],
        locations @ rotation.T + translation,
        10,
    )

    assert transformed.point == pytest.approx(
        rotation @ original.point + translation,
        abs=1e-12,
    )

    assert transformed.tangent == pytest.approx(
        rotation @ original.tangent,
        abs=1e-12,
    )

    assert transformed.first_normal == pytest.approx(
        rotation @ original.first_normal,
        abs=1e-12,
    )

    assert transformed.second_normal == pytest.approx(
        rotation @ original.second_normal,
        abs=1e-12,
    )


def test_boundary_turning_angle_is_zero_for_smooth_continuation():
    angle = boundary_turning_angle(
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    )

    assert angle == pytest.approx(0.0, abs=1e-14)


def test_boundary_turning_angle_is_pi_over_two_for_right_angle():
    angle = boundary_turning_angle(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    )

    assert angle == pytest.approx(np.pi/2, abs=1e-14)


def test_boundary_turning_angle_recovers_known_angle():
    expected = np.pi/3

    first = np.array([1.0, 0.0, 0.0])
    second = np.array([
        np.cos(expected),
        np.sin(expected),
        0.0,
    ])

    assert boundary_turning_angle(
        first,
        second,
    ) == pytest.approx(expected, abs=1e-14)


def test_boundary_turning_angle_is_independent_of_surface_order():
    first = np.array([0.3, 0.4, 0.866025403784])
    first /= np.linalg.norm(first)

    second = np.array([-0.2, 0.9, 0.38])
    second /= np.linalg.norm(second)

    forward = boundary_turning_angle(first, second)
    reverse = boundary_turning_angle(second, first)

    assert forward == pytest.approx(reverse, abs=1e-14)


def test_aw_cell_turning_angle_matches_local_incident_normals():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    indices = [10, 20, 30]

    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    angle = aw_cell_turning_angle(
        edge_geometry=edge,
        t=0.5,
        generator_indices=indices,
        generator_locations=locations,
        cell_index=10,
    )

    assert angle == pytest.approx(np.pi/2, abs=1e-12)


def test_aw_turning_angle_along_curved_edge():
    locations = np.array([
        [0., 0., 0.],
        [4., 0., 0.],
        [0., 5., 0.],
    ])
    radii = np.array([1., 1.5, 2.])

    edge = AdditivelyWeightedTrisectorConic(
        locations,
        radii,
        -0.7,
        1.1,
        component=1,
    )

    indices = [10, 20, 30]

    parameters = np.linspace(edge.t_min, edge.t_max, 9)

    angles = []

    for t in parameters:
        local = aw_edge_point_geometry(
            edge_geometry=edge,
            t=t,
            generator_indices=indices,
            generator_locations=locations,
            cell_index=10,
        )

        angle = aw_cell_turning_angle(
            edge_geometry=edge,
            t=t,
            generator_indices=indices,
            generator_locations=locations,
            cell_index=10,
        )

        direct = np.arccos(np.clip(
            np.dot(local.first_normal, local.second_normal),
            -1.0,
            1.0,
        ))

        assert angle == pytest.approx(direct, abs=1e-12)
        assert np.isfinite(angle)

        # The analytic edge tangent must remain tangent to both incident
        # surfaces everywhere along the curved edge.
        assert np.dot(
            local.tangent,
            local.first_normal,
        ) == pytest.approx(0.0, abs=1e-12)

        assert np.dot(
            local.tangent,
            local.second_normal,
        ) == pytest.approx(0.0, abs=1e-12)

        angles.append(angle)

    # This particular unequal-radius hyperbolic edge has a genuinely
    # varying incident angle along its length.
    assert np.ptp(angles) > 1e-8


def test_aw_cell_edge_mean_curvature_right_angle_has_exact_value():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    value = aw_cell_edge_mean_curvature(
        edge_geometry=edge,
        generator_indices=[10, 20, 30],
        generator_locations=locations,
        cell_index=10,
    )

    assert value == pytest.approx(np.pi, abs=1e-12)


def test_aw_cell_edge_mean_curvature_scales_linearly():
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    original = aw_cell_edge_mean_curvature(
        edge,
        [10, 20, 30],
        locations,
        10,
    )

    scale = 3.7

    scaled = aw_cell_edge_mean_curvature(
        LineEdgeGeometry(
            scale * np.array([1.0, 1.0, -2.0]),
            scale * np.array([1.0, 1.0, 2.0]),
        ),
        [10, 20, 30],
        scale * locations,
        10,
    )

    assert scaled == pytest.approx(
        scale * original,
        rel=1e-12,
    )


def test_aw_cell_edge_mean_curvature_curved_edge_converges():
    locations = np.array([
        [0., 0., 0.],
        [4., 0., 0.],
        [0., 5., 0.],
    ])

    radii = np.array([1., 1.5, 2.])

    edge = AdditivelyWeightedTrisectorConic(
        locations,
        radii,
        -0.7,
        1.1,
        component=1,
    )

    kwargs = dict(
        edge_geometry=edge,
        generator_indices=[10, 20, 30],
        generator_locations=locations,
        cell_index=10,
    )

    value_16 = aw_cell_edge_mean_curvature(
        **kwargs,
        order=16,
    )

    value_32 = aw_cell_edge_mean_curvature(
        **kwargs,
        order=32,
    )

    value_64 = aw_cell_edge_mean_curvature(
        **kwargs,
        order=64,
    )

    assert np.isfinite(value_16)
    assert np.isfinite(value_32)
    assert np.isfinite(value_64)

    assert value_32 == pytest.approx(
        value_64,
        rel=1e-10,
        abs=1e-12,
    )


def test_aw_cell_edge_mean_curvature_is_rigid_transform_invariant():
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    start = np.array([1.0, 1.0, -2.0])
    end = np.array([1.0, 1.0, 2.0])

    original = aw_cell_edge_mean_curvature(
        edge_geometry=LineEdgeGeometry(start, end),
        generator_indices=[10, 20, 30],
        generator_locations=locations,
        cell_index=10,
        order=32,
    )

    theta = 0.63
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta),  np.cos(theta), 0.0],
        [0.0,            0.0,           1.0],
    ])
    translation = np.array([4.0, -7.0, 2.5])

    transformed_locations = (
        locations @ rotation.T + translation
    )

    transformed_start = rotation @ start + translation
    transformed_end = rotation @ end + translation

    transformed = aw_cell_edge_mean_curvature(
        edge_geometry=LineEdgeGeometry(
            transformed_start,
            transformed_end,
        ),
        generator_indices=[10, 20, 30],
        generator_locations=transformed_locations,
        cell_index=10,
        order=32,
    )

    assert transformed == pytest.approx(
        original,
        rel=1e-12,
        abs=1e-12,
    )


def test_aw_cell_edge_mean_curvature_vs_midpoint_diagnostic():
    locations = np.array([
        [0., 0., 0.],
        [4., 0., 0.],
        [0., 5., 0.],
    ])

    radii = np.array([1., 1.5, 2.])

    edge = AdditivelyWeightedTrisectorConic(
        locations,
        radii,
        -0.7,
        1.1,
        component=1,
    )

    indices = [10, 20, 30]

    integrated = aw_cell_edge_mean_curvature(
        edge_geometry=edge,
        generator_indices=indices,
        generator_locations=locations,
        cell_index=10,
        order=64,
    )

    midpoint_t = 0.5 * (edge.t_min + edge.t_max)

    midpoint_angle = aw_cell_turning_angle(
        edge_geometry=edge,
        t=midpoint_t,
        generator_indices=indices,
        generator_locations=locations,
        cell_index=10,
    )

    midpoint_approximation = (
        0.5
        * midpoint_angle
        * edge.arc_length(order=64)
    )

    difference = integrated - midpoint_approximation

    print()
    print("Curved AW edge mean-curvature diagnostic")
    print("----------------------------------------")
    print(f"Integrated value:       {integrated:.15g}")
    print(f"Midpoint approximation: {midpoint_approximation:.15g}")
    print(f"Difference:             {difference:.15g}")

    assert np.isfinite(integrated)
    assert np.isfinite(midpoint_approximation)
    assert np.isfinite(difference)


def test_cube_integrated_mean_curvature_is_three_pi_a():
    """A cube has zero face curvature and all of M concentrated on its edges."""
    side = 3.4

    # Adjacent outward face normals of a cube are perpendicular.
    theta = boundary_turning_angle(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    )

    assert theta == pytest.approx(np.pi / 2.0, abs=1e-14)

    number_of_edges = 12

    surface_contribution = 0.0
    edge_contribution = (
        0.5
        * number_of_edges
        * theta
        * side
    )

    total = surface_contribution + edge_contribution

    assert edge_contribution == pytest.approx(
        3.0 * np.pi * side,
        rel=1e-13,
        abs=1e-13,
    )

    assert total == pytest.approx(
        3.0 * np.pi * side,
        rel=1e-13,
        abs=1e-13,
    )


@pytest.mark.parametrize("side", [0.5, 1.0, 2.0, 7.3])
def test_cube_mean_curvature_scales_linearly(side):
    theta = np.pi / 2.0
    total = 0.5 * 12.0 * theta * side

    assert total == pytest.approx(
        3.0 * np.pi * side,
        rel=1e-13,
        abs=1e-13,
    )


def test_closed_equal_radius_aw_cubic_cell_has_exact_mean_curvature():
    """Six equal-radius AW neighbors produce an exact cubic central cell.

    The central generator is at the origin. Neighbors at +/-a on the
    coordinate axes generate the six planes x,y,z = +/-a/2.

    The resulting central AW cell is a cube of side a, so

        M_surface = 0
        M_edge = 3*pi*a
        M_total = 3*pi*a.
    """
    side = 4.0
    half = side / 2.0

    central_index = 0
    central = np.array([0.0, 0.0, 0.0])

    neighbors = {
        "xp": (1, np.array([ side, 0.0, 0.0])),
        "xm": (2, np.array([-side, 0.0, 0.0])),
        "yp": (3, np.array([0.0,  side, 0.0])),
        "ym": (4, np.array([0.0, -side, 0.0])),
        "zp": (5, np.array([0.0, 0.0,  side])),
        "zm": (6, np.array([0.0, 0.0, -side])),
    }

    # Every cube edge is the intersection of two of the six pairwise
    # bisector planes. Each tuple is:
    #
    #   (neighbor 1, neighbor 2, edge start, edge end)
    #
    # There are four x-directed, four y-directed, and four z-directed edges.
    edges = [
        # x-directed edges: y = +/-half, z = +/-half
        ("yp", "zp", [-half,  half,  half], [ half,  half,  half]),
        ("yp", "zm", [-half,  half, -half], [ half,  half, -half]),
        ("ym", "zp", [-half, -half,  half], [ half, -half,  half]),
        ("ym", "zm", [-half, -half, -half], [ half, -half, -half]),

        # y-directed edges: x = +/-half, z = +/-half
        ("xp", "zp", [ half, -half,  half], [ half,  half,  half]),
        ("xp", "zm", [ half, -half, -half], [ half,  half, -half]),
        ("xm", "zp", [-half, -half,  half], [-half,  half,  half]),
        ("xm", "zm", [-half, -half, -half], [-half,  half, -half]),

        # z-directed edges: x = +/-half, y = +/-half
        ("xp", "yp", [ half,  half, -half], [ half,  half,  half]),
        ("xp", "ym", [ half, -half, -half], [ half, -half,  half]),
        ("xm", "yp", [-half,  half, -half], [-half,  half,  half]),
        ("xm", "ym", [-half, -half, -half], [-half, -half,  half]),
    ]

    edge_contributions = []

    for first_name, second_name, start, end in edges:
        first_index, first_location = neighbors[first_name]
        second_index, second_location = neighbors[second_name]

        edge = LineEdgeGeometry(start, end)

        value = aw_cell_edge_mean_curvature(
            edge_geometry=edge,
            generator_indices=[
                central_index,
                first_index,
                second_index,
            ],
            generator_locations=np.array([
                central,
                first_location,
                second_location,
            ]),
            cell_index=central_index,
            order=32,
        )

        # Each cube edge has length a and turning angle pi/2:
        #
        #   M_edge,e = 1/2 * (pi/2) * a = pi*a/4.
        assert value == pytest.approx(
            np.pi * side / 4.0,
            rel=1e-12,
            abs=1e-12,
        )

        edge_contributions.append(value)

    assert len(edge_contributions) == 12

    surface_contribution = 0.0
    edge_contribution = sum(edge_contributions)
    total = surface_contribution + edge_contribution

    assert edge_contribution == pytest.approx(
        3.0 * np.pi * side,
        rel=1e-12,
        abs=1e-12,
    )

    assert total == pytest.approx(
        3.0 * np.pi * side,
        rel=1e-12,
        abs=1e-12,
    )


def test_aw_edge_mean_curvatures_matches_independent_cell_integrals():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    indices = [10, 20, 30]
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    combined = aw_edge_mean_curvatures(
        edge,
        indices,
        locations,
        order=32,
    )

    for cell_index in indices:
        reference = aw_cell_edge_mean_curvature(
            edge_geometry=edge,
            generator_indices=indices,
            generator_locations=locations,
            cell_index=cell_index,
            order=32,
        )

        assert combined[cell_index] == pytest.approx(
            reference,
            rel=1e-12,
            abs=1e-12,
        )

def test_straight_line_has_zero_geodesic_curvature_on_plane():
    edge = LineEdgeGeometry(
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    )

    value = signed_geodesic_curvature(
        edge_geometry=edge,
        t=0.5,
        surface_normal=[0.0, 0.0, 1.0],
    )

    assert value == pytest.approx(0.0, abs=1e-14)


def test_planar_circle_geodesic_curvature_is_inverse_radius():
    radius = 2.5

    edge = ParametricEdgeGeometry(
        point=lambda t: np.array([
            radius * np.cos(t),
            radius * np.sin(t),
            0.0,
        ]),
        tangent=lambda t: np.array([
            -radius * np.sin(t),
            radius * np.cos(t),
            0.0,
        ]),
        second_derivative=lambda t: np.array([
            -radius * np.cos(t),
            -radius * np.sin(t),
            0.0,
        ]),
        t_min=0.0,
        t_max=2.0 * np.pi,
    )

    for t in np.linspace(0.0, 2.0 * np.pi, 9):
        value = signed_geodesic_curvature(
            edge_geometry=edge,
            t=t,
            surface_normal=[0.0, 0.0, 1.0],
        )

        assert value == pytest.approx(
            1.0 / radius,
            rel=1e-12,
            abs=1e-12,
        )


def test_planar_circle_integrated_geodesic_curvature_is_two_pi():
    radius = 2.5

    edge = ParametricEdgeGeometry(
        point=lambda t: np.array([
            radius * np.cos(t),
            radius * np.sin(t),
            0.0,
        ]),
        tangent=lambda t: np.array([
            -radius * np.sin(t),
            radius * np.cos(t),
            0.0,
        ]),
        second_derivative=lambda t: np.array([
            -radius * np.cos(t),
            -radius * np.sin(t),
            0.0,
        ]),
        t_min=0.0,
        t_max=2.0 * np.pi,
    )

    value = integrated_geodesic_curvature(
        edge_geometry=edge,
        surface_normal=lambda t: np.array([0.0, 0.0, 1.0]),
        order=64,
    )

    assert value == pytest.approx(
        2.0 * np.pi,
        rel=1e-12,
        abs=1e-12,
    )


def test_great_circle_has_zero_geodesic_curvature():
    radius = 3.0

    edge = ParametricEdgeGeometry(
        point=lambda t: np.array([
            radius * np.cos(t),
            radius * np.sin(t),
            0.0,
        ]),
        tangent=lambda t: np.array([
            -radius * np.sin(t),
            radius * np.cos(t),
            0.0,
        ]),
        second_derivative=lambda t: np.array([
            -radius * np.cos(t),
            -radius * np.sin(t),
            0.0,
        ]),
        t_min=0.0,
        t_max=2.0 * np.pi,
    )

    for t in np.linspace(0.0, 2.0 * np.pi, 9):
        point = edge.point(t)
        sphere_normal = point / np.linalg.norm(point)

        value = signed_geodesic_curvature(
            edge_geometry=edge,
            t=t,
            surface_normal=sphere_normal,
        )

        assert value == pytest.approx(
            0.0,
            abs=1e-12,
        )


def test_planar_circle_geodesic_curvature_reverses_with_surface_normal():
    radius = 2.5

    edge = ParametricEdgeGeometry(
        point=lambda t: np.array([
            radius * np.cos(t),
            radius * np.sin(t),
            0.0,
        ]),
        tangent=lambda t: np.array([
            -radius * np.sin(t),
            radius * np.cos(t),
            0.0,
        ]),
        second_derivative=lambda t: np.array([
            -radius * np.cos(t),
            -radius * np.sin(t),
            0.0,
        ]),
        t_min=0.0,
        t_max=2.0 * np.pi,
    )

    positive = signed_geodesic_curvature(
        edge,
        0.7,
        [0.0, 0.0, 1.0],
    )

    negative = signed_geodesic_curvature(
        edge,
        0.7,
        [0.0, 0.0, -1.0],
    )

    assert negative == pytest.approx(-positive, abs=1e-12)


def test_aw_planar_face_straight_edge_has_zero_geodesic_curvature():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    indices = [10, 20, 30]
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    value = aw_face_edge_geodesic_curvature(
        edge_geometry=edge,
        t=0.5,
        generator_indices=indices,
        generator_locations=locations,
        cell_index=10,
        other_index=20,
    )

    assert value == pytest.approx(0.0, abs=1e-14)


def test_aw_planar_face_straight_edge_integral_is_zero():
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    indices = [10, 20, 30]
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    value = integrated_aw_face_edge_geodesic_curvature(
        edge_geometry=edge,
        generator_indices=indices,
        generator_locations=locations,
        cell_index=10,
        other_index=20,
        order=32,
    )

    assert value == pytest.approx(0.0, abs=1e-14)


def test_aw_curved_face_edge_geodesic_curvature_is_finite():
    locations = np.array([
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
    ])

    radii = np.array([1.0, 1.5, 2.0])

    edge = AdditivelyWeightedTrisectorConic(
        locations,
        radii,
        -0.7,
        1.1,
        component=1,
    )

    indices = [10, 20, 30]
    values = []

    for t in np.linspace(edge.t_min, edge.t_max, 9):
        value = aw_face_edge_geodesic_curvature(
            edge_geometry=edge,
            t=t,
            generator_indices=indices,
            generator_locations=locations,
            cell_index=10,
            other_index=20,
        )

        assert np.isfinite(value)
        values.append(value)

    # This unequal-radius curved AW edge should have a nonzero geodesic
    # curvature contribution on the selected pairwise face.
    assert np.max(np.abs(values)) > 1e-10


def test_network_surface_boundary_cycle_orders_scrambled_edges():
    from types import SimpleNamespace
    import pandas as pd

    net = SimpleNamespace(
        surfs=pd.DataFrame({
            "edges": [[12, 3, 7, 5]],
        }),
        edges=pd.DataFrame(
            {
                "verts": [
                    None, None, None,
                    [9, 14],       # edge 3
                    None,
                    [2, 6],        # edge 5
                    None,
                    [2, 9],        # edge 7
                    None, None, None, None,
                    [14, 6],       # edge 12
                ]
            }
        ),
    )

    cycle = network_surface_boundary_cycle(net, 0)

    assert len(cycle) == 4
    assert {edge for edge, _, _ in cycle} == {3, 5, 7, 12}

    for i, (_, _, end_vertex) in enumerate(cycle):
        _, next_start_vertex, _ = cycle[(i + 1) % len(cycle)]
        assert end_vertex == next_start_vertex


def test_network_surface_boundary_cycle_is_independent_of_surface_edge_order():
    from types import SimpleNamespace
    import pandas as pd

    edge_table = pd.DataFrame(
        {
            "verts": [
                None, None, None,
                [9, 14],
                None,
                [2, 6],
                None,
                [2, 9],
                None, None, None, None,
                [14, 6],
            ]
        }
    )

    net_a = SimpleNamespace(
        surfs=pd.DataFrame({"edges": [[12, 3, 7, 5]]}),
        edges=edge_table.copy(),
    )

    net_b = SimpleNamespace(
        surfs=pd.DataFrame({"edges": [[5, 7, 12, 3]]}),
        edges=edge_table.copy(),
    )

    assert network_surface_boundary_cycle(net_a, 0) == \
           network_surface_boundary_cycle(net_b, 0)


def test_aw_face_edge_boundary_orientation_known_case():
    locations = np.array([
        [0.0, 0.0, 0.0],  # cell 10
        [2.0, 0.0, 0.0],  # +x neighbor 20
        [0.0, 2.0, 0.0],  # +y neighbor 30
    ])

    # x=1, y=1 edge traversed from -z to +z.
    edge = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    sign = aw_face_edge_boundary_orientation(
        edge_geometry=edge,
        generator_indices=[10, 20, 30],
        generator_locations=locations,
        cell_index=10,
        face_other_index=20,
    )

    assert sign == 1


def test_aw_face_edge_boundary_orientation_reverses_with_edge_direction():
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    forward = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    reverse = LineEdgeGeometry(
        [1.0, 1.0, 2.0],
        [1.0, 1.0, -2.0],
    )

    kwargs = dict(
        generator_indices=[10, 20, 30],
        generator_locations=locations,
        cell_index=10,
        face_other_index=20,
    )

    sign_forward = aw_face_edge_boundary_orientation(
        edge_geometry=forward,
        **kwargs,
    )

    sign_reverse = aw_face_edge_boundary_orientation(
        edge_geometry=reverse,
        **kwargs,
    )

    assert sign_forward == 1
    assert sign_reverse == -1


def test_oriented_aw_face_edge_geodesic_curvature_is_direction_invariant():
    locations = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])

    forward = LineEdgeGeometry(
        [1.0, 1.0, -2.0],
        [1.0, 1.0, 2.0],
    )

    reverse = LineEdgeGeometry(
        [1.0, 1.0, 2.0],
        [1.0, 1.0, -2.0],
    )

    kwargs = dict(
        generator_indices=[10, 20, 30],
        generator_locations=locations,
        cell_index=10,
        face_other_index=20,
        order=32,
    )

    value_forward = oriented_aw_face_edge_geodesic_curvature(
        edge_geometry=forward,
        **kwargs,
    )

    value_reverse = oriented_aw_face_edge_geodesic_curvature(
        edge_geometry=reverse,
        **kwargs,
    )

    assert value_forward == pytest.approx(value_reverse, abs=1e-14)

@pytest.mark.parametrize("curved", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_fused_edge_curvature_matches_independent_integrals(curved, reverse):
    from vorpy.src.calculations.edge_geometry import aw_edge_curvature_measures

    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    if curved:
        edge = AdditivelyWeightedTrisectorConic(
            locations, np.array([1., 1.5, 2.]), -0.7, 1.1, component=1,
        )
    else:
        edge = LineEdgeGeometry([2., 2.5, -2.], [2., 2.5, 2.])
    if reverse:
        from vorpy.src.calculations.edge_geometry import AffineEdgeGeometry
        edge = AffineEdgeGeometry(edge, edge.t_max, edge.t_min)
    indices = [10, 20, 30]
    pairs = [(cell, other) for cell in indices for other in indices if cell != other]
    result = aw_edge_curvature_measures(edge, indices, locations, pairs)
    for cell in indices:
        reference = aw_cell_edge_mean_curvature(edge, indices, locations, cell)
        assert result['mean'][cell] == pytest.approx(reference, rel=1e-12, abs=1e-12)
    for cell, other in pairs:
        reference = oriented_aw_face_edge_geodesic_curvature(
            edge, indices, locations, cell, other,
        )
        assert result['gaussian'][cell, other] == pytest.approx(
            reference, rel=1e-12, abs=1e-12,
        )


def test_fused_edge_timing_accumulates_per_call(monkeypatch):
    import itertools
    import vorpy.src.calculations.edge_geometry as geometry

    ticks = itertools.count()
    monkeypatch.setattr(geometry, 'perf_counter', lambda: float(next(ticks)))
    edge = LineEdgeGeometry([1., 1., -2.], [1., 1., 2.])
    locations = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]])
    timing = {}
    kwargs = dict(edge_geometry=edge, generator_indices=[10, 20, 30],
                  generator_locations=locations, timing=timing)
    first = geometry.aw_edge_curvature_measures(**kwargs)
    initial = timing.copy()
    second = geometry.aw_edge_curvature_measures(**kwargs)
    assert first == second
    assert initial['other'] > 0
    assert timing == {key: 2 * value for key, value in initial.items()}
