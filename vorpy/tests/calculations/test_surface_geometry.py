import numpy as np
import pytest

from vorpy.src.calculations import gaussian_curvature
from vorpy.src.calculations import mean_curvature
from vorpy.src.calculations import calc_surf_func

from vorpy.src.calculations import PlaneSurfaceGeometry
from vorpy.src.calculations import QuadraticSurfaceGeometry
from vorpy.src.calculations import aw_clearance_difference
from vorpy.src.calculations import aw_outward_normal


def test_quadratic_api_matches_existing_curvature_functions():
    coefficients = np.array([
        1.2, 0.8, 1.5, 0.1, -0.2, 0.3,
        -0.5, 0.7, 0.2, -4.0, 0.0, 0.0, 0.0, 0.0,
    ])
    point = np.array([0.6, -0.2, 1.1])
    geometry = QuadraticSurfaceGeometry(coefficients)

    assert geometry.mean_curvature(point) == pytest.approx(
        mean_curvature(coefficients, point), rel=1e-13, abs=1e-13
    )
    assert geometry.gaussian_curvature(point) == pytest.approx(
        gaussian_curvature(coefficients, point), rel=1e-13, abs=1e-13
    )


def test_sphere_has_expected_curvature_and_unit_normal():
    radius = 2.0
    geometry = QuadraticSurfaceGeometry([
        1, 1, 1, 0, 0, 0, 0, 0, 0, -radius**2, 0, 0, 0, 0
    ])
    point = np.array([radius, 0.0, 0.0])

    assert geometry.implicit_value(point) == pytest.approx(0.0)
    assert geometry.normal(point) == pytest.approx([1.0, 0.0, 0.0])
    assert geometry.normal(point, orientation=-1) == pytest.approx([-1.0, 0.0, 0.0])
    assert geometry.mean_curvature(point) == pytest.approx(1.0/radius)
    assert geometry.gaussian_curvature(point) == pytest.approx(1.0/radius**2)


def test_plane_geometry_is_flat():
    geometry = PlaneSurfaceGeometry([0, 0, 2], [0, 0, 5])
    assert geometry.implicit_value([3, -1, 2]) == pytest.approx(0.0)
    assert geometry.normal([3, -1, 2]) == pytest.approx([0, 0, 1])
    assert geometry.mean_curvature([3, -1, 2]) == 0.0
    assert geometry.gaussian_curvature([3, -1, 2]) == 0.0


def test_undefined_quadratic_normal_is_explicit_error():
    geometry = QuadraticSurfaceGeometry(np.zeros(14))
    with pytest.raises(ValueError, match="undefined"):
        geometry.normal([0, 0, 0])


def test_aw_clearance_difference_identifies_pairwise_boundary():
    point = np.array([1.5, 0.0, 0.0])

    value = aw_clearance_difference(
        point=point,
        cell_location=[0.0, 0.0, 0.0],
        cell_radius=1.0,
        neighbor_location=[4.0, 0.0, 0.0],
        neighbor_radius=2.0,
    )

    assert value == pytest.approx(0.0, abs=1e-12)


def test_aw_outward_normal_points_from_cell_to_neighbor_for_planar_case():
    normal = aw_outward_normal(
        point=[1.0, 0.0, 0.0],
        cell_location=[0.0, 0.0, 0.0],
        neighbor_location=[2.0, 0.0, 0.0],
    )

    assert normal == pytest.approx([1.0, 0.0, 0.0])


def test_swapping_aw_cell_and_neighbor_reverses_normal():
    point = np.array([1.5, 0.3, -0.2])
    first = np.array([0.0, 0.0, 0.0])
    second = np.array([4.0, 1.0, 0.0])

    forward = aw_outward_normal(point, first, second)
    reverse = aw_outward_normal(point, second, first)

    assert reverse == pytest.approx(-forward, abs=1e-13)


def test_aw_outward_normal_is_unit_length():
    normal = aw_outward_normal(
        point=[1.2, 0.7, -0.1],
        cell_location=[0.0, 0.0, 0.0],
        neighbor_location=[3.0, 1.0, 0.4],
    )

    assert np.linalg.norm(normal) == pytest.approx(1.0, abs=1e-13)


def test_aw_outward_normal_is_rigid_transform_covariant():
    point = np.array([1.4, 0.6, 0.2])
    first = np.array([0.0, 0.0, 0.0])
    second = np.array([3.0, 1.0, -0.4])

    theta = 0.71
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta),  np.cos(theta), 0.0],
        [0.0,            0.0,           1.0],
    ])
    translation = np.array([5.0, -2.0, 4.0])

    original = aw_outward_normal(point, first, second)

    transformed = aw_outward_normal(
        rotation @ point + translation,
        rotation @ first + translation,
        rotation @ second + translation,
    )

    assert transformed == pytest.approx(
        rotation @ original,
        abs=1e-12,
    )


def test_aw_quadratic_and_clearance_normals_define_same_tangent_plane():
    first = np.array([0.0, 0.0, 0.0])
    second = np.array([4.0, 0.0, 0.0])
    first_radius = 1.0
    second_radius = 2.0

    # Exact point on the AW pairwise surface along the center-center axis.
    point = np.array([1.5, 0.0, 0.0])

    coefficients = calc_surf_func(
        first,
        first_radius,
        second,
        second_radius,
    )
    quadratic = QuadraticSurfaceGeometry(coefficients)

    quadratic_normal = quadratic.normal(point)
    clearance_normal = aw_outward_normal(
        point,
        first,
        second,
    )

    alignment = abs(
        float(np.dot(quadratic_normal, clearance_normal))
    )

    assert alignment == pytest.approx(1.0, abs=1e-12)


def test_aw_outward_normal_is_undefined_at_generator_center():
    with pytest.raises(ValueError, match="generator center"):
        aw_outward_normal(
            point=[0.0, 0.0, 0.0],
            cell_location=[0.0, 0.0, 0.0],
            neighbor_location=[2.0, 0.0, 0.0],
        )


def test_aw_normal_does_not_depend_on_smaller_ball_construction_order():
    small_location = np.array([0.0, 0.0, 0.0])
    large_location = np.array([4.0, 0.0, 0.0])

    small_radius = 1.0
    large_radius = 2.0

    # AW equality:
    # x - 1 = (4 - x) - 2
    # x = 1.5
    point = np.array([1.5, 0.0, 0.0])

    outward_from_large = aw_outward_normal(
        point,
        cell_location=large_location,
        neighbor_location=small_location,
    )

    outward_from_small = aw_outward_normal(
        point,
        cell_location=small_location,
        neighbor_location=large_location,
    )

    assert outward_from_large == pytest.approx([-1.0, 0.0, 0.0])
    assert outward_from_small == pytest.approx([1.0, 0.0, 0.0])


def test_sphere_integrated_mean_curvature_is_four_pi_r():
    """A smooth sphere carries all of M on its surface."""
    radius = 2.7

    geometry = QuadraticSurfaceGeometry([
        1.0, 1.0, 1.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        -radius**2,
        0.0, 0.0, 0.0, 0.0,
    ])

    # Any point on a sphere has H = 1/R for the outward orientation.
    sample_points = np.array([
        [ radius, 0.0, 0.0],
        [-radius, 0.0, 0.0],
        [0.0,  radius, 0.0],
        [0.0, -radius, 0.0],
        [0.0, 0.0,  radius],
        [0.0, 0.0, -radius],
    ])

    for point in sample_points:
        assert geometry.mean_curvature(point) == pytest.approx(
            1.0 / radius,
            rel=1e-13,
            abs=1e-13,
        )

    surface_area = 4.0 * np.pi * radius**2
    integrated_mean_curvature = (1.0 / radius) * surface_area

    assert integrated_mean_curvature == pytest.approx(
        4.0 * np.pi * radius,
        rel=1e-13,
        abs=1e-13,
    )
    