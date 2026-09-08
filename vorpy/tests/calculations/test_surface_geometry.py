import numpy as np
import pytest

from vorpy.src.calculations.curvature import gaussian_curvature, mean_curvature
from vorpy.src.calculations.surface_geometry import (
    PlaneSurfaceGeometry,
    QuadraticSurfaceGeometry,
)


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
