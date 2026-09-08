import numpy as np
import pytest

from vorpy.src.calculations.edge_geometry import (
    LineEdgeGeometry,
    ParametricEdgeGeometry,
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
