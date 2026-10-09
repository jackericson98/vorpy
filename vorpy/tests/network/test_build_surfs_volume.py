import numpy as np
import pytest

from vorpy.src.calculations import calc_tetra_vol
from vorpy.src.network.build_surfs import _surface_tetra_volumes


def test_batched_surface_volumes_match_scalar_tetrahedron_sums():
    references = np.array([
        [0.25, -0.5, 1.75],
        [2.0, 1.5, -0.75],
    ])
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.5],
    ])
    triangles = np.array([[0, 1, 2], [0, 3, 1], [1, 3, 4], [1, 4, 2]])
    expected = np.array([
        sum(
            calc_tetra_vol(reference, points[triangle[0]], points[triangle[1]], points[triangle[2]])
            for triangle in triangles
        )
        for reference in references
    ])

    actual = _surface_tetra_volumes(references, points, triangles)

    assert actual == pytest.approx(expected, rel=1.0e-14, abs=1.0e-14)


def test_batched_surface_volume_preserves_empty_surface_zero():
    actual = _surface_tetra_volumes(
        [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        np.empty((0, 3)),
        np.empty((0, 3), dtype=int),
    )

    assert np.array_equal(actual, np.zeros(2))


def test_batched_surface_volume_is_invariant_to_triangle_orientation():
    references = np.array([[0.25, -0.5, 1.75]])
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])

    forward = _surface_tetra_volumes(references, points, [[0, 1, 2]])
    reverse = _surface_tetra_volumes(references, points, [[0, 2, 1]])

    assert reverse == pytest.approx(forward, rel=0.0, abs=0.0)
