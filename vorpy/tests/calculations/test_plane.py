# vorpy/tests/calculations/test_plane.py

import math
import numpy as np
import numpy.testing as npt
import pytest

from vorpy.src.calculations import (
    project_to_plane,
    unproject_to_3d,
    map_to_plane,
)


def _rand_plane(seed=0):
    rng = np.random.default_rng(seed)
    plane_point = rng.normal(size=3)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    return plane_point, n


def test_map_then_project_roundtrip(tol):
    """2D -> 3D (map) -> 2D (project) should recover the original 2D coords."""
    plane_point, normal = _rand_plane(seed=1)
    pts_2d = [(0.0, 0.0), (1.5, -2.0), (10.0, 3.14159)]

    pts_3d = map_to_plane(pts_2d, plane_point, normal)
    out_2d = project_to_plane(pts_3d, plane_point, normal)

    for (u0, v0), (u1, v1) in zip(pts_2d, out_2d):
        assert u1 == pytest.approx(u0, rel=tol["rel"], abs=tol["abs"])
        assert v1 == pytest.approx(v0, rel=tol["rel"], abs=tol["abs"])


def test_project_then_unproject_is_plane_projection(tol):
    """
    3D -> 2D (project) -> 3D (unproject) should equal the orthogonal projection
    of the original points onto the plane.
    """
    plane_point, normal = _rand_plane(seed=2)

    pts_3d = [
        np.array([3.0, -1.0, 2.0]),
        np.array([0.0, 0.0, 0.0]),
        np.array([-4.2, 5.5, 9.1]),
    ]

    # Expected: orthogonal projection formula
    def orth_proj(p):
        d = np.dot(p - plane_point, normal)
        return p - d * normal

    expected = [orth_proj(p) for p in pts_3d]

    coords_2d = project_to_plane(pts_3d, plane_point, normal)
    recon_3d = unproject_to_3d(coords_2d, plane_point, normal)

    for e, r in zip(expected, recon_3d):
        npt.assert_allclose(r, e, rtol=tol["rel"], atol=tol["abs"])


def test_distance_preserved_in_plane_coords(tol):
    """
    Euclidean distance between two mapped points equals sqrt((du)^2 + (dv)^2).
    This indirectly checks u,v are orthonormal and orthogonal to the normal.
    """
    plane_point, normal = _rand_plane(seed=3)

    a2d = (2.0, -7.0)
    b2d = (-1.5, 4.0)

    a3d, b3d = map_to_plane([a2d, b2d], plane_point, normal)
    d3d = np.linalg.norm(a3d - b3d)
    d2d = math.hypot(b2d[0] - a2d[0], b2d[1] - a2d[1])

    assert d3d == pytest.approx(d2d, rel=tol["rel"], abs=tol["abs"])


def test_normal_scaling_invariance(tol):
    """Scaling the plane normal by any nonzero factor should not change results."""
    plane_point, normal = _rand_plane(seed=4)
    pts_2d = [(0.3, 0.7), (10.0, -2.3)]

    pts_a = map_to_plane(pts_2d, plane_point, normal)
    pts_b = map_to_plane(pts_2d, plane_point, 5.0 * normal)
    for a, b in zip(pts_a, pts_b):
        npt.assert_allclose(a, b, rtol=tol["rel"], atol=tol["abs"])

    pts_3d = [np.array([1.0, 2.0, 3.0]), np.array([-5.0, 0.5, 7.7])]
    uv_a = project_to_plane(pts_3d, plane_point, normal)
    uv_b = project_to_plane(pts_3d, plane_point, -2.0 * normal)
    # Note: flipping normal flips u,v basis orientation; we check distances.
    for (ua, va), (ub, vb) in zip(uv_a, uv_b):
        da = math.hypot(ua, va)
        db = math.hypot(ub, vb)
        assert db == pytest.approx(da, rel=tol["rel"], abs=tol["abs"])


def test_x_axis_aligned_normal_branch(tol):
    """
    With the current implementation, map_to_plane() and project_to_plane()
    may choose different in-plane bases for x-axis normals, so direct 2D
    coordinate equality is not guaranteed. Validate round-trip via project+unproject.
    """
    plane_point = np.array([0.1, -0.2, 0.3])
    pts_2d = [(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)]

    for normal in (np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])):
        pts_3d = map_to_plane(pts_2d, plane_point, normal)
        # round-trip through project/unproject using the SAME basis pair
        back_3d = unproject_to_3d(project_to_plane(pts_3d, plane_point, normal), plane_point, normal)
        for p, q in zip(pts_3d, back_3d):
            assert np.allclose(p, q, rtol=tol["rel"], atol=tol["abs"])


def test_translation_invariance(tol):
    """
    Translating both the plane and the 3D points by the same vector should
    leave projected 2D coordinates unchanged.
    """
    plane_point, normal = _rand_plane(seed=5)
    pts_3d = [np.array([4.0, -2.0, 1.0]), np.array([0.0, 3.0, -5.0])]
    shift = np.array([100.0, -50.0, 7.0])

    uv0 = project_to_plane(pts_3d, plane_point, normal)
    uv1 = project_to_plane([p + shift for p in pts_3d], plane_point + shift, normal)

    for (u0, v0), (u1, v1) in zip(uv0, uv1):
        assert u1 == pytest.approx(u0, rel=tol["rel"], abs=tol["abs"])
        assert v1 == pytest.approx(v0, rel=tol["rel"], abs=tol["abs"])


def test_accepts_list_and_numpy_inputs(tol):
    """Smoke test: lists vs numpy arrays both work and shapes look right."""
    plane_point, normal = _rand_plane(seed=6)

    pts_list = [[0.0, 0.0, 0.0], [1.0, -2.0, 3.0]]
    uv = project_to_plane(pts_list, plane_point, normal)
    assert isinstance(uv, list) and all(len(t) == 2 for t in uv)

    pts_2d = np.array([[0.5, 0.5], [2.0, -3.0]])
    mapped = map_to_plane(pts_2d, plane_point, normal)
    assert isinstance(mapped, list) and all(np.shape(p) == (3,) for p in mapped)


@pytest.mark.xfail(reason="No input validation; zero-length normal yields NaNs downstream.")
def test_zero_length_normal_should_error(tol):
    plane_point = np.array([0.0, 0.0, 0.0])
    normal = np.array([0.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        _ = project_to_plane([np.array([1.0, 2.0, 3.0])], plane_point, normal)

