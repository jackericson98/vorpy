import numpy as np
from numba import njit


@njit(cache=True)
def _norm3(v):
    return np.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


@njit(cache=True)
def _dot3(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


@njit(cache=True)
def _cross3(a, b):
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ], dtype=np.float64)


@njit(cache=True)
def _draw_line_numba(pts, radius):
    """
    Numba-compiled triangular tube construction.

    Returns
    -------
    draw_points : ndarray, shape (3 * N, 3)
    draw_tris : ndarray, shape (6 * (N - 1), 3)
    error_code : int
        0 = success
        1 = degenerate tangent
        2 = invalid initial normal
        3 = unstable transported normal
        4 = invalid binormal
    error_index : int
        Point index associated with an error.
    """
    n_pts = len(pts)
    eps = 1e-12

    if n_pts < 2:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64), 0, -1

    tangents = np.empty((n_pts, 3), dtype=np.float64)
    normals = np.empty((n_pts, 3), dtype=np.float64)
    binormals = np.empty((n_pts, 3), dtype=np.float64)

    # Tangents
    tangents[0] = pts[1] - pts[0]
    tangents[n_pts - 1] = pts[n_pts - 1] - pts[n_pts - 2]

    for i in range(1, n_pts - 1):
        tangents[i] = pts[i + 1] - pts[i - 1]

    for i in range(n_pts):
        norm_t = _norm3(tangents[i])

        if norm_t < eps:
            return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64), 1, i

        tangents[i, 0] /= norm_t
        tangents[i, 1] /= norm_t
        tangents[i, 2] /= norm_t

    # Initial frame
    t0 = tangents[0]

    ref = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    if abs(_dot3(ref, t0)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)

    n0 = _cross3(t0, ref)
    norm_n0 = _norm3(n0)

    if norm_n0 < eps:
        ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        n0 = _cross3(t0, ref)
        norm_n0 = _norm3(n0)

        if norm_n0 < eps:
            return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64), 2, 0

    n0 /= norm_n0

    b0 = _cross3(t0, n0)
    norm_b0 = _norm3(b0)

    if norm_b0 < eps:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64), 4, 0

    b0 /= norm_b0

    normals[0] = n0
    binormals[0] = b0

    # Parallel transport
    for i in range(1, n_pts):
        t_prev = tangents[i - 1]
        t_curr = tangents[i]
        n_prev = normals[i - 1]

        axis = _cross3(t_prev, t_curr)
        axis_norm = _norm3(axis)

        if axis_norm < eps:
            n_curr = n_prev.copy()

        else:
            axis /= axis_norm

            dot_tt = _dot3(t_prev, t_curr)

            if dot_tt < -1.0:
                dot_tt = -1.0
            elif dot_tt > 1.0:
                dot_tt = 1.0

            angle = np.arccos(dot_tt)
            sin_angle = np.sin(angle)
            cos_angle = np.cos(angle)

            axis_cross_n = _cross3(axis, n_prev)
            axis_dot_n = _dot3(axis, n_prev)

            n_curr = (
                n_prev * cos_angle
                + axis_cross_n * sin_angle
                + axis * axis_dot_n * (1.0 - cos_angle)
            )

        proj = _dot3(n_curr, t_curr)

        n_curr[0] -= proj * t_curr[0]
        n_curr[1] -= proj * t_curr[1]
        n_curr[2] -= proj * t_curr[2]

        norm_n_curr = _norm3(n_curr)

        if norm_n_curr < eps:
            n_curr = normals[i - 1].copy()
            proj = _dot3(n_curr, t_curr)

            n_curr[0] -= proj * t_curr[0]
            n_curr[1] -= proj * t_curr[1]
            n_curr[2] -= proj * t_curr[2]

            norm_n_curr = _norm3(n_curr)

            if norm_n_curr < eps:
                return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64), 3, i

        n_curr /= norm_n_curr

        b_curr = _cross3(t_curr, n_curr)
        norm_b_curr = _norm3(b_curr)

        if norm_b_curr < eps:
            return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64), 4, i

        b_curr /= norm_b_curr

        normals[i] = n_curr
        binormals[i] = b_curr

    # Tube vertices
    draw_points = np.empty((n_pts * 3, 3), dtype=np.float64)
    sqrt3_half = 0.5 * np.sqrt(3.0)

    for i in range(n_pts):
        p = pts[i]
        n = normals[i]
        b = binormals[i]

        j = 3 * i

        draw_points[j, 0] = p[0] + radius * n[0]
        draw_points[j, 1] = p[1] + radius * n[1]
        draw_points[j, 2] = p[2] + radius * n[2]

        draw_points[j + 1, 0] = p[0] + radius * (-0.5 * n[0] + sqrt3_half * b[0])
        draw_points[j + 1, 1] = p[1] + radius * (-0.5 * n[1] + sqrt3_half * b[1])
        draw_points[j + 1, 2] = p[2] + radius * (-0.5 * n[2] + sqrt3_half * b[2])

        draw_points[j + 2, 0] = p[0] + radius * (-0.5 * n[0] - sqrt3_half * b[0])
        draw_points[j + 2, 1] = p[1] + radius * (-0.5 * n[1] - sqrt3_half * b[1])
        draw_points[j + 2, 2] = p[2] + radius * (-0.5 * n[2] - sqrt3_half * b[2])

    # Tube triangles
    draw_tris = np.empty(((n_pts - 1) * 6, 3), dtype=np.int64)

    for i in range(n_pts - 1):
        a0 = 3 * i
        a1 = a0 + 1
        a2 = a0 + 2

        b0i = a0 + 3
        b1 = a0 + 4
        b2 = a0 + 5

        j = 6 * i

        draw_tris[j] = (a0, a1, b0i)
        draw_tris[j + 1] = (b0i, b1, a1)
        draw_tris[j + 2] = (a1, a2, b1)
        draw_tris[j + 3] = (b1, b2, a2)
        draw_tris[j + 4] = (a2, a0, b2)
        draw_tris[j + 5] = (b2, b0i, a0)

    return draw_points, draw_tris, 0, -1


def draw_line(points, radius=0.05, color="#000000", edge_org=None):
    """
    Create triangular tube geometry around a 3D polyline.

    edge_org is retained for compatibility with the previous API.
    """
    if len(points) < 2:
        return [], []

    pts = np.asarray(points, dtype=np.float64)

    draw_points, draw_tris, error_code, error_index = _draw_line_numba(pts, float(radius))

    if error_code == 1:
        raise ValueError(f"Degenerate segment near point index {error_index}.")
    elif error_code == 2:
        raise ValueError("Could not construct an initial normal vector.")
    elif error_code == 3:
        raise ValueError(f"Could not stabilize frame at point index {error_index}.")
    elif error_code == 4:
        raise ValueError(f"Could not construct binormal at point index {error_index}.")

    # Keep the old output type so the rest of the exporter behaves identically.
    return draw_points.tolist(), draw_tris.tolist()


def draw_edge(edge, radius=0.05, color=None):
    """Generate drawable triangular tube geometry for an edge."""
    return draw_line(edge.points, radius=radius, color=color)