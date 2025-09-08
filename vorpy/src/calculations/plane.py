import numpy as np


__all__ = [
    "project_to_plane",
    "unproject_to_3d",
    "map_to_plane",
]


def _normalize_normal(plane_normal: np.ndarray) -> np.ndarray:
    """
    Normalize a plane normal and validate it is finite and non-zero.
    """
    plane_normal = np.asarray(plane_normal, dtype=float)
    if not np.isfinite(plane_normal).all():
        raise ValueError("plane_normal must be finite.")
    norm = np.linalg.norm(plane_normal)
    if norm == 0.0:
        raise ValueError("plane_normal must be non-zero.")
    return plane_normal / norm


def _plane_basis(n_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Build an orthonormal basis (u, v) for the plane whose unit normal is n_hat.
    Chooses a stable axis to cross with to avoid near-collinearity.
    """
    n_hat = np.asarray(n_hat, dtype=float)

    # Choose the axis least aligned with n_hat to keep cross product stable.
    ax = np.argmax(np.abs(n_hat))
    if ax == 0:
        # Normal mostly along x -> use y-axis to build u
        a = np.array([0.0, 1.0, 0.0])
    elif ax == 1:
        # Normal mostly along y -> use z-axis
        a = np.array([0.0, 0.0, 1.0])
    else:
        # Normal mostly along z -> use x-axis
        a = np.array([1.0, 0.0, 0.0])

    u = np.cross(n_hat, a)
    u_norm = np.linalg.norm(u)
    if u_norm == 0.0:
        # Extremely pathological, but guard anyway
        a = np.array([1.0, 0.0, 0.0])
        u = np.cross(n_hat, a)
        u_norm = np.linalg.norm(u)
        if u_norm == 0.0:
            raise ValueError("Failed to construct plane basis from normal.")

    u /= u_norm
    v = np.cross(n_hat, u)
    v /= np.linalg.norm(v)
    return u, v


def project_to_plane(points, plane_point, plane_normal):
    """
    Project 3D points onto a plane and return 2D (u, v) coordinates in that plane.

    Parameters
    ----------
    points : Iterable[ArrayLike]
        Sequence of 3D points to project. Each point must be length-3.
    plane_point : ArrayLike
        A 3D point lying on the plane.
    plane_normal : ArrayLike
        The (non-zero) normal vector of the plane.

    Returns
    -------
    list[tuple[float, float]]
        (u, v) coordinates of each projected point in an orthonormal in-plane basis.

    Notes
    -----
    - The basis (u, v) is orthonormal and constructed deterministically from the normal.
    - Input validity: normal must be finite and non-zero.
    """
    n_hat = _normalize_normal(plane_normal)
    plane_point = np.asarray(plane_point, dtype=float)
    u, v = _plane_basis(n_hat)

    projected_points = []
    for point in points:
        point = np.asarray(point, dtype=float)
        pv = point - plane_point
        projected_points.append((float(np.dot(pv, u)), float(np.dot(pv, v))))

    return projected_points


def unproject_to_3d(projected_points, plane_point, plane_normal):
    """
    Reconstruct 3D points from their 2D (u, v) coordinates on a plane.

    Parameters
    ----------
    projected_points : Iterable[tuple[float, float] | ArrayLike]
        Sequence of (u, v) coordinates previously obtained by projecting onto the plane.
    plane_point : ArrayLike
        A 3D point lying on the plane.
    plane_normal : ArrayLike
        The (non-zero) normal vector of the plane.

    Returns
    -------
    list[np.ndarray]
        Reconstructed 3D points (shape (3,)) that lie on the plane.

    Notes
    -----
    - This is the inverse of `project_to_plane` (up to floating point).
    - The same deterministic orthonormal basis is used.
    """
    n_hat = _normalize_normal(plane_normal)
    plane_point = np.asarray(plane_point, dtype=float)
    u, v = _plane_basis(n_hat)

    reconstructed = []
    for uv in projected_points:
        u_coord, v_coord = np.asarray(uv, dtype=float)
        p3 = plane_point + u_coord * u + v_coord * v
        reconstructed.append(p3)

    return reconstructed


def map_to_plane(points_2d, plane_point, plane_normal):
    """
    Map arbitrary 2D (u, v) coordinates into 3D points on a plane.

    Parameters
    ----------
    points_2d : Iterable[tuple[float, float] | ArrayLike]
        Sequence of (u, v) coordinates to place on the plane.
    plane_point : ArrayLike
        A 3D point lying on the plane.
    plane_normal : ArrayLike
        The (non-zero) normal vector of the plane.

    Returns
    -------
    list[np.ndarray]
        3D points (shape (3,)) corresponding to the given (u, v) coordinates.

    Notes
    -----
    - Uses the same deterministic orthonormal basis (u, v) defined by the plane normal.
    """
    n_hat = _normalize_normal(plane_normal)
    plane_point = np.asarray(plane_point, dtype=float)
    u, v = _plane_basis(n_hat)

    mapped = []
    for uv in points_2d:
        u_coord, v_coord = np.asarray(uv, dtype=float)
        p3 = plane_point + u_coord * u + v_coord * v
        mapped.append(p3)

    return mapped
