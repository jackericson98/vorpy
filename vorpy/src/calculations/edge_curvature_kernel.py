"""Compiled quadrature accumulation shared by AW edge M and G."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def integrate_edge_samples(points, firsts, seconds, locations, pairs, weights, tol):
    """Accumulate in node order without relaxing floating-point precision."""
    mean = np.zeros(3)
    gaussian = np.zeros(len(pairs))
    for node in range(len(weights)):
        first = firsts[node]
        speed_sq = np.dot(first, first)
        if not np.isfinite(speed_sq) or speed_sq < tol * tol:
            raise ValueError("Edge tangent is undefined during AW curvature integration.")
        radial = points[node] - locations
        for cell in range(3):
            distance = np.sqrt(np.dot(radial[cell], radial[cell]))
            if not np.isfinite(distance) or distance < tol:
                raise ValueError("AW edge point coincides with a generator.")
            radial[cell] /= distance
        normals = np.empty((3, 3, 3))
        for cell in range(3):
            for other in range(cell + 1, 3):
                normal = radial[cell] - radial[other]
                length = np.sqrt(np.dot(normal, normal))
                if length < tol:
                    raise ValueError("AW pairwise surface normal is undefined on edge.")
                normals[cell, other] = normal / length
                normals[other, cell] = -normals[cell, other]
        speed = np.sqrt(speed_sq)
        for cell in range(3):
            other1 = (cell + 1) % 3
            other2 = (cell + 2) % 3
            cosine = np.dot(normals[cell, other1], normals[cell, other2])
            mean[cell] += weights[node] * np.arccos(min(1.0, max(-1.0, cosine))) * speed
        if len(pairs):
            second = seconds[node]
            if not np.all(np.isfinite(second)):
                raise ValueError("Edge second derivative is non-finite.")
            for face in range(len(pairs)):
                cell, other = pairs[face]
                value = np.dot(second, np.cross(normals[cell, other], first)) / speed_sq
                if not np.isfinite(value):
                    raise ValueError("AW geodesic-curvature integrand is non-finite.")
                gaussian[face] += weights[node] * value
    return mean, gaussian
