"""Convex vertex bounds derived from PDB ``CRYST1`` lattice metadata."""

from itertools import product

import numpy as np
from scipy.spatial import HalfspaceIntersection, QhullError


def vertex_region_from_cryst1(line):
    """Return the centered Wigner-Seitz cell encoded by one CRYST1 record.

    Lengths are in Angstroms and angles are in degrees, as required by PDB.
    The returned plane equations use ``normal . (point - center) <= offset``.
    """
    try:
        lengths = np.array([
            float(line[6:15]), float(line[15:24]), float(line[24:33]),
        ])
        alpha, beta, gamma = np.deg2rad([
            float(line[33:40]), float(line[40:47]), float(line[47:54]),
        ])
    except (TypeError, ValueError) as error:
        raise ValueError("CRYST1 must contain three lengths and three angles") from error
    if np.any(lengths <= 0) or not all(0 < angle < np.pi for angle in (alpha, beta, gamma)):
        raise ValueError("CRYST1 lengths must be positive and angles must be between 0 and 180 degrees")

    sin_gamma = np.sin(gamma)
    if abs(sin_gamma) < 1e-12:
        raise ValueError("CRYST1 gamma angle produces a singular lattice")
    a, b, c = lengths
    c_y = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sin_gamma
    c_z_squared = c * c - (c * np.cos(beta)) ** 2 - c_y ** 2
    if c_z_squared <= 0:
        raise ValueError("CRYST1 parameters do not define a three-dimensional cell")
    vectors = np.array([
        [a, 0.0, 0.0],
        [b * np.cos(gamma), b * sin_gamma, 0.0],
        [c * np.cos(beta), c_y, np.sqrt(c_z_squared)],
    ])

    # shortcut: two lattice shells cover ordinary PDB cells; expand the search
    # if highly skewed, non-reduced crystallographic bases must be supported.
    lattice = np.array([
        coefficients @ vectors
        for coefficients in product(range(-2, 3), repeat=3)
        if coefficients != (0, 0, 0)
    ])
    offsets = np.einsum("ij,ij->i", lattice, lattice) / 2.0
    try:
        vertices = HalfspaceIntersection(
            np.column_stack((lattice, -offsets)), np.zeros(3)
        ).intersections
    except QhullError as error:
        raise ValueError("CRYST1 lattice does not produce a bounded unit cell") from error

    facet_indices = []
    for index, (normal, offset) in enumerate(zip(lattice, offsets)):
        on_plane = vertices[
            np.abs(vertices @ normal - offset) <= 1e-7 * max(1.0, offset)
        ]
        if len(on_plane) >= 3 and np.linalg.matrix_rank(
            on_plane[1:] - on_plane[0], tol=1e-7
        ) >= 2:
            facet_indices.append(index)

    # shortcut: CRYST1 is treated as a conventional zero-origin cell; add
    # SCALE/ORIGX support if transformed crystallographic inputs require it.
    return {
        "kind": "cryst1_wigner_seitz",
        "center": (vectors.sum(axis=0) / 2.0).tolist(),
        "normals": lattice[facet_indices].tolist(),
        "offsets": offsets[facet_indices].tolist(),
        "cell_vectors": vectors.tolist(),
        "face_count": len(facet_indices),
    }


def scaled_vertex_region(region, multiplier):
    """Scale a CRYST1 region about its center by the box multiplier."""
    if not region:
        return None
    scaled = dict(region)
    scaled["offsets"] = [float(multiplier) * value for value in region["offsets"]]
    scaled["multiplier"] = float(multiplier)
    return scaled


def point_in_vertex_box(point, box):
    """Test a CRYST1 convex region, or the legacy AABB when absent."""
    if point is None:
        return False
    if box is not None and len(box) >= 3 and box[2]:
        region = box[2]
        delta = np.asarray(point, dtype=float) - np.asarray(region["center"], dtype=float)
        normals = np.asarray(region["normals"], dtype=float)
        offsets = np.asarray(region["offsets"], dtype=float)
        tolerance = 1e-9 * max(1.0, float(np.max(offsets)))
        return bool(np.all(normals @ delta <= offsets + tolerance))
    return box is None or not any(
        box[0][axis] > point[axis] or point[axis] > box[1][axis]
        for axis in range(3)
    )
