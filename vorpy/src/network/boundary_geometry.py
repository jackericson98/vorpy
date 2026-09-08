

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass(frozen=True)
class BoundarySurfaceView:
    """One outward-facing surface perspective adjoining a barrier edge."""

    surface_index: int
    group_ball_index: int
    external_ball_index: int


@dataclass(frozen=True)
class BarrierEdgeTopology:
    """Topology needed before any signed edge geometry is evaluated."""

    classification: str
    generator_indices: tuple
    group_generators: tuple
    external_generators: tuple
    boundary_surfaces: tuple


def outward_surface_normal(surface_geometry, point, group_ball_location,
                           external_ball_location, tol=1e-12):
    """Orient an implicit-surface normal from group ball toward external ball."""
    native = np.asarray(surface_geometry.normal(point), dtype=float)
    outward = (
        np.asarray(external_ball_location, dtype=float)
        - np.asarray(group_ball_location, dtype=float)
    )
    outward_norm = float(np.linalg.norm(outward))
    if not np.isfinite(outward_norm) or outward_norm < tol:
        raise ValueError("Group and external ball locations do not define an orientation.")
    dot = float(native @ outward)
    if abs(dot) <= tol*outward_norm:
        raise ValueError("Surface orientation is ambiguous at this point.")
    return native if dot > 0.0 else -native


def signed_surface_intersection_angle(tangent, normal_1, normal_2, tol=1e-12):
    """Return the signed angle between two surface normals about an edge tangent."""
    tangent = np.asarray(tangent, dtype=float)
    normal_1 = np.asarray(normal_1, dtype=float)
    normal_2 = np.asarray(normal_2, dtype=float)
    vectors = (tangent, normal_1, normal_2)
    magnitudes = [float(np.linalg.norm(vector)) for vector in vectors]
    if any(not np.isfinite(value) or value < tol for value in magnitudes):
        raise ValueError("Tangent and normals must be finite, nonzero vectors.")
    tangent = tangent/magnitudes[0]
    normal_1 = normal_1/magnitudes[1]
    normal_2 = normal_2/magnitudes[2]
    sine = float(tangent @ np.cross(normal_1, normal_2))
    cosine = float(np.clip(normal_1 @ normal_2, -1.0, 1.0))
    return float(np.arctan2(sine, cosine))


def _surface_balls(surfaces, surface_index):
    """Read a surface's generator pair from a DataFrame or mapping/sequence."""
    if hasattr(surfaces, "loc"):
        try:
            row = surfaces.loc[surface_index]
        except KeyError:
            row = surfaces.iloc[int(surface_index)]
    else:
        row = surfaces[surface_index]
    balls = row["balls"] if hasattr(row, "__getitem__") else row.balls
    return tuple(sorted(int(ball) for ball in balls))


def get_boundary_surfaces(edge_balls, edge_surface_indices, surfaces,
                          group_indices, strict=True):
    """Classify an edge and select its two group/external boundary surfaces.

    The result is invariant to generator ordering, incident-surface ordering,
    and group-index ordering. Edges with zero or three group generators are not
    barrier edges and return ``None``.
    """
    generators = tuple(sorted(int(ball) for ball in edge_balls))
    if len(generators) != 3 or len(set(generators)) != 3:
        raise ValueError("An edge must have exactly three distinct generators.")

    group = {int(ball) for ball in group_indices}
    group_generators = tuple(ball for ball in generators if ball in group)
    external_generators = tuple(ball for ball in generators if ball not in group)

    if len(group_generators) in (0, 3):
        return None
    classification = (
        "1_group_2_external"
        if len(group_generators) == 1
        else "2_group_1_external"
    )

    expected_pairs = {
        tuple(sorted((group_ball, external_ball))):
            (group_ball, external_ball)
        for group_ball, external_ball
        in product(group_generators, external_generators)
    }
    matches = {pair: [] for pair in expected_pairs}
    for surface_index in sorted(int(index) for index in edge_surface_indices):
        pair = _surface_balls(surfaces, surface_index)
        if pair in matches:
            matches[pair].append(surface_index)

    invalid = {pair: indices for pair, indices in matches.items() if len(indices) != 1}
    if invalid and strict:
        raise ValueError(
            "Barrier edge does not have exactly one incident surface for each "
            f"group/external pair: {invalid}"
        )

    views = []
    for pair in sorted(expected_pairs):
        if len(matches[pair]) != 1:
            continue
        group_ball, external_ball = expected_pairs[pair]
        views.append(BoundarySurfaceView(
            surface_index=matches[pair][0],
            group_ball_index=group_ball,
            external_ball_index=external_ball,
        ))

    return BarrierEdgeTopology(
        classification=classification,
        generator_indices=generators,
        group_generators=group_generators,
        external_generators=external_generators,
        boundary_surfaces=tuple(views),
    )
