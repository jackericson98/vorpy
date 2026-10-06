"""Finite explicit-solvent hull-distance audit for sampled interface geometry."""
from __future__ import annotations

import numpy as np
from scipy.spatial import ConvexHull, QhullError


def solvent_hull(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 4 or not np.isfinite(points).all():
        raise ValueError("At least four finite 3D solvent oxygen coordinates are required")
    try:
        hull = ConvexHull(points)
    except QhullError as exc:
        raise ValueError(f"Solvent coordinates do not define a 3D convex hull: {exc}") from exc
    equations = np.asarray(hull.equations, dtype=float)
    return points[hull.vertices], equations


def hull_clearances(sample_points, equations):
    """Signed shortest supporting-plane clearance; negative means outside hull."""
    points = np.asarray(sample_points, dtype=float)
    equations = np.asarray(equations, dtype=float)
    if points.size == 0:
        return np.empty(0)
    normals = equations[:, :3]
    norms = np.linalg.norm(normals, axis=1)
    if not np.isfinite(norms).all() or (norms <= 0).any():
        raise ValueError("Invalid solvent hull planes")
    # ConvexHull uses outward unit normals; normalizing also guards stored hulls.
    signed = -(points @ normals.T + equations[:, 3]) / norms[None, :]
    return np.min(signed, axis=1)


def audit_feature_samples(samples, hull_equations, influence_radii, *, uncertainty_A=0.2):
    """Return per-feature boundary distances and conservative local margins.

    Influence radii are derived from the existing tessellation at sampled
    interface vertices. The 0.2 A allowance is the frozen AW surface resolution.
    Positive margins classify sampled features as interior to the supplied
    solvent hull; a negative margin flags possible crop-edge truncation.
    """
    rows = []
    for feature_id, points in samples.items():
        clearance = hull_clearances(points, hull_equations)
        reach = float(influence_radii[feature_id])
        if not np.isfinite(reach) or reach < 0:
            rows.append({"feature_id": feature_id, "sample_count": len(points),
                         "minimum_hull_distance_A": None,
                         "maximum_influence_radius_A": "Infinity (unresolved/unbounded)" if np.isinf(reach) else "NaN (invalid)",
                         "minimum_boundary_margin_A": None, "boundary_status": "unresolved influence radius"})
            continue
        distance = float(np.min(clearance)) if len(clearance) else None
        margin = distance - reach - uncertainty_A if distance is not None else None
        rows.append({"feature_id": feature_id, "sample_count": len(points),
                     "minimum_hull_distance_A": distance, "maximum_influence_radius_A": reach,
                     "uncertainty_allowance_A": uncertainty_A,
                     "minimum_boundary_margin_A": margin,
                     "boundary_status": "interior_at_sampled_geometry" if margin is not None and margin >= 0 else "boundary_sensitive_or_outside_solvent_hull"})
    return rows
