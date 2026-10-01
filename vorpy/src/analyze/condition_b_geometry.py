"""Independent mathematical diagnostics for Cazals condition beta.

For two probe-expanded balls, an orthogonal sphere centered at x exists only
when their common power value is nonnegative.  This helper documents/tests
that geometry; the production facet filter remains in
``power_interface_curvature.py``.
"""
from __future__ import annotations

import math
import numpy as np


def common_power_squared_radius(point, center_i, radius_i, center_j, radius_j):
    """Return the common squared orthogonal-ball radius at a bisector point.

    The input point must lie on the radical plane.  A positive result gives a
    real orthogonal sphere of radius sqrt(q); zero gives the degenerate
    radius-zero sphere; a negative result gives no real orthogonal sphere at
    that center.
    """
    x = np.asarray(point, dtype=float)
    pi, pj = np.asarray(center_i, dtype=float), np.asarray(center_j, dtype=float)
    ri, rj = float(radius_i), float(radius_j)
    qi = float(np.dot(x-pi, x-pi) - ri*ri)
    qj = float(np.dot(x-pj, x-pj) - rj*rj)
    if not np.isclose(qi, qj, rtol=1e-10, atol=1e-10):
        raise ValueError("point is not on the two-ball Power bisector")
    return (qi + qj) / 2.0


def condition_b_decision(candidate_m_squared, smaller_radius, threshold=5.0,
                         *, tolerance=1e-12):
    """Apply the documented rejection rule to candidate facet vertices.

    On a bounded convex power facet, the common power is convex, so its
    maximum is attained at a polygon vertex. If that maximum is negative,
    there is no real orthogonal-ball candidate on the facet. That is distinct
    from a radius-zero sphere (which occurs only at exactly q=0). Since
    condition beta rejects only when a real candidate has m/r > M, absence of
    such a candidate cannot trigger rejection. ``m=0`` is therefore only a
    reporting placeholder used by the production implementation, not the
    radius of a real orthogonal sphere.
    """
    values = tuple(float(q) for q in candidate_m_squared)
    if not values or not all(math.isfinite(q) for q in values):
        raise ValueError("finite candidate squared radii are required")
    r, M = float(smaller_radius), float(threshold)
    if r <= 0 or M <= 0:
        raise ValueError("smaller radius and threshold must be positive")
    qmax = max(values)
    if qmax < -tolerance:
        return {"max_m_squared": qmax, "m": None, "m_over_r": None,
                "has_real_candidate": False, "status": "no_real_orthogonal_ball",
                "accepted": True}
    m = math.sqrt(max(0.0, qmax))
    ratio = m / r
    return {"max_m_squared": qmax, "m": m, "m_over_r": ratio,
            "has_real_candidate": True,
            "status": "zero_radius_candidate" if abs(qmax) <= tolerance else "positive_radius_candidate",
            "accepted": ratio <= M}
