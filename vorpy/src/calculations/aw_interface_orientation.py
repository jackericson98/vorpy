"""Orientation conventions for signed AW interface curvature.

The oriented bicolor interface normal points from group A to group B.  At a
three-generator AW edge, the selected bicolor patches have a common singleton
generator ``s`` and the remaining generator ``q`` on the other side of the
group partition.  Each patch is restricted to the side
``d_s <= d_q``.  Its outward co-normal is therefore the projection of
``grad(d_s-d_q)`` into that patch's tangent plane.

For an oriented patch with unit normal ``n`` and outward co-normal ``m``, the
induced boundary tangent is ``t = n x m``.  This fixes edge direction from the
surface orientation and local domain, rather than from an arbitrary analytic
curve parameter or table ordering.
"""

from __future__ import annotations

import numpy as np

from vorpy.src.calculations.surface_geometry import aw_clearance_gradient


def _unit(vector, label, tolerance):
    vector = np.asarray(vector, dtype=float)
    magnitude = float(np.linalg.norm(vector))
    if vector.shape != (3,) or not np.all(np.isfinite(vector)) or magnitude <= tolerance:
        raise ValueError(f"{label} is undefined or non-finite.")
    return vector / magnitude


def induced_boundary_tangent(normal, outward_conormal, tolerance=1e-10):
    """Orient a patch boundary using ``t = n x m_out``."""
    normal = _unit(normal, "surface normal", tolerance)
    outward_conormal = _unit(outward_conormal, "outward co-normal", tolerance)
    tangent = np.cross(normal, outward_conormal)
    tangent -= normal * float(tangent @ normal)
    return _unit(tangent, "induced boundary tangent", tolerance)


def signed_dihedral(tangent, normal_1, normal_2, tolerance=1e-10):
    """Signed turning angle around an already oriented common edge."""
    tangent = _unit(tangent, "edge tangent", tolerance)
    normal_1 = _unit(normal_1, "first surface normal", tolerance)
    normal_2 = _unit(normal_2, "second surface normal", tolerance)
    return float(np.arctan2(
        tangent @ np.cross(normal_1, normal_2),
        np.clip(normal_1 @ normal_2, -1.0, 1.0),
    ))


def aw_interface_edge_orientation(
        point,
        surface_pairs,
        generator_ids,
        group_a,
        group_b,
        locations,
        *,
        tolerance=1e-10,
        tangent_tolerance=1e-6):
    """Return the signed dihedral and induced tangent for an AW interface edge.

    ``surface_pairs`` must contain the two bicolor pair surfaces incident to
    the edge.  The three edge generators must have a 1/2 group split; the
    singleton is the cell whose local patch domain is used to orient its
    boundary. ``locations`` maps generator ID to a three-vector.

    The geometry tangent is used only as an unoriented consistency check.  Its
    parameter direction has no influence on the returned sign.
    """
    point = np.asarray(point, dtype=float)
    ids = tuple(int(value) for value in generator_ids)
    group_a = {int(value) for value in group_a}
    group_b = {int(value) for value in group_b}
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError("The AW edge point must be a finite three-vector.")
    if len(ids) != 3 or len(set(ids)) != 3:
        raise ValueError("An AW interface edge requires three unique generators.")
    if any((value in group_a) == (value in group_b) for value in ids):
        raise ValueError("Every edge generator must belong to exactly one group.")
    membership = [value in group_a for value in ids]
    if sum(membership) not in (1, 2):
        raise ValueError("A bicolor interface edge requires a 1/2 group split.")
    singleton = ids[membership.index(True)] if sum(membership) == 1 else ids[membership.index(False)]

    pairs = [tuple(int(value) for value in pair) for pair in surface_pairs]
    if len(pairs) != 2 or any(len(pair) != 2 or set(pair) - set(ids) for pair in pairs):
        raise ValueError("Exactly two incident bicolor pair surfaces are required.")
    if any(singleton not in pair for pair in pairs):
        raise ValueError("The two bicolor surfaces must contain the singleton generator.")

    def surface_data(pair):
        other = pair[0] if pair[1] == singleton else pair[1]
        if (other in group_a) == (singleton in group_a):
            raise ValueError("Incident selected surfaces must be bicolor.")
        competitor = next(value for value in ids if value not in pair)
        a_id = singleton if singleton in group_a else other
        b_id = other if singleton in group_a else singleton
        n_ab = _unit(
            aw_clearance_gradient(point, locations[a_id], locations[b_id]),
            "A-to-B interface normal", tolerance,
        )
        # The pair boundary also satisfies d_singleton == d_other.  The patch
        # remains on the singleton-cell side of the third generator.
        grad_domain = aw_clearance_gradient(
            point, locations[singleton], locations[competitor]
        )
        conormal = grad_domain - n_ab * float(grad_domain @ n_ab)
        m_out = _unit(conormal, "AW interface-patch outward co-normal", tolerance)
        tangent = induced_boundary_tangent(n_ab, m_out, tolerance)
        return other, n_ab, m_out, tangent, competitor

    # Use physical generator geometry for deterministic ordering; table IDs
    # and surface-row ordering must not select the sign convention.
    def physical_key(pair):
        other = pair[0] if pair[1] == singleton else pair[1]
        loc = tuple(float(value) for value in np.asarray(locations[other], dtype=float))
        return loc, other in group_a

    ordered_pairs = sorted(pairs, key=physical_key)
    first = surface_data(ordered_pairs[0])
    second = surface_data(ordered_pairs[1])
    alignment = float(first[3] @ second[3])
    if alignment > -1.0 + tangent_tolerance:
        raise ValueError(
            "Incident AW patch boundary tangents are not oppositely oriented "
            f"(dot={alignment:.9g}); edge incidence/orientation is unresolved."
        )

    n1, n2 = first[1], second[1]
    tangent = first[3]
    beta = signed_dihedral(tangent, n1, n2, tolerance)
    return {
        "beta": beta,
        "tangent": tangent,
        "normal_1": n1,
        "normal_2": n2,
        "conormal_1": first[2],
        "conormal_2": second[2],
        "singleton_generator": singleton,
        "competitor_generators": (first[4], second[4]),
        }
