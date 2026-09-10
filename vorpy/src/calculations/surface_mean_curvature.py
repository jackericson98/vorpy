"""Cell-oriented smooth-surface contributions to integrated mean curvature."""

import numpy as np

from vorpy.src.calculations.surface_geometry import (
    QuadraticSurfaceGeometry,
    aw_clearance_difference,
    aw_outward_normal,
)


def calculate_aw_network_surface_mean_curvatures(
        net,
        clearance_tolerance=1e-6,
        normal_tolerance=1e-6,
        zero_tolerance=1e-14):
    """Calculate cell-oriented smooth-surface mean-curvature contributions.

    VorPy stores one native ``int_mean_curv`` value per pairwise surface.
    That value is oriented according to the native implicit-function gradient.

    A pairwise surface belongs to two AW cells, whose outward normals are
    opposite. Therefore the physical integrated mean-curvature contribution
    must be stored separately for each participating cell.

    This function returns one dictionary per surface:

        {
            ball_a: oriented_integral_for_a,
            ball_b: oriented_integral_for_b,
        }

    The existing ``net.surfs['int_mean_curv']`` field is not modified.
    """
    settings = getattr(net, "settings", None) or {}
    if settings.get("net_type", "aw") != "aw":
        raise ValueError(
            "AW surface mean-curvature orientation requires an AW network."
        )

    if getattr(net, "surfs", None) is None:
        raise ValueError("The network has no surface table.")

    required_columns = {
        "balls",
        "func",
        "com",
        "flat",
        "int_mean_curv",
    }
    missing = required_columns.difference(net.surfs.columns)

    if missing:
        raise ValueError(
            "Surface mean-curvature orientation requires columns: "
            + ", ".join(sorted(missing))
        )

    contributions = []

    for surf_index, surf in net.surfs.iterrows():
        balls = tuple(int(value) for value in surf["balls"])

        if len(balls) != 2:
            raise ValueError(
                f"Surface {surf_index} has {len(balls)} generators; "
                "a regular pairwise AW surface requires exactly two."
            )

        ball_a, ball_b = balls
        native_integral = float(surf["int_mean_curv"])

        if not np.isfinite(native_integral):
            raise ValueError(
                f"Surface {surf_index} has non-finite integrated "
                f"mean curvature: {native_integral}"
            )

        # Flat pairwise surfaces have H = 0 regardless of orientation.
        if bool(surf["flat"]) or abs(native_integral) <= zero_tolerance:
            contributions.append({
                ball_a: 0.0,
                ball_b: 0.0,
            })
            continue

        point = np.asarray(surf["com"], dtype=float)

        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError(
                f"Surface {surf_index} has an invalid representative point."
            )

        loc_a = np.asarray(net.balls.loc[ball_a, "loc"], dtype=float)
        loc_b = np.asarray(net.balls.loc[ball_b, "loc"], dtype=float)

        rad_a = float(net.balls.loc[ball_a, "rad"])
        rad_b = float(net.balls.loc[ball_b, "rad"])

        # The representative point must actually lie on the pairwise AW
        # boundary before we use it to determine orientation.
        clearance_error = abs(
            aw_clearance_difference(
                point=point,
                cell_location=loc_a,
                cell_radius=rad_a,
                neighbor_location=loc_b,
                neighbor_radius=rad_b,
            )
        )

        if (
            not np.isfinite(clearance_error)
            or clearance_error > clearance_tolerance
        ):
            raise ValueError(
                f"Surface {surf_index} representative point is not on the "
                f"AW boundary: clearance error = {clearance_error:.6g}"
            )

        geometry = QuadraticSurfaceGeometry(
            np.asarray(surf["func"], dtype=float)
        )

        native_normal = geometry.normal(point)

        outward_a = aw_outward_normal(
            point=point,
            cell_location=loc_a,
            neighbor_location=loc_b,
        )

        alignment = float(np.dot(native_normal, outward_a))

        # Both constructions describe the same tangent plane, so their unit
        # normals must be parallel or antiparallel.
        if (
            not np.isfinite(alignment)
            or abs(abs(alignment) - 1.0) > normal_tolerance
        ):
            raise ValueError(
                f"Surface {surf_index} native and AW normals disagree: "
                f"dot = {alignment:.12g}"
            )

        orientation_a = 1.0 if alignment >= 0.0 else -1.0

        value_a = orientation_a * native_integral
        value_b = -value_a

        contributions.append({
            ball_a: float(value_a),
            ball_b: float(value_b),
        })

    if len(contributions) != len(net.surfs):
        raise ValueError(
            "Surface mean-curvature orientation did not produce one "
            "result for every surface."
        )

    return contributions