import os
import numpy as np
from vorpy.src.output import write_off_verts
from vorpy.src.output import write_edges
from vorpy.src.output import write_surfs
from vorpy.src.output import write_pdb
from vorpy.src.output import write_interface_logs


def get_interface_atoms(iface):
    """Return the unique atom indices participating in the interface definition."""
    return sorted(
        set(getattr(iface, "group1_indices", []) or [])
        | set(getattr(iface, "group2_indices", []) or [])
    )

def _get_column(dataframe, *candidate_names):
    """
    Return the first matching DataFrame column from candidate_names.

    Matching is case-insensitive and ignores spaces and underscores so this
    works with both internal Network column names and exported-log headings.
    """
    if dataframe is None:
        return None

    normalized_columns = {
        str(column).lower().replace(" ", "").replace("_", ""): column
        for column in dataframe.columns
    }

    for candidate in candidate_names:
        normalized_candidate = (
            str(candidate)
            .lower()
            .replace(" ", "")
            .replace("_", "")
        )

        if normalized_candidate in normalized_columns:
            return normalized_columns[normalized_candidate]

    return None


def _numeric_series(dataframe, *candidate_names):
    """
    Return a numeric Series for the first matching column, or None.
    """
    column = _get_column(dataframe, *candidate_names)

    if column is None:
        return None

    return dataframe[column].dropna()


def _safe_sum(series, empty_value=0.0):
    if series is None:
        return None

    if len(series) == 0:
        return empty_value

    return float(series.sum())


def _safe_mean(series):
    return None if series is None or len(series) == 0 else float(series.mean())


def _safe_min(series):
    return None if series is None or len(series) == 0 else float(series.min())


def _safe_max(series):
    return None if series is None or len(series) == 0 else float(series.max())


def _weighted_mean(values, weights):
    """
    Return the area-weighted mean of values.

    Rows with missing values, missing weights, or non-positive weights are
    excluded.
    """
    if values is None or weights is None:
        return None

    valid_indices = values.index.intersection(weights.index)

    if len(valid_indices) == 0:
        return None

    valid_values = values.loc[valid_indices]
    valid_weights = weights.loc[valid_indices]

    valid_mask = (
        valid_values.notna()
        & valid_weights.notna()
        & (valid_weights > 0)
    )

    valid_values = valid_values[valid_mask]
    valid_weights = valid_weights[valid_mask]

    if len(valid_values) == 0 or float(valid_weights.sum()) == 0:
        return None

    return float(
        (valid_values * valid_weights).sum()
        / valid_weights.sum()
    )


def _format_metric(value, decimals=6):
    if value is None:
        return "Not available"

    return f"{value:.{decimals}f}"



def _build_interface_orientation_context(iface):
    """Build topology -> system/location/radius maps for an Interface network."""
    net_balls = iface.net.balls
    columns = list(net_balls.columns)
    num_idx = columns.index("num") if "num" in columns else None
    system_idx = columns.index("system_num") if "system_num" in columns else None
    loc_idx = columns.index("loc") if "loc" in columns else None
    rad_idx = columns.index("rad") if "rad" in columns else None

    topology_to_system = {}
    topology_to_loc = {}
    topology_to_rad = {}

    for row_pos, row in enumerate(net_balls.itertuples(index=False, name=None)):
        try:
            topology_id = int(row[num_idx]) if num_idx is not None else int(row_pos)
        except (TypeError, ValueError):
            topology_id = int(row_pos)

        system_id = topology_id
        if system_idx is not None:
            value = row[system_idx]
            try:
                if value is not None and not (isinstance(value, float) and np.isnan(value)):
                    system_id = int(value)
            except (TypeError, ValueError):
                pass

        topology_to_system[topology_id] = system_id

        if loc_idx is not None:
            try:
                topology_to_loc[topology_id] = np.asarray(row[loc_idx], dtype=float)
            except (TypeError, ValueError):
                pass

        if rad_idx is not None:
            try:
                topology_to_rad[topology_id] = float(row[rad_idx])
            except (TypeError, ValueError):
                pass

    return {
        "group1_system_ids": set(int(i) for i in iface.group1_indices),
        "group2_system_ids": set(int(i) for i in iface.group2_indices),
        "topology_to_system": topology_to_system,
        "topology_to_loc": topology_to_loc,
        "topology_to_rad": topology_to_rad,
    }


def _surface_interface_orientation(surf, context, tol=1e-12):
    """Return orientation multiplier for the Group 1 -> Group 2 normal."""
    try:
        ball0, ball1 = [int(v) for v in surf["balls"]]
    except (KeyError, TypeError, ValueError):
        return 0, None, None

    topo_to_system = context["topology_to_system"]
    system0 = topo_to_system.get(ball0)
    system1 = topo_to_system.get(ball1)
    g1 = context["group1_system_ids"]
    g2 = context["group2_system_ids"]

    if system0 in g1 and system1 in g2:
        ball_g1, ball_g2 = ball0, ball1
    elif system1 in g1 and system0 in g2:
        ball_g1, ball_g2 = ball1, ball0
    else:
        return 0, None, None

    loc1 = context["topology_to_loc"].get(ball_g1)
    loc2 = context["topology_to_loc"].get(ball_g2)
    if loc1 is None or loc2 is None:
        return 0, ball_g1, ball_g2

    try:
        point = np.asarray(surf.get("com"), dtype=float)
        func = np.asarray(surf.get("func"), dtype=float)
    except (TypeError, ValueError):
        return 0, ball_g1, ball_g2

    if point.shape != (3,) or func.size < 9:
        return 0, ball_g1, ball_g2

    A, B, C, D, E, F, G, Hc, Ic = func[:9]
    x, y, z = point
    grad = np.array([
        2.0 * A * x + D * y + F * z + G,
        2.0 * B * y + D * x + E * z + Hc,
        2.0 * C * z + F * x + E * y + Ic,
    ], dtype=float)

    outward = loc2 - loc1
    grad_norm = np.linalg.norm(grad)
    outward_norm = np.linalg.norm(outward)
    if (
        not np.isfinite(grad_norm)
        or not np.isfinite(outward_norm)
        or grad_norm < tol
        or outward_norm < tol
    ):
        return 0, ball_g1, ball_g2

    dot = float(np.dot(grad, outward))
    if abs(dot) <= tol * grad_norm * outward_norm:
        return 0, ball_g1, ball_g2

    return (1 if dot > 0.0 else -1), ball_g1, ball_g2


def _triangle_areas(points, tris):
    try:
        pts = np.asarray(points, dtype=float)
        tri_idx = np.asarray(tris, dtype=int)
        if tri_idx.ndim != 2 or tri_idx.shape[1] != 3 or len(tri_idx) == 0:
            return None
        tri_pts = pts[tri_idx]
        return 0.5 * np.linalg.norm(
            np.cross(tri_pts[:, 1] - tri_pts[:, 0], tri_pts[:, 2] - tri_pts[:, 0]),
            axis=1,
        )
    except (TypeError, ValueError, IndexError):
        return None


def summarize_direct_interface_orientation(iface, direct_surfaces, flat_tol=1e-10):
    """
    Summarize direct interface geometry from both side perspectives.

    Group 1 -> Group 2 is the primary orientation. Group 2 -> Group 1 is its
    exact normal reversal, so C changes sign while A, Q, and X are invariant.
    """
    summary = {
        "surface_count": 0,
        "area": 0.0,
        "c12": 0.0,
        "q": 0.0,
        "x": 0.0,
        "surf_energy": 0.0,
        "g1_convex_sa": 0.0,
        "g1_concave_sa": 0.0,
        "flat_sa": 0.0,
        "failed_sa": 0.0,
        "curved_tested": 0,
        "radius_sign_agree": 0,
        "radius_sign_disagree": 0,
        "orientation_failed": 0,
        "equal_radius": 0,
        "equal_radius_flat": 0,
        "equal_radius_nonflat": 0,
        "equal_radius_nonzero_h": 0,
        "missing_radius": 0,
    }

    if iface.net is None or direct_surfaces is None or len(direct_surfaces) == 0:
        return summary

    context = _build_interface_orientation_context(iface)
    tol = 1e-10

    for _, surf in direct_surfaces.iterrows():
        sa = float(surf.get("sa", 0.0) or 0.0)
        raw_c = float(surf.get("int_mean_curv", 0.0) or 0.0)
        q = float(surf.get("int_mean_curv_sq", 0.0) or 0.0)
        x = float(surf.get("int_gauss_curv", 0.0) or 0.0)
        energy = float(surf.get("surf_energy", 0.0) or 0.0)

        summary["surface_count"] += 1
        summary["area"] += sa
        summary["q"] += q
        summary["x"] += x
        summary["surf_energy"] += energy

        orientation, ball_g1, ball_g2 = _surface_interface_orientation(surf, context)
        is_flat = bool(surf.get("flat", False))

        rin = context["topology_to_rad"].get(ball_g1, np.nan) if ball_g1 is not None else np.nan
        rout = context["topology_to_rad"].get(ball_g2, np.nan) if ball_g2 is not None else np.nan

        if is_flat:
            summary["flat_sa"] += sa
        elif orientation == 0:
            summary["failed_sa"] += sa
            summary["orientation_failed"] += 1
        else:
            oriented_c = orientation * raw_c
            summary["c12"] += oriented_c

            tri_h = surf.get("mean_tri_curvs", None)
            tri_sa = _triangle_areas(surf.get("points", None), surf.get("tris", None))
            classified = False
            if tri_h is not None and tri_sa is not None:
                try:
                    h = orientation * np.asarray(tri_h, dtype=float)
                    if len(h) == len(tri_sa):
                        valid = np.isfinite(h) & np.isfinite(tri_sa) & (tri_sa > 0.0)
                        h = h[valid]
                        areas = tri_sa[valid]
                        convex = h > flat_tol
                        concave = h < -flat_tol
                        near_flat = ~(convex | concave)
                        summary["g1_convex_sa"] += float(np.sum(areas[convex]))
                        summary["g1_concave_sa"] += float(np.sum(areas[concave]))
                        summary["flat_sa"] += float(np.sum(areas[near_flat]))
                        classified = True
                except (TypeError, ValueError):
                    pass

            if not classified:
                avg_h = oriented_c / sa if sa > 0.0 else 0.0
                if avg_h > flat_tol:
                    summary["g1_convex_sa"] += sa
                elif avg_h < -flat_tol:
                    summary["g1_concave_sa"] += sa
                else:
                    summary["flat_sa"] += sa

        if not np.isfinite(rin) or not np.isfinite(rout):
            summary["missing_radius"] += 1
            continue

        if abs(rin - rout) <= tol:
            summary["equal_radius"] += 1
            if is_flat:
                summary["equal_radius_flat"] += 1
            else:
                summary["equal_radius_nonflat"] += 1
            if abs(raw_c) > tol:
                summary["equal_radius_nonzero_h"] += 1
        else:
            summary["curved_tested"] += 1
            if orientation == 0:
                summary["radius_sign_disagree"] += 1
                continue
            oriented_avg_h = orientation * raw_c / sa if sa > 0.0 else 0.0
            expected_sign = 1 if rin < rout else -1
            if abs(oriented_avg_h) <= tol or np.sign(oriented_avg_h) == expected_sign:
                summary["radius_sign_agree"] += 1
            else:
                summary["radius_sign_disagree"] += 1

    summary["c21"] = -summary["c12"]
    summary["g2_convex_sa"] = summary["g1_concave_sa"]
    summary["g2_concave_sa"] = summary["g1_convex_sa"]
    summary["validation_pass"] = (
        summary["radius_sign_disagree"] == 0
        and summary["orientation_failed"] == 0
        and summary["equal_radius_nonflat"] == 0
        and summary["equal_radius_nonzero_h"] == 0
        and summary["missing_radius"] == 0
    )
    return summary


def write_interface_orientation_statistics(info, iface, direct_surfaces):
    summary = summarize_direct_interface_orientation(iface, direct_surfaces)
    group1_name = getattr(iface.group1, "name", "group1")
    group2_name = getattr(iface, "group2_name", None) or getattr(iface.group2, "name", "group2")
    area = summary["area"]

    info.write("Direct interface oriented curvature:\n")
    info.write(f"  Primary orientation: {group1_name} -> {group2_name}\n")
    info.write(f"  Direct surfaces: {summary['surface_count']}\n")
    info.write(f"  Direct surface area: {_format_metric(area)} Å²\n")
    info.write(f"  Integrated H dA ({group1_name} -> {group2_name}): {_format_metric(summary['c12'])} Å\n")
    info.write(f"  Integrated H dA ({group2_name} -> {group1_name}): {_format_metric(summary['c21'])} Å\n")
    info.write(f"  Integrated H² dA (orientation invariant): {_format_metric(summary['q'])}\n")
    info.write(f"  Integrated K dA (orientation invariant): {_format_metric(summary['x'])}\n")
    info.write(f"  Representative surface energy (orientation invariant): {_format_metric(summary['surf_energy'])} kBT\n\n")

    info.write(f"  {group1_name} perspective:\n")
    info.write(f"    Convex area: {_format_metric(summary['g1_convex_sa'])} Å²\n")
    info.write(f"    Concave area: {_format_metric(summary['g1_concave_sa'])} Å²\n")
    info.write(f"    Flat/near-flat area: {_format_metric(summary['flat_sa'])} Å²\n")
    info.write(f"    Failed curved area: {_format_metric(summary['failed_sa'])} Å²\n")

    info.write(f"  {group2_name} perspective:\n")
    info.write(f"    Convex area: {_format_metric(summary['g2_convex_sa'])} Å²\n")
    info.write(f"    Concave area: {_format_metric(summary['g2_concave_sa'])} Å²\n")
    info.write(f"    Flat/near-flat area: {_format_metric(summary['flat_sa'])} Å²\n")
    info.write(f"    Failed curved area: {_format_metric(summary['failed_sa'])} Å²\n\n")

    tested = summary["curved_tested"]
    pct = 100.0 * summary["radius_sign_agree"] / tested if tested else 100.0
    info.write("  Orientation validation:\n")
    info.write(f"    Curved surfaces tested: {tested}\n")
    info.write(f"    Radius/sign agreements: {summary['radius_sign_agree']}\n")
    info.write(f"    Radius/sign disagreements: {summary['radius_sign_disagree']}\n")
    info.write(f"    Curved orientation failures: {summary['orientation_failed']}\n")
    info.write(f"    Radius/sign agreement: {pct:.3f} %\n")
    info.write(f"    Equal-radius surfaces: {summary['equal_radius']}\n")
    info.write(f"    Equal-radius marked flat: {summary['equal_radius_flat']}\n")
    info.write(f"    Equal-radius non-flat: {summary['equal_radius_nonflat']}\n")
    info.write(f"    Equal-radius nonzero H: {summary['equal_radius_nonzero_h']}\n")
    info.write(f"    Missing radius metadata: {summary['missing_radius']}\n")
    info.write(f"    Result: {'PASS' if summary['validation_pass'] else 'FAIL'}\n\n")


def get_interface_surface_sets(iface):
    """
    Separate interface-network surfaces into:

    direct_surfaces
        One defining ball belongs to Group 1 and the other belongs to Group 2.

    group1_surfaces
        Both defining balls belong to Group 1.

    group2_surfaces
        Both defining balls belong to Group 2.

    support_surfaces
        At least one defining ball is outside the two explicit interface groups.

    The support category may contain geometry required to construct or close
    the dedicated interface network without representing direct inter-group
    contact.
    """
    if iface.net is None or iface.net.surfs is None:
        return None, None, None, None

    surfaces = iface.net.surfs

    ball1_column = _get_column(
        surfaces,
        "Ball 1",
        "ball 1",
        "ball1",
        "b1",
    )
    ball2_column = _get_column(
        surfaces,
        "Ball 2",
        "ball 2",
        "ball2",
        "b2",
    )

    # Some internal surface tables store both balls in one iterable column.
    balls_column = _get_column(
        surfaces,
        "balls",
        "surface balls",
    )

    group1_indices = set(int(index) for index in iface.group1_indices)
    group2_indices = set(int(index) for index in iface.group2_indices)

    if ball1_column is not None and ball2_column is not None:
        ball_pairs = [
            (int(ball1), int(ball2))
            for ball1, ball2 in zip(
                surfaces[ball1_column],
                surfaces[ball2_column],
            )
        ]

    elif balls_column is not None:
        ball_pairs = [
            (int(balls[0]), int(balls[1]))
            for balls in surfaces[balls_column]
        ]

    else:
        # The topology can still be summarized, but surface membership cannot
        # be classified without defining-ball information.
        empty = surfaces.iloc[0:0]
        return empty, empty, empty, surfaces

    direct_mask = []
    group1_mask = []
    group2_mask = []
    support_mask = []

    for ball1, ball2 in ball_pairs:
        ball1_in_group1 = ball1 in group1_indices
        ball2_in_group1 = ball2 in group1_indices
        ball1_in_group2 = ball1 in group2_indices
        ball2_in_group2 = ball2 in group2_indices

        is_direct = (
            (ball1_in_group1 and ball2_in_group2)
            or
            (ball1_in_group2 and ball2_in_group1)
        )

        is_group1 = ball1_in_group1 and ball2_in_group1
        is_group2 = ball1_in_group2 and ball2_in_group2

        is_support = not (
            is_direct
            or is_group1
            or is_group2
        )

        direct_mask.append(is_direct)
        group1_mask.append(is_group1)
        group2_mask.append(is_group2)
        support_mask.append(is_support)

    return (
        surfaces.loc[direct_mask],
        surfaces.loc[group1_mask],
        surfaces.loc[group2_mask],
        surfaces.loc[support_mask],
    )


def write_surface_statistics(info, title, surfaces):
    """
    Write area, curvature, contact, and overlap statistics for a collection
    of surfaces.
    """
    if surfaces is None:
        info.write(f"{title}:\n")
        info.write("  Not available\n\n")
        return

    surface_areas = _numeric_series(
        surfaces,
        "sa",
    )

    mean_curvatures = _numeric_series(
        surfaces,
        "mean_curv",
    )

    average_mean_curvatures = _numeric_series(
        surfaces,
        "avg_mean_curv",
    )

    gaussian_curvatures = _numeric_series(
        surfaces,
        "gauss_curv",
    )

    average_gaussian_curvatures = _numeric_series(
        surfaces,
        "avg_gauss_curv",
    )

    # Integrated curvature geometry
    integrated_mean_curvatures = _numeric_series(
        surfaces,
        "int_mean_curv",
        "Integrated Mean Curvature",
    )

    integrated_mean_curvature_squared = _numeric_series(
        surfaces,
        "int_mean_curv_sq",
        "Integrated Mean Curvature Squared",
    )

    integrated_gaussian_curvatures = _numeric_series(
        surfaces,
        "int_gauss_curv",
        "Integrated Gaussian Curvature",
    )

    contact_areas = _numeric_series(
        surfaces,
        "contact_area",
    )

    contact_areas = _numeric_series(
        surfaces,
        "contact_area",
    )

    overlaps = _numeric_series(
        surfaces,
        "overlap",
    )

    info.write(f"{title}:\n")
    info.write(f"  Number of surfaces: {len(surfaces)}\n")
    info.write(
        f"  Total surface area: "
        f"{_format_metric(_safe_sum(surface_areas))} Å²\n"
    )
    info.write(
        f"  Mean surface area: "
        f"{_format_metric(_safe_mean(surface_areas))} Å²\n"
    )
    info.write(
        f"  Minimum surface area: "
        f"{_format_metric(_safe_min(surface_areas))} Å²\n"
    )
    info.write(
        f"  Maximum surface area: "
        f"{_format_metric(_safe_max(surface_areas))} Å²\n"
    )

    info.write("\n  Mean curvature:\n")
    info.write(
        f"    Mean of surface-average mean curvature: "
        f"{_format_metric(_safe_mean(average_mean_curvatures))} Å⁻¹\n"
    )
    info.write(f"    Area-weighted surface-average mean curvature: {_format_metric(_weighted_mean(average_mean_curvatures, surface_areas))} Å⁻¹\n")
    info.write(
        f"    Minimum surface-average mean curvature: "
        f"{_format_metric(_safe_min(average_mean_curvatures))} Å⁻¹\n"
    )
    info.write(
        f"    Maximum surface-average mean curvature: "
        f"{_format_metric(_safe_max(average_mean_curvatures))} Å⁻¹\n"
    )
    info.write(
        f"    Maximum local mean curvature: "
        f"{_format_metric(_safe_max(mean_curvatures))} Å⁻¹\n"
    )

    info.write("\n  Gaussian curvature:\n")
    info.write(f"    Mean of surface-average Gaussian curvature: "
        f"{_format_metric(_safe_mean(average_gaussian_curvatures))} Å⁻²\n")
    info.write(
        f"    Area-weighted surface-average Gaussian curvature: "
        f"{_format_metric(_weighted_mean(average_gaussian_curvatures, surface_areas))} Å⁻²\n")
    info.write(
        f"    Minimum surface-average Gaussian curvature: "
        f"{_format_metric(_safe_min(average_gaussian_curvatures))} Å⁻²\n"
    )
    info.write(
        f"    Maximum surface-average Gaussian curvature: "
        f"{_format_metric(_safe_max(average_gaussian_curvatures))} Å⁻²\n")
    info.write(
        f"    Maximum local Gaussian curvature: "
        f"{_format_metric(_safe_max(gaussian_curvatures))} Å⁻²\n"
    )

    # ------------------------------------------------------------------
    # Integrated curvature geometry
    # ------------------------------------------------------------------

    info.write("\n  Integrated curvature geometry:\n")

    total_int_mean_curv = _safe_sum(
        integrated_mean_curvatures
    )

    total_int_mean_curv_sq = _safe_sum(
        integrated_mean_curvature_squared
    )

    total_int_gauss_curv = _safe_sum(
        integrated_gaussian_curvatures
    )

    info.write(
        f"    Integrated mean curvature (∫H dA): "
        f"{_format_metric(total_int_mean_curv)} Å\n"
    )

    info.write(
        f"    Integrated squared mean curvature (∫H² dA): "
        f"{_format_metric(total_int_mean_curv_sq)}\n"
    )

    info.write(
        f"    Integrated Gaussian curvature (∫K dA): "
        f"{_format_metric(total_int_gauss_curv)}\n"
    )

    # Derive area-normalized curvature from the integrated quantities
    total_surface_area = _safe_sum(surface_areas)

    if total_surface_area is not None and total_surface_area > 0.0:

        if total_int_mean_curv is not None:
            info.write(
                f"    Area-normalized mean curvature: "
                f"{_format_metric(total_int_mean_curv / total_surface_area)} Å⁻¹\n"
            )

        if total_int_gauss_curv is not None:
            info.write(
                f"    Area-normalized Gaussian curvature: "
                f"{_format_metric(total_int_gauss_curv / total_surface_area)} Å⁻²\n"
            )

    info.write("\n  Contact and overlap:\n")
    info.write(
        f"    Total contact area: "
        f"{_format_metric(_safe_sum(contact_areas))} Å²\n"
    )
    info.write(
        f"    Surfaces with positive contact area: "
        f"{0 if contact_areas is None else int((contact_areas > 0).sum())}\n"
    )
    info.write(
        f"    Total overlap measure: "
        f"{_format_metric(_safe_sum(overlaps))} Å\n"
    )
    info.write(
        f"    Surfaces with positive overlap: "
        f"{0 if overlaps is None else int((overlaps > 0).sum())}\n"
    )

    info.write("\n")



def write_interface_water_topology_statistics(info, iface):
    """Write a compact production summary of interface-water classification."""
    topology = getattr(iface, "water_topology", None)

    info.write("Interface waters:\n")
    if not topology:
        info.write("  No interface-water topology analysis available.\n\n")
        return

    waters = topology.get("waters", {})
    cycle_sets = topology.get("cycle_sets", [])
    buried_sets = topology.get("buried_cycle_sets", [])

    counts = {"buried": 0, "semi_buried": 0, "peripheral": 0}
    for water in waters.values():
        burial_class = water.get("burial_class")
        if burial_class in counts:
            counts[burial_class] += 1

    non_buried_count = counts["peripheral"] + counts["semi_buried"]

    info.write(f"  Total interface waters: {len(waters)}\n")
    info.write(f"  Non-buried interface waters: {non_buried_count}\n")
    info.write(f"  Buried interface waters: {counts['buried']}\n")
    info.write(f"  Buried water groups: {len(buried_sets)}\n")

    qc = topology.get("discovery_qc", {})
    mapping_failures = (
        int(qc.get("unmapped_ball_references", 0))
        + int(qc.get("edges_with_mapping_failure", 0))
    )
    info.write(f"  Topology mapping failures: {mapping_failures}\n")

    if buried_sets:
        info.write("\n  Buried water groups:\n")
        for group_index, cycle_set in enumerate(buried_sets, start=1):
            labels = cycle_set.get("water_labels", [])
            solved_name = cycle_set.get("solved_group_name", f"buried_{group_index}")
            info.write(
                f"    {group_index}: {', '.join(labels)} "
                f"[{len(cycle_set.get('edge_indices', []))} cycle edges] "
                f"-> waters/{solved_name}/\n"
            )

    info.write("\n  Detailed buried-water geometry: waters/\n\n")


def _water_class_atom_indices(iface):
    """Return complete atom selections for non-buried and buried waters.

    ``peripheral`` and ``semi_buried`` remain separate internal topology labels
    because they are useful diagnostics, but neither defines a closed buried
    water volume.  Production exports therefore combine both as ``non_buried``.
    """
    topology = getattr(iface, "water_topology", None) or {}
    class_atoms = {"non_buried": set(), "buried": set()}

    for water in topology.get("waters", {}).values():
        burial_class = water.get("burial_class")
        if burial_class == "buried":
            export_class = "buried"
        elif burial_class in {"peripheral", "semi_buried"}:
            export_class = "non_buried"
        else:
            continue

        residue = water.get("residue")
        if residue is None:
            continue
        for atom_index in getattr(residue, "atoms", []) or []:
            try:
                class_atoms[export_class].add(int(atom_index))
            except (TypeError, ValueError):
                continue

    return {key: sorted(values) for key, values in class_atoms.items()}


def export_water_class_pdbs(iface):
    """Export non-buried and buried interface waters as complete-residue PDBs."""
    waters_dir = os.path.join(iface.dir, "waters")
    os.makedirs(waters_dir, exist_ok=True)

    selections = _water_class_atom_indices(iface)
    filenames = {
        "non_buried": "non_buried_waters",
        "buried": "buried_waters",
    }

    cwd = os.getcwd()
    try:
        os.chdir(waters_dir)
        for burial_class, atoms in selections.items():
            if not atoms:
                continue
            write_pdb(
                atoms=atoms,
                file_name=filenames[burial_class],
                sys=iface.sys,
            )
    finally:
        os.chdir(cwd)

    return selections


def _write_buried_water_summary(iface):
    """Write aggregate geometry for solved buried-water Groups."""
    waters_dir = os.path.join(iface.dir, "waters")
    os.makedirs(waters_dir, exist_ok=True)
    path = os.path.join(waters_dir, "info.txt")
    groups = list(getattr(iface, "buried_water_groups", []) or [])

    with open(path, "w", encoding="utf-8") as info:
        info.write(f"Buried interface-water analysis - {iface.name}\n\n")
        info.write(f"Buried water groups: {len(groups)}\n")
        info.write(f"Total waters in buried groups: {sum(len(getattr(g, 'buried_water_metadata', {}).get('residues', [])) for g in groups)}\n\n")

        if not groups:
            info.write("No closed buried-water groups were detected.\n")
            return path

        info.write(
            "Group | Waters | Shell Color | Volume (A^3) | Surface Area (A^2) | "
            "Oriented int H dA (A) | int H^2 dA | int K dA | Representative Energy (kBT)\n"
        )
        info.write("-" * 140 + "\n")

        total_volume = 0.0
        total_sa = 0.0
        total_c = 0.0
        total_q = 0.0
        total_x = 0.0

        for group in groups:
            metadata = getattr(group, "buried_water_metadata", {})
            labels = ", ".join(metadata.get("water_labels", []))
            shell_color = metadata.get("shell_color_map", "default")
            volume = float(getattr(group, "vol", 0.0) or 0.0)
            sa = float(getattr(group, "sa", 0.0) or 0.0)
            c = float(getattr(group, "oriented_int_mean_curv", 0.0) or 0.0)
            q = float(getattr(group, "int_mean_curv_sq", 0.0) or 0.0)
            x = float(getattr(group, "int_gauss_curv", 0.0) or 0.0)
            energy = 2.0 * q

            total_volume += volume
            total_sa += sa
            total_c += c
            total_q += q
            total_x += x

            info.write(
                f"{group.name} | {labels} | {shell_color} | {volume:.6f} | {sa:.6f} | "
                f"{c:.6f} | {q:.6f} | {x:.6f} | {energy:.6f}\n"
            )

        info.write("\nAggregate sums across buried groups:\n")
        info.write(f"  Volume: {total_volume:.6f} A^3\n")
        info.write(f"  Surface area: {total_sa:.6f} A^2\n")
        info.write(f"  Oriented int H dA: {total_c:.6f} A\n")
        info.write(f"  int H^2 dA: {total_q:.6f}\n")
        info.write(f"  int K dA: {total_x:.6f}\n")
        info.write(f"  Representative surface energy: {2.0 * total_q:.6f} kBT\n")

    return path


def _update_water_export_progress(iface, process, progress=None):
    """Relabel the active interface-export progress without starting a nested percentage scale."""
    sys = getattr(iface, "sys", None)
    if sys is None:
        return

    updater = getattr(sys, "update_progress", None)
    if updater is None:
        return

    # Preserve the enclosing export percentage.  The previous implementation
    # ran a second 0->100 counter inside one interface export step, which caused
    # percentage jumps and stale console text.
    current_progress = getattr(sys, "run_progress", None)
    if current_progress is None:
        current_progress = 0.0 if progress is None else progress

    try:
        updater(
            process=process,
            progress=current_progress,
            network=iface.name,
        )
    except TypeError:
        try:
            updater(process, current_progress, name=iface.name)
        except TypeError:
            updater(process=process, progress=current_progress)


def export_buried_water_groups(iface):
    """Export a compact companion set for every solved buried-water group.

    Buried-water geometry is intentionally a subordinate interface product, so
    only the files needed to identify and visualize the closed water cells are
    emitted: class-level PDBs, aggregate info, per-group info/PDB, and the
    external shell surfaces/edges/vertices.  Full Group network tables and logs
    remain available in memory but are not written by default.
    """
    if getattr(iface, "_buried_water_exports_complete", False):
        return

    groups = list(getattr(iface, "buried_water_groups", []) or [])
    total_steps = max(2 + len(groups), 1)
    step = 0

    _update_water_export_progress(
        iface,
        "Exporting interface waters: PDBs",
        100.0 * step / total_steps,
    )
    export_water_class_pdbs(iface)
    step += 1

    _update_water_export_progress(
        iface,
        "Exporting interface waters: summary",
        100.0 * step / total_steps,
    )
    _write_buried_water_summary(iface)
    step += 1

    for group_index, group in enumerate(groups, start=1):
        _update_water_export_progress(
            iface,
            f"Exporting interface water {group_index}/{len(groups)}: {group.name}",
            100.0 * step / total_steps,
        )

        os.makedirs(group.dir, exist_ok=True)

        # Compact buried-water product: detailed Group info plus the external
        # cell shell.  Do not emit full surfs/edges/verts or a Group log here;
        # this is a small subordinate object of the parent Interface.
        group.exports(
            info=True,
            shell_surfs=True,
            shell_edges=True,
            shell_verts=True,
        )

        # One plainly named PDB per buried group; avoid the duplicate
        # group_atoms.pdb generated by the normal Group ``atoms=True`` export.
        cwd = os.getcwd()
        try:
            os.chdir(group.dir)
            write_pdb(
                atoms=list(group.ball_ndxs),
                file_name=group.name,
                sys=group.sys,
            )
        finally:
            os.chdir(cwd)

        step += 1

    _update_water_export_progress(
        iface,
        "Exporting interface waters",
        100.0,
    )
    iface._buried_water_exports_complete = True


def export_info(iface, directory=None):
    """
    Export a detailed geometric and topological summary of an Interface.
    """
    if directory is None:
        directory = iface.dir

    os.makedirs(directory, exist_ok=True)

    net = iface.net

    num_verts = (
        0
        if net is None or net.verts is None
        else len(net.verts)
    )
    num_edges = (
        0
        if net is None or net.edges is None
        else len(net.edges)
    )
    num_surfs = (
        0
        if net is None or net.surfs is None
        else len(net.surfs)
    )

    group1_atoms = sorted(set(iface.group1_indices))
    group2_atoms = sorted(set(iface.group2_indices))
    interface_atoms = sorted(
        set(group1_atoms) | set(group2_atoms)
    )

    (
        direct_surfaces,
        group1_surfaces,
        group2_surfaces,
        support_surfaces,
    ) = get_interface_surface_sets(iface)

    # One-sided or partially populated interfaces may not have explicit
    # surface-classification sets. Normalize them to empty DataFrames so
    # topology/info export remains valid.
    empty_surfaces = None if net is None or net.surfs is None else net.surfs.iloc[0:0].copy()
    if direct_surfaces is None:
        direct_surfaces = empty_surfaces
    if group1_surfaces is None:
        group1_surfaces = empty_surfaces
    if group2_surfaces is None:
        group2_surfaces = empty_surfaces
    if support_surfaces is None:
        support_surfaces = empty_surfaces

    group1_name = getattr(getattr(iface, "group1", None), "name", "group1")
    group2_name = getattr(iface, "group2_name", None) or getattr(getattr(iface, "group2", None), "name", "surrounding")
    interface_id = getattr(iface, "interface_id", getattr(iface, "name", "interface"))

    file_path = os.path.join(directory, "info.txt")

    with open(file_path, "w", encoding="utf-8") as info:
        info.write(f"{iface.name} - {iface.sys.name}\n\n")

        info.write("Interface definition:\n")
        info.write(f"  Interface ID: {interface_id}\n")
        info.write(f"  Group 1: {group1_name}\n")
        info.write(f"  Group 2: {group2_name}\n\n")

        info.write("Interface atom membership:\n")
        info.write(f"  Group 1 atoms: {len(group1_atoms)}\n")
        info.write(f"  Group 2 atoms: {len(group2_atoms)}\n")
        info.write(
            f"  Unique interface-network atoms: "
            f"{len(interface_atoms)}\n"
        )

        shared_atoms = set(group1_atoms) & set(group2_atoms)

        info.write(
            f"  Atoms shared by both definitions: "
            f"{len(shared_atoms)}\n\n"
        )

        info.write("Interface network topology:\n")
        info.write(f"  Vertices: {num_verts}\n")
        info.write(f"  Edges: {num_edges}\n")
        info.write(f"  Surfaces: {num_surfs}\n")

        if num_verts > 0:
            info.write(
                f"  Edges per vertex: "
                f"{num_edges / num_verts:.6f}\n"
            )
            info.write(
                f"  Surfaces per vertex: "
                f"{num_surfs / num_verts:.6f}\n"
            )

        info.write("\n")

        info.write("Surface classification:\n")
        info.write(
            f"  Direct Group 1–Group 2 surfaces: "
            f"{0 if direct_surfaces is None else len(direct_surfaces)}\n"
        )
        info.write(
            f"  Group 1 internal surfaces: "
            f"{0 if group1_surfaces is None else len(group1_surfaces)}\n"
        )
        info.write(
            f"  Group 2 internal surfaces: "
            f"{0 if group2_surfaces is None else len(group2_surfaces)}\n"
        )
        info.write(
            f"  Supporting/unclassified surfaces: "
            f"{0 if support_surfaces is None else len(support_surfaces)}\n\n"
        )

        write_surface_statistics(
            info,
            "Full interface-network surface statistics",
            net.surfs if net is not None else None,
        )

        write_surface_statistics(
            info,
            "Direct inter-group surface statistics",
            direct_surfaces,
        )

        write_interface_orientation_statistics(
            info=info,
            iface=iface,
            direct_surfaces=direct_surfaces,
        )

        write_interface_water_topology_statistics(
            info=info,
            iface=iface,
        )

        # Zero-length classes add nearly a page of "Not available" values and
        # carry no information. Only report classes that are actually present.
        if group1_surfaces is not None and len(group1_surfaces) > 0:
            write_surface_statistics(
                info,
                f"{group1_name} internal surface statistics",
                group1_surfaces,
            )

        if group2_surfaces is not None and len(group2_surfaces) > 0:
            write_surface_statistics(
                info,
                f"{group2_name} internal surface statistics",
                group2_surfaces,
            )

        if support_surfaces is not None and len(support_surfaces) > 0:
            write_surface_statistics(
                info,
                "Supporting surface statistics",
                support_surfaces,
            )

        if net is not None and net.edges is not None:
            edge_lengths = _numeric_series(net.edges, "length")

            info.write("Edge geometry:\n")
            info.write(
                f"  Total edge length: "
                f"{_format_metric(_safe_sum(edge_lengths))} Å\n"
            )
            info.write(
                f"  Average edge length: "
                f"{_format_metric(_safe_mean(edge_lengths))} Å\n"
            )
            info.write(
                f"  Minimum edge length: "
                f"{_format_metric(_safe_min(edge_lengths))} Å\n"
            )
            info.write(
                f"  Maximum edge length: "
                f"{_format_metric(_safe_max(edge_lengths))} Å\n\n"
            )

        if net is not None and net.verts is not None:
            vertex_radii = _numeric_series(net.verts, "rad")

            info.write("Vertex geometry:\n")
            info.write(
                f"  Average vertex radius: "
                f"{_format_metric(_safe_mean(vertex_radii))} Å\n"
            )
            info.write(
                f"  Minimum vertex radius: "
                f"{_format_metric(_safe_min(vertex_radii))} Å\n"
            )
            info.write(
                f"  Maximum vertex radius: "
                f"{_format_metric(_safe_max(vertex_radii))} Å\n"
            )

            negative_radii = (
                0
                if vertex_radii is None
                else int((vertex_radii < 0).sum())
            )

            info.write(
                f"  Vertices with negative radius: "
                f"{negative_radii}\n\n"
            )


def interface_exports(iface, all_=False, atoms=False, surfs=False, edges=False, verts=False, logs=False, info=False,
                      group_info=False, round_to=3):
    """
    Export data belonging to an Interface and its dedicated Network.
    """
    if iface.net is None:
        print(f'Interface "{iface.name}" has no network to export.')
        return

    if iface.dir is None:
        iface.dir = os.path.join(
            iface.sys.files["dir"],
            iface.name,
        )

    os.makedirs(iface.dir, exist_ok=True)

    # Buried-water networks are subordinate interface products. Export them
    # once with the interface information pass rather than during an arbitrary
    # geometry export (surfaces/edges/verts). The helper retains its own guard
    # for explicit repeated info exports.
    if info or all_:
        export_buried_water_groups(iface)

    if group_info or all_:
        exported_group_ids = set()

        for group in (iface.group1, iface.group2):
            if group is None:
                continue

            group_id = str(getattr(group, "group_id", None) or group.name)

            if group_id in exported_group_ids:
                continue

            exported_group_ids.add(group_id)

            group_directory = os.path.join(
                iface.sys.files["dir"],
                group.name,
            )

            export_interface_group_info(
                group=group,
                directory=group_directory,
            )

    if atoms or all_:
        interface_atoms = get_interface_atoms(iface)

        if iface.sys.files["base_file"][-3:].lower() != "txt":
            write_pdb(
                atoms=interface_atoms,
                file_name="interface_atoms",
                directory=iface.dir,
                sys=iface.sys,
            )

    if info or all_:
        export_info(iface)

    if logs or all_:
        write_interface_logs(
            iface,
            round_to=round_to,
        )

    if verts or all_:
        if iface.net.verts is not None and len(iface.net.verts) > 0:
            write_off_verts(
                iface.net,
                list(range(len(iface.net.verts))),
                directory=iface.dir,
                file_name="verts",
                color=iface.net.settings["vert_col"],
            )

    if edges or all_:
        if iface.net.edges is not None and len(iface.net.edges) > 0:
            write_edges(
                iface.net,
                list(range(len(iface.net.edges))),
                directory=iface.dir,
                file_name="edges",
                color=iface.net.settings["edge_col"],
            )

    if surfs or all_:
        if iface.net.surfs is not None and len(iface.net.surfs) > 0:
            write_surfs(
                iface.net,
                list(range(len(iface.net.surfs))),
                directory=iface.dir,
                file_name="surfs",
            )

def export_interface_group_info(group, directory):
    """
    Export information for a Group that participated in an interface-only
    calculation.

    The group does not need a complete Network. This exports the group's
    selected atoms and interface relationship metadata only.
    """
    os.makedirs(directory, exist_ok=True)

    group_indices = sorted(
        set(int(index) for index in (getattr(group, "ball_ndxs", []) or []))
    )

    info_path = os.path.join(directory, "info.txt")

    with open(info_path, "w", encoding="utf-8") as info:
        info.write(f"{group.name} - {group.sys.name}\n\n")

        info.write("Group definition:\n")
        info.write(f"  Group ID: {getattr(group, 'group_id', group.name)}\n")
        info.write(f"  Selected atoms: {len(group_indices)}\n")
        info.write(f"  Full group network built: {group.net is not None}\n\n")

        info.write("Selections:\n")
        info.write(f"  Atom selections: {len(group.atms or [])}\n")
        info.write(f"  Residue selections: {len(group.rsds or [])}\n")
        info.write(f"  Chain selections: {len(group.chns or [])}\n")
        info.write(f"  Molecule selections: {len(group.mols or [])}\n\n")

        metadata = getattr(group, "interface_metadata", {}) or {}

        info.write("Interfaces:\n")

        if not metadata:
            info.write("  None\n")
        else:
            for interface_id, record in metadata.items():
                info.write(f"  Interface ID: {interface_id}\n")
                info.write(
                    f"    Name: "
                    f"{record.get('interface_name', interface_id)}\n"
                )
                info.write(
                    f"    Other group: "
                    f"{record.get('other_group_name', 'surrounding')}\n"
                )
                info.write(
                    f"    Side: {record.get('side')}\n"
                )
                info.write(
                    f"    Status: {record.get('status', 'defined')}\n"
                )
                info.write(
                    f"    Built: {record.get('built', False)}\n"
                )

    base_file = group.sys.files.get("base_file")

    if (
            base_file is not None
            and not str(base_file).lower().endswith(".txt")
    ):
        write_pdb(
            atoms=group_indices,
            file_name="group_atoms",
            directory=directory,
            sys=group.sys,
        )