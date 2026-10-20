import csv
import os
import shutil
import time
from copy import copy
from dataclasses import replace
import numpy as np
from vorpy.src.output import write_off_verts
from vorpy.src.output import write_edges
from vorpy.src.output import write_surfs
from vorpy.src.output import write_pdb
from vorpy.src.output import write_interface_logs
from vorpy.src.output.draw import DEFAULT_EDGE_RADIUS, DEFAULT_VERTEX_RADIUS
from vorpy.src.output.export_scope import export_scope


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


def _safe_sum(series, empty_value=None):
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


def _write_empty_off(directory, file_name):
    """Write a valid empty OFF file for an explicitly requested layer."""
    with open(os.path.join(directory, f"{file_name}.off"), "w", encoding="utf-8") as stream:
        stream.write("OFF\n0 0 0\n")


def _interface_export_profile_enabled():
    value = os.environ.get(
        "VORPY_INTERFACE_EXPORT_PROFILE",
        os.environ.get("VORPY_EXPORT_TIMING", "0"),
    )
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _profile_stage(profile, name, callback):
    stage_callback = profile.get('_cli_stage_callback') if profile is not None else None
    stage_label = {
        'scope_selection': 'Export preparation',
        'mesh_network': 'Export preparation',
        'contact_atom_pdb': 'Contact/summary exports',
        'group_info': 'Contact/summary exports',
        'interface_info': 'Contact/summary exports',
        'legacy_interface_logs': 'Contact/summary exports',
        'physical_interface_tables': 'Contact/summary exports',
        'vertex_mesh': 'Exporting physical geometry',
        'edge_mesh': 'Exporting physical geometry',
        'surface_mesh': 'Exporting physical geometry',
        'dual_visualization': 'Exporting visualization',
        'buried_water_exports': 'Contact/summary exports',
    }.get(name, 'Contact/summary exports')
    if callable(stage_callback) and profile.get('_cli_stage_current') != stage_label:
        current = profile.get('_cli_stage_current')
        if current is not None:
            stage_callback(current, 'finish')
        stage_callback(stage_label, 'start')
        profile['_cli_stage_current'] = stage_label
    if profile is None:
        return callback()
    started = time.perf_counter()
    try:
        return callback()
    finally:
        profile["timings"][name] = profile["timings"].get(name, 0.0) + (
            time.perf_counter() - started
        )


def _interface_export_counts(scope):
    """Count selected objects without copying source geometry."""
    network = scope.network
    def table_length(name):
        table = getattr(network, name, None)
        return 0 if table is None else len(table)

    counts = {
        "source_generators": table_length("balls"),
        "source_surfaces": table_length("surfs"),
        "source_edges": table_length("edges"),
        "source_vertices": table_length("verts"),
        "selected_surfaces": len(scope.surfaces),
        "selected_edges": len(scope.edges),
        "selected_vertices": len(scope.vertices),
        "selected_surface_points": 0,
        "selected_surface_triangles": 0,
        "selected_edge_points": 0,
        "selected_edge_segments": 0,
    }
    for position in scope.surfaces:
        row = network.surfs.iloc[int(position)]
        points = row.get("points", ())
        tris = row.get("tris", ())
        counts["selected_surface_points"] += 0 if points is None else len(points)
        counts["selected_surface_triangles"] += 0 if tris is None else len(tris)
    for position in scope.edges:
        points = network.edges.iloc[int(position)].get("points", ())
        if points is None:
            points = ()
        counts["selected_edge_points"] += len(points)
        counts["selected_edge_segments"] += max(len(points) - 1, 0)
    return counts


def _interface_export_profile_summary(profile):
    if profile is None:
        return
    print("INTERFACE EXPORT PROFILE")
    for name, elapsed in sorted(
        profile["timings"].items(), key=lambda item: item[1], reverse=True
    ):
        print(f"  {name}: {elapsed:.6f} s")
    print("  object counts:")
    for name, count in profile["counts"].items():
        print(f"    {name}: {count}")



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


def summarize_direct_interface_orientation(
    iface, direct_surfaces, flat_tol=1e-10, *, include_surface_energy=False
):
    """
    Summarize direct interface geometry from both side perspectives.

    Group 1 -> Group 2 is the primary orientation. Group 2 -> Group 1 is its
    exact normal reversal, so C changes sign while A, Q, and X are invariant.
    """
    summary = {
        "surface_count": 0,
        "area": 0.0,
        "c12": 0.0,
        "c21": 0.0,
        "q": 0.0,
        "x": 0.0,
        "surf_energy": 0.0 if include_surface_energy else None,
        "g1_convex_sa": 0.0,
        "g1_concave_sa": 0.0,
        "g2_convex_sa": 0.0,
        "g2_concave_sa": 0.0,
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
        "validation_pass": False,
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

        summary["surface_count"] += 1
        summary["area"] += sa
        summary["q"] += q
        summary["x"] += x
        if include_surface_energy:
            summary["surf_energy"] += float(surf.get("surf_energy", 0.0) or 0.0)

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


def write_interface_orientation_statistics(
    info, iface, direct_surfaces, *, include_surface_energy=False
):
    summary = summarize_direct_interface_orientation(
        iface, direct_surfaces, include_surface_energy=include_surface_energy
    )
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
    if include_surface_energy:
        info.write(
            f"  Representative surface energy (orientation invariant): "
            f"{_format_metric(summary['surf_energy'])} kBT\n\n"
        )
    else:
        info.write(
            "  Representative surface energy: NOT_CALCULATED "
            "(experimental; opt-in only)\n\n"
        )

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


def get_interface_surface_sets(iface, network=None):
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
    net = iface.net if network is None else network
    if net is None or net.surfs is None:
        return None, None, None, None

    surfaces = net.surfs

    topology_to_system = {}
    balls = getattr(net, "balls", None)
    if balls is not None:
        for position, (row_index, ball) in enumerate(balls.iterrows()):
            try:
                topology_id = int(ball.get("num", position))
            except (TypeError, ValueError):
                topology_id = position
            system_value = ball.get("system_num", None)
            try:
                if system_value is None or (isinstance(system_value, float) and np.isnan(system_value)):
                    system_id = int(row_index)
                else:
                    system_id = int(system_value)
            except (TypeError, ValueError):
                system_id = topology_id
            for identifier in (position, row_index, topology_id):
                try:
                    topology_to_system[int(identifier)] = system_id
                except (TypeError, ValueError):
                    continue

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

    group1_indices = set(int(index) for index in (
        getattr(iface, "group1_indices", None)
        or getattr(getattr(iface, "group1", None), "ball_ndxs", ())
        or ()
    ))
    group2_indices = set(int(index) for index in (
        getattr(iface, "group2_indices", None)
        or getattr(getattr(iface, "group2", None), "ball_ndxs", ())
        or ()
    ))

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
        ball1_system = topology_to_system.get(ball1, ball1)
        ball2_system = topology_to_system.get(ball2, ball2)
        ball1_in_group1 = ball1_system in group1_indices
        ball2_in_group1 = ball2_system in group1_indices
        ball1_in_group2 = ball1_system in group2_indices
        ball2_in_group2 = ball2_system in group2_indices

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


def summarize_interface_surface_areas(iface, network=None, *, supported_surface_ids=None):
    """Return named area ledgers for one constructed interface network.

    Direct A/B area is deliberately independent from the all-surface network
    total.  ``supported_surface_ids`` is supplied by a materialized partial
    view; for a direct solve, the row's local ``complete`` flag is the only
    available support evidence.
    """
    network = iface.net if network is None else network
    surfaces = getattr(network, "surfs", None) if network is not None else None
    if surfaces is None:
        return {
            "status": "NOT_CALCULATED",
            "direct_facet_count": 0,
            "supported_direct_facet_count": 0,
            "unresolved_direct_facet_count": 0,
            "direct_interface_area_A2": None,
            "known_supported_partial_area_A2": None,
            "unresolved_direct_area_A2": None,
            "internal_supporting_network_area_A2": None,
            "total_constructed_network_area_A2": None,
            "known_constructed_network_area_A2": None,
        }

    direct, group1, group2, supporting = get_interface_surface_sets(
        iface, network=network
    )

    def rows(table):
        return [] if table is None else list(table.iterrows())

    def area(row):
        try:
            value = float(row.get("sa"))
        except (TypeError, ValueError):
            return None
        return value if np.isfinite(value) and value >= 0.0 else None

    def area_stats(table):
        entries = rows(table)
        values = [value for _, row in entries if (value := area(row)) is not None]
        missing = len(entries) - len(values)
        known = float(sum(values)) if values else None
        exact = known if not missing else None
        return entries, exact, known, missing

    direct_rows, direct_area, _, direct_missing = area_stats(direct)
    group1_rows, group1_area, group1_known, group1_missing = area_stats(group1)
    group2_rows, group2_area, group2_known, group2_missing = area_stats(group2)
    support_rows, support_area, support_known, support_missing = area_stats(supporting)
    all_rows, total_area, total_known, total_missing = area_stats(surfaces)

    supported_ids = None if supported_surface_ids is None else set(supported_surface_ids)
    supported_rows = []
    unresolved_rows = []
    for source_id, row in direct_rows:
        finite_area = area(row) is not None
        if supported_ids is not None:
            is_supported = source_id in supported_ids and finite_area
        else:
            complete = row.get("complete", True)
            is_supported = bool(complete) and finite_area
        (supported_rows if is_supported else unresolved_rows).append((source_id, row))

    supported_values = [area(row) for _, row in supported_rows]
    unresolved_values = [area(row) for _, row in unresolved_rows]
    supported_area = float(sum(supported_values)) if supported_values else None
    unresolved_area = (
        float(sum(unresolved_values))
        if all(value is not None for value in unresolved_values)
        else None
    )
    direct_count = len(direct_rows)
    supported_count = len(supported_rows)
    unresolved_count = len(unresolved_rows)
    if not direct_count:
        status = "NOT_CALCULATED"
    elif unresolved_count == 0 and direct_missing == 0:
        status = "VERIFIED"
    elif supported_count:
        status = "PARTIAL"
    else:
        status = "UNRESOLVED"

    non_direct_exact = (
        float(sum(value for value in (group1_area, group2_area, support_area)
                 if value is not None))
        if not any(value is None and count for value, count in (
            (group1_area, group1_missing),
            (group2_area, group2_missing),
            (support_area, support_missing),
        ))
        else None
    )
    non_direct_known = float(sum(value for value in (group1_known, group2_known, support_known)
                                 if value is not None))
    return {
        "status": status,
        "direct_facet_count": direct_count,
        "supported_direct_facet_count": supported_count,
        "unresolved_direct_facet_count": unresolved_count,
        "direct_interface_area_A2": direct_area,
        "known_supported_partial_area_A2": supported_area,
        "unresolved_direct_area_A2": unresolved_area,
        "internal_supporting_network_area_A2": non_direct_exact,
        "known_internal_supporting_network_area_A2": non_direct_known,
        "total_constructed_network_area_A2": total_area,
        "known_constructed_network_area_A2": total_known,
        "constructed_facet_count": len(all_rows),
        "constructed_facet_area_missing_count": total_missing,
    }


def _interface_physical_scope(iface):
    """Return the direct A/B surface closure used by every physical writer."""
    scope = export_scope(iface)

    def topology_position(table, identifier):
        try:
            position = table.index.get_loc(identifier)
        except (AttributeError, KeyError, TypeError, ValueError):
            position = int(identifier)
        if isinstance(position, slice):
            position = position.start
        return int(position)

    network = scope.network
    direct_surfaces, _, _, _ = get_interface_surface_sets(iface, network=network)
    direct_positions = set()
    for source_id in getattr(direct_surfaces, "index", ()):
        try:
            position = topology_position(network.surfs, source_id)
        except (AttributeError, KeyError, TypeError, ValueError):
            continue
        if position is not None and 0 <= int(position) < len(network.surfs):
            direct_positions.add(int(position))

    # A materialized view is still intersected with the physical A/B
    # classification. This prevents an overly broad view from exporting
    # intragroup or solvent geometry as an interface.
    surface_positions = [
        int(position) for position in scope.surfaces
        if int(position) in direct_positions
    ]

    edge_ids = {
        edge_id
        for position in surface_positions
        for edge_id in (network.surfs.iloc[position].get("edges", ()) or ())
    }
    vertex_ids = {
        vertex_id
        for position in surface_positions
        for vertex_id in (network.surfs.iloc[position].get("verts", ()) or ())
    }

    def positions(table, identifiers):
        result = []
        for identifier in sorted(identifiers):
            try:
                position = topology_position(table, identifier)
            except (AttributeError, KeyError, TypeError, ValueError):
                continue
            if position is not None and 0 <= int(position) < len(table):
                result.append(int(position))
        return tuple(result)

    edge_positions = positions(network.edges, edge_ids)
    vertex_positions = positions(network.verts, vertex_ids)
    support = set()
    for position in surface_positions:
        support.update(
            topology_position(network.balls, value)
            for value in (network.surfs.iloc[position].get("balls", ()) or ())
        )
    for position in edge_positions:
        support.update(
            topology_position(network.balls, value)
            for value in (network.edges.iloc[position].get("balls", ()) or ())
        )
    for position in vertex_positions:
        support.update(
            topology_position(network.balls, value)
            for value in (network.verts.iloc[position].get("balls", ()) or ())
        )
    support = tuple(sorted(value for value in support if 0 <= value < len(network.balls)))
    return replace(
        scope,
        cells=support,
        surfaces=tuple(sorted(set(surface_positions))),
        edges=edge_positions,
        vertices=vertex_positions,
        support_generators=support,
    )


def _physical_scope_closure_diagnostics(physical_scope):
    """Return missing selected incidence records without hiding them."""
    network = physical_scope.network

    def topology_position(table, identifier):
        try:
            position = table.index.get_loc(identifier)
        except (AttributeError, KeyError, TypeError, ValueError):
            position = int(identifier)
        if isinstance(position, slice):
            position = position.start
        position = int(position)
        if position < 0 or position >= len(table):
            raise IndexError(position)
        return position

    missing = []
    for surface_position in physical_scope.surfaces:
        surface = network.surfs.iloc[int(surface_position)]
        for field, table, label in (
            ("edges", network.edges, "edge"),
            ("verts", network.verts, "vertex"),
        ):
            for identifier in surface.get(field, ()) or ():
                try:
                    topology_position(table, identifier)
                except (AttributeError, KeyError, TypeError, ValueError, IndexError):
                    missing.append(
                        f"surface {surface_position} references missing {label} {identifier}"
                    )
    for edge_position in physical_scope.edges:
        edge = network.edges.iloc[int(edge_position)]
        for identifier in edge.get("verts", ()) or ():
            try:
                topology_position(network.verts, identifier)
            except (AttributeError, KeyError, TypeError, ValueError, IndexError):
                missing.append(
                    f"edge {edge_position} references missing vertex {identifier}"
                )
    return tuple(dict.fromkeys(missing))


def _csv_value(value):
    if value is None:
        return ""
    try:
        if not np.isfinite(float(value)):
            return ""
        return value
    except (TypeError, ValueError):
        return value


def _completeness(value):
    if value is None:
        return ""
    try:
        if isinstance(value, float) and np.isnan(value):
            return ""
    except TypeError:
        pass
    return "COMPLETE" if bool(value) else "INCOMPLETE"


def _interface_generator_metadata(network):
    metadata = {}
    balls = getattr(network, "balls", None)
    if balls is None:
        return metadata
    boundary = set(int(value) for value in (getattr(network, "boundary_indices", ()) or ()))
    for position, (row_index, row) in enumerate(balls.iterrows()):
        try:
            topology_id = int(row.get("num", position))
        except (TypeError, ValueError):
            topology_id = position
        system_id = row.get("system_num", topology_id)
        try:
            if system_id is None or (isinstance(system_id, float) and np.isnan(system_id)):
                system_id = None
            else:
                system_id = int(system_id)
        except (TypeError, ValueError):
            system_id = None
        if system_id is None:
            try:
                system_id = int(row_index)
            except (TypeError, ValueError):
                system_id = position
        stable_id = row.get("stable_id") or row.get("generator_stable_id")
        if not stable_id:
            stable_id = f"boundary:{topology_id}" if topology_id in boundary else f"generator:{topology_id}"
        record = {
            "position": position,
            "topology_id": topology_id,
            "system_id": system_id,
            "stable_id": str(stable_id),
            "is_boundary": position in boundary or topology_id in boundary,
        }
        for identifier in (position, row_index, topology_id):
            try:
                metadata[int(identifier)] = record
            except (TypeError, ValueError):
                continue
    return metadata


def _system_atom_metadata(system, system_id):
    if system_id is None or getattr(system, "balls", None) is None:
        return {"atom_id": "", "name": "", "residue": "", "residue_number": "", "chain": ""}
    row = None
    try:
        row = system.balls.loc[int(system_id)]
    except (IndexError, KeyError, TypeError, ValueError):
        try:
            matches = system.balls[system.balls["system_num"] == int(system_id)]
            if len(matches):
                row = matches.iloc[0]
        except (KeyError, TypeError, ValueError):
            pass
    if row is None:
        try:
            row = system.balls.iloc[int(system_id)]
        except (IndexError, KeyError, TypeError, ValueError):
            return {"atom_id": int(system_id), "name": "", "residue": "", "residue_number": "", "chain": ""}
    residue = row.get("res_name", row.get("res", ""))
    residue = getattr(residue, "name", residue)
    chain = row.get("chain_name", row.get("chain", ""))
    chain = getattr(chain, "name", chain)
    return {
        "atom_id": int(system_id),
        "name": row.get("name", ""),
        "residue": residue,
        "residue_number": row.get("res_seq", ""),
        "chain": chain,
    }


def _interface_scope_atom_sets(iface, physical_scope):
    network = physical_scope.network
    generator_metadata = _interface_generator_metadata(network)
    group_a = set(int(value) for value in (
        getattr(iface, "group1_indices", None)
        or getattr(getattr(iface, "group1", None), "ball_ndxs", ())
        or ()
    ))
    group_b = set(int(value) for value in (
        getattr(iface, "group2_indices", None)
        or getattr(getattr(iface, "group2", None), "ball_ndxs", ())
        or ()
    ))
    surface_selected = set()
    for position in physical_scope.surfaces:
        surface_selected.update(
            int(value) for value in (network.surfs.iloc[position].get("balls", ()) or ())
        )
    selected = set(surface_selected)
    for position in physical_scope.edges:
        selected.update(int(value) for value in (network.edges.iloc[position].get("balls", ()) or ()))
    for position in physical_scope.vertices:
        selected.update(int(value) for value in (network.verts.iloc[position].get("balls", ()) or ()))
    surface_system_ids = {
        generator_metadata[topology_id]["system_id"]
        for topology_id in surface_selected
        if topology_id in generator_metadata
        and generator_metadata[topology_id]["system_id"] is not None
    }
    closure_system_ids = {
        generator_metadata[topology_id]["system_id"]
        for topology_id in selected
        if topology_id in generator_metadata
        and generator_metadata[topology_id]["system_id"] is not None
    }
    return (
        sorted(system_id for system_id in surface_system_ids if system_id in group_a),
        sorted(system_id for system_id in surface_system_ids if system_id in group_b),
        sorted(
            system_id for system_id in closure_system_ids
            if system_id not in (group_a | group_b)
            and any(
                record["system_id"] == system_id and not record["is_boundary"]
                for record in generator_metadata.values()
            )
        ),
    )


def _validated_surrounding_atoms(iface, fallback):
    """Use Interface/GAUSS surrounding selection when it is available."""
    selector = getattr(iface, "get_surrounding_indices", None)
    if not callable(selector):
        return fallback
    group_a = getattr(iface, "group1_indices", None)
    if group_a is None:
        group_a = getattr(getattr(iface, "group1", None), "ball_ndxs", ())
    try:
        return sorted(set(int(value) for value in selector(group_a or ())))
    except (AttributeError, IndexError, KeyError, TypeError, ValueError):
        return fallback


def _residue_expanded_atoms(system, atoms):
    table = getattr(system, "balls", None)
    if table is None:
        return []
    keys = set()
    for atom in atoms:
        try:
            row = table.iloc[int(atom)]
        except (IndexError, KeyError, TypeError, ValueError):
            continue
        keys.add((row.get("chain_name", row.get("chain", "")),
                  row.get("res_seq", ""), row.get("res_name", row.get("res", ""))))
    selected = []
    for position, row in table.iterrows():
        key = (row.get("chain_name", row.get("chain", "")),
               row.get("res_seq", ""), row.get("res_name", row.get("res", "")))
        if key in keys:
            selected.append(int(position))
    return selected


def _write_interface_atom_exports(iface, physical_scope, include_surrounding=True):
    group_a, group_b, surrounding = _interface_scope_atom_sets(iface, physical_scope)
    surrounding = _validated_surrounding_atoms(iface, surrounding)
    system = iface.sys
    base_file = str((getattr(system, "files", {}) or {}).get("base_file") or "")
    canonical = getattr(system, "_canonical_output_layout", False)
    if base_file.lower().endswith(".txt") and not canonical:
        return {"group_a": group_a, "group_b": group_b, "surrounding": surrounding}
    if canonical:
        from vorpy.src.output.layout import canonical_groups_directory

        group_a_dir = canonical_groups_directory(system) / "A"
        group_b_dir = canonical_groups_directory(system) / "B"
        group_a_dir.mkdir(parents=True, exist_ok=True)
        group_b_dir.mkdir(parents=True, exist_ok=True)
        selection_a = sorted(set(int(value) for value in (getattr(iface.group1, "ball_ndxs", ()) or ())))
        selection_b = sorted(set(int(value) for value in (getattr(iface.group2, "ball_ndxs", ()) or ())))
        write_pdb(selection_a, "A", system, directory=group_a_dir)
        write_pdb(selection_b, "B", system, directory=group_b_dir)
        contact_atoms = sorted(set(group_a) | set(group_b))
        write_pdb(contact_atoms, "interface_atoms", system, directory=iface.dir)
        contact_path = os.path.join(iface.dir, "interface_atoms.pdb")
        if not os.path.exists(contact_path):
            open(contact_path, "w", encoding="utf-8").close()
        if include_surrounding:
            write_pdb(surrounding, "surr_atoms", system, directory=iface.dir)
            write_pdb(_residue_expanded_atoms(system, surrounding), "surr_residues", system, directory=iface.dir)
        return {"group_a": group_a, "group_b": group_b, "surrounding": surrounding}
    write_pdb(group_a, "group_a_atoms", system, directory=iface.dir)
    write_pdb(group_b, "group_b_atoms", system, directory=iface.dir)
    if include_surrounding:
        write_pdb(surrounding, "surr_atoms", system, directory=iface.dir)
        write_pdb(_residue_expanded_atoms(system, surrounding), "surr_residues", system, directory=iface.dir)
    return {"group_a": group_a, "group_b": group_b, "surrounding": surrounding}


def _write_interface_surrounding_exports(iface, physical_scope, *, atoms=True, residues=True):
    """Write only the explicitly requested validated surrounding selections."""
    _, _, surrounding = _interface_scope_atom_sets(iface, physical_scope)
    surrounding = _validated_surrounding_atoms(iface, surrounding)
    system = iface.sys
    if atoms:
        write_pdb(surrounding, "surr_atoms", system, directory=iface.dir)
    if residues:
        write_pdb(
            _residue_expanded_atoms(system, surrounding),
            "surr_residues", system, directory=iface.dir,
        )


def _generator_text(network, topology_ids, metadata):
    values = []
    for topology_id in topology_ids or ():
        record = metadata.get(int(topology_id))
        if record is not None:
            values.append(record["stable_id"])
    return ";".join(values)


def _surface_generator_records(network, balls, metadata, group_ids):
    records = []
    for topology_id in balls or ():
        record = metadata.get(int(topology_id))
        if record is None:
            continue
        if record["system_id"] in group_ids:
            records.append(record)
    return records


def _write_interface_physical_tables(iface, physical_scope):
    """Write readable physical A/B tables alongside the legacy group log."""
    network = physical_scope.network
    metadata = _interface_generator_metadata(network)
    group_a_ids = set(int(value) for value in (
        getattr(iface, "group1_indices", None)
        or getattr(getattr(iface, "group1", None), "ball_ndxs", ()) or ()
    ))
    group_b_ids = set(int(value) for value in (
        getattr(iface, "group2_indices", None)
        or getattr(getattr(iface, "group2", None), "ball_ndxs", ()) or ()
    ))
    source_reference = getattr(network, "source_network_id", None) or getattr(
        getattr(physical_scope, "view", None), "source_network_id", None
    ) or ""
    system = iface.sys
    edge_surfaces = {int(position): [] for position in physical_scope.edges}
    for surface_position in physical_scope.surfaces:
        surface_id = network.surfs.index[int(surface_position)]
        for edge_id in (network.surfs.iloc[int(surface_position)].get("edges", ()) or ()):
            edge_position = int(edge_id)
            if edge_position in edge_surfaces:
                edge_surfaces[edge_position].append(surface_id)
    vertex_surfaces = {int(position): [] for position in physical_scope.vertices}
    for surface_position in physical_scope.surfaces:
        surface_id = network.surfs.index[int(surface_position)]
        for vertex_id in (network.surfs.iloc[int(surface_position)].get("verts", ()) or ()):
            vertex_position = int(vertex_id)
            if vertex_position in vertex_surfaces:
                vertex_surfaces[vertex_position].append(surface_id)
    vertex_edges = {int(position): [] for position in physical_scope.vertices}
    for edge_position in physical_scope.edges:
        edge_id = network.edges.index[int(edge_position)]
        for vertex_id in (network.edges.iloc[int(edge_position)].get("verts", ()) or ()):
            vertex_position = int(vertex_id)
            if vertex_position in vertex_edges:
                vertex_edges[vertex_position].append(edge_id)

    surface_fields = [
        "Surface ID", "Group A Generator ID", "Group B Generator ID",
        "Group A Atom ID", "Group A Atom Name", "Group A Residue", "Group A Residue Number", "Group A Chain",
        "Group B Atom ID", "Group B Atom Name", "Group B Residue", "Group B Residue Number", "Group B Chain",
        "Surface Area (Å²)", "Mean Curvature (Å⁻¹)", "Gaussian Curvature (Å⁻²)",
        "Integrated Mean Curvature (Å)", "Integrated Mean Curvature Squared",
        "Integrated Gaussian Curvature", "Physical Completeness", "Source Network",
    ]
    surface_rows = []
    for position in physical_scope.surfaces:
        row = network.surfs.iloc[int(position)]
        balls = [int(value) for value in (row.get("balls", ()) or ())]
        a_records = _surface_generator_records(network, balls, metadata, group_a_ids)
        b_records = _surface_generator_records(network, balls, metadata, group_b_ids)
        a_atom = _system_atom_metadata(system, a_records[0]["system_id"] if a_records else None)
        b_atom = _system_atom_metadata(system, b_records[0]["system_id"] if b_records else None)
        surface_rows.append({
            "Surface ID": network.surfs.index[int(position)],
            "Group A Generator ID": ";".join(record["stable_id"] for record in a_records),
            "Group B Generator ID": ";".join(record["stable_id"] for record in b_records),
            "Group A Atom ID": _csv_value(a_atom["atom_id"]),
            "Group A Atom Name": a_atom["name"], "Group A Residue": a_atom["residue"],
            "Group A Residue Number": a_atom["residue_number"], "Group A Chain": a_atom["chain"],
            "Group B Atom ID": _csv_value(b_atom["atom_id"]),
            "Group B Atom Name": b_atom["name"], "Group B Residue": b_atom["residue"],
            "Group B Residue Number": b_atom["residue_number"], "Group B Chain": b_atom["chain"],
            "Surface Area (Å²)": _csv_value(row.get("sa")),
            "Mean Curvature (Å⁻¹)": _csv_value(row.get("mean_curv")),
            "Gaussian Curvature (Å⁻²)": _csv_value(row.get("gauss_curv")),
            "Integrated Mean Curvature (Å)": _csv_value(row.get("int_mean_curv")),
            "Integrated Mean Curvature Squared": _csv_value(row.get("int_mean_curv_sq")),
            "Integrated Gaussian Curvature": _csv_value(row.get("int_gauss_curv")),
            "Physical Completeness": _completeness(row.get("complete")),
            "Source Network": source_reference,
        })

    edge_fields = ["Edge ID", "Incident Selected Interface Surfaces", "Endpoint Vertex IDs",
                   "Supporting Generator IDs", "Edge Length (Å)", "Mean-Curvature Contribution (Å)",
                   "Validity/Completeness", "Source Network"]
    edge_rows = []
    for position in physical_scope.edges:
        row = network.edges.iloc[int(position)]
        curvature = row.get("int_mean_curv_by_ball")
        curvature_text = ""
        if isinstance(curvature, dict):
            values = []
            for topology_id, value in curvature.items():
                number = _csv_value(value)
                if number != "":
                    values.append(f"{metadata.get(int(topology_id), {}).get('stable_id', topology_id)}={number}")
            curvature_text = ";".join(values)
        edge_rows.append({
            "Edge ID": network.edges.index[int(position)],
            "Incident Selected Interface Surfaces": ";".join(str(value) for value in edge_surfaces[int(position)]),
            "Endpoint Vertex IDs": ";".join(str(network.verts.index[int(value)]) for value in (row.get("verts", ()) or ())),
            "Supporting Generator IDs": _generator_text(network, row.get("balls", ()), metadata),
            "Edge Length (Å)": _csv_value(row.get("length")),
            "Mean-Curvature Contribution (Å)": curvature_text,
            "Validity/Completeness": _completeness(row.get("complete")),
            "Source Network": source_reference,
        })

    vertex_fields = ["Vertex ID", "X (Å)", "Y (Å)", "Z (Å)", "Supporting Generator IDs",
                     "Incident Selected Surfaces", "Incident Selected Edges", "Gaussian-Curvature Contribution",
                     "Validity/Completeness", "Source Network"]
    vertex_rows = []
    for position in physical_scope.vertices:
        row = network.verts.iloc[int(position)]
        location = list(row.get("loc", ()))
        location += [None] * max(0, 3 - len(location))
        gaussian = row.get("int_gauss_curv")
        if gaussian is None:
            gaussian = row.get("int_gauss_curv_by_ball")
        if isinstance(gaussian, dict):
            gaussian = ";".join(
                f"{metadata.get(int(key), {}).get('stable_id', key)}={_csv_value(value)}"
                for key, value in gaussian.items() if _csv_value(value) != ""
            )
        vertex_rows.append({
            "Vertex ID": network.verts.index[int(position)],
            "X (Å)": _csv_value(location[0]), "Y (Å)": _csv_value(location[1]), "Z (Å)": _csv_value(location[2]),
            "Supporting Generator IDs": _generator_text(network, row.get("balls", ()), metadata),
            "Incident Selected Surfaces": ";".join(str(value) for value in vertex_surfaces[int(position)]),
            "Incident Selected Edges": ";".join(str(value) for value in vertex_edges[int(position)]),
            "Gaussian-Curvature Contribution": _csv_value(gaussian),
            "Validity/Completeness": _completeness(row.get("complete")),
            "Source Network": source_reference,
        })

    for filename, fields, rows in (
        ("interface_surfaces.csv", surface_fields, surface_rows),
        ("interface_edges.csv", edge_fields, edge_rows),
        ("interface_vertices.csv", vertex_fields, vertex_rows),
    ):
        with open(os.path.join(iface.dir, filename), "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def write_surface_statistics(info, title, surfaces):
    """
    Write area, curvature, contact, and overlap statistics for a collection
    of surfaces.
    """
    if surfaces is None:
        info.write(f"{title}:\n")
        info.write("  Not available\n\n")
        return

    if len(surfaces) == 0:
        info.write(f"{title}:\n")
        info.write("  Number of surfaces: 0\n")
        info.write("  Surface area and curvature: Not available\n\n")
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
    selections = _water_class_atom_indices(iface)
    if not any(selections.values()):
        return selections

    waters_dir = os.path.join(iface.dir, "waters")
    os.makedirs(waters_dir, exist_ok=True)
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
    groups = list(getattr(iface, "buried_water_groups", []) or [])
    if not groups:
        return None

    waters_dir = os.path.join(iface.dir, "waters")
    os.makedirs(waters_dir, exist_ok=True)
    path = os.path.join(waters_dir, "info.txt")

    with open(path, "w", encoding="utf-8") as info:
        info.write(f"Buried interface-water analysis - {iface.name}\n\n")
        info.write(f"Buried water groups: {len(groups)}\n")
        info.write(f"Total waters in buried groups: {sum(len(getattr(g, 'buried_water_metadata', {}).get('residues', [])) for g in groups)}\n\n")

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


def export_info(iface, directory=None, physical_scope=None):
    """
    Export a detailed geometric and topological summary of an Interface.
    """
    if physical_scope is None and getattr(iface, "net", None) is not None:
        physical_scope = _interface_physical_scope(iface)

    selected_atom_sets = (
        _interface_scope_atom_sets(iface, physical_scope)
        if physical_scope is not None else ([], [], [])
    )
    if physical_scope is not None:
        selected_atom_sets = (
            selected_atom_sets[0],
            selected_atom_sets[1],
            _validated_surrounding_atoms(iface, selected_atom_sets[2]),
        )
    source_network = physical_scope.network if physical_scope is not None else getattr(iface, "net", None)
    if physical_scope is not None:
        iface = copy(iface)
        iface.net = physical_scope.scoped_network()

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

    group1_atoms = sorted(set(
        getattr(iface, "group1_indices", None)
        or getattr(getattr(iface, "group1", None), "ball_ndxs", ())
        or ()
    ))
    group2_atoms = sorted(set(
        getattr(iface, "group2_indices", None)
        or getattr(getattr(iface, "group2", None), "ball_ndxs", ())
        or ()
    ))
    interface_atoms = sorted(
        set(group1_atoms) | set(group2_atoms)
    )

    (
        direct_surfaces,
        group1_surfaces,
        group2_surfaces,
        support_surfaces,
    ) = get_interface_surface_sets(iface)
    supported_surface_ids = None
    if physical_scope is not None and physical_scope.view is not None:
        supported_surface_ids = set(physical_scope.source_ids("surface"))
    area_ledger = summarize_interface_surface_areas(
        iface,
        network=source_network,
        supported_surface_ids=supported_surface_ids,
    )

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

    file_path = os.path.join(directory, "interface_info.txt")

    with open(file_path, "w", encoding="utf-8") as info:
        info.write(
            f"{getattr(iface, 'name', 'interface')} - "
            f"{getattr(getattr(iface, 'sys', None), 'name', 'system')}\n\n"
        )

        info.write("Interface definition:\n")
        info.write(f"  Interface ID: {interface_id}\n")
        info.write(f"  Group 1: {group1_name}\n")
        info.write(f"  Group 2: {group2_name}\n\n")

        settings = getattr(source_network, "settings", {}) or {}
        source_system = getattr(iface, "sys", None)
        source_files = getattr(source_system, "files", {}) or {}
        info.write("Input and network settings:\n")
        info.write(f"  Input structure: {source_files.get('base_file') or 'Not available'}\n")
        info.write(f"  Frame: {getattr(source_system, 'frame_index', 'Not available')}\n")
        info.write(f"  Scheme: {settings.get('net_type', 'Not available')}\n")
        info.write(f"  Surface resolution: {settings.get('surf_res', 'Not available')}\n")
        info.write(f"  Box size: {settings.get('box_size', 'Not available')}\n")
        info.write(f"  Maximum vertex radius: {settings.get('max_vert', 'Not available')}\n")
        info.write(f"  Solvent mode: {getattr(source_network, 'solvent_mode', None) or settings.get('solvent_mode', 'Not available')}\n")
        source_balls = getattr(source_network, "balls", None)
        info.write(f"  Source generators: {len(source_balls) if source_balls is not None else 'Not available'}\n\n")

        component_status = getattr(iface, "export_component_status", {}) or {}
        if component_status:
            info.write("Requested component status:\n")
            for component, status in sorted(component_status.items()):
                info.write(f"  {component}: {status}\n")
            info.write("\n")

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
        info.write(f"  Selected Group A contact atoms: {len(selected_atom_sets[0])}\n")
        info.write(f"  Selected Group B contact atoms: {len(selected_atom_sets[1])}\n")
        info.write(f"  Selected surrounding atoms: {len(selected_atom_sets[2])}\n\n")

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
        info.write("Scientific area accounting (constructed source network):\n")
        info.write(
            f"  Direct Group 1-Group 2 facet area (A²): "
            f"{_format_metric(area_ledger['direct_interface_area_A2'])}\n"
        )
        info.write(
            f"  Internal/supporting network area (A²): "
            f"{_format_metric(area_ledger['internal_supporting_network_area_A2'])}\n"
        )
        info.write(
            f"  Total constructed network area (A²): "
            f"{_format_metric(area_ledger['total_constructed_network_area_A2'])}\n"
        )
        info.write(
            f"  Known supported direct area (A²): "
            f"{_format_metric(area_ledger['known_supported_partial_area_A2'])}\n"
        )
        info.write(
            f"  Unresolved direct area (A²): "
            f"{_format_metric(area_ledger['unresolved_direct_area_A2'])}\n"
        )
        info.write(
            f"  Direct geometry status: {area_ledger['status']} "
            f"({area_ledger['supported_direct_facet_count']}/"
            f"{area_ledger['direct_facet_count']} supported facets)\n\n"
        )
        direct_surface_count = len(direct_surfaces) if direct_surfaces is not None else None
        direct_area = (
            None
            if not direct_surface_count
            else _safe_sum(_numeric_series(direct_surfaces, 'sa'))
        )
        info.write(f"  Shared physical surface count: {direct_surface_count if direct_surface_count is not None else 'Not available'}\n")
        info.write(f"  Interface area (Å²): {_format_metric(direct_area)}\n")
        info.write(f"  Selected edges: {num_edges if num_edges is not None else 'Not available'}\n")
        info.write(f"  Selected vertices: {num_verts if num_verts is not None else 'Not available'}\n")
        if physical_scope is not None:
            info.write("  Exported physical subset:\n")
            info.write(f"    Surfaces: {len(physical_scope.surfaces)}\n")
            info.write(f"    Edges: {len(physical_scope.edges)}\n")
            info.write(f"    Vertices: {len(physical_scope.vertices)}\n")

        view_completeness = (
            getattr(physical_scope.view, "completeness", None)
            if physical_scope is not None and physical_scope.view is not None
            else None
        )
        unresolved_surface_ids = tuple(
            getattr(view_completeness, "unresolved_surface_ids", ()) or ()
        )
        if view_completeness is not None and not view_completeness.proven:
            completeness = "PARTIAL" if direct_surface_count else "UNRESOLVED"
            diagnostics = "; ".join(view_completeness.reasons) or "Unresolved physical coverage"
        elif not direct_surface_count:
            completeness = "NOT_CALCULATED"
            diagnostics = "No direct A/B physical facets selected"
        else:
            surface_table = getattr(net, "surfs", None) if net is not None else None
            if surface_table is None or "complete" not in surface_table.columns:
                completeness = "Not available"
                diagnostics = "Not available"
            else:
                selected_positions = (
                    physical_scope.surfaces
                    if physical_scope is not None
                    else tuple(range(len(surface_table)))
                )
                values = [surface_table.iloc[int(position)].get("complete")
                          for position in selected_positions]
                complete = all(
                    value is not None
                    and not (isinstance(value, float) and np.isnan(value))
                    and bool(value)
                    for value in values
                )
                closure = _physical_scope_closure_diagnostics(physical_scope)
                completeness = "COMPLETE" if complete and not closure else "INCOMPLETE"
                diagnostics = "None" if completeness == "COMPLETE" else (
                    "Selected physical facets include incomplete records"
                    if not closure else "Selected physical closure is incomplete: " + "; ".join(closure)
                )
        if physical_scope is not None and physical_scope.view is not None:
            info.write(f"  Missing physical surfaces: {len(unresolved_surface_ids)}\n")
            if unresolved_surface_ids:
                info.write(
                    "  Missing surface IDs: "
                    + ";".join(str(value) for value in unresolved_surface_ids)
                    + "\n"
                )
            unresolved_area = getattr(view_completeness, "unresolved_surface_area", None)
            if unresolved_area is not None:
                info.write(f"  Missing physical area (Å²): {_format_metric(unresolved_area)}\n")
        info.write(f"  Completeness: {completeness}\n")
        info.write(f"  Unresolved diagnostics: {diagnostics}\n")
        info.write("\n")
        representation = getattr(iface, "_geometry_representation", None)
        selection = (
            "selected_source_view" if physical_scope is not None and physical_scope.view is not None
            else "direct_inter_group" if physical_scope is not None
            else "full_network (all cached network surfaces)"
        )
        info.write(
            f"  Physical export selection: {selection}\n"
        )
        info.write(
            f"  Alpha metadata selection: "
            f"{getattr(representation, 'selection_mode', 'not recorded')}\n\n"
        )

        write_surface_statistics(
            info,
            "Selected interface-network surface statistics"
            if physical_scope is not None else "Full interface-network surface statistics",
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

        from vorpy.src.interface.geometry_analysis import write_geometry_info
        try:
            write_geometry_info(info, getattr(iface, 'geometry_analysis', None))
        except (AttributeError, KeyError, TypeError, ValueError):
            info.write("Geometry analysis: unavailable for this cached interface.\n\n")
        # Physical interface logs/info must remain usable when an optional
        # cached dual has missing simplex records.
        if physical_scope is None or not physical_scope.is_materialized_view:
            from vorpy.src.geometry.visualization.interface import write_interface_dual_summary
            try:
                write_interface_dual_summary(info, iface)
            except (AttributeError, KeyError, IndexError, TypeError, ValueError) as error:
                detail = error.args[0] if isinstance(error, KeyError) and error.args else str(error)
                info.write(f"\nApollonius dual visualization: unavailable ({detail})\n")

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

    # Keep the established group-style filename only for legacy directories.
    if not getattr(getattr(iface, "sys", None), "_canonical_output_layout", False):
        shutil.copyfile(file_path, os.path.join(directory, "network_info.txt"))


def interface_exports(iface, all_=False, atoms=False, interface_atoms=False, surr_atoms=False,
                      surr_resids=False, surfs=False, edges=False, verts=False, logs=False, info=False,
                      group_info=False, round_to=3, dual=False, buried_water=None, contact_tables=False):
    """
    Export data belonging to an Interface and its dedicated Network.
    """
    if iface.net is None:
        print(f'Interface "{iface.name}" has no network to export.')
        return
    if getattr(iface, 'dir', None) is None and hasattr(iface, 'set_dir'):
        iface.set_dir()

    profile = {"timings": {}, "counts": {}, "measurements": {}}
    stage_callback = getattr(
        getattr(iface, 'sys', None), '_cli_export_stage_callback', None
    )
    if not _interface_export_profile_enabled() and not callable(stage_callback):
        profile = None
    elif callable(stage_callback):
        if profile is None:
            profile = {"timings": {}, "counts": {}, "measurements": {}}
            total_started = time.perf_counter()
        profile['_cli_stage_callback'] = stage_callback
    total_started = time.perf_counter() if profile is not None else None

    scope = _profile_stage(
        profile,
        "scope_selection",
        lambda: _interface_physical_scope(iface),
    )
    if profile is not None:
        profile["counts"].update(_interface_export_counts(scope))
        area_values = []
        for position in scope.surfaces:
            value = scope.network.surfs.iloc[int(position)].get("sa")
            try:
                value = float(value)
            except (TypeError, ValueError):
                value = None
            if value is None or not np.isfinite(value):
                area_values = None
                break
            area_values.append(value)
        profile["measurements"]["selected_surface_area_A2"] = (
            None if not area_values else float(sum(area_values))
        )

    mesh_net = None
    if surfs or edges or verts or all_:
        mesh_net = _profile_stage(
            profile,
            "mesh_network",
            scope.mesh_network,
        )
    scoped_iface = None
    if info or logs or all_:
        scoped_iface = copy(iface)
        scoped_iface.net = scope.scoped_network()

    if iface.dir is None:
        set_dir = getattr(iface, "set_dir", None)
        if set_dir is not None:
            set_dir()
        else:
            iface.dir = os.path.join(iface.sys.files["dir"], iface.name)

    os.makedirs(iface.dir, exist_ok=True)
    if (dual or (all_ and iface.net.settings.get('net_type', 'aw') == 'aw')) and not scope.is_materialized_view:
        def export_dual():
            from vorpy.src.geometry.visualization.interface import export_interface_dual_visualization
            export_interface_dual_visualization(iface)

        _profile_stage(profile, "dual_visualization", export_dual)

    # Buried-water networks are subordinate interface products. Export them
    # once with the interface information pass rather than during an arbitrary
    # geometry export (surfaces/edges/verts). The helper retains its own guard
    # for explicit repeated info exports.
    if (info or all_) and buried_water is not False and not scope.is_materialized_view:
        _profile_stage(
            profile,
            "buried_water_exports",
            lambda: export_buried_water_groups(iface),
        )

    if group_info or all_:
        def export_group_info():
            exported_group_ids = set()

            for position, group in enumerate((iface.group1, iface.group2), start=1):
                if group is None:
                    continue

                group_id = str(getattr(group, "group_id", None) or group.name)

                if group_id in exported_group_ids:
                    continue

                exported_group_ids.add(group_id)

                export_interface_group_info(
                    group=group,
                    directory=iface.dir,
                    file_prefix=f"group_{position}",
                )

        _profile_stage(profile, "group_info", export_group_info)

    contact_atoms_requested = atoms or interface_atoms or all_
    if contact_atoms_requested:
        atom_counts = _profile_stage(
            profile,
            "contact_atom_pdb",
            lambda: _write_interface_atom_exports(
                iface, scope,
                include_surrounding=(all_ or (atoms and not interface_atoms)),
            ),
        )
        if profile is not None and atom_counts is not None:
            profile["counts"].update({
                "group_a_contact_atoms": len(atom_counts["group_a"]),
                "group_b_contact_atoms": len(atom_counts["group_b"]),
                "surrounding_atoms": len(atom_counts["surrounding"]),
            })

    if surr_atoms or surr_resids:
        _profile_stage(
            profile,
            "surrounding_atom_pdb",
            lambda: _write_interface_surrounding_exports(
                iface, scope, atoms=surr_atoms, residues=surr_resids,
            ),
        )

    if info or all_:
        _profile_stage(
            profile,
            "interface_info",
            lambda: export_info(scoped_iface, physical_scope=scope),
        )

    if logs or all_:
        _profile_stage(
            profile,
            "legacy_interface_logs",
            lambda: write_interface_logs(scoped_iface, round_to=round_to),
        )
    if logs or contact_tables or all_:
        _profile_stage(
            profile,
            "physical_interface_tables",
            lambda: _write_interface_physical_tables(iface, scope),
        )

    if verts or all_:
        def export_vertices():
            if scope.vertices:
                write_off_verts(
                    mesh_net,
                    scope.vertices,
                    directory=iface.dir,
                    file_name="interface_verts",
                    color=mesh_net.settings.get("vert_col", "red"),
                    vert_rad=mesh_net.settings.get("vertex_size", DEFAULT_VERTEX_RADIUS),
                )
            else:
                _write_empty_off(iface.dir, "interface_verts")

        _profile_stage(profile, "vertex_mesh", export_vertices)

    if edges or all_:
        def export_edges():
            if scope.edges:
                write_edges(
                    mesh_net,
                    scope.edges,
                    directory=iface.dir,
                    file_name="interface_edges",
                    color=mesh_net.settings.get("edge_col", "grey"),
                    radius=mesh_net.settings.get("edge_width", DEFAULT_EDGE_RADIUS),
                )
            else:
                _write_empty_off(iface.dir, "interface_edges")

        _profile_stage(profile, "edge_mesh", export_edges)

    if surfs or all_:
        def export_surfaces():
            if scope.surfaces:
                write_surfs(
                    mesh_net,
                    scope.surfaces,
                    directory=iface.dir,
                    file_name="interface_surfs",
                    universal_max=False,
                )
            else:
                _write_empty_off(iface.dir, "interface_surfs")

        _profile_stage(profile, "surface_mesh", export_surfaces)

    if profile is not None:
        stage_callback = profile.get('_cli_stage_callback')
        current_stage = profile.get('_cli_stage_current')
        if callable(stage_callback) and current_stage is not None:
            stage_callback(current_stage, 'finish')
        profile.pop('_cli_stage_callback', None)
        profile.pop('_cli_stage_current', None)
        profile["timings"]["total"] = time.perf_counter() - total_started
        iface.interface_export_profile = profile
        if _interface_export_profile_enabled():
            _interface_export_profile_summary(profile)

    return profile

def export_interface_group_info(group, directory, file_prefix="group"):
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

    info_path = os.path.join(directory, f"{file_prefix}_info.txt")

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
            file_name=f"{file_prefix}_atoms",
            directory=directory,
            sys=group.sys,
        )
