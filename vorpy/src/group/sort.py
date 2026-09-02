import os
import numpy as np
from vorpy.src.calculations import combine_inertia_tensors
from vorpy.src.calculations import calc_total_inertia_tensor
from vorpy.src.calculations import ndx_search
from vorpy.src.calculations import calc_surf_sa


def _group_net_row_positions(group):
    """
    Map Group SYSTEM atom indices to positional rows in group.net.balls.

    This prevents group.ball_ndxs from being used directly with .iloc when
    net.balls has a distinct topology/system_num index mapping.
    """
    net_balls = group.net.balls
    system_to_row = {}

    for row_pos, (_, atom) in enumerate(net_balls.iterrows()):
        try:
            topology_id = int(atom.get("num", row_pos))
        except (TypeError, ValueError):
            topology_id = int(row_pos)

        system_id = None
        if "system_num" in net_balls.columns:
            value = atom.get("system_num", None)
            try:
                if value is not None and not (isinstance(value, float) and np.isnan(value)):
                    system_id = int(value)
            except (TypeError, ValueError):
                system_id = None

        if system_id is None:
            system_id = topology_id

        system_to_row.setdefault(system_id, row_pos)

    rows = []
    missing = []

    for system_id in (int(_) for _ in group.ball_ndxs):
        row_pos = system_to_row.get(system_id)
        if row_pos is None:
            missing.append(system_id)
        else:
            rows.append(row_pos)

    return rows



def _build_orientation_context(group):
    """
    Build all Group/surface-orientation lookup data once per get_info() call.

    Surface topology uses network topology IDs, while Group membership uses
    parent-System IDs.  This function resolves that mapping once, rather than
    rebuilding it for every boundary surface.
    """
    net_balls = group.net.balls
    topology_to_system = {}
    topology_to_loc = {}
    topology_to_rad = {}

    has_system_num = "system_num" in net_balls.columns

    # itertuples is substantially cheaper than iterrows for this one-time pass.
    columns = list(net_balls.columns)
    num_idx = columns.index("num") if "num" in columns else None
    system_idx = columns.index("system_num") if has_system_num else None
    loc_idx = columns.index("loc") if "loc" in columns else None
    rad_idx = columns.index("rad") if "rad" in columns else None

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
        "group_system_ids": set(int(_) for _ in group.ball_ndxs),
        "topology_to_system": topology_to_system,
        "topology_to_loc": topology_to_loc,
        "topology_to_rad": topology_to_rad,
    }


def _surface_group_orientation(surf, context, tol=1e-12):
    """
    Return +1 when the implicit-surface gradient points from the Group interior
    toward the exterior atom, -1 when it points inward, and 0 when orientation
    cannot be determined robustly.

    The sign applies only to mean curvature. Gaussian curvature and H^2 are
    invariant under normal reversal.
    """
    try:
        ball0, ball1 = [int(_) for _ in surf["balls"]]
    except (KeyError, TypeError, ValueError):
        return 0

    topology_to_system = context["topology_to_system"]
    group_system_ids = context["group_system_ids"]

    ball0_inside = topology_to_system.get(ball0) in group_system_ids
    ball1_inside = topology_to_system.get(ball1) in group_system_ids

    if ball0_inside == ball1_inside:
        return 0

    inside_topology = ball0 if ball0_inside else ball1
    outside_topology = ball1 if ball0_inside else ball0

    topology_to_loc = context["topology_to_loc"]
    inside_loc = topology_to_loc.get(inside_topology)
    outside_loc = topology_to_loc.get(outside_topology)
    if inside_loc is None or outside_loc is None:
        return 0

    try:
        point = np.asarray(surf.get("com"), dtype=float)
        func = np.asarray(surf.get("func"), dtype=float)
    except (TypeError, ValueError):
        return 0

    if point.shape != (3,) or func.size < 9:
        return 0

    A, B, C, D, E, F, G, Hc, Ic = func[:9]
    x, y, z = point

    grad = np.array([
        2.0 * A * x + D * y + F * z + G,
        2.0 * B * y + D * x + E * z + Hc,
        2.0 * C * z + F * x + E * y + Ic,
    ], dtype=float)

    outward = outside_loc - inside_loc
    grad_norm = np.linalg.norm(grad)
    outward_norm = np.linalg.norm(outward)

    if (
        not np.isfinite(grad_norm)
        or not np.isfinite(outward_norm)
        or grad_norm < tol
        or outward_norm < tol
    ):
        return 0

    dot = float(np.dot(grad, outward))
    if abs(dot) <= tol * grad_norm * outward_norm:
        return 0

    return 1 if dot > 0.0 else -1



def _surface_group_orientation_details(surf, context, tol=1e-12):
    """
    Return the same orientation decision as _surface_group_orientation together
    with the geometry used to make that decision.  This is for validation/debug
    output only and does not alter the production calculation.
    """
    details = {
        "orientation": 0,
        "inside_topology": None,
        "outside_topology": None,
        "inside_system": None,
        "outside_system": None,
        "inside_rad": np.nan,
        "outside_rad": np.nan,
        "gradient": None,
        "outward": None,
        "dot": np.nan,
        "reason": "",
    }

    try:
        ball0, ball1 = [int(_) for _ in surf["balls"]]
    except (KeyError, TypeError, ValueError):
        details["reason"] = "invalid ball pair"
        return details

    topology_to_system = context["topology_to_system"]
    group_system_ids = context["group_system_ids"]

    ball0_system = topology_to_system.get(ball0)
    ball1_system = topology_to_system.get(ball1)
    ball0_inside = ball0_system in group_system_ids
    ball1_inside = ball1_system in group_system_ids

    if ball0_inside == ball1_inside:
        details["reason"] = "both endpoints have same group membership"
        return details

    inside_topology = ball0 if ball0_inside else ball1
    outside_topology = ball1 if ball0_inside else ball0
    details["inside_topology"] = inside_topology
    details["outside_topology"] = outside_topology
    details["inside_system"] = topology_to_system.get(inside_topology)
    details["outside_system"] = topology_to_system.get(outside_topology)

    topology_to_rad = context.get("topology_to_rad", {})
    details["inside_rad"] = topology_to_rad.get(inside_topology, np.nan)
    details["outside_rad"] = topology_to_rad.get(outside_topology, np.nan)

    topology_to_loc = context["topology_to_loc"]
    inside_loc = topology_to_loc.get(inside_topology)
    outside_loc = topology_to_loc.get(outside_topology)
    if inside_loc is None or outside_loc is None:
        details["reason"] = "missing endpoint coordinates"
        return details

    try:
        point = np.asarray(surf.get("com"), dtype=float)
        func = np.asarray(surf.get("func"), dtype=float)
    except (TypeError, ValueError):
        details["reason"] = "invalid COM/function"
        return details

    if point.shape != (3,) or func.size < 9:
        details["reason"] = "invalid COM/function shape"
        return details

    A, B, C, D, E, F, G, Hc, Ic = func[:9]
    x, y, z = point
    grad = np.array([
        2.0 * A * x + D * y + F * z + G,
        2.0 * B * y + D * x + E * z + Hc,
        2.0 * C * z + F * x + E * y + Ic,
    ], dtype=float)
    outward = outside_loc - inside_loc

    details["gradient"] = grad
    details["outward"] = outward

    grad_norm = np.linalg.norm(grad)
    outward_norm = np.linalg.norm(outward)

    if (
        not np.isfinite(grad_norm)
        or not np.isfinite(outward_norm)
        or grad_norm < tol
        or outward_norm < tol
    ):
        details["reason"] = "degenerate gradient/outward vector"
        return details

    dot = float(np.dot(grad, outward))
    details["dot"] = dot

    if abs(dot) <= tol * grad_norm * outward_norm:
        details["reason"] = "ambiguous gradient/outward dot product"
        return details

    details["orientation"] = 1 if dot > 0.0 else -1
    details["reason"] = "ok"
    return details


def _print_orientation_debug(surf_ndx, surf, surf_sa, details, classification):
    """Print one human-readable surface-orientation validation record."""
    raw_c = surf.get("int_mean_curv", 0.0)
    raw_c = 0.0 if raw_c is None else float(raw_c)
    raw_avg_h = raw_c / surf_sa if surf_sa > 0.0 else 0.0
    orientation = int(details.get("orientation", 0) or 0)
    oriented_c = orientation * raw_c
    oriented_avg_h = orientation * raw_avg_h

    print()
    print("=" * 72)
    print(f"ORIENTATION DEBUG - {classification.upper()} SURFACE")
    print("=" * 72)
    print(f"Surface index:              {surf_ndx}")
    print(f"VorPy flat flag:            {bool(surf.get('flat', False))}")
    print(f"Surface area:               {surf_sa:.8f} A^2")
    print(f"Inside topology atom:       {details.get('inside_topology')}")
    print(f"Outside topology atom:      {details.get('outside_topology')}")
    print(f"Inside system atom:         {details.get('inside_system')}")
    print(f"Outside system atom:        {details.get('outside_system')}")
    print(f"Inside radius:              {details.get('inside_rad', np.nan):.8f} A")
    print(f"Outside radius:             {details.get('outside_rad', np.nan):.8f} A")
    print(f"Raw integrated H dA:        {raw_c:.8f} A")
    print(f"Raw area-weighted H:        {raw_avg_h:.8f} A^-1")
    print(f"Orientation multiplier:     {orientation:+d}")
    print(f"Oriented integrated H dA:   {oriented_c:.8f} A")
    print(f"Oriented area-weighted H:   {oriented_avg_h:.8f} A^-1")
    print(f"grad(F) dot outward:        {details.get('dot', np.nan):.12g}")
    print(f"Orientation status:         {details.get('reason', '')}")
    grad = details.get("gradient")
    outward = details.get("outward")
    if grad is not None:
        print("grad(F) at surface COM:     "
              f"[{grad[0]:.8f}, {grad[1]:.8f}, {grad[2]:.8f}]")
    if outward is not None:
        print("inside -> outside vector:   "
              f"[{outward[0]:.8f}, {outward[1]:.8f}, {outward[2]:.8f}]")
    print("=" * 72)


def _accumulate_oriented_surface_curvature(group, surf, surf_sa, orientation, flat_tol=1e-10):
    """
    Accumulate Group-oriented mean curvature plus convex/concave/flat area.

    True flat surfaces are classified before orientation because H = K = 0 and
    normal reversal has no curvature consequence.  Only non-flat surfaces that
    cannot be oriented are counted as failed/unoriented QC cases.
    """
    raw_int_mean = surf.get("int_mean_curv", 0.0)
    raw_int_mean = 0.0 if raw_int_mean is None else float(raw_int_mean)

    # VorPy marks AW equal-radius bisectors (and prm/pow surfaces) as flat.
    # These belong in the flat bucket directly, not in orientation failures.
    if bool(surf.get("flat", False)):
        group.flat_sa += surf_sa
        group.flat_surface_count += 1
        return

    if orientation == 0:
        group.unoriented_sa += surf_sa
        group.unoriented_surface_count += 1
        group.curved_unoriented_sa += surf_sa
        group.curved_unoriented_surface_count += 1
        return

    oriented_int_mean = orientation * raw_int_mean
    group.oriented_int_mean_curv += oriented_int_mean
    group.oriented_surface_count += 1

    points = surf.get("points", None)
    tris = surf.get("tris", None)
    tri_curvs = surf.get("mean_tri_curvs", None)

    try:
        if points is not None and tris is not None and tri_curvs is not None:
            pts = np.asarray(points, dtype=float)
            tri_idx = np.asarray(tris, dtype=int)
            h = orientation * np.asarray(tri_curvs, dtype=float)

            if (
                tri_idx.ndim == 2
                and tri_idx.shape[1] == 3
                and len(tri_idx) == len(h)
                and len(tri_idx) > 0
            ):
                tri_pts = pts[tri_idx]
                ab = tri_pts[:, 1, :] - tri_pts[:, 0, :]
                ac = tri_pts[:, 2, :] - tri_pts[:, 0, :]
                tri_sa = 0.5 * np.linalg.norm(np.cross(ab, ac), axis=1)

                valid = np.isfinite(tri_sa) & np.isfinite(h) & (tri_sa > 0.0)
                tri_sa = tri_sa[valid]
                h = h[valid]

                convex = h > flat_tol
                concave = h < -flat_tol
                near_flat = ~(convex | concave)

                if np.any(convex):
                    group.convex_sa += float(np.sum(tri_sa[convex]))
                    group.convex_int_mean_curv += float(
                        np.sum(h[convex] * tri_sa[convex])
                    )

                if np.any(concave):
                    group.concave_sa += float(np.sum(tri_sa[concave]))
                    group.concave_int_mean_curv += float(
                        np.sum(h[concave] * tri_sa[concave])
                    )

                # This is numerical near-flatness on otherwise curved surfaces.
                if np.any(near_flat):
                    group.flat_sa += float(np.sum(tri_sa[near_flat]))

                return
    except (TypeError, ValueError, IndexError):
        pass

    # Backward-compatible fallback for old/checkpoint surfaces lacking
    # triangle-level curvature.
    avg_h = oriented_int_mean / surf_sa if surf_sa > 0.0 else 0.0
    if avg_h > flat_tol:
        group.convex_sa += surf_sa
        group.convex_int_mean_curv += oriented_int_mean
    elif avg_h < -flat_tol:
        group.concave_sa += surf_sa
        group.concave_int_mean_curv += oriented_int_mean
    else:
        group.flat_sa += surf_sa

def get_info(group):
    """
    Gathers and calculates comprehensive information about a molecular group, including geometric, physical, and structural properties.

    This function performs a detailed analysis of a molecular group, calculating:
    - Surface area and volume metrics
    - Center of mass (both geometric and van der Waals)
    - Mass and density properties
    - Moment of inertia tensors
    - Layer-based surface properties

    Parameters:
        group (Group): The Group object to analyze. Must have a valid network (group.net) and ball indices (group.ball_ndxs).

    Returns:
        None. Results are stored in the Group object's attributes:
        - sa (float): Total surface area of the group
        - vol (float): Total volume of the group
        - vdw_vol (float): Van der Waals volume
        - density (float): Ratio of van der Waals volume to total volume
        - mass (float): Total mass of atoms
        - com (numpy.ndarray): Center of mass coordinates
        - vdw_com (list): Van der Waals center of mass
        - spatial_moment (list): Spatial moment tensor
        - moi (list): Moment of inertia tensor

    Examples:
        >>> from vorpy.src.group import Group
        >>> # Create and analyze a group
        >>> group = Group(sys=my_system, name='protein_A')
        >>> group.add_balls(atom_indices=range(100))
        >>> group.build()
        >>> # Calculate group properties
        >>> get_info(group)
        >>> # Access calculated properties
        >>> print(f"Surface area: {group.sa:.2f} Å²")
        >>> print(f"Volume: {group.vol:.2f} Å³")
        >>> print(f"Center of mass: {group.com}")
        >>> print(f"Moment of inertia: {group.moi}")
    """
    # Reset the group's data attributes
    group.sa, group.vol, group.vdw_vol, group.density, group.mass = 0, 0, 0, 0, 0
    # Reset the group's curvature geometry
    group.avg_mean_curv = 0.0
    group.avg_gauss_curv = 0.0

    group.int_mean_curv = 0.0
    group.int_mean_curv_sq = 0.0
    group.int_gauss_curv = 0.0
    # Group-oriented curvature. Raw surface/network values above remain unchanged.
    group.oriented_int_mean_curv = 0.0
    group.oriented_avg_mean_curv = 0.0
    group.convex_sa = 0.0
    group.concave_sa = 0.0
    group.flat_sa = 0.0
    group.unoriented_sa = 0.0
    group.flat_surface_count = 0
    group.curved_unoriented_sa = 0.0
    group.curved_unoriented_surface_count = 0
    group.convex_sa_fraction = 0.0
    group.concave_sa_fraction = 0.0
    group.flat_sa_fraction = 0.0
    group.convex_int_mean_curv = 0.0
    group.concave_int_mean_curv = 0.0
    group.abs_int_mean_curv = 0.0
    group.convex_avg_mean_curv = 0.0
    group.concave_avg_mean_curv = 0.0
    group.oriented_surface_count = 0
    group.unoriented_surface_count = 0
    # Center of masses
    com, vdw_com = [0, 0, 0], [0, 0, 0]
    # Get the balls in the group
    group_net_rows = _group_net_row_positions(group)
    group_balls = group.net.balls.iloc[group_net_rows].to_dict(orient='records')
    # Get the volume of the group
    for i, ball in enumerate(group_balls):
        # Check for the ball to be complets
        if not ball['complete']:
            continue
        # Add the volume to that of the group
        group.vol += ball['vol']
        # Add the vdw volume to that of the group
        group.vdw_vol += ball['vdw_vol']
        # Add the mass to that of the group
        group.mass += ball['mass']
        # Add to the coms
        com = [com[j] + ball['com'][j] * ball['vol'] for j in range(3)]
        vdw_com = [vdw_com[j] + ball['loc'][j] * ball['mass'] for j in range(3)]
    # Check to see if the volume is greater than 0
    if group.vol > 0:
        # Calculate the density
        group.density = group.vdw_vol / group.vol
        # Calculate the center of mass
        group.com = np.array([com[j] / group.vol for j in range(3)])
        # Calculate the vdw center of mass
        group.vdw_com = [vdw_com[j] / group.vdw_vol for j in range(3)]
    # Check to see if the moi has been calculated
    if group_net_rows and 'moi' in group.net.balls.iloc[group_net_rows[0]]:
        # Calculate the spatial moment
        group.spatial_moment = combine_inertia_tensors([_['moi'] for _ in group_balls], [_['com'] for _ in group_balls],
                                                       group.com, [_['vol'] for _ in group_balls])
    if group.vdw_vol > 0:
        group.moi = calc_total_inertia_tensor(group_balls, group.vdw_com)
    # Check to see if the first layer has been calculated
    if group.layer_surfs is None or len(group.layer_surfs) == 0:
        group.get_layers(max_layers=1)

    # Build the topology/system/location orientation maps once.  The previous
    # implementation rebuilt these maps once per surface, which was O(S * N).
    orientation_context = _build_orientation_context(group)

    # Optional full AW orientation validation. Enable with:
    #   PowerShell: $env:VORPY_ORIENTATION_DEBUG="1"
    #
    # Every Group-boundary surface is checked against the AW radius/sign rule:
    #   r_inside < r_outside  -> oriented H >= 0  (convex)
    #   r_inside > r_outside  -> oriented H <= 0  (concave)
    #   r_inside = r_outside  -> flat=True and H ~= 0
    #
    # Detailed records are printed only for disagreements.
    orientation_debug = os.environ.get("VORPY_ORIENTATION_DEBUG", "0").strip().lower() in {
        "1", "true", "yes", "on"
    }

    orientation_validation = {
        "curved_tested": 0,
        "curved_agree": 0,
        "curved_disagree": 0,
        "equal_radius_tested": 0,
        "equal_radius_flat": 0,
        "equal_radius_nonflat": 0,
        "equal_radius_nonzero_h": 0,
        "orientation_failed": 0,
        "missing_radius": 0,
    }
    orientation_validation_tol = 1e-10

    # Check to see if there are any layers
    if group.layer_surfs is not None and len(group.layer_surfs) > 0:

        # The first surface layer defines the exposed group boundary.
        for surf_ndx in group.layer_surfs[0]:

            # Get the surface
            surf = group.net.surfs.iloc[surf_ndx]

            # --------------------------------------------------------------
            # Surface area
            # --------------------------------------------------------------

            if surf['sa'] is None or surf['sa'] == 0:
                surf_sa = calc_surf_sa(
                    tris=surf['tris'],
                    points=surf['points']
                )
            else:
                surf_sa = float(surf['sa'])

            group.sa += surf_sa

            # --------------------------------------------------------------
            # Integrated curvature
            # --------------------------------------------------------------

            int_mean_curv = surf.get('int_mean_curv', 0.0)
            int_mean_curv_sq = surf.get('int_mean_curv_sq', 0.0)
            int_gauss_curv = surf.get('int_gauss_curv', 0.0)

            # Protect against missing/None values
            if int_mean_curv is not None:
                group.int_mean_curv += float(int_mean_curv)

            if int_mean_curv_sq is not None:
                group.int_mean_curv_sq += float(int_mean_curv_sq)

            if int_gauss_curv is not None:
                group.int_gauss_curv += float(int_gauss_curv)

            # --------------------------------------------------------------
            # Group-oriented mean curvature and convex/concave decomposition
            # --------------------------------------------------------------
            if bool(surf.get("flat", False)):
                orientation = 0
                orientation_details = None
            elif orientation_debug:
                orientation_details = _surface_group_orientation_details(
                    surf, orientation_context
                )
                orientation = orientation_details["orientation"]
            else:
                orientation_details = None
                orientation = _surface_group_orientation(
                    surf, orientation_context
                )

            _accumulate_oriented_surface_curvature(
                group=group,
                surf=surf,
                surf_sa=surf_sa,
                orientation=orientation,
            )

            if orientation_debug:
                if orientation_details is None:
                    orientation_details = _surface_group_orientation_details(
                        surf, orientation_context
                    )

                rin = orientation_details.get("inside_rad", np.nan)
                rout = orientation_details.get("outside_rad", np.nan)
                raw_c = surf.get("int_mean_curv", 0.0)
                raw_c = 0.0 if raw_c is None else float(raw_c)
                oriented_avg_h = (
                    orientation * raw_c / surf_sa
                    if surf_sa > 0.0 else 0.0
                )
                is_flat = bool(surf.get("flat", False))

                if not np.isfinite(rin) or not np.isfinite(rout):
                    orientation_validation["missing_radius"] += 1

                elif abs(rin - rout) <= orientation_validation_tol:
                    orientation_validation["equal_radius_tested"] += 1

                    if is_flat:
                        orientation_validation["equal_radius_flat"] += 1
                    else:
                        orientation_validation["equal_radius_nonflat"] += 1
                        _print_orientation_debug(
                            surf_ndx, surf, surf_sa,
                            orientation_details,
                            "EQUAL-RADIUS NONFLAT DISAGREEMENT"
                        )

                    if abs(raw_c) > orientation_validation_tol:
                        orientation_validation["equal_radius_nonzero_h"] += 1
                        _print_orientation_debug(
                            surf_ndx, surf, surf_sa,
                            orientation_details,
                            "EQUAL-RADIUS NONZERO-H DISAGREEMENT"
                        )

                else:
                    orientation_validation["curved_tested"] += 1

                    if orientation == 0:
                        orientation_validation["orientation_failed"] += 1
                        orientation_validation["curved_disagree"] += 1
                        _print_orientation_debug(
                            surf_ndx, surf, surf_sa,
                            orientation_details,
                            "CURVED ORIENTATION FAILURE"
                        )
                    else:
                        expected_sign = 1 if rin < rout else -1

                        if (
                            abs(oriented_avg_h) <= orientation_validation_tol
                            or np.sign(oriented_avg_h) == expected_sign
                        ):
                            orientation_validation["curved_agree"] += 1
                        else:
                            orientation_validation["curved_disagree"] += 1
                            _print_orientation_debug(
                                surf_ndx, surf, surf_sa,
                                orientation_details,
                                "RADIUS/SIGN DISAGREEMENT"
                            )

    # ----------------------------------------------------------------------
    # Area-weighted group curvature
    # ----------------------------------------------------------------------

    if group.sa > 0.0:
        group.avg_mean_curv = (
                group.int_mean_curv / group.sa
        )

        group.avg_gauss_curv = (
                group.int_gauss_curv / group.sa
        )


        group.oriented_avg_mean_curv = (
                group.oriented_int_mean_curv / group.sa
        )

        group.convex_sa_fraction = group.convex_sa / group.sa
        group.concave_sa_fraction = group.concave_sa / group.sa
        group.flat_sa_fraction = group.flat_sa / group.sa

    if group.convex_sa > 0.0:
        group.convex_avg_mean_curv = (
                group.convex_int_mean_curv / group.convex_sa
        )

    if group.concave_sa > 0.0:
        group.concave_avg_mean_curv = (
                group.concave_int_mean_curv / group.concave_sa
        )

    group.abs_int_mean_curv = (
            group.convex_int_mean_curv - group.concave_int_mean_curv
    )


    if orientation_debug:
        tested_total = (
            orientation_validation["curved_tested"]
            + orientation_validation["equal_radius_tested"]
        )
        curved_tested = orientation_validation["curved_tested"]
        curved_agree = orientation_validation["curved_agree"]
        curved_disagree = orientation_validation["curved_disagree"]

        agreement_pct = (
            100.0 * curved_agree / curved_tested
            if curved_tested > 0 else 100.0
        )

        print()
        print("=" * 72)
        print("AW ORIENTATION VALIDATION SUMMARY")
        print("=" * 72)
        print(f"Boundary surfaces tested:        {tested_total}")
        print()
        print(f"Curved surfaces tested:          {curved_tested}")
        print(f"Radius/sign agreements:          {curved_agree}")
        print(f"Radius/sign disagreements:       {curved_disagree}")
        print(f"Curved orientation failures:     {orientation_validation['orientation_failed']}")
        print(f"Radius/sign agreement:           {agreement_pct:.3f} %")
        print()
        print(f"Equal-radius surfaces:           {orientation_validation['equal_radius_tested']}")
        print(f"Equal-radius marked flat:        {orientation_validation['equal_radius_flat']}")
        print(f"Equal-radius non-flat:           {orientation_validation['equal_radius_nonflat']}")
        print(f"Equal-radius nonzero H:          {orientation_validation['equal_radius_nonzero_h']}")
        print(f"Missing radius metadata:         {orientation_validation['missing_radius']}")
        print("-" * 72)

        passed = (
            curved_disagree == 0
            and orientation_validation["orientation_failed"] == 0
            and orientation_validation["equal_radius_nonflat"] == 0
            and orientation_validation["equal_radius_nonzero_h"] == 0
            and orientation_validation["missing_radius"] == 0
        )

        if passed:
            print("RESULT: PASS")
            print("All tested AW Group-boundary surfaces agree with the radius/sign rule.")
        else:
            print("RESULT: FAIL")
            print("One or more AW Group-boundary surfaces require inspection.")
        print("=" * 72)
        print()


def add_balls(grp, ball_list):
    """
    Adds atoms to a group while maintaining sorted order and preventing duplicates.

    This function efficiently integrates new atoms into a group's existing atom list. It uses binary search
    to maintain a sorted order of atom indices and ensures no duplicate atoms are added. The function
    handles various input types including lists of atoms from molecules, residues, or direct atom selections.

    Parameters:
        grp (Group): The Group object to which atoms will be added
        ball_list (list): List of atom indices to be added to the group. Can be from various sources
                         (e.g., molecule.atoms, residue.atoms, or direct atom selections)

    Returns:
        None. The group's ball_ndxs and atms attributes are updated with the new atoms.
    """
    # Check to see if the index list has been instantiated
    if grp.ball_ndxs is None:
        grp.ball_ndxs = []
    # Go through the atom_list
    for sphere in ball_list:
        # Get the atom's location
        sphere_ndx = ndx_search(np.array(grp.ball_ndxs), sphere)
        # Check to see if we have found this atom before
        if sphere_ndx >= len(grp.ball_ndxs) or grp.ball_ndxs[sphere_ndx] != sphere:
            grp.ball_ndxs.insert(sphere_ndx, sphere)
            grp.atms.insert(sphere_ndx, sphere)
