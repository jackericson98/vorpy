import time
import numpy as np
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import ndx_search
from vorpy.src.calculations import verify_site
from vorpy.src.calculations import calc_flat_vert
from vorpy.src.calculations import calc_vert
from vorpy.src.calculations import verify_aw_cached


def _verify_seed_aw(loc, rad, vert_balls, locs, rads):
    """Check a seed against all generators without rebuilding neighborhoods."""
    if box_search(loc) is None:
        return False
    return verify_aw_cached(np.asarray(loc), rad, locs, rads, *vert_balls)


def find_site_container_slow(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, net_type, box=None,
                             group_ndxs=None, iface_grps=None, metrics=None, printing=False,
                             full_max_vert=None):
    """
    Search thoroughly for a valid vertex associated with a three-ball edge.

    The search begins locally and expands geometrically until a valid site is
    found or the allowable vertex radius is reached. Verification arrays are
    prepared once and rejected fourth-ball candidates are retained
    between expansion steps so they are not recalculated.

    Parameters
    ----------
    edge_balls : list of int
        Three ball indices defining the edge.
    locs : array-like
        Ball-center coordinates.
    rads : array-like
        Ball radii.
    b_verts : list of list
        Ball-to-vertex adjacency.
    vert_ndxs : list
        Defining balls for already discovered vertices.
    max_vert : float
        Initial search extent and accepted absolute vertex radius. For AW seed
        searches this is intentionally the fast initial limit (normally the
        user's limit divided by ten).
    full_max_vert : float, optional
        User's full maximum allowable weighted Voronoi vertex radius. AW seed
        searches progressively retry up to this value while reusing the
        verification arrays. Other network types ignore this argument.
    net_type : {'aw', 'pow', 'prm'}
        Network geometry being solved.
    box : list, optional
        Geometric bounds for accepted primary vertex locations.
    group_ndxs : collection, optional
        Group constraint for the fourth defining ball.
    iface_grps : tuple, optional
        Two collections defining an interface.
    metrics : dict, optional
        Retained for API compatibility.
    printing : bool, optional
        Retained for API compatibility.

    Returns
    -------
    object or None
        The verified vertex returned by ``find_site`` or ``None``.
    """
    initial_max_vert = float(max_vert)
    if net_type == 'aw':
        locs = np.asarray(locs, dtype=float)
        rads = np.asarray(rads, dtype=float)
    final_max_vert = float(full_max_vert if full_max_vert is not None else max_vert)
    if final_max_vert < initial_max_vert:
        final_max_vert = initial_max_vert

    # AW uses a deliberately cheap first seed search. If it finds no valid
    # seed, expand the limit geometrically rather than treating the fast
    # limit as a correctness bound. The final value is always the user's
    # exact limit.
    if net_type == 'aw' and full_max_vert is not None and initial_max_vert < final_max_vert:
        search_limits = []
        limit = initial_max_vert
        while True:
            limit = min(limit, final_max_vert)
            search_limits.append(limit)
            if limit >= final_max_vert:
                break
            limit = min(limit * 2.0, final_max_vert)
    else:
        search_limits = [initial_max_vert]

    # The fourth ball only needs group filtering when none of the edge balls
    # already belongs to the requested group.
    check_ndxs = group_ndxs is not None and not any(ball in group_ndxs for ball in edge_balls)

    my_boxes = [box_search(loc=locs[ball]) for ball in edge_balls]
    # AW verifies against the prepared full arrays. Other geometries retain
    # the shared spatial neighborhood used by their verification routines.
    surr_balls = [] if net_type == 'aw' else get_balls(cells=my_boxes, dist=final_max_vert)

    for attempt, attempt_max_vert in enumerate(search_limits, start=1):
        attempt_start = time.perf_counter()
        invalid_ndxs = []
        vert = None

        # Always test the current maximum extent rather than jumping over it.
        mv_inc = min(0.45, attempt_max_vert)
        while vert is None:
            vert, invalid_ndxs = find_site(
                edge_balls=edge_balls, locs=locs, rads=rads, b_verts=b_verts, vert_ndxs=vert_ndxs,
                max_vert=attempt_max_vert, mv_inc=mv_inc, net_type=net_type, invalid_ndxs=invalid_ndxs,
                check_balls=check_ndxs, surr_balls=surr_balls, my_boxes=my_boxes, group_ndxs=group_ndxs,
                iface_grps=iface_grps, metrics=metrics, box=box
            )

            if vert is not None or mv_inc >= attempt_max_vert:
                break

            mv_inc = min(mv_inc * 10, attempt_max_vert)

        if metrics is not None:
            metrics['seed_search_attempts'] = metrics.get('seed_search_attempts', 0) + 1
            metrics.setdefault('seed_search_limits', []).append(attempt_max_vert)
            metrics.setdefault('seed_search_times', []).append(
                time.perf_counter() - attempt_start
            )

        if vert is not None:
            if metrics is not None:
                metrics['seed_success_max_vert'] = attempt_max_vert
            return vert

    return None


def find_site(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, net_type, invalid_ndxs=None,
              check_balls=True, surr_balls=None, vn_1=None, vn_1_loc=None, group_ndxs=None, iface_grps=None,
              metrics=None, my_boxes=None, box=None):
    """
    Find a valid fourth-ball vertex for an existing three-ball edge.

    Candidate fourth balls are gathered from the current spatial search range,
    filtered by group/interface constraints and previously rejected candidates,
    solved geometrically, and verified against the full surrounding search set.

    AW candidates use full-array verification without neighborhood setup, while
    POW and PRM retain full surrounding-ball verification. Rejected fourth-ball
    indices are stored as integers so expanding searches do not reconsider them.

    Returns
    -------
    tuple
        ``(vertex, invalid_ndxs)`` where ``vertex`` is the selected verified
        vertex or ``None``.
    """
    if invalid_ndxs is None:
        invalid_ndxs = []

    edge_ndxs = edge_balls[:]
    if net_type == 'aw':
        locs = np.asarray(locs, dtype=float)
        rads = np.asarray(rads, dtype=float)
    invalid_set = set(invalid_ndxs)

    # Determine whether the edge is missing either side of an interface.
    missing_iface_grps = []
    if iface_grps is not None:
        edge_ball_set = set(edge_balls)
        missing_iface_grps = [set(group_indices) for group_indices in iface_grps
                              if not edge_ball_set.intersection(group_indices)]

    vert_ball_ndxs = edge_ndxs if vn_1 is None else vn_1

    if my_boxes is None:
        my_boxes = [box_search(loc=locs[ball]) for ball in edge_balls]

    test_balls = [ball for ball in get_balls(cells=my_boxes, dist=mv_inc) if ball not in invalid_set]

    # Reuse the full surrounding-ball set supplied by the container.
    if surr_balls is None:
        surr_balls = get_balls(cells=my_boxes, dist=max_vert)

    new_test_balls = []

    for ball in test_balls:
        if ball in vert_ball_ndxs:
            continue
        if check_balls and ball not in group_ndxs:
            continue
        if missing_iface_grps and not all(ball in group_indices for group_indices in missing_iface_grps):
            continue

        ball_ndxs = sorted(edge_ndxs + [ball])

        # If this four-ball vertex already exists, this edge does not need a
        # second geometric search from the current traversal.
        if vert_ndxs is not None and len(vert_ndxs) > 0:
            check_verts = [vert_ndxs[vert_ndx] for vert_ndx in b_verts[ball_ndxs[0]]]
            my_vert_ndx = ndx_search(check_verts, ball_ndxs)
            if my_vert_ndx < len(check_verts) and ball_ndxs == check_verts[my_vert_ndx]:
                return None, invalid_ndxs

        new_test_balls.append(ball)

    verts = []

    for ball in new_test_balls:
        vert_balls = sorted(edge_ndxs + [ball])
        vert_loc2, vert_rad2 = None, None

        if net_type == 'pow':
            vert_loc, vert_rad = calc_flat_vert(locs=[locs[ndx] for ndx in vert_balls],
                                                rads=[rads[ndx] for ndx in vert_balls], power=True)
        elif net_type == 'prm':
            vert_loc, vert_rad = calc_flat_vert(locs=[locs[ndx] for ndx in vert_balls],
                                                rads=[rads[ndx] for ndx in vert_balls], power=False)
        else:
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(
                locs=[locs[ndx] for ndx in vert_balls],
                rads=[rads[ndx] for ndx in vert_balls]
            )

        if vert_loc is None:
            invalid_ndxs.append(ball)
            continue

        if box is not None and any(box[0][axis] > vert_loc[axis] or vert_loc[axis] > box[1][axis] for axis in range(3)):
            invalid_ndxs.append(ball)
            continue

        # Verify the primary and optional secondary solutions. AW scans the
        # prepared arrays and exits on the first blocker; POW and PRM
        # retain their existing full surrounding-ball verification behavior.
        if net_type == 'aw':
            primary_valid = abs(vert_rad) < max_vert and _verify_seed_aw(
                vert_loc, vert_rad, vert_balls, locs, rads
            )
            secondary_valid = (
                vert_loc2 is not None
                and abs(vert_rad2) < max_vert
                and _verify_seed_aw(vert_loc2, vert_rad2, vert_balls, locs, rads)
            )
        else:
            filtered_test_balls = [test_ball for test_ball in surr_balls if test_ball not in vert_balls]
            test_locs = np.asarray([locs[test_ball] for test_ball in filtered_test_balls])
            test_rads = np.asarray([rads[test_ball] for test_ball in filtered_test_balls])

            primary_valid = abs(vert_rad) < max_vert and verify_site(
                loc=np.asarray(vert_loc), rad=vert_rad, test_locs=test_locs, test_rads=test_rads, net_type=net_type
            )
            secondary_valid = (
                vert_loc2 is not None
                and abs(vert_rad2) < max_vert
                and verify_site(
                    loc=np.asarray(vert_loc2), rad=vert_rad2,
                    test_locs=test_locs, test_rads=test_rads, net_type=net_type
                )
            )

        if primary_valid:
            if len(verts) > 0 and verts[0]['rad'] < vert_rad:
                return [verts[0], metrics], invalid_ndxs

            verts.append({'balls': vert_balls, 'loc': vert_loc, 'rad': vert_rad, 'loc2': None, 'rad2': None})

            if secondary_valid:
                verts[-1]['loc2'], verts[-1]['rad2'] = vert_loc2, vert_rad2

        elif secondary_valid:
            verts.append({'balls': vert_balls, 'loc': vert_loc2, 'rad': vert_rad2, 'loc2': None, 'rad2': None})

        if ball not in invalid_set:
            invalid_ndxs.append(ball)
            invalid_set.add(ball)

    if len(verts) == 0:
        return None, invalid_ndxs

    if len(verts) == 1 or verts[0]['rad'] < verts[1]['rad']:
        return [verts[0], metrics], invalid_ndxs

    return [verts[1], metrics], invalid_ndxs
