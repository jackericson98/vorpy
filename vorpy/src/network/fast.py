import time
import bisect
import numpy as np
import warnings
from numba import jit
from vorpy.src.calculations import calc_flat_vert
from vorpy.src.calculations import calc_vert
from vorpy.src.calculations import verify_aw
from vorpy.src.calculations import verify_pow
from vorpy.src.calculations import verify_prm
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import calc_circ

warnings.filterwarnings("error", category=RuntimeWarning)

KNOWN_EDGE = object()
POW_PRM_METRICS = {
    'container_calls': 0,
    'surrounding': 0.0,
    'candidate_gather': 0.0,
    'candidate_filter': 0.0,
    'calc_vert': 0.0,
    'verify_arrays': 0.0,
    'verify': 0.0,
    'candidates': 0,
    'verify_balls': 0,
}


@jit(nopython=True, cache=True)
def verify_pow_cached(loc, rad, test_locs, test_rads, skip0=-1, skip1=-1, skip2=-1, skip3=-1):
    """Verify Power geometry against cached edge arrays while skipping defining balls."""
    for i in range(len(test_locs)):
        if i == skip0 or i == skip1 or i == skip2 or i == skip3:
            continue
        dx = test_locs[i, 0] - loc[0]
        dy = test_locs[i, 1] - loc[1]
        dz = test_locs[i, 2] - loc[2]
        if dx * dx + dy * dy + dz * dz - test_rads[i] * test_rads[i] < rad:
            return False
    return True


@jit(nopython=True, cache=True)
def verify_prm_cached(loc, rad, test_locs, skip0=-1, skip1=-1, skip2=-1, skip3=-1):
    """Verify primitive geometry against cached edge arrays while skipping defining balls."""
    rad2 = rad * rad
    for i in range(len(test_locs)):
        if i == skip0 or i == skip1 or i == skip2 or i == skip3:
            continue
        dx = test_locs[i, 0] - loc[0]
        dy = test_locs[i, 1] - loc[1]
        dz = test_locs[i, 2] - loc[2]
        if dx * dx + dy * dy + dz * dz < rad2:
            return False
    return True


def find_site_container(edge_balls, locs, rads, b_verts, vert_ndxs,
                        max_vert, net_type, box=None, vn_1=None, vn_1_loc=None,
                        group_ndxs=None, metrics=None, printing=False, max_ball_rad=None):
    """
    Find a neighboring vertex associated with a three-ball edge.

    The search begins in a small spatial neighborhood around the defining
    edge and expands until a valid neighboring vertex is found or
    ``max_vert`` is reached.

    For additively weighted (AW) networks, candidate balls are gathered
    progressively around the edge and candidate vertices are verified using
    a local neighborhood around the calculated vertex. Power and primitive
    networks retain the full surrounding-ball verification path.

    Group and interface searches may require the fourth defining ball to
    belong to a particular selection.

    Parameters
    ----------
    edge_balls : list of int
        Three ball indices defining the edge.
    locs : array-like
        Ball-center coordinates for the full system.
    rads : array-like
        Ball radii for the full system.
    b_verts : list of list
        Vertex indices associated with each ball.
    vert_ndxs : list
        Defining ball indices for previously discovered vertices.
    max_vert : float
        Maximum permitted vertex radius/search extent.
    net_type : {'aw', 'pow', 'prm'}
        Network geometry being solved.
    box : list, optional
        Geometric bounds for accepted vertex locations.
    vn_1 : list of int, optional
        Defining balls of the previous vertex.
    vn_1_loc : array-like, optional
        Location of the previous vertex.
    group_ndxs : collection or tuple, optional
        Group constraint for the fourth defining ball. A tuple represents
        the two sides of an interface.
    metrics : dict, optional
        Performance metrics populated when profiling is enabled.
    printing : bool, optional
        Enable detailed diagnostic printing.
    max_ball_rad : float, optional
        Maximum radius in the system. Used by local AW verification.

    Returns
    -------
    object or None
        The neighboring vertex returned by the selected site-search routine,
        or ``None`` if no valid neighbor is found.
    """
    # Initialize the search state.
    invalid_ndxs, vert = [], None

    if max_ball_rad is None:
        max_ball_rad = max(rads)

    # Without a previous vertex, treat the edge itself as the known definition.
    if vn_1 is None:
        vn_1 = edge_balls
    # If this edge already belongs to another known vertex, its neighboring
    # topology has already been traversed and no geometric search is required.
    # Skip geometric searching when this edge already connects to another
    # previously discovered vertex.
    if vert_ndxs is not None and b_verts is not None:
        common_verts = set(b_verts[edge_balls[0]]).intersection(
            b_verts[edge_balls[1]],
            b_verts[edge_balls[2]]
        )

        current_balls = set(vn_1)

        for vert_ndx in common_verts:
            if set(vert_ndxs[vert_ndx]) != current_balls:
                return None
    # Determine whether the fourth ball must come from a particular group.
    required_group = None

    if isinstance(group_ndxs, tuple):

        satisfied = [
            any(ball in grp for ball in edge_balls)
            for grp in group_ndxs
        ]

        if all(satisfied):
            required_group = None

        elif satisfied.count(True) == 1:
            required_group = group_ndxs[satisfied.index(False)]

        else:
            return None

    elif group_ndxs is not None:

        required_group = group_ndxs

        for ball in edge_balls:
            if ball in group_ndxs:
                required_group = None
                break

    my_boxes = [box_search(loc=locs[edge_balls[_]]) for _ in range(3)]
    # AW verifies locally. Power caches one verification neighborhood per edge.
    # Primitive retains its existing surrounding-ball verification path.
    surr_balls = None
    surr_locs = None
    surr_rads = None
    surr_lookup = None

    if net_type == 'pow':
        POW_PRM_METRICS['container_calls'] += 1
        metric_start = time.perf_counter()

        surr_balls = get_balls(cells=my_boxes, dist=max_vert)
        surr_locs = np.asarray([locs[ball] for ball in surr_balls], dtype=float)
        surr_rads = np.asarray([rads[ball] for ball in surr_balls], dtype=float)
        surr_lookup = {ball: i for i, ball in enumerate(surr_balls)}

        POW_PRM_METRICS['surrounding'] += time.perf_counter() - metric_start

    elif net_type == 'prm':
        POW_PRM_METRICS['container_calls'] += 1
        metric_start = time.perf_counter()

        # Primitive verification also reuses one edge-level neighborhood.
        # Verification depends only on center distances, not ball ordering.
        surr_balls = get_balls(cells=my_boxes, dist=max_vert)
        surr_locs = np.asarray([locs[ball] for ball in surr_balls], dtype=float)
        surr_lookup = {ball: i for i, ball in enumerate(surr_balls)}

        POW_PRM_METRICS['surrounding'] += time.perf_counter() - metric_start

    mv_inc = min(0.45, max_vert)

    while vert is None:
        # Search for the neighboring vertex within the current range.
        if net_type == 'aw':
            vert, invalid_ndxs = find_site_aw(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                                              required_group is not None, surr_balls, my_boxes, invalid_ndxs, vn_1,
                                              vn_1_loc, box=box, group_balls=required_group, metrics=metrics,
                                              printing=printing, max_ball_rad=max_ball_rad)
            if vert is KNOWN_EDGE:
                return None
        elif net_type == 'pow':
            vert, invalid_ndxs = find_site_pow(
                edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                required_group is not None, surr_balls, my_boxes, invalid_ndxs, vn_1,
                box, vn_1_loc, group_ndxs=required_group, metrics=metrics,
                surr_locs=surr_locs, surr_rads=surr_rads, surr_lookup=surr_lookup
            )
        elif net_type == 'prm':
            vert, invalid_ndxs = find_site_del(
                edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                required_group is not None, surr_balls, my_boxes, invalid_ndxs, vn_1,
                box, vn_1_loc, group_ndxs=required_group, metrics=metrics,
                surr_locs=surr_locs, surr_lookup=surr_lookup
            )

        if vert is not None or mv_inc >= max_vert:
            break

        mv_inc = min(mv_inc * 10, max_vert)
    # Return the vertex if found
    return vert


def find_site_del(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, check_ndxs, surr_balls,
                  my_boxes, invalid_ndxs, vn_1, box=None, vn_1_loc=None, group_ndxs=None, metrics=None,
                  surr_locs=None, surr_lookup=None):
    """
    Finds a new vertex in a Delaunay network by searching for valid ball combinations.

    This function searches for a valid vertex that can be formed by combining the edge balls with
    a new ball from the surrounding region. It verifies that the potential vertex:
    - Is not already part of the previous vertex
    - Is not in the list of invalid balls
    - Has not been previously found
    - Meets any group membership criteria if specified

    The search is performed by:
    1. Getting surrounding balls within the specified range
    2. Sorting them by distance from the previous vertex
    3. Checking each ball for potential vertex formation
    4. Returning None if no valid vertex is found

    Parameters
    ----------
    edge_balls : list
        List of ball indices that form the edge
    locs : list
        List of ball locations
    rads : list
        List of ball radii
    b_verts : dict
        Dictionary mapping ball indices to their vertices
    vert_ndxs : list
        List of vertex indices
    max_vert : float
        Maximum distance to search for vertices
    mv_inc : float
        Current search increment
    check_ndxs : bool
        Whether to check group membership
    surr_balls : list
        List of surrounding ball indices
    my_boxes : list
        List of search boxes
    invalid_ndxs : list
        List of invalid ball indices
    vn_1 : list
        List of balls in previous vertex
    box : list, optional
        Search box boundaries
    vn_1_loc : numpy.ndarray, optional
        Location of previous vertex
    group_ndxs : list, optional
        List of ball indices in the group
    metrics : dict, optional
        Dictionary for storing performance metrics

    Returns
    -------
    tuple
        - The new vertex if found, None otherwise
        - Updated list of invalid ball indices
    """
    # Get the balls that should not be a part of the new vertex
    edge_ndxs = edge_balls[:]
    invalid_ndxs_set = set(invalid_ndxs)
    metric_start = time.perf_counter()
    # Get the balls not in the invalid balls that are within the range specified
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs_set]
    POW_PRM_METRICS['candidate_gather'] += time.perf_counter() - metric_start
    # Sort the test balls to be in order by distance from the previous vert location
    if vn_1_loc is None:
        vn_1_loc = calc_com([locs[_] for _ in edge_ndxs])

    dists = [calc_dist(np.array(locs[_]), np.array(vn_1_loc)) for _ in test_balls]
    test_balls = [_ for x, _ in sorted(zip(dists, test_balls))]
    metric_start = time.perf_counter()
    # Instantiate the list for test vertices to be calculated later. This saves us from sorting the vertices balls twice
    test_verts = []
    # Go through the surrounding balls to look for vertices that have been found before and filter out edge balls
    for ball in test_balls:
        # If the ball is in the previous vertex move on
        if ball in vn_1:
            continue
        # Check if we need to check and if so check for the ball in the list
        if check_ndxs and ball not in group_ndxs:
            continue
        # If we have found the vertex before it is not the previous vertex return
        ball_ndxs = edge_ndxs + [ball]
        ball_ndxs.sort()
        # Get the vertex's index/insert index
        check_verts = [vert_ndxs[_] for _ in b_verts[ball_ndxs[0]]]
        # Take the potential index of the ball in group
        my_vert_ndx = bisect.bisect_left(check_verts, ball_ndxs)
        # If the index returned is larger than the list or the vertex at the index is not equal to the ball_ndxs were ok
        if my_vert_ndx < len(check_verts) and ball_ndxs == check_verts[my_vert_ndx]:
            return None, invalid_ndxs
        # Add the vertex indices to the test_vertices for calculation
        test_verts.append((ball_ndxs, ball))
    POW_PRM_METRICS['candidate_filter'] += time.perf_counter() - metric_start
    POW_PRM_METRICS['candidates'] += len(test_verts)
    # Go through each ball in the given test balls. Extremely optimized
    for vert in test_verts:

        # Add the vertex ball to the
        vert_balls, ball = vert
        metric_start = time.perf_counter()
        try:
            v_loc, vert_rad = calc_flat_vert(
                locs=[locs[_] for _ in vert_balls],
                rads=[rads[_] for _ in vert_balls],
                power=False
            )
        finally:
            POW_PRM_METRICS['calc_vert'] += time.perf_counter() - metric_start

        # Catch the none location case
        if v_loc is None:
            invalid_ndxs.append(ball)
            continue
        # Check if the vert is outside the box
        if box is not None and any([box[0][k] > v_loc[k] or v_loc[k] > box[1][k] for k in range(3)]):
            continue
        # Reject oversized Primitive vertices before verification.
        if vert_rad >= max_vert:
            invalid_ndxs.append(ball)
            continue

        # Reuse the edge-level verification locations and skip the four defining
        # balls in the compiled empty-sphere test.
        metric_start = time.perf_counter()
        skip0 = surr_lookup.get(vert_balls[0], -1)
        skip1 = surr_lookup.get(vert_balls[1], -1)
        skip2 = surr_lookup.get(vert_balls[2], -1)
        skip3 = surr_lookup.get(vert_balls[3], -1)
        POW_PRM_METRICS['verify_arrays'] += time.perf_counter() - metric_start
        POW_PRM_METRICS['verify_balls'] += len(surr_locs)

        metric_start = time.perf_counter()
        valid = verify_prm_cached(np.asarray(v_loc), vert_rad, surr_locs, skip0, skip1, skip2, skip3)
        POW_PRM_METRICS['verify'] += time.perf_counter() - metric_start

        if valid:
            # Return the validated ball and the invalidated list
            return [{'balls': vert_balls, 'loc': v_loc, 'rad': vert_rad}, metrics], invalid_ndxs
        else:
            # Add the ball to the invalid balls list if it isn't verified
            invalid_ndxs.append(ball)
    # Return the non-vertex and invalid balls
    return None, invalid_ndxs


def find_site_pow(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, check_ndxs, surr_balls,
                  my_boxes, invalid_ndxs, vn_1, box=None, vn_1_loc=None, group_ndxs=None, metrics=None,
                  surr_locs=None, surr_rads=None, surr_lookup=None):
    """
    Finds a new vertex in a power network by searching for valid ball combinations.

    This function searches for a valid vertex that can be formed by combining the edge balls with
    a new ball from the surrounding region. It verifies that the potential vertex:
    - Is not already part of the previous vertex
    - Is not in the list of invalid balls
    - Has not been previously found
    - Meets any group membership criteria if specified

    The search is performed by:
    1. Getting surrounding balls within the specified range
    2. Sorting them by distance from the previous vertex
    3. Checking each ball for potential vertex formation
    4. Returning None if no valid vertex is found
    """
    # Get the balls that should not be a part of the new vertex
    edge_ndxs = edge_balls[:]

    # Get the balls not in the invalid balls that are within the range specified
    invalid_ndxs_set = set(invalid_ndxs)
    metric_start = time.perf_counter()
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs_set]
    POW_PRM_METRICS['candidate_gather'] += time.perf_counter() - metric_start
    # Sort the test balls to be in order by distance from the previous vert location
    if vn_1_loc is None:
        vn_1_loc = calc_com([locs[_] for _ in edge_ndxs])

    vn_1_loc_array = np.array(vn_1_loc)
    dists = [calc_dist(np.array(locs[_]), vn_1_loc_array) for _ in test_balls]
    metric_start = time.perf_counter()
    test_balls = [_ for x, _ in sorted(zip(dists, test_balls))]
    metric_start = time.perf_counter()
    # Instantiate the list for test vertices to be calculated later. This saves us from sorting the vertices balls twice
    test_verts = []
    # Go through the surrounding balls to look for vertices that have been found before and filter out edge balls
    for ball in test_balls:
        # If the ball is in the previous vertex move on
        if ball in vn_1:
            continue
        # Check if we need to check and if so check for the ball in the list
        if check_ndxs and ball not in group_ndxs:
            continue
        # If we have found the vertex before it is not the previous vertex return
        ball_ndxs = edge_ndxs + [ball]
        ball_ndxs.sort()
        # Get the vertex's index/insert index
        check_verts = [vert_ndxs[_] for _ in b_verts[ball_ndxs[0]]]
        # Take the potential index of the ball in group
        my_vert_ndx = bisect.bisect_left(check_verts, ball_ndxs)
        # If the index returned is larger than the list or the vertex at the index is not equal to the ball_ndxs were ok
        if my_vert_ndx < len(check_verts) and ball_ndxs == check_verts[my_vert_ndx]:
            return None, invalid_ndxs
        # Add the vertex indices to the test_vertices for calculation
        test_verts.append((ball_ndxs, ball))
    POW_PRM_METRICS['candidate_filter'] += time.perf_counter() - metric_start
    POW_PRM_METRICS['candidates'] += len(test_verts)
    # Go through each ball in the given test balls. Extremely optimized
    for vert in test_verts:

        # Add the vertex ball to the
        vert_balls, ball = vert
        metric_start = time.perf_counter()
        try:
            v_loc, vert_rad = calc_flat_vert(
                locs=[locs[_] for _ in vert_balls],
                rads=[rads[_] for _ in vert_balls],
                power=True
            )
        except RuntimeWarning:
            invalid_ndxs.append(ball)
            continue
        finally:
            POW_PRM_METRICS['calc_vert'] += time.perf_counter() - metric_start

        # Catch the none location case
        if v_loc is None:
            invalid_ndxs.append(ball)
            continue

        # Check if the vert is outside the box
        if box is not None and any([box[0][k] > v_loc[k] or v_loc[k] > box[1][k] for k in range(3)]):
            continue
        # Reject oversized Power vertices before verification.
        max_power = max_vert ** 2 - min(rads[_] for _ in vert_balls) ** 2
        if vert_rad >= max_power:
            invalid_ndxs.append(ball)
            continue

        # Reuse the edge-level Power verification arrays and skip only this
        # candidate's four defining balls. No candidate-level spatial query or
        # NumPy verification-array rebuild is required.
        metric_start = time.perf_counter()
        skip0 = surr_lookup.get(vert_balls[0], -1)
        skip1 = surr_lookup.get(vert_balls[1], -1)
        skip2 = surr_lookup.get(vert_balls[2], -1)
        skip3 = surr_lookup.get(vert_balls[3], -1)
        POW_PRM_METRICS['verify_arrays'] += time.perf_counter() - metric_start
        POW_PRM_METRICS['verify_balls'] += len(surr_locs)

        metric_start = time.perf_counter()
        valid = verify_pow_cached(np.asarray(v_loc), vert_rad, surr_locs, surr_rads, skip0, skip1, skip2, skip3)
        POW_PRM_METRICS['verify'] += time.perf_counter() - metric_start

        if valid:
            # Return the validated ball and the invalidated list
            return [{'balls': vert_balls, 'loc': v_loc, 'rad': vert_rad}, metrics], invalid_ndxs
        else:
            # Add the ball to the invalid balls list if it isn't verified
            invalid_ndxs.append(ball)
    # Return the non-vertex and invalid balls
    return None, invalid_ndxs


def find_site_aw(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, check_ndxs, surr_balls,
                 my_boxes, invalid_ndxs, vn_1, vn_1_loc, box=None, group_balls=None, metrics=None, printing=False, max_ball_rad=None):
    """
    Find the neighboring additively weighted Voronoi vertex for an edge.

    Candidate fourth balls are gathered from progressively larger grid
    neighborhoods around the three defining edge balls. Previously discovered
    vertices and balls excluded by group/interface constraints are removed before
    the four-sphere AW calculation is performed.

    When multiple mathematical solutions remain, candidates are separated by
    their position relative to the edge and ordered relative to the previous
    vertex. Candidate sites are then verified against only the spatial
    neighborhood capable of invalidating the AW vertex.

    Parameters
    ----------
    edge_balls : list of int
        Three defining ball indices for the current edge.
    locs : array-like
        Ball-center coordinates.
    rads : array-like
        Ball radii.
    b_verts : list of list
        Existing vertex indices associated with each ball.
    vert_ndxs : list
        Defining ball indices of previously discovered vertices.
    max_vert : float
        Maximum permitted AW vertex radius.
    mv_inc : float
        Current candidate-search distance.
    check_ndxs : bool
        Whether the fourth ball is constrained to ``group_balls``.
    surr_balls : optional
        Retained for compatibility with the other network search functions.
        AW verification no longer uses a global surrounding-ball list.
    my_boxes : list
        Grid cells containing the defining edge balls.
    invalid_ndxs : list
        Candidate fourth balls previously rejected for this edge.
    vn_1 : list
        Defining balls of the previous vertex.
    vn_1_loc : array-like
        Location of the previous vertex.
    box : list, optional
        Geometric bounds for accepted vertex locations.
    group_balls : collection, optional
        Required group for the candidate fourth ball.
    metrics : dict, optional
        Performance metrics populated when profiling is enabled.
    printing : bool, optional
        Enable detailed diagnostic printing.
    max_ball_rad : float, optional
        Maximum ball radius used to bound local verification.

    Returns
    -------
    tuple
        ``(vertex, invalid_ndxs)`` where ``vertex`` is the verified neighboring
        vertex or ``None``.
    """

    # Get the balls that should not be a part of the new vertex
    edge_ndxs = edge_balls[:]

    # Get the balls not in the invalid balls that are within the range specified
    invalid_ndxs_set = set(invalid_ndxs)
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs_set]

    # Instantiate the list for test vertices to be calculated later. This saves us from sorting the vertices balls twice
    new_test_balls = []
    # Go through the surrounding balls to look for vertices that have been found before and filter out edge balls
    for ball in test_balls:
        # If the ball is in the previous vertex move on
        if ball in vn_1:
            continue
        # Check if we need to check and if so check for the ball in the list
        if check_ndxs and ball not in group_balls:
            continue
        # If we have found the vertex before it is not the previous vertex return
        ball_ndxs = edge_ndxs + [ball]
        ball_ndxs.sort()
        # Get the vertices for the first ball. All balls will contain the vertex so only one ball needs to be checked
        check_verts = [vert_ndxs[_] for _ in b_verts[ball_ndxs[0]]]
        # Use the ndx_search function to quickly search the list of sorted vertices
        my_vert_ndx = bisect.bisect_left(check_verts, ball_ndxs)
        # If the index returned is larger than the list or the vertex at the index is not equal to the ball_ndxs were ok
        if my_vert_ndx < len(check_verts) and ball_ndxs == check_verts[my_vert_ndx]:
            return KNOWN_EDGE, invalid_ndxs
        # Add the vertex indices to the test_vertices for calculation
        new_test_balls.append(ball)

    # Instantiate the calculated vertices list
    calc_verts = []
    # Go through each ball in the given test balls. Extremely optimized
    for ball in new_test_balls:

        # Combine the new ball with the edge balls and sort
        vert_balls = edge_balls + [ball]
        vert_balls.sort()
        # Calculate the Voronoi vertex values
        v_loc, v_rad, v_loc2, v_rad2 = calc_vert([locs[_] for _ in vert_balls], [rads[_] for _ in vert_balls])

        min_allowed_rad = -min(rads[_] for _ in vert_balls)

        if v_loc is None or v_rad < min_allowed_rad or v_rad > max_vert:
            continue

        # Delete the second location for the vertex if it is too large
        if v_rad2 is not None and (v_rad2 < min_allowed_rad or v_rad2 > max_vert):
            v_loc2, v_rad2 = None, None

        # Check if the vert is outside the box
        if box is not None and any(box[0][k] > v_loc[k] or v_loc[k] > box[1][k] for k in range(3)):
            continue

        # Add the vertex to the list of calculated vertices
        calc_verts.append({'balls': vert_balls, 'loc': np.array(v_loc), 'rad': v_rad, 'loc2': v_loc2, 'rad2': v_rad2})

    # If no vertices survived return
    if  not calc_verts:
        return None, invalid_ndxs
    # If there is only one vertex left, no need to sort. Just verify it
    if len(calc_verts) == 1:
        return choose_vert(calc_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, max_ball_rad=max_ball_rad)[0], invalid_ndxs

    # Instantiate the left and right vertex lists
    left_verts, right_verts = [], []
    # Get the centers of the edge balls
    c0, c1, c2 = [locs[_] for _ in edge_ndxs]
    # Get the center of the inscribed circle
    edge_center, edge_radius = calc_circ(locs[edge_ndxs[0]], locs[edge_ndxs[1]], locs[edge_ndxs[2]],
                                         rads[edge_ndxs[0]], rads[edge_ndxs[1]], rads[edge_ndxs[2]])

    # Calculate the edge normal  direction - take cross product of vector centers of edge balls - a0 a1 X a1, a2
    edge_direction = np.cross(c0 - c1, c0 - c2)
    edge_normal = edge_direction / np.linalg.norm(edge_direction)

    # Calculate the projection of the previous vertex onto the edge normal (value) or edge_normal dot prev vert center
    pv_dist = np.dot(edge_normal, edge_center - vn_1_loc)
    # Go through the calculated vertices made by the edge balls and the surrounding balls - filtering process
    for vert in calc_verts:
        # Get the vertex's projected distance
        vert_proj_dist = np.dot(edge_normal, edge_center - vert['loc'])
        # Calculate the distance to the previous vertex and assign it as a value in the vertex dictionary
        vert['d2pv'] = abs(pv_dist - vert_proj_dist)

        # If the other balls projection (value1) is less than the previous vertex's projection (value)
        if pv_dist < vert_proj_dist:
            # Add the vertex to the list of filtered vertices
            left_verts.append(vert)
        else:
            # Add the vertex to the list of filtered vertices
            right_verts.append(vert)

        vert['d2pv2'] = None
        if vert['loc2'] is not None:
            vert_proj_dist = np.dot(edge_normal, edge_center - vert['loc2'])
            # Treat the secondary solution as an independent candidate while preserving
            # the primary location as its alternate doublet solution.
            flipped_vert = {'balls': vert['balls'], 'loc': vert['loc2'], 'rad': vert['rad2'], 'd2pv': abs(pv_dist - vert_proj_dist), 'loc2': vert['loc'], 'rad2': vert['rad']}
            # If the other balls projection (value1) is less than the previous vertex's projection (value)
            if pv_dist < vert_proj_dist:
                # Add the vertex to the list of filtered vertices
                left_verts.append(flipped_vert)
            else:
                # Add the vertex to the list of filtered vertices
                right_verts.append(flipped_vert)

    # Sort the filtered vertices by distance to the previous vertex
    left_verts.sort(key=lambda my_vert: my_vert['d2pv'])
    right_verts.sort(key=lambda my_vert: my_vert['d2pv'])

    # Set up the left neighbor and the right neighbor variables for assignment
    left_neighbor, right_neighbor = None, None
    # If all vertices lie on the left side of the previous vertex
    if len(right_verts) == 0:
        # Get the leftmost vertex and the rightmost vertex
        vl, vr = left_verts[-1]['loc'], vn_1_loc
        # Counter variable
        i = 0
        # Loop through the vertices looking for the left and right neighbor
        while (left_neighbor is None or right_neighbor is None) and i < len(left_verts) - 1:
            # Grab the current vertex in the loop
            vi = left_verts[i]
            # Calculate the determinant of the vertex and the left most and right most vertices
            my_det = np.linalg.det([vl, vr, vi['loc']])
            # If the edge is straight, verify/return the leftmost vertex on the right
            if my_det == 0:
                # Verification
                return choose_vert(left_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, max_ball_rad=max_ball_rad)[0], invalid_ndxs
            # If the vertex falls in the lower hull it is the left neighbor
            elif my_det > 0 and left_neighbor is None:
                left_neighbor = vi
            # If the vertex falls in the upper hull it is the right neighbor
            elif my_det < 0 and right_neighbor is None:
                right_neighbor = vi
            # Increment the counter
            i += 1
        if left_neighbor is None:
            left_neighbor = left_verts[-1]
        elif right_neighbor is None:
            right_neighbor = left_verts[-1]
    # If all vertices lie on the right side of the previous vertex
    elif len(left_verts) == 0:
        # Get the leftmost vertex and the rightmost vertex
        vr, vl = right_verts[-1]['loc'], vn_1_loc
        # Counter variable
        i = 0
        # Loop through the vertices looking for the left and right neighbor
        while (left_neighbor is None or right_neighbor is None) and i < len(right_verts) - 1:
            # Grab the current vertex in the loop
            vi = right_verts[i]
            # Calculate the determinant of the vertex and the left most and right most vertices
            my_det = np.linalg.det([vl, vr, vi['loc']])
            # If the edge is straight, verify/return the leftmost vertex on the right
            if my_det == 0:
                # Verification
                return choose_vert(right_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, max_ball_rad=max_ball_rad)[0], invalid_ndxs
            # If the vertex falls in the upper hull it is the left neighbor
            elif my_det < 0 and left_neighbor is None:
                left_neighbor = vi
            # If the vertex falls in the lower hull it is the right neighbor
            elif my_det > 0 and right_neighbor is None:
                right_neighbor = vi
            # Increment the counter
            i += 1

        if left_neighbor is None:
            left_neighbor = right_verts[-1]
        elif right_neighbor is None:
            right_neighbor = right_verts[-1]
    # If there are vertices on either side
    else:
        # Find the left most and right most vertices
        vl, vr = left_verts[-1], right_verts[-1]
        vert_det = np.linalg.det([vl['loc'], vr['loc'], vn_1_loc])
        # Assign the left and right neighbor variables
        left_neighbor, right_neighbor = None, None
        # Counter variable
        i = 0
        # Go through the vertices on the left og the vertex
        while left_neighbor is None and i < len(left_verts):
            # Get the current vertex in the loop
            vi = left_verts[i]
            # Calculate the determinant of the left most, right most and current vertex
            my_det = np.linalg.det([vl['loc'], vr['loc'], vi['loc']])
            # If they share a sign, we have found the vertex
            if my_det <= 0 and vert_det <= 0 or my_det >= 0 and vert_det >= 0:
                left_neighbor = vi
            # Increment the counter
            i += 1
        # Reset the counter variable
        i = 0
        # Go through the vertices on the right of the previous vertex
        while right_neighbor is None and i < len(right_verts):
            # Get the current vertex in the loop
            vi = right_verts[i]
            # Calculate the determinant of the left most, right most and current vertex
            my_det = np.linalg.det([vl['loc'], vr['loc'], vi['loc']])
            # If they share a sign, we have found the vertex
            if my_det <= 0 and vert_det <= 0 or my_det >= 0 and vert_det >= 0:
                right_neighbor = vi
            # Increment the counter
            i += 1

    # Check the left neighbor vertex
    if left_neighbor is not None:
        my_vert, extra_ball = choose_vert(left_neighbor, edge_ndxs, surr_balls, locs, rads, metrics, max_ball_rad=max_ball_rad)

        if my_vert is not None:
            return my_vert, invalid_ndxs
        invalid_ndxs.append(extra_ball)
    # Check the right neighbor vertex
    if right_neighbor is not None:
        my_vert, extra_ball = choose_vert(right_neighbor, edge_ndxs, surr_balls, locs, rads, metrics, max_ball_rad=max_ball_rad)
        if my_vert is not None:
            return my_vert, invalid_ndxs
        invalid_ndxs.append(extra_ball)
    return None, invalid_ndxs


def verify_aw_local(loc, rad, vert_balls, b_locs, b_rads, max_ball_rad):
    """
    Verify an AW vertex using only balls that can geometrically invalidate it.

    For a candidate AW vertex with radius ``rad``, a ball can invalidate the
    vertex only if its center lies within ``rad + ball_radius`` of the vertex.
    Searching to ``rad + max_ball_rad`` therefore provides a conservative
    spatial bound that contains every possible invalidating ball while avoiding
    a system-wide verification search.

    Parameters
    ----------
    loc : array-like
        Candidate vertex location.
    rad : float
        Candidate AW vertex radius.
    vert_balls : collection of int
        Four balls defining the candidate vertex.
    b_locs : array-like
        Locations of all balls.
    b_rads : array-like
        Radii of all balls.
    max_ball_rad : float
        Maximum ball radius in the system.

    Returns
    -------
    bool
        True when no nearby non-defining ball invalidates the candidate.
    """
    verify_dist = max(0.0, rad + max_ball_rad)
    verify_box = box_search(loc)

    if verify_box is None:
        return False

    check_balls = [_ for _ in get_balls([verify_box], dist=verify_dist) if _ not in vert_balls]
    test_locs = np.asarray([b_locs[_] for _ in check_balls])
    test_rads = np.asarray([b_rads[_] for _ in check_balls])

    return verify_aw(np.asarray(loc), rad, test_locs, test_rads)


def choose_vert(my_vert, edge_ndxs, test_balls, b_locs, b_rads, metrics, max_ball_rad=None):
    """
    Verify a candidate AW neighbor and resolve optional doublet solutions.

    The primary candidate is verified first using ``verify_aw_local``. If the
    candidate has a secondary geometric solution, that location is independently
    verified. If the primary fails but the secondary succeeds, the secondary
    solution becomes the returned vertex.

    Parameters
    ----------
    my_vert : dict
        Candidate vertex containing ``balls``, ``loc``, ``rad`` and optional
        ``loc2``/``rad2`` values.
    edge_ndxs : list of int
        Three balls defining the current edge.
    test_balls : optional
        Retained for compatibility with the previous verification API. AW local
        verification does not use this argument.
    b_locs : array-like
        Ball locations for the full system.
    b_rads : array-like
        Ball radii for the full system.
    metrics : dict, optional
        Performance metrics populated when profiling is enabled.
    max_ball_rad : float, optional
        Maximum system ball radius.

    Returns
    -------
    tuple
        ``([vertex, metrics], None)`` when a valid candidate is found, otherwise
        ``(None, extra_ball)`` where ``extra_ball`` is the rejected fourth ball.
    """

    if max_ball_rad is None:
        max_ball_rad = max(b_rads)

    # Verify the primary vertex location
    primary_valid = verify_aw_local(
        my_vert['loc'], my_vert['rad'], my_vert['balls'],
        b_locs, b_rads, max_ball_rad
    )

    if primary_valid:
        # If a secondary/doublet location exists, independently verify it
        if my_vert['loc2'] is not None:
            secondary_valid = verify_aw_local(
                my_vert['loc2'], my_vert['rad2'], my_vert['balls'],
                b_locs, b_rads, max_ball_rad
            )

            if not secondary_valid:
                my_vert['loc2'], my_vert['rad2'] = None, None

        return [my_vert, metrics], None

    # If the primary location failed, try the secondary location
    if my_vert['loc2'] is not None:
        secondary_valid = verify_aw_local(
            my_vert['loc2'], my_vert['rad2'], my_vert['balls'],
            b_locs, b_rads, max_ball_rad
        )

        if secondary_valid:
            my_vert = {
                'balls': my_vert['balls'],
                'loc': my_vert['loc2'],
                'rad': my_vert['rad2'],
                'loc2': None,
                'rad2': None
            }

            return [my_vert, metrics], None

    # Neither location is valid
    extra_ball = [_ for _ in my_vert['balls'] if _ not in edge_ndxs][0]
    return None, extra_ball