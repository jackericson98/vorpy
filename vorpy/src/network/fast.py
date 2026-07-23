import bisect
import numpy as np
import time
import warnings
from vorpy.src.calculations import calc_flat_vert
from vorpy.src.calculations import calc_vert
from vorpy.src.calculations import verify_aw
from vorpy.src.calculations import verify_pow
from vorpy.src.calculations import verify_prm
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import sort_lists
from vorpy.src.calculations import calc_circ

warnings.filterwarnings("error", category=RuntimeWarning)


TARGET_VERTS = {
    (27, 334, 448, 582),
    (27, 334, 582, 583),
    (27, 448, 582, 583),
}

TARGET_EDGES = {
    (27, 334, 448),
    (27, 334, 583),
    (27, 448, 583),
}


def debug_target_edge(edge_balls):
    return tuple(sorted(edge_balls)) in TARGET_EDGES


def debug_target_vert(vert_balls):
    return tuple(sorted(vert_balls)) in TARGET_VERTS


def find_site_container(edge_balls, locs, rads, b_verts, vert_ndxs,
                        max_vert, net_type, box=None, vn_1=None, vn_1_loc=None,
                        group_ndxs=None, metrics=None, printing=False):
    """
    Searches for a valid vertex site by iteratively expanding the search area around edge balls.

    This function implements a progressive search strategy that starts with a small area around
    the edge balls and gradually increases the search radius until either a valid vertex is found
    or the maximum search distance is reached. The search can be constrained to specific groups
    of balls and supports different network types (aw, pow, prm).

    Parameters
    ----------
    edge_balls : list of int
        Indices of the balls forming the edge
    locs : list of numpy.ndarray
        List of ball locations in 3D space
    rads : list of float
        List of ball radii
    b_verts : list of list
        List of vertices associated with each ball
    vert_ndxs : list of int
        List of vertex indices
    max_vert : float
        Maximum search distance for finding vertices
    net_type : str
        Type of network ('aw', 'pow', or 'prm')
    box : list, optional
        Bounding box for spatial partitioning
    vn_1 : list of int, optional
        List of balls to consider for vertex formation
    vn_1_loc : numpy.ndarray, optional
        Location of the vn_1 balls' center of mass
    group_ndxs : list of int, optional
        Indices of balls belonging to a specific group
    metrics : dict, optional
        Dictionary for storing performance metrics
    printing : bool, optional
        Flag to enable progress printing

    Returns
    -------
    tuple
        A tuple containing:
        - vert : list or None
            List of balls forming the vertex if found, None otherwise
        - invalid_ndxs : list of int
            List of invalid vertex indices encountered during search
    """
    # Set up the vert and invalid indices parameters
    invalid_ndxs, vert = [], None

    # If no vn_1 is provided set it to the edge_balls
    if vn_1 is None:
        vn_1 = edge_balls

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

    # Find the 3 boxes the edge balls are in
    my_boxes = [box_search(loc=locs[edge_balls[_]]) for _ in range(3)]
    # Gather the surrounding balls or the entire list of balls we could be comparing to
    surr_balls = get_balls(cells=my_boxes, dist=max_vert)
    # Calculate the edge balls' com
    edge_com = calc_com([locs[_] for _ in edge_balls])
    # Get the distance from the edge com for each of the surr balls
    dists = [calc_dist(edge_com, locs[_]) for _ in surr_balls]
    # Sort the distances and surr balls
    dists, surr_balls = sort_lists(dists, surr_balls)
    # Se the initial vert size
    mv_inc = 0.45
    # Look for the vert and keep increasing box size until the vert is found
    while vert is None and mv_inc < max_vert:
        # Search for the vertx in the current range
        if net_type == 'aw':
            vert, invalid_ndxs = find_site_aw(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                                              required_group is not None, surr_balls, my_boxes, invalid_ndxs, vn_1,
                                              vn_1_loc, box=box, group_balls=required_group, metrics=metrics,
                                              printing=printing)
        elif net_type == 'pow':
            vert, invalid_ndxs = find_site_pow(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                                               required_group is not None, surr_balls, my_boxes, invalid_ndxs, vn_1,
                                               box, vn_1_loc, group_ndxs=required_group, metrics=metrics)
        elif net_type == 'prm':
            vert, invalid_ndxs = find_site_del(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                                               required_group is not None, surr_balls, my_boxes, invalid_ndxs, vn_1,
                                               box, vn_1_loc, group_ndxs=required_group, metrics=metrics)
        # If a vertex is found exit the loop
        if vert is not None:
            break
        # Increment the range for the search
        mv_inc *= 10
    # Return the vertex if found
    return vert


def find_site_del(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, check_ndxs, surr_balls,
                  my_boxes, invalid_ndxs, vn_1, box=None, vn_1_loc=None, group_ndxs=None, metrics=None):
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
    # Get the balls that should not ba a part of the new vertex
    edge_ndxs = edge_balls[:]

    # Time printing metrics <-- Delete later
    start = time.perf_counter()

    # Time printing metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()
    # Get the balls not in the invalid balls that are within the range specified
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs]
    # Sort the test balls to be in order by distance from the previous vert location
    if vn_1_loc is None:
        vn_1_loc = calc_com([locs[_] for _ in edge_ndxs])

    dists = [calc_dist(np.array(locs[_]), np.array(vn_1_loc)) for _ in test_balls]
    test_balls = [_ for x, _ in sorted(zip(dists, test_balls))]

    # Gather balls metrics <-- Delete later
    if metrics is not None:
        metrics['gather_balls'] += time.perf_counter() - start
        start = time.perf_counter()

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
    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start

    # Go through each ball in the given test balls. Extremely optimized
    for i, vert in enumerate(test_verts):

        # Add the vertex ball to the
        vert_balls, ball = vert
        # Calculate the 181L vertex values
        start = time.perf_counter()
        v_loc, vert_rad = calc_flat_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls], power=False)

        # Record the calculate vertex metrics
        if metrics is not None:
            metrics['calc_vert'] += time.perf_counter() - start

        # Catch the none location case
        if v_loc is None:
            invalid_ndxs.append(ball)
            continue
        # Check if the vert is outside the box
        if box is not None and any([box[0][k] > v_loc[k] or v_loc[k] > box[1][k] for k in range(3)]):
            continue

        # Restart ste start time to only record verify site time to the verify site metrics
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_balls = [_ for _ in surr_balls if _ not in vert_balls]
        # Get the locations from the test balls
        test_locs = np.array([locs[_] for _ in filtered_test_balls])
        # Compare the vertex to the maximum allowed vertex and verify it
        if vert_rad < max_vert and verify_prm(loc=np.array(v_loc), rad=vert_rad, test_locs=test_locs):
            # Add the time for verification to the verify_site metrics
            if metrics is not None:
                metrics['verify_site'] += time.perf_counter() - start
            # Return the validated ball and the invalidated ist
            return [{'balls': vert_balls, 'loc': v_loc, 'rad': vert_rad}, metrics], invalid_ndxs
        else:
            # Add the ball to the invalid balls list if it isn't verified
            invalid_ndxs.append(ball)
    # Return the non-vertex and invalid balls
    return None, invalid_ndxs


def find_site_pow(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, check_ndxs, surr_balls,
                  my_boxes, invalid_ndxs, vn_1, box=None, vn_1_loc=None, group_ndxs=None, metrics=None):
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
    # Get the balls that should not ba a part of the new vertex
    edge_ndxs = edge_balls[:]

    # Time printing metrics <-- Delete later
    start = time.perf_counter()

    # Time printing metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()
    # Get the balls not in the invalid balls that are within the range specified
    invalid_ndxs_set = set(invalid_ndxs)
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs_set]
    # Sort the test balls to be in order by distance from the previous vert location
    if vn_1_loc is None:
        vn_1_loc = calc_com([locs[_] for _ in edge_ndxs])

    dists = [calc_dist(np.array(locs[_]), np.array(vn_1_loc)) for _ in test_balls]
    test_balls = [_ for x, _ in sorted(zip(dists, test_balls))]

    # Gather balls metrics <-- Delete later
    if metrics is not None:
        metrics['gather_balls'] += time.perf_counter() - start
        start = time.perf_counter()

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
    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start

    # Go through each ball in the given test balls. Extremely optimized
    for i, vert in enumerate(test_verts):

        # Add the vertex ball to the
        vert_balls, ball = vert
        # Calculate the 181L vertex values
        start = time.perf_counter()
        try:
            v_loc, vert_rad = calc_flat_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls], power=True)
        except RuntimeWarning:
            invalid_ndxs.append(ball)
            continue
        # Record the calculate vertex metrics
        if metrics is not None:
            metrics['calc_vert'] += time.perf_counter() - start

        # Catch the none location case
        if v_loc is None:
            invalid_ndxs.append(ball)
            continue

        # Check if the vert is outside the box
        if box is not None and any([box[0][k] > v_loc[k] or v_loc[k] > box[1][k] for k in range(3)]):
            continue

        # Restart ste start time to only record verify site time to the verify site metrics
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_balls = [_ for _ in surr_balls if _ not in vert_balls]
        # Get the locations from the test balls
        test_locs = np.array([locs[_] for _ in filtered_test_balls])
        test_rads = np.array([rads[_] for _ in filtered_test_balls])
        # Compare the vertex to the maximum allowed vertex and verify it
        if (vert_rad < max_vert ** 2 - min([rads[_] for _ in vert_balls]) ** 2 and
                verify_pow(loc=np.array(v_loc), rad=vert_rad, test_locs=test_locs, test_rads=test_rads)):
            # Add the time for verification to the verify_site metrics
            if metrics is not None:
                metrics['verify_site'] += time.perf_counter() - start
            # Return the validated ball and the invalidated ist
            return [{'balls': vert_balls, 'loc': v_loc, 'rad': vert_rad}, metrics], invalid_ndxs
        else:
            # Add the ball to the invalid balls list if it isn't verified
            invalid_ndxs.append(ball)
    # Return the non-vertex and invalid balls
    return None, invalid_ndxs


def find_site_aw(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, check_ndxs, surr_balls,
                 my_boxes, invalid_ndxs, vn_1, vn_1_loc, box=None, group_balls=None, metrics=None, printing=False):
    """
    Finds a new vertex in an AW (Additive Weighted) network by searching for valid ball combinations.

    This function searches for a valid vertex that can be formed by combining the edge balls with
    a new ball from the surrounding region. It verifies that the potential vertex:
    - Is not already part of the previous vertex
    - Is not in the list of invalid balls
    - Has not been previously found
    - Meets any group membership criteria if specified

    Parameters
    ----------
    edge_balls : list of int
        List of ball indices that form the edge
    locs : list of numpy.ndarray
        List of ball locations in 3D space
    rads : list of float
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
    vn_1_loc : numpy.ndarray
        Location of previous vertex
    box : list, optional
        Search box boundaries
    group_balls : list, optional
        List of ball indices in the group
    metrics : dict, optional
        Dictionary for storing performance metrics
    printing : bool, optional
        Flag to enable progress printing

    Returns
    -------
    tuple
        A tuple containing:
        - vert : dict or None
            Dictionary containing vertex information if found, None otherwise
        - invalid_ndxs : list
            Updated list of invalid ball indices

    Notes
    -----
    - The function uses a progressive search strategy that expands outward from the edge balls
    - Performance metrics can be tracked if a metrics dictionary is provided
    - Group membership checks can be skipped by setting check_ndxs to False
    - The search can be constrained to a specific region using the box parameter
    """
    # Get the balls that should not ba a part of the new vertex
    edge_ndxs = edge_balls[:]
    # Time printing metrics <-- Delete later
    start = time.perf_counter()

    # Box search metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()

    # Get the balls not in the invalid balls that are within the range specified
    invalid_ndxs_set = set(invalid_ndxs)
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs_set]

    # Gather balls metrics <-- Delete later
    if metrics is not None:
        metrics['gather_balls'] += time.perf_counter() - start
        start - time.perf_counter()

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
            return None, invalid_ndxs
        # Add the vertex indices to the test_vertices for calculation
        new_test_balls.append(ball)

    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
        start = time.perf_counter()

    # Instantiate the calculated vertices list
    calc_verts = []
    # Go through each ball in the given test balls. Extremely optimized
    for i, ball in enumerate(new_test_balls):

        # Combine the new ball with the edge balls and sort
        vert_balls = edge_balls + [ball]
        vert_balls.sort()
        # Calculate the Voronoi vertex values
        v_loc, v_rad, v_loc2, v_rad2 = calc_vert([locs[_] for _ in vert_balls], [rads[_] for _ in vert_balls])
        # Catch the none location and the too large vertex cases
        if v_loc is None or v_rad > max_vert:
            continue

        # Delete the second location for the vertex if it is too large
        if v_rad2 is not None and v_rad2 > max_vert:
            v_loc2, v_rad2 = None, None

        # Check if the vert is outside the box
        if box is not None and any([box[0][k] > v_loc[k] > box[1][k] for k in range(3)]):
            continue

        # Add the vertex to the list of calculated vertices
        calc_verts.append({'balls': vert_balls, 'loc': np.array(v_loc), 'rad': v_rad, 'loc2': v_loc2, 'rad2': v_rad2})

    # Calculate vertices metrics <-- Delete later
    if metrics is not None:
        metrics['calc_vert'] += time.perf_counter() - start
        start = time.perf_counter()

    # If no vertices survived return
    if len(calc_verts) == 0:
        return None, invalid_ndxs
    # If there is only one vertex left, no need to sort. Just verify it
    elif len(calc_verts) == 1:
        return choose_vert(calc_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, start)[0], invalid_ndxs

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
                return choose_vert(left_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, start)[0], invalid_ndxs
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
                return choose_vert(right_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, start)[0], invalid_ndxs
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
        my_vert, extra_ball = choose_vert(left_neighbor, edge_ndxs, surr_balls, locs, rads, metrics, start)

        if my_vert is not None:
            return my_vert, invalid_ndxs
        invalid_ndxs.append(extra_ball)
    # Check the right neighbor vertex
    if right_neighbor is not None:
        my_vert, extra_ball = choose_vert(right_neighbor, edge_ndxs, surr_balls, locs, rads, metrics, start)
        if my_vert is not None:
            return my_vert, invalid_ndxs
        invalid_ndxs.append(extra_ball)
    return None, invalid_ndxs


def choose_vert(my_vert, edge_ndxs, test_balls, b_locs, b_rads, metrics, start):
    """
    Chooses a valid vertex from a set of potential vertices by verifying their validity.

    This function takes a potential vertex and verifies if it can form a valid vertex with the given edge balls.
    It checks both possible vertex locations (if they exist) and returns the first valid one found.
    If no valid vertex is found, it returns None along with an extra ball that was part of the invalid vertex.

    Parameters
    ----------
    my_vert : dict
        Dictionary containing vertex information including:
        - balls: List of ball indices forming the vertex
        - loc: Primary vertex location
        - rad: Primary vertex radius
        - loc2: Secondary vertex location (optional)
        - rad2: Secondary vertex radius (optional)
    edge_ndxs : list
        List of ball indices forming the edge
    test_balls : list
        List of ball indices to test against
    b_locs : list
        List of ball locations
    b_rads : list
        List of ball radii
    metrics : dict, optional
        Dictionary for storing performance metrics
    start : float
        Start time for performance measurement

    Returns
    -------
    tuple
        A tuple containing:
        - my_vert : dict or None
            The valid vertex if found, None otherwise
        - extra_ball : int or None
            An extra ball from the invalid vertex if no valid vertex is found
    """
    # Create the extra ball variable
    extra_ball = None
    # Get the balls surrounding the vertex, not including the vertex balls
    my_check_balls = [_ for _ in test_balls if _ not in my_vert['balls']]
    # Gather the locations and radii of the balls
    test_locs = np.array([b_locs[_] for _ in my_check_balls])
    test_rads = np.array([b_rads[_] for _ in my_check_balls])
    # Check the first location for the vertex
    if verify_aw(np.array(my_vert['loc']), my_vert['rad'], test_locs, test_rads):
        # Check the second location if it exists, if it is within the allowed size range and if it is verified
        if my_vert['rad2'] is None or not verify_aw(np.array(my_vert['loc2']), my_vert['rad2'], test_locs, test_rads):
            my_vert['loc2'], my_vert['rad2'] = None, None

        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        # Return what is left of the left vertex
        return [my_vert, metrics], extra_ball

    # If the first site is unverified try the other vertex site
    elif my_vert['loc2'] is not None and verify_aw(loc=np.array(my_vert['loc2']), rad=my_vert['rad2'],
                                                     test_locs=test_locs, test_rads=test_rads):
        # Reset the left_vert variable with the other location and return it
        my_vert = {'balls': my_vert['balls'], 'loc': my_vert['loc2'], 'rad': my_vert['rad2'], 'loc2': None,
                      'rad2': None}
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        return [my_vert, metrics], extra_ball
    # We still need to return invalid balls if they are not included
    return None, [_ for _ in my_vert['balls'] if _ not in edge_ndxs][0]

