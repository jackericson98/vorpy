from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
from System.Network.verts.verify_site import verify_site
from System.sys_funcs.calcs.calcs import calc_dist, calc_com
from System.sys_funcs.calcs.sorting import box_search, get_balls
from System.sys_funcs.calcs.circle import calc_circ
import bisect
import numpy as np
import time


def find_site_container(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, net_type, vn_1=None, vn_1_loc=None,
                        group_ndxs=None, metrics=None, printing=False):
    """
    Cycles through larger and larger areas searching for
    """
    # Set up the vert and invalid indices parameters
    invalid_ndxs, vert = [], None

    # If no vn_1 is provided set it to the edge_balls
    if vn_1 is None:
        vn_1 = edge_balls

    # Check if the edge contains a group ball, to see if the next ball needs to be checked or not
    # Start with check balls as false if no group is defined
    check_ndxs = False
    if group_ndxs is not None:
        # If a group exists default to checking each ball
        check_ndxs = True
        # Go through the edge balls checking if they are in the group --> any vert found from another ball is included
        for ball in edge_balls:
            # Take the potential index of the ball in group
            my_index = bisect.bisect_left(group_ndxs, ball)
            # If the index is in the list check if the ball matches the index's element
            if my_index != len(group_ndxs) and group_ndxs[my_index] == ball:
                # If the element is found no need to check the balls and break the for loop
                check_ndxs = False
                break

    # Find the 3 boxes the edge balls are in
    my_boxes = [box_search(loc=locs[edge_balls[_]]) for _ in range(3)]
    # Gather the surrounding balls or the entire list of balls we could be comparing to
    surr_balls = get_balls(cells=my_boxes, dist=max_vert)
    # Se the initial vert size
    mv_inc = 0.45
    # Look for the vert and keep increasing box size until the vert is found
    while vert is None and mv_inc < max_vert:
        # Search for the vertx in the current range
        if net_type == 'aw':
            vert, invalid_ndxs = find_site_aw(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, net_type,
                                              check_ndxs, surr_balls, my_boxes, invalid_ndxs, vn_1, vn_1_loc,
                                              group_balls=group_ndxs, metrics=metrics, printing=printing)
        else:
            vert, invalid_ndxs = find_site_pd(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc,
                                              net_type, check_ndxs, surr_balls, my_boxes, invalid_ndxs, vn_1,
                                              vn_1_loc, group_ndxs=group_ndxs, metrics=metrics)
        # If a vertex is found exit the loop
        if vert is not None:
            break
        # Increment the range for the search
        mv_inc *= 10
    # Return the vertex if found
    return vert


def find_site_pd(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, net_type, check_ndxs, surr_balls,
                 my_boxes, invalid_ndxs, vn_1, vn_1_loc=None, group_ndxs=None, metrics=None):
    """
    Used a vertex and a combination of it's edge balls to find the connecting vertex
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
        if net_type == 'pow':
            vert_loc, vert_rad = calc_flat_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls], power=True)
        elif net_type == 'del':
            vert_loc, vert_rad = calc_flat_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls], power=False)
        else:
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls])

        # Record the calculate vertex metrics
        if metrics is not None:
            metrics['calc_vert'] += time.perf_counter() - start

        # Catch the none location case
        if vert_loc is None:
            invalid_ndxs.append(ball)
            continue

        # Restart ste start time to only record verify site time to the verify site metrics
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_balls = [_ for _ in surr_balls if _ not in vert_balls]
        # Get the locations from the test balls
        test_locs = np.array([locs[_] for _ in filtered_test_balls])
        test_rads = np.array([rads[_] for _ in filtered_test_balls])
        # Compare the vertex to the maximum allowed vertex and verify it
        if vert_rad < max_vert and verify_site(loc=np.array(vert_loc), rad=vert_rad, test_locs=test_locs,
                                               test_rads=test_rads, net_type=net_type):
            # Add the time for verification to the verify_site metrics
            if metrics is not None:
                metrics['verify_site'] += time.perf_counter() - start
            # Return the validated ball and the invalidated ist
            return [{'balls': vert_balls, 'loc': vert_loc, 'rad': vert_rad}, metrics], invalid_ndxs
        else:
            # Add the ball to the invalid balls list if it isn't verified
            invalid_ndxs.append(ball)
    # Return the non-vertex and invalid balls
    return None, invalid_ndxs


def find_site_aw(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, net_type, check_ndxs, surr_balls,
                 my_boxes, invalid_ndxs, vn_1, vn_1_loc, group_balls=None, metrics=None, printing=False):
    """
    Used a vertex and a combination of it's edge balls to find the connecting vertex
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
    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs]

    # Sort the test balls to be in order by distance from the previous vert location
    if net_type != 'aw' and vn_1_loc is not None:
        dists = [calc_dist(np.array(locs[_]), np.array(vn_1_loc)) for _ in test_balls]
        test_balls = [_ for x, _ in sorted(zip(dists, test_balls))]

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
        # Make sure the vertex values are defined
        v_loc, v_rad, v_loc2, v_rad2 = None, None, None, None
        # If the network type is power
        if net_type == 'pow':
            # Calculate the power vertex values
            v_loc, v_rad = calc_flat_vert([locs[_] for _ in vert_balls], [rads[_] for _ in vert_balls], True)
        # If the network type is Delaunay
        elif net_type == 'del':
            # Calculate the Delaunay vertex values
            v_loc, v_rad = calc_flat_vert([locs[_] for _ in vert_balls], [rads[_] for _ in vert_balls], False)
        # If the network type is Voronoi
        elif net_type == 'aw':
            # Calculate the Voronoi vertex values
            v_loc, v_rad, v_loc2, v_rad2 = calc_vert([locs[_] for _ in vert_balls], [rads[_] for _ in vert_balls])
        # Catch the none location and the too large vertex cases
        if v_loc is None or v_rad > max_vert:
            continue

        # Delete the second location for the vertex if it is too large
        if v_rad2 is not None and v_rad2 > max_vert:
            v_loc2, v_rad2 = None, None

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
        return choose_vert(calc_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, start, net_type)[0], invalid_ndxs

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
                return choose_vert(left_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, start, net_type)[0], invalid_ndxs
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
                return choose_vert(right_verts[0], edge_ndxs, surr_balls, locs, rads, metrics, start, net_type)[0], invalid_ndxs
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
        my_vert, extra_ball = choose_vert(left_neighbor, edge_ndxs, surr_balls, locs, rads, metrics, start, net_type)

        if my_vert is not None:
            return my_vert, invalid_ndxs
        invalid_ndxs.append(extra_ball)
    # Check the right neighbor vertex
    if right_neighbor is not None:
        my_vert, extra_ball = choose_vert(right_neighbor, edge_ndxs, surr_balls, locs, rads, metrics, start, net_type)
        if my_vert is not None:
            return my_vert, invalid_ndxs
        invalid_ndxs.append(extra_ball)
    return None, invalid_ndxs


def choose_vert(my_vert, edge_ndxs, test_balls, b_locs, b_rads, metrics, start, net_type):
    # Create the extra ball variable
    extra_ball = None
    # Get the balls surrounding the vertex, not including the vertex balls
    my_check_balls = [_ for _ in test_balls if _ not in my_vert['balls']]
    # Gather the locations and radii of the balls
    test_locs = np.array([b_locs[_] for _ in my_check_balls])
    test_rads = np.array([b_rads[_] for _ in my_check_balls])
    # Check the first location for the vertex
    if verify_site(np.array(my_vert['loc']), my_vert['rad'], test_locs, test_rads):
        # Check the second location if it exists, if it is within the allowed size range and if it is verified
        if my_vert['rad2'] is None or not verify_site(np.array(my_vert['loc2']), my_vert['rad2'], test_locs, test_rads):
            my_vert['loc2'], my_vert['rad2'] = None, None

        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        # Return what is left of the left vertex
        return [my_vert, metrics], extra_ball

    # If the first site is unverified try the other vertex site
    elif my_vert['loc2'] is not None and verify_site(loc=np.array(my_vert['loc2']), rad=my_vert['rad2'],
                                                     test_locs=test_locs, test_rads=test_rads, net_type=net_type):
        # Reset the left_vert variable with the other location and return it
        my_vert = {'balls': my_vert['balls'], 'loc': my_vert['loc2'], 'rad': my_vert['rad2'], 'loc2': None,
                      'rad2': None}
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        return [my_vert, metrics], extra_ball
    # We still need to return invalid balls if they are not included
    return None, [_ for _ in my_vert['balls'] if _ not in edge_ndxs][0]

