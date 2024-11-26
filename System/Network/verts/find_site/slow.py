from System.sys_funcs.calcs.sorting import box_search, get_balls, ndx_search
from System.sys_funcs.calcs.vert import verify_site
from System.sys_funcs.calcs.vert import calc_flat_vert, calc_vert
import bisect
import numpy as np
import time


def find_site_container_slow(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, net_type, group_ndxs=None, metrics=None,
                             printing=False):
    """
    Cycles through larger and larger areas searching for
    """
    # Set up the vert and invalid indices parameters
    invalid_ndxs, vert = [], None

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
        vert, invalid_ndxs = find_site(edge_balls=edge_balls, locs=locs, rads=rads, b_verts=b_verts,
                                       vert_ndxs=vert_ndxs, max_vert=max_vert, mv_inc=mv_inc, net_type=net_type,
                                       invalid_ndxs=invalid_ndxs, check_balls=check_ndxs, surr_balls=surr_balls,
                                       my_boxes=my_boxes, group_ndxs=group_ndxs, metrics=metrics)
        # If a vertex is found exit the loop
        if vert is not None:
            break
        # Increment the range for the search
        mv_inc *= 10
    # Return the vertex if found
    return vert


def find_site(edge_balls, locs, rads, b_verts, vert_ndxs, max_vert, mv_inc, net_type, invalid_ndxs=None,
              check_balls=True, surr_balls=None, vn_1=None, vn_1_loc=None, group_ndxs=None, metrics=None, my_boxes=None):
    """
    Used a vertex and a combination of it's edge balls to find the connecting vertex
    """
    if invalid_ndxs is None:
        invalid_ndxs = []
    # Get the balls that should not ba a part of the new vertex
    edge_ndxs = edge_balls[:]

    # If the previous vertex has been provided, add the other  to the not allowed balls
    vert_ball_ndxs = vn_1
    if vn_1 is None:
        vert_ball_ndxs = edge_ndxs

    # Time printing metrics <-- Delete later
    start = time.perf_counter()
    # Grab the balls we want to test against
    if my_boxes is None:
        my_boxes = [box_search(loc=locs[edge_balls[_]]) for _ in range(3)]
    # Time printing metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()

    test_balls = [_ for _ in get_balls(cells=my_boxes, dist=mv_inc) if _ not in invalid_ndxs]
    if surr_balls is not None:
        surr_balls = get_balls(cells=my_boxes, dist=max_vert)

    if metrics is not None:
        metrics['gather_balls'] += time.perf_counter() - start
    # First look for vertices that have been found before
    new_test_balls = []
    start = time.perf_counter()
    for ball in test_balls:
        # If the ball is in the previous vertex move on
        if ball in vert_ball_ndxs:
            continue
        # Check if we need to check and if so check for the ball in the list
        if check_balls and ball not in group_ndxs:
            continue
        # If we have found the vertex before it is not the previous vertex return
        ball_ndxs = edge_ndxs + [ball]
        ball_ndxs.sort()
        # Get the vertex's index/insert index
        if vert_ndxs is not None and len(vert_ndxs) > 0:
            check_verts = [vert_ndxs[_] for _ in b_verts[ball_ndxs[0]]]
            my_vert_ndx = ndx_search(check_verts, ball_ndxs)
            if my_vert_ndx < len(check_verts) and ball_ndxs == check_verts[my_vert_ndx]:
                return None, invalid_ndxs
        new_test_balls.append(ball)
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
    # Instantiate the vertex list and the size limit for vertices found
    verts = []
    new_start = time.perf_counter()
    # Go through each ball in the given test balls
    for i, ball in enumerate(new_test_balls):
        ball_new_star = time.perf_counter()
        # Create the vertex and calculate its value
        vert_balls = edge_balls + [ball]
        vert_balls.sort()
        vert_loc2, vert_rad2 = None, None
        # Calculate the 181L vertex values
        start = time.perf_counter()
        if net_type == 'pow':
            vert_loc, vert_rad = calc_flat_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls], power=True)
        elif net_type == 'prm':
            vert_loc, vert_rad = calc_flat_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls], power=False)
        else:
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(locs=[locs[_] for _ in vert_balls], rads=[rads[_] for _ in vert_balls])
        if metrics is not None:
            metrics['calc_vert'] += time.perf_counter() - start
        # Catch the none location case
        if vert_loc is None:
            invalid_ndxs.append([_ for _ in vert_balls if _ not in edge_ndxs])
            continue
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_balls = [_ for _ in surr_balls if _ not in vert_balls]
        test_locs = np.array([locs[_] for _ in filtered_test_balls])
        test_rads = np.array([rads[_] for _ in filtered_test_balls])
        if abs(vert_rad) < max_vert and verify_site(loc=np.array(vert_loc), rad=vert_rad, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            if len(verts) > 0 and verts[0]['rad'] < vert_rad:
                return [verts[0], metrics], invalid_ndxs
            verts.append({'balls': vert_balls, 'loc': vert_loc, 'rad': vert_rad, 'loc2': None, 'rad2': None})
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if vert_loc2 is not None and abs(vert_rad2) < max_vert and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
                verts[-1]['loc2'], verts[-1]['rad2'] = vert_loc2, vert_rad2
        # Check to see if the doublet's site is verified
        elif vert_loc2 is not None and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            verts.append({'balls': vert_balls, 'loc': vert_loc2, 'rad': vert_rad2, 'loc2': None, 'rad2': None})
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        else:
            invalid_ndxs.append([_ for _ in vert_balls if _ not in edge_ndxs])
    # If no verts have been found return
    if len(verts) == 0:
        return None, invalid_ndxs
    # If we find only 1 vertex, return it
    elif len(verts) == 1 or verts[0]['rad'] < verts[1]['rad']:
        return [verts[0], metrics], invalid_ndxs
    return [verts[1], metrics], invalid_ndxs
