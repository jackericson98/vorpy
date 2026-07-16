from vorpy.src.calculations import calc_vert
from vorpy.src.network.find_v0 import find_v0
from vorpy.src.network.fast import find_site_container
from vorpy.src.calculations import get_time
from vorpy.src.calculations import ndx_search
import time
from numpy import sqrt


# Find network function. Keeps searching the network until all verts are found
def find_verts(locs, rads, max_vert, net_type, check_ndxs, b0=None, my_group=None,
               b_verts=None, vert_ndxs=None, vlocs=None, vrads=None, vloc2s=None, vrad2s=None, start_time=0,
               print_metrics=False, box=None, vert_box=None, group_box=None, tot_ball_num=None, printing=False,
               start_vert=0, split=False):
    """
    Finds vertices in a network by searching through combinations of balls and verifying their validity.

    This function serves as the main entry point for vertex finding in the network. It can:
    1. Find the initial vertex (v0) using various strategies
    2. Search for additional vertices based on existing ones
    3. Handle different network types and constraints

    Parameters
    ----------
    locs : list of numpy.ndarray
        List of ball locations in 3D space
    rads : list of float
        List of ball radii
    max_vert : float
        Maximum distance to search for vertices
    net_type : str
        Type of network ('aw', 'pow', or 'prm')
    check_ndxs : list
        List of ball indices to check for vertices
    b0 : int, optional
        Index of initial ball to start search from
    my_group : list, optional
        List of ball indices in the group to constrain search
    b_verts : list of lists, optional
        List mapping ball indices to their vertices
    vert_ndxs : list, optional
        List of existing vertex indices
    vlocs : list, optional
        List of vertex locations
    vrads : list, optional
        List of vertex radii
    vloc2s : list, optional
        List of secondary vertex locations (for doublets)
    vrad2s : list, optional
        List of secondary vertex radii (for doublets)
    start_time : float, optional
        Start time for performance measurement
    print_metrics : bool, optional
        Whether to print performance metrics
    box : list, optional
        Overall search box boundaries
    vert_box : list, optional
        Bounding box for vertex search
    group_box : list, optional
        Bounding box for group search
    tot_ball_num : int, optional
        Total number of balls in the network
    printing : bool, optional
        Whether to print progress information
    start_vert : int, optional
        Starting vertex index
    split : bool, optional
        Whether to split the search

    Returns
    -------
    tuple or None
        If vertices are found, returns a tuple containing:
        - List of vertex indices
        - List of vertex locations
        - List of vertex radii
        - List of secondary vertex locations
        - List of secondary vertex radii
        - List of remaining unvisited balls
        - Dictionary mapping balls to their vertices
        Returns None if no valid vertices are found

    Notes
    -----
    - The function uses different strategies based on the number of balls in the group
    - For single ball groups, it directly finds v0
    - For groups of 4 balls, it calculates the vertex directly
    - For other cases, it uses a progressive search strategy
    - The function handles different network types and can be constrained to specific groups
    - Performance metrics can be tracked if enabled
    """
    # Metrics measuring < -- Deleting later
    metrics = None
    start = time.perf_counter()
    if print_metrics:
        metrics = {'ndx_search': 0, 'box_search': 0, 'gather_balls': 0, 'verify_site': 0, 'calc_vert': 0, 'other': 0}

    # Get the group balls from which to check vertices against
    if my_group is None or len(my_group) == len(locs):
        # Set the group balls to just the integers in up to the number of balls
        my_group = [_ for _ in range(len(locs))]
        # Calculate the rough number of vertices
        if tot_ball_num is None:
            tot_verts = 7 * len(locs)
    # If a group was provided make sure to get its indices
    elif my_group is not None:
        # Calculate the number of vertices
        if tot_ball_num is None:
            tot_verts = 7 * len(my_group) + int(60 * sqrt(len(my_group)))
    else:
        return

    if tot_ball_num is not None:
        tot_verts = 7 * tot_ball_num + int(60 * sqrt(tot_ball_num))
    if b_verts is None:
        b_verts = [[] for _ in range(len(locs))]
    # Find the first verified vertex
    if len(my_group) == 1:
        v0 = find_v0(locs=locs, rads=rads, b_verts=b_verts, max_vert=max_vert, net_type=net_type, b0=my_group[0],
                     group_ndxs=my_group, metrics=metrics, vert_ndxs=vert_ndxs, group_box=group_box)

    elif len(my_group) == 4:
        v0_loc, v0_rad, v0_loc2, v0_rad2 = calc_vert(locs=[locs[_] for _ in my_group],
                                                     rads=[rads[_] for _ in my_group])
        v0 = {'balls': my_group, 'loc': v0_loc, 'rad': v0_rad, 'loc2': v0_loc2, 'rad2': v0_rad2}
    else:
        v0 = find_v0(locs=locs, rads=rads, b_verts=b_verts, max_vert=max_vert, net_type=net_type, b0=b0,
                     group_ndxs=my_group, metrics=metrics, vert_ndxs=vert_ndxs, group_box=group_box, box=box)
        j = 1
        while v0 is None and j < len(check_ndxs):

            v0 = find_v0(locs=locs, rads=rads, b_verts=b_verts, max_vert=max_vert, net_type=net_type, b0=check_ndxs[j],
                         group_ndxs=my_group, metrics=metrics, vert_ndxs=vert_ndxs, group_box=group_box, box=box)

            j += 1
    # If no v0 is possible (e.g., a lone ball) return
    if v0 is None:
        return
    # Check if this is the first go around
    if vert_ndxs is None:
        for ball in v0['balls']:
            # noinspection PyTypeChecker
            b_verts[ball].append(0)
        vert_ndxs = [v0['balls']]
        vlocs = [v0['loc']]
        vrads = [v0['rad']]
        if 'loc2' in v0:
            vloc2s = [v0['loc2']]
            vrad2s = [v0['rad2']]
        else:
            vloc2s = [[None, None, None]]
            vrad2s = [None]
    else:
        for ball in v0['balls']:
            b_vert_ndxs = [vert_ndxs[_] for _ in b_verts[ball]]
            # noinspection PyTypeChecker
            b_verts[ball].insert(ndx_search(b_vert_ndxs, v0['balls']), len(vert_ndxs))
        vert_ndxs.append(v0['balls'])
        vlocs.append(v0['loc'])
        vrads.append(v0['rad'])
        if 'loc2' in v0:
            vloc2s.append(v0['loc2'])
            vrad2s.append(v0['rad2'])
        else:
            vloc2s.append([None, None, None])
            vrad2s.append(None)
    # Set up the vertex stack
    vert_stack = [v0]
    # While the verts stack is not empty
    while vert_stack:
        # Get the vertex from the bottom of the stack
        vert = vert_stack.pop()
        # Set up the edge stack
        e_stack = [[[vert['balls'][i], vert['balls'][(i + 1) % 4], vert['balls'][(i + 2) % 4]], vert] for i in range(4)]
        # While the edge stack is not empty
        while e_stack:
            # Get the percentage and print it
            percentage = min((len(vlocs) / tot_verts) * 100, 100)
            my_time = time.perf_counter() - start_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{:02d}:{:2.2f} - Process: finding vertices: {} verts - {:.2f} %"
                  .format(int(h), int(m), round(s, 2), len(vert_ndxs) + start_vert, percentage), end="")
            # Get the edge from the top of the stack
            edge_balls, vert = e_stack.pop()
            # Find the next site in the network
            vert_ndx_pr = find_site_container(edge_balls=edge_balls, locs=locs, rads=rads, b_verts=b_verts,
                                              vert_ndxs=vert_ndxs, max_vert=max_vert, net_type=net_type,
                                              vn_1=vert['balls'], box=box, vn_1_loc=vert['loc'],
                                              group_ndxs=my_group, metrics=metrics, printing=printing)
            # If the vertex is none continue
            if vert_ndx_pr is None:
                continue
            # Set the vertex and its index
            my_vert, metrics = vert_ndx_pr
            # # Check if there is a retaining box for the vertices and if the vertex is outside the box
            # if vert_box is not None and (any([vert_box[0][i] > my_vert['loc'][i] for i in range(3)]) or
            #                              any([vert_box[1][i] < my_vert['loc'][i] for i in range(3)])):
            #     continue
            if my_vert['loc'] is None:
                continue
            # print(my_vert['balls'], box, my_vert['loc'], [box[0][k] > my_vert['loc'][k] or my_vert['loc'][k] > box[1][k] for k in range(3)])
            if box is not None and any([box[0][k] > my_vert['loc'][k] or my_vert['loc'][k] > box[1][k] for k in range(3)]):
                continue
            if box is not None and 'loc2' in my_vert and my_vert['loc2'] is not None and any([box[0][k] > my_vert['loc2'][k] or my_vert['loc2'][k] > box[1][k] for k in range(3)]):
                my_vert['loc2'], my_vert['rad2'] = None, None
            # Add the vertex to the stack and the network
            vert_stack.append(my_vert)
            # Insert the vertices in order of increasing ball indices
            vert_ndxs.append(my_vert['balls'])
            vlocs.append(my_vert['loc'])
            vrads.append(my_vert['rad'])
            if 'loc2' in my_vert:
                vloc2s.append(my_vert['loc2'])
                vrad2s.append(my_vert['rad2'])
            else:
                vloc2s.append([None, None, None])
                vrad2s.append(None)
            # Remove the balls from the
            for ball in my_vert['balls']:
                # noinspection PyTypeChecker
                b_vert_ndxs = [vert_ndxs[_] for _ in b_verts[ball]]
                # noinspection PyTypeChecker
                b_verts[ball].insert(ndx_search(b_vert_ndxs, my_vert['balls']), len(vert_ndxs) - 1)
                if ball in check_ndxs:
                    check_ndxs.remove(ball)
    # Printing out metrics < --- Delete later
    # if print_metrics:
    #     metrics['total'] = time.perf_counter() - start
    #     metrics['other'] = metrics['total'] - (metrics['ndx_search'] + metrics['box_search'] + metrics['gather_balls'] +
    #                                            metrics['verify_site'] + metrics['calc_vert'])
    #     print('\n\nVertex Finding Time Metrics: \n'
    #           '  Index Search      = {:.3f} s\n'
    #           '  Box Search        = {:.3f} s\n'
    #           '  Get balls         = {:.3f} s\n'
    #           '  Verify Site       = {:.3f} s\n'
    #           '  Calculate Vertex  = {:.3f} s\n'
    #           '  Other             = {:.3f} s\n'
    #           '  Total             = {:.3f} s\n\n'
    #           .format(metrics['ndx_search'], metrics['box_search'], metrics['gather_balls'], metrics['verify_site'],
    #                   metrics['calc_vert'], metrics['other'], metrics['total']))
    # Return the values of the vertices
    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, check_ndxs, b_verts
