import os
import numpy as np
import pandas as pd
import time
from vorpy.src.calculations import divide_box
from vorpy.src.calculations import global_vars
from vorpy.src.calculations import ndx_search
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import calc_dist
from vorpy.src.network.find_verts import find_verts
from vorpy.src.output.net import add_metrics


def split_net(sys, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None, my_group=None,
              print_actions=None, num_balls_sub_net=100, add_net_metrics=True, min_ball_split=30):
    """
    Sorts the balls in the main network by their locations to optimize spatial organization.

    This function sorts the balls in the network to improve spatial organization and subsequent operations.
    The sorting is based on the ball locations and is essential for efficient network construction and analysis.

    Parameters
    ----------
    sys : Network
        The main network system containing the balls to be sorted

    Notes
    -----
    - Sorting is performed on the ball locations to optimize spatial organization
    - This step is crucial for efficient network construction and analysis
    - The sorting is done in-place on the network's ball DataFrame
    """
    sys.net.sort_balls()
    # Calculate the group box
    group_box = sys.net.calc_box([sys.balls['loc'][_] for _ in my_group.group_ndxs],
                                 [sys.balls['rad'][_] for _ in my_group.group_ndxs], return_val=True, box_size=1.1)
    # Get the sub boxes
    sub_boxes = divide_box(group_box, round(len(my_group.group_ndxs) / num_balls_sub_net), c=1.2)
    print('num splits', len(sub_boxes))
    # Check for a max_vert that isn't defined
    if max_vert is None:
        max_vert = sys.net.settings['max_vert']

    # Sort the balls into their sub_boxes
    balls_lists = [[] for _ in range(len(sub_boxes))]
    ball_locs = sys.balls['loc']
    # Loop through the ball locations and sort the balls
    for ball in my_group.group_ndxs:
        loc = ball_locs[ball]
        # Loop through the sub boxes to find the placement of the ball
        for j, sub_box in enumerate(sub_boxes):
            if [sub_box[0][k] <= loc[k] <= sub_box[1][k] for k in range(3)] == [True, True, True]:
                balls_lists[j].append(ball)
    # If a list of balls is too small add it to another
    skip_boxes = []
    for i, balls_list in enumerate(balls_lists):
        # If no balls exist nothing to deal with
        if len(balls_list) == 0:
            skip_boxes.append(i)
            continue
        if len(balls_list) < min_ball_split:
            # Get the com of the balls to find the closes sub_box to add to
            balls_com = calc_com([ball_locs[_] for _ in balls_list])
            min_dist = np.inf
            closest_sub_box = None
            for j, sub_box in enumerate(sub_boxes):
                # Make sure we aren't adding to a sub_box scheduled for deletion
                if j in skip_boxes or j == i:
                    continue
                # Calculate the distance of the com of the sub_box from the balls_com
                my_dist = calc_dist(calc_com(sub_box), balls_com)
                # Replace the variables if they are closer
                if my_dist < min_dist:
                    closest_sub_box, min_dist = j, my_dist
            # Add the balls to the new sub_box
            balls_lists[closest_sub_box] += balls_list
            skip_boxes.append(i)
    for i, ball_list in enumerate(balls_lists):
        # Skip the boxes to be skipped
        if i in skip_boxes:
            continue
    # Instantiate the global variables
    global_vars(sys.net.sub_boxes, sys.net.box, sys.net.num_splits, sys.max_ball_rad, sys.net.sub_box_size)
    vert_ndxs, v_locs, v_rads, v_loc2s, v_rad2s, ball_nums, b_verts = None, None, None, None, None, None, None
    # Create the subnetworks
    for i, ball_list in enumerate(balls_lists):
        # Skip the boxes to be skipped
        if i in skip_boxes:
            continue
        # Get the balls we are tying to find vertices for
        check_balls = [_ for _ in ball_list if _ in my_group.group_ndxs]
        ball_nums = check_balls[:]
        # Find the initial vertices for the vertex group
        init_verts = find_verts(locs=sys.balls['loc'].to_numpy(), rads=sys.balls['rad'].to_numpy(),
                                max_vert=max_vert, net_type=net_type, check_ndxs=check_balls,
                                my_group=ball_nums, start_time=sys.net.start_time,
                                vert_box=sys.foam_box, group_box=sub_boxes[i], vert_ndxs=vert_ndxs, vlocs=v_locs,
                                vrads=v_rads, vloc2s=v_loc2s, vrad2s=v_rad2s, b_verts=b_verts,
                                tot_ball_num=len(my_group.group_ndxs))
        # Check to see if find_verts fails
        if init_verts is not None:
            vert_ndxs, v_locs, v_rads, v_loc2s, v_rad2s, ball_nums, b_verts = init_verts

        # Check for disconnects in the network
        while len(ball_nums) > 0:
            # Grab the initial ball for the next search
            b0 = ball_nums.pop()

            # Find verts again
            more_verts = find_verts(b0=b0, locs=sys.balls['loc'].to_numpy(), rads=sys.balls['rad'].to_numpy(),
                                    max_vert=max_vert, net_type=net_type, check_ndxs=ball_nums,
                                    my_group=check_balls, vert_ndxs=vert_ndxs, vlocs=v_locs, vrads=v_rads,
                                    vloc2s=v_loc2s, vrad2s=v_rad2s, start_time=sys.net.start_time,
                                    vert_box=sys.foam_box, b_verts=b_verts, group_box=sub_boxes[i],
                                    tot_ball_num=len(my_group.balls))
            # Check to see if find_verts fails
            if more_verts is not None:
                vert_ndxs, v_locs, v_rads, v_loc2s, v_rad2s, ball_nums, b_verts = more_verts
            # Every sphere needs a vert
            if sys.type == 'foam' and len(ball_nums) <= 0.25 * len(sys.balls['loc']):
                break
    # Create the doublets list
    doublets = [0 for _ in range(len(vert_ndxs))]
    # Incorporate the doublets into the vlocs, vballs, vrads lists and lose the vloc2s and vrad2s
    i = 0
    while i < len(v_locs):
        # Check for doubletness
        if v_rad2s[i] is not None:
            # Insert the relevant information into their respective lists
            vert_ndxs.insert(i + 1, vert_ndxs[i])
            v_locs.insert(i + 1, v_loc2s[i])
            v_rads.insert(i + 1, v_rad2s[i])
            doublets.insert(i + 1, 1)
            # Preserve the relational aspects of vrad2s and vloc2s
            v_rad2s.insert(i + 1, None)
            v_loc2s.insert(i + 1, [None, None, None])
        i += 1
    # Make the dataframe
    sys.net.verts = pd.DataFrame({"vballs": vert_ndxs, 'vloc': v_locs, 'vrad': v_rads, 'vdub': doublets})
    # Clear the print statement
    if sys.print_actions:
        print("\r                                                                  ", end="")
    sys.net.metrics['vert'] = time.perf_counter() - sys.net.start_time
    sys.net.build(surf_res=surf_res, max_vert=max_vert, box_size=box_size, build_surfs=build_surfs,
                  calc_verts=False, net_type=net_type, my_group=my_group, print_actions=print_actions)
    sys.net.metrics['splits'] = len(sub_boxes)
    if add_net_metrics:
        add_metrics(sys.net)


def split_net_slow(sys, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None, my_group=None,
                   print_actions=None, num_balls_sub_net=50, add_net_metrics=True, min_ball_split=30):
    """
    Splits a network into smaller sub-networks for parallel processing using a slow but thorough approach.

    This function divides a network into smaller sub-networks by:
    1. Sorting balls into spatial sub-boxes
    2. Merging small sub-boxes with nearby larger ones
    3. Processing each sub-network independently
    4. Recombining results into the main network

    Parameters
    ----------
    sys : Network
        The main network system to be split
    surf_res : float, optional
        Resolution for surface generation
    max_vert : float, optional
        Maximum distance to search for vertices
    box_size : float, optional
        Size of spatial partitioning boxes
    build_surfs : bool, optional
        Whether to build surfaces
    net_type : str, optional
        Type of network ('aw', 'pow', or 'prm')
    my_group : Group, optional
        Group of balls to process
    print_actions : bool, optional
        Whether to print progress information
    num_balls_sub_net : int, optional
        Target number of balls per sub-network
    add_net_metrics : bool, optional
        Whether to add network metrics
    min_ball_split : int, optional
        Minimum number of balls required for a sub-network

    Returns
    -------
    None

    Notes
    -----
    - Uses spatial partitioning to divide the network
    - Merges small sub-networks with nearby larger ones
    - Maintains network connectivity across sub-networks
    - Can be used for parallel processing of large networks
    """
    # Check to see if the pdb directory is suitable
    if sys.dir is None:
        sys.set_output_directory()
        os.mkdir(sys.dir + '/verts')
    # Sort the balls in the main network
    sys.net.sort_balls()
    # Calculate the group box
    group_box = sys.net.calc_box([sys.balls['loc'][_] for _ in my_group.balls],
                                 [sys.balls['rad'][_] for _ in my_group.balls], return_val=True, box_size=1.1)
    # Get the sub boxes
    sub_boxes = divide_box(group_box, round(len(my_group.balls) / num_balls_sub_net), c=3.0)
    print('num splits', len(sub_boxes))
    # Check for a max_vert that isn't defined
    if max_vert is None:
        max_vert = sys.net.settings['max_vert']

    # Sort the balls into their sub_boxes
    balls_lists = [[] for _ in range(len(sub_boxes))]
    ball_locs = sys.balls['loc']
    # Loop through the ball locations and sort the balls
    for ball in my_group.balls:
        loc = ball_locs[ball]
        # Loop through the sub boxes to find the placement of the ball
        for j, sub_box in enumerate(sub_boxes):
            if [sub_box[0][k] <= loc[k] <= sub_box[1][k] for k in range(3)] == [True, True, True]:
                balls_lists[j].append(ball)
    # If a list of balls is too small add it to another
    skip_boxes = []
    for i, balls_list in enumerate(balls_lists):
        # If no balls exist nothing to deal with
        if len(balls_list) == 0:
            skip_boxes.append(i)
            continue
        if len(balls_list) < min_ball_split:
            # Get the com of the balls to find the closes sub_box to add to
            balls_com = calc_com([ball_locs[_] for _ in balls_list])
            min_dist = np.inf
            closest_sub_box = None
            for j, sub_box in enumerate(sub_boxes):
                # Make sure we aren't adding to a sub_box scheduled for deletion
                if j in skip_boxes or j == i:
                    continue
                # Calculate the distance of the com of the sub_box from the balls_com
                my_dist = calc_dist(calc_com(sub_box), balls_com)
                # Replace the variables if they are closer
                if my_dist < min_dist:
                    closest_sub_box, min_dist = j, my_dist
            # Add the balls to the new sub_box
            balls_lists[closest_sub_box] += balls_list
            skip_boxes.append(i)
    # Instantiate the global variables
    global_vars(sys.net.sub_boxes, sys.net.box, sys.net.num_splits, sys.max_ball_rad, sys.net.sub_box_size)
    # Create the vertices file
    vert_file_name = sys.dir + '/verts/verts.txt'
    with open(vert_file_name, 'w') as verts_file:
        # Header
        verts_file.write(sys.name + " vertices")
        count = 0
        # Vertices
        for i, ball_list in enumerate(balls_lists):
            # Reset the variables
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, ball_nums, averts = None, None, None, None, None, None, None
            # Write the box info
            verts_file.write('\nBox {} - {}'.format(i, sub_boxes[i]))
            # Skip the boxes to be skipped
            if i in skip_boxes:
                continue
            # Get the balls we are tying to find vertices for
            check_balls = [_ for _ in ball_list if _ in my_group.balls]
            ball_nums = check_balls[:]
            # Find the initial vertices for the vertex group
            init_verts = find_verts(locs=sys.balls['loc'].to_numpy(), rads=sys.balls['rad'].to_numpy(),
                                    max_vert=max_vert, net_type=net_type, check_ndxs=check_balls,
                                    my_group=ball_nums, start_time=sys.net.start_time,
                                    vert_box=sys.foam_box, group_box=sub_boxes[i], vert_ndxs=vert_ndxs, vlocs=vlocs,
                                    vrads=vrads, vloc2s=vloc2s, vrad2s=vrad2s, b_verts=averts,
                                    tot_ball_num=len(my_group.balls), start_vert=count, split=True)
            # Check to see if find_verts fails
            if init_verts is not None:
                vert_ndxs, vlocs, vrads, vloc2s, vrad2s, ball_nums, averts = init_verts
                count += len(vert_ndxs)
            # Check for disconnects in the network
            while len(ball_nums) > 0:
                # Grab the initial ball for the next search
                a0 = ball_nums.pop()

                # Find verts again
                more_verts = find_verts(b0=a0, locs=sys.balls['loc'].to_numpy(), rads=sys.balls['rad'].to_numpy(),
                                        max_vert=max_vert, net_type=net_type, check_ndxs=ball_nums,
                                        my_group=check_balls, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                                        vloc2s=vloc2s, vrad2s=vrad2s, start_time=sys.net.start_time,
                                        vert_box=sys.foam_box, b_verts=averts, group_box=sub_boxes[i],
                                        tot_ball_num=len(my_group.balls), printing=True if i == 5 else False,
                                        start_vert=count, split=True)
                # Check to see if find_verts fails
                if more_verts is not None:
                    vert_ndxs, vlocs, vrads, vloc2s, vrad2s, ball_nums, averts = more_verts
                    count += len(vert_ndxs)
                # Every sphere needs a vert
                if sys.type == 'foam' and len(ball_nums) <= 0.25 * len(sys.balls['loc']):
                    break

            # Write the vertex information into the vert file
            for j, vert in enumerate(vert_ndxs):
                # Write the line
                loc2 = 'None'
                if vloc2s[j] is not None:
                    loc2 = str([_ for _ in vloc2s[j]])
                verts_file.write('\n' + str(j) + ';' + str(vert) + ';' + str([_ for _ in vlocs[j]]) + ';' + str(vrads[j]) + ';' +
                                 loc2 + ';' + str(vrad2s[j]))

    # Combine the subnetworks
    vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts = combine_nets(vert_file_name, len(sys.balls['loc']))
    # Create the doublets list
    doublets = [0 for _ in range(len(vert_ndxs))]
    # Incorporate the doublets into the vlocs, vballs, vrads lists and lose the vloc2s and vrad2s
    i = 0
    while i < len(vlocs):
        # Check for doubletness
        if vrad2s[i] is not None:
            # Insert the relevant information into their respective lists
            vert_ndxs.insert(i + 1, vert_ndxs[i])
            vlocs.insert(i + 1, vloc2s[i])
            vrads.insert(i + 1, vrad2s[i])
            doublets.insert(i + 1, 1)
            # Preserve the relational aspects of vrad2s and vloc2s
            vrad2s.insert(i + 1, None)
            vloc2s.insert(i + 1, [None, None, None])
        i += 1
    # Make the dataframe
    sys.net.verts = pd.DataFrame({"vballs": vert_ndxs, 'vloc': vlocs, 'vrad': vrads, 'vdub': doublets})
    # Clear the print statement
    if sys.print_actions:
        print("\r                                                                  ", end="")
    sys.net.metrics['vert'] = time.perf_counter() - sys.net.start_time
    sys.net.build(surf_res=surf_res, max_vert=max_vert, box_size=box_size, build_surfs=build_surfs,
                  calc_verts=False, net_type=net_type, my_group=my_group, print_actions=print_actions)
    sys.net.metrics['splits'] = len(sub_boxes)
    if add_net_metrics:
        add_metrics(sys.net)


def combine_nets(verts, num_balls):
    """
    Combines multiple subnetworks into a single network by merging vertices and maintaining ball-vertex relationships.

    This function:
    1. Reads vertex information from a file containing multiple subnetworks
    2. Combines vertices while avoiding duplicates
    3. Maintains the relationship between balls and their vertices
    4. Handles doublet vertices (vertices sharing the same set of balls)

    Parameters
    ----------
    verts : str
        Path to the file containing vertex information from subnetworks
    num_balls : int
        Total number of balls in the system

    Returns
    -------
    tuple
        A tuple containing:
        - vert_ndxs : list of lists
            List of ball indices for each vertex
        - vlocs : list of numpy.ndarray
            List of vertex locations
        - vrads : list of float
            List of vertex radii
        - vloc2s : list of numpy.ndarray or None
            List of secondary vertex locations for doublets
        - vrad2s : list of float or None
            List of secondary vertex radii for doublets
        - averts : list of lists
            List mapping ball indices to their vertices

    Notes
    -----
    - Handles both regular vertices and doublet vertices
    - Maintains spatial relationships between balls and vertices
    - Uses binary search for efficient vertex insertion
    - Preserves the order of vertices while combining networks
    """
    # Create the dictionary based on the sub_boxes
    file_verts = []
    # Go through the regular vertices
    with open(verts, 'r') as verts_file:
        # Go through the vertices
        for i, line in enumerate(verts_file):
            if i == 0 or line[:3] == 'Box':
                continue
            # Split the line
            line = line.split(';')
            balls = line[1][1:-1].split(',')
            loc = line[2][1:-1].split(',')
            if line[4][0] == 'N':
                loc2 = None
            else:
                loc2 = line[4][1:-1].split(',')
                loc2 = [float(_) for _ in loc2]
            if line[5][0] == 'N':
                rad2 = None
            else:
                rad2 = float(line[5])
            my_vert = {'balls': [int(_) for _ in balls if _ != ''], 'loc': [float(_) for _ in loc if _ != ''], 'rad': float(line[3]),
                               'loc2': loc2, 'rad2': rad2}
            file_verts.append(my_vert)
    # Go through the sub_boxes one by one
    averts = [[] for _ in range(num_balls)]
    verts = []

    for vert in file_verts:
        # Create a tracking variable for whether the vertex has been found or not
        vert_found = False
        # Go through the ball vertices for the first ball in the vertex
        for my_vert in averts[vert['balls'][0]]: ## This is where ndx search comes in
            if my_vert < len(verts) and verts[my_vert]['balls'] == vert['balls']:
                vert_found = True
                break
        if vert_found:
            continue
        else:
            # Add the vertex to the list of balls
            verts.append(vert)
            vert_ndx = len(verts) - 1
            # Add the vertex to each ball in the vertex:
            for ndx in vert['balls']:
                # Get the list of vertex indices

                ball_vert_ndxs = [verts[_]['balls'] for _ in averts[ndx]]
                # find the place where to insert the vertex
                my_ndx = ndx_search(ball_vert_ndxs, vert['balls'])
                averts[ndx].insert(my_ndx, vert_ndx)
    vert_ndxs, vlocs, vrads, vloc2s, vrad2s = [], [], [], [], []
    for i in range(len(verts)):
        vert_ndxs.append(verts[i]['balls'])
        vlocs.append(verts[i]['loc'])
        vrads.append(verts[i]['rad'])
        vloc2s.append(verts[i]['loc2'])
        vrad2s.append(verts[i]['rad2'])
    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts
