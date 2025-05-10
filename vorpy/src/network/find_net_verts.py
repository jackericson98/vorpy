import time
import pandas as pd
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import calc_dist
from vorpy.src.network.find_verts import find_verts
from vorpy.src.output import write_verts


def find_net_verts(net):

    """
    Finds vertices in a network by iteratively searching for valid vertex configurations.

    This function implements the main vertex finding algorithm for different network types (aw, pow, prm).
    It starts by finding an initial set of vertices and then continues to find additional vertices
    until all balls in the network are either part of a vertex or determined to be encapsulated.

    Parameters
    ----------
    net : Network
        Network object containing:
        - balls : pandas.DataFrame
            DataFrame containing ball information including locations and radii
        - settings : dict
            Dictionary of network settings including:
            - max_vert : float
                Maximum vertex radius
            - net_type : str
                Type of network ('aw', 'pow', or 'prm')
            - print_metrics : bool
                Flag to enable progress printing
            - foam_box : list
                Bounding box for foam vertices
            - ball_type : str
                Type of balls ('foam' or other)
        - group : list, optional
            List of ball indices in the group
        - metrics : dict
            Dictionary for storing performance metrics
        - box : dict
            Dictionary containing bounding boxes for different components

    Returns
    -------
    tuple
        A tuple containing:
        - vert_ndxs : list
            List of vertex indices
        - vlocs : list
            List of vertex locations
        - vrads : list
            List of vertex radii
        - vloc2s : list
            List of secondary vertex locations (for doublets)
        - vrad2s : list
            List of secondary vertex radii (for doublets)
        - sphere_check_list : list
            List of remaining unvisited balls
        - averts : dict
            Dictionary mapping balls to their vertices

    Notes
    -----
    - The function handles encapsulated balls by checking if any ball is fully contained within another
    - For foam networks, the search stops when less than 25% of balls remain unvisited
    - Doublets (vertices with two possible locations) are handled by keeping track of both locations
    """
    # Create the group indices
    if net.group is None:
        net.group = [_['num'] for i, _ in net.balls.iterrows()]
    # Create the sphere check list
    sphere_check_list = net.group.copy()
    # Get the indices of the balls in the network to keep track of the balls that haven't been visited
    my_guuy = find_verts(locs=net.balls['loc'].to_numpy(), rads=net.balls['rad'].to_numpy(),
                         max_vert=net.settings['max_vert'], net_type=net.settings['net_type'], check_ndxs=sphere_check_list,
                         my_group=net.group, start_time=net.metrics['start'], print_metrics=net.settings['print_metrics'],
                         vert_box=net.settings['foam_box'], box=net.box['verts'])
    # If the function returns a valid vertex, set the variables
    if my_guuy is not None:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
    # Check to see if any of the balls are encapsulated
    if len(sphere_check_list) > 0:
        # Create the skip numbers list
        skip_nums = []
        # Iterate through the sphere check list
        for sphere in sphere_check_list:
            # Get the radius and location of the sphere
            sphere_rad, sphere_loc = net.balls['rad'][sphere], net.balls['loc'][sphere]
            # Create the sphere box
            sphere_box = box_search(sphere_loc)
            # Get the balls within the sphere box
            close_spheres = get_balls(sphere_box, dist=max(net.balls['rad']) - sphere_rad)
            # Iterate through the close spheres
            for sphere2 in close_spheres:
                # Check if the sphere is fully encapsulated by another sphere
                if calc_dist(sphere_loc, net.balls['loc'][sphere2]) < abs(net.balls['rad'][sphere2] - sphere_rad):
                    print("\nUh oh! Ball # {} is fully encapsulated by ball # {}! Skipping {}"
                          .format(sphere, sphere2, sphere))
                    skip_nums.append(sphere)
                    break
        # Iterate through the skip numbers
        for _ in skip_nums:
            sphere_check_list.pop(sphere_check_list.index(_))
    # Check for disconnects in the network
    while len(sphere_check_list) > 0:
        # Get the next sphere to check
        a0 = sphere_check_list.pop()
        # Find the vertices
        my_guuy = find_verts(b0=a0, locs=net.balls['loc'].to_numpy(), rads=net.balls['rad'].to_numpy(),
                             max_vert=net.settings['max_vert'], net_type=net.settings['net_type'], check_ndxs=sphere_check_list,
                             my_group=net.group, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                             vloc2s=vloc2s, vrad2s=vrad2s, start_time=net.metrics['start'],
                             vert_box=net.settings['foam_box'], b_verts=averts, box=net.box['verts'])
        # If the function returns a valid vertex, set the variables
        if my_guuy is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
        # If the network is a foam network and less than 25% of the balls remain unvisited, break
        if net.settings['ball_type'] == 'foam' and len(sphere_check_list) <= 0.25 * len(net.balls['loc']):
            print(f'Missing Ball Indices:\n{sphere_check_list}\n')
            break
    # Create the doublets list
    doublets = [0 for _ in range(len(vert_ndxs))]
    # Incorporate the doublets into the v_locs, balls, v_rads lists and lose the v_loc2s and v_rad2s
    i = 0
    while i < len(vlocs):
        # Check for doubletness
        if vrad2s[i] is not None:
            # Insert the relevant information into their respective lists
            vert_ndxs.insert(i + 1, vert_ndxs[i])
            vlocs.insert(i + 1, vloc2s[i])
            vrads.insert(i + 1, vrad2s[i])
            doublets[i] = 2
            doublets.insert(i + 1, 1)
            # Preserve the relational aspects of vrad2s and vloc2s
            vrad2s.insert(i + 1, None)
            vloc2s.insert(i + 1, [None, None, None])
        i += 1

    # Make the dataframe
    net.verts = pd.DataFrame({"balls": vert_ndxs, 'loc': vlocs, 'rad': vrads, 'dub': doublets})
    # Clear the print statement
    print("\r                                                                  ", end="")
    net.metrics['vert'] = time.perf_counter() - net.metrics['start']
    write_verts(net)
