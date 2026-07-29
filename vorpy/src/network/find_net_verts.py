import time
import pandas as pd
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import calc_dist
from vorpy.src.network.find_verts import find_verts
from vorpy.src.output import write_verts


def vertex_is_allowed(net, vertex_balls):
    """
    Determine whether a candidate vertex should be retained.

    Normal network
    --------------
    The vertex must contain at least one ball from net.group.

    Interface network
    -----------------
    The vertex must contain at least one ball from every selection
    stored in net.iface_grps.
    """
    vertex_balls = set(vertex_balls)

    if net.iface_grps is not None:
        return all(
            bool(vertex_balls.intersection(group_indices))
            for group_indices in net.iface_grps
        )

    if net.group is None:
        return True

    return bool(vertex_balls.intersection(net.group))


def _load_cached_vertex_state(net):
    """Return cached vertices as independent native ``find_verts`` arrays."""
    records = list(getattr(net, 'cached_interface_vertex_state', None) or [])
    if not records:
        return None

    records.sort(key=lambda record: record.get('source_index', 0))
    vert_ndxs = []
    vlocs = []
    vrads = []
    vloc2s = []
    vrad2s = []
    b_verts = [[] for _ in range(len(net.balls))]

    for record in records:
        row = record['data']
        balls = [int(ball) for ball in row.get('balls', [])]
        if not vertex_is_allowed(net, balls):
            continue

        vertex_index = len(vert_ndxs)
        vert_ndxs.append(balls)
        vlocs.append([float(value) for value in row.get('loc', [])])
        vrads.append(float(row.get('rad', 0.0)))

        loc2 = row.get('loc2', None)
        if loc2 is None or all(value is None for value in loc2):
            vloc2s.append([None, None, None])
            vrad2s.append(None)
        else:
            vloc2s.append([float(value) for value in loc2])
            rad2 = row.get('rad2', None)
            vrad2s.append(None if rad2 is None else float(rad2))

        for ball in balls:
            b_verts[ball].append(vertex_index)

    if not vert_ndxs:
        return None

    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, b_verts


def _store_native_vertex_state(net, vert_ndxs, vlocs, vrads, vloc2s, vrad2s):
    """Persist the completed native search state for System-level caching."""
    net.interface_vertex_state = []
    for balls, loc, rad, loc2, rad2 in zip(
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s):
        net.interface_vertex_state.append({
            'balls': [int(ball) for ball in balls],
            'loc': [float(value) for value in loc],
            'rad': float(rad),
            'loc2': (
                [None, None, None]
                if loc2 is None or all(value is None for value in loc2)
                else [float(value) for value in loc2]
            ),
            'rad2': None if rad2 is None else float(rad2),
        })


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
    # print("\n[VERTEX RUN SETTINGS]")
    # print(f"  network mode = {getattr(net, 'network_mode', None)}")
    # print(f"  net type     = {net.settings['net_type']}")
    # print(f"  max vert     = {net.settings['max_vert']}")
    # print(f"  group size   = {len(net.group) if net.group is not None else None}")
    # print(f"  group        = {net.group}")
    # print(f"  iface grps   = {net.iface_grps}")
    # print(f"  vert box     = {net.box['verts']}")
    # print(f"  foam box     = {net.settings['foam_box']}")
    print("  group size   =", len(net.group))
    print(
        "  group range  =",
        (min(net.group), max(net.group)) if net.group else None
    )

    if net.iface_grps is not None:
        print(
            "  iface sizes  =",
            tuple(len(group_indices) for group_indices in net.iface_grps)
        )
    # Create the group indices
    if net.group is None:
        net.group = [_['num'] for i, _ in net.balls.iterrows()]
    # Create the sphere check list
    sphere_check_list = net.group.copy()

    cached_state = _load_cached_vertex_state(net)
    if cached_state is None:
        vert_ndxs = vlocs = vrads = vloc2s = vrad2s = averts = None
        cached_count = 0
    else:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts = cached_state
        cached_count = len(vert_ndxs)
        print(
            f"[INTERFACE CACHE LOAD] verts={cached_count} "
            f"covered_balls={sum(bool(vertices) for vertices in averts)}",
            flush=True,
        )

    # Continue normal discovery with cached vertices available for duplicate
    # detection and traversal adjacency.
    my_guuy = find_verts(
        locs=net.balls['loc'].to_numpy(),
        rads=net.balls['rad'].to_numpy(),
        max_vert=net.settings['max_vert'],
        net_type=net.settings['net_type'],
        check_ndxs=sphere_check_list,

        # The complete interface network search selection.
        my_group=net.group,

        # The two original interface sides. These are distinct from
        # my_group and will later control candidate-vertex acceptance.
        iface_grps=net.iface_grps,
        vert_ndxs=vert_ndxs,
        vlocs=vlocs,
        vrads=vrads,
        vloc2s=vloc2s,
        vrad2s=vrad2s,
        b_verts=averts,

        start_time=net.metrics['start'],
        print_metrics=net.settings['print_metrics'],
        vert_box=net.settings['foam_box'],
        box=net.box['verts'],
    )
    # If the function returns a valid vertex, set the variables.
    if my_guuy is not None:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
    elif cached_state is None:
        return

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
        # if net.iface_grps is not None:
        #     interface_side_1 = set(net.iface_grps[0])
        #
        #     sphere_check_list = [
        #         ball
        #         for ball in sphere_check_list
        #         if ball in interface_side_1
        #     ]
        # Interface vertex traversal begins from a verified interface seed
        # and follows all connected interface vertices. Do not perform the
        # standard group-mode disconnected-component seed search afterward.
        if net.iface_grps is not None:
            break
        # Get the next sphere to check
        a0 = sphere_check_list.pop()
        # Find the vertices
        my_guuy = find_verts(
            b0=a0,
            locs=net.balls['loc'].to_numpy(),
            rads=net.balls['rad'].to_numpy(),
            max_vert=net.settings['max_vert'],
            net_type=net.settings['net_type'],
            check_ndxs=sphere_check_list,
            my_group=net.group,
            iface_grps=net.iface_grps,
            vert_ndxs=vert_ndxs,
            vlocs=vlocs,
            vrads=vrads,
            vloc2s=vloc2s,
            vrad2s=vrad2s,
            start_time=net.metrics['start'],
            vert_box=net.settings['foam_box'],
            b_verts=averts,
            box=net.box['verts'],
        )
        # If the function returns a valid vertex, set the variables
        if my_guuy is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
        # If the network is a foam network and less than 25% of the balls remain unvisited, break
        if net.settings['ball_type'] == 'foam' and len(sphere_check_list) <= 0.25 * len(net.balls['loc']):
            print(f'Missing Ball Indices:\n{sphere_check_list}\n')
            break
    # Save the lossless native representation before doublets are expanded
    # into separate dataframe rows.
    if net.iface_grps is not None:
        _store_native_vertex_state(
            net, vert_ndxs, vlocs, vrads, vloc2s, vrad2s
        )

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