import time
import pandas as pd
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import calc_dist
from vorpy.src.network.find_verts import find_verts
from vorpy.src.output import write_verts


def vertex_is_allowed(net, vertex_balls):
    """
    Determine whether a vertex belongs to the requested network scope.

    Normal networks retain a vertex when at least one defining ball belongs
    to ``net.group``. Interface networks retain a vertex only when its four
    defining balls include at least one ball from each interface side.

    Parameters
    ----------
    net : Network
        Network providing ``group`` and optional ``iface_grps`` constraints.
    vertex_balls : collection of int
        Defining ball indices for the candidate vertex.

    Returns
    -------
    bool
        True when the vertex belongs to the requested network scope.
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
    """
    Reconstruct native vertex-search arrays from cached interface state.

    Cached records are filtered against the current network scope and converted
    back into the parallel arrays expected by ``find_verts``. Ball-to-vertex
    adjacency is rebuilt so cached vertices participate normally in duplicate
    detection and traversal.

    Returns
    -------
    tuple or None
        ``(vert_ndxs, vlocs, vrads, vloc2s, vrad2s, b_verts)`` when usable cached
        state exists, otherwise ``None``.
    """
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
    """
Store the native vertex-search representation for interface reuse.

The unexpanded primary/secondary vertex representation is preserved so future
interface calculations can resume discovery without reconstructing doublets
from exported dataframe rows.
"""
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
    Find and store all vertices belonging to a network.

    The search begins from the network's requested group and traverses connected
    vertices using ``find_verts``. Previously cached interface vertices may be
    loaded as an initial native search state so duplicate detection and adjacency
    remain available during continued discovery.

    Any balls left uncovered after the initial traversal are checked for complete
    encapsulation. Remaining balls are then used as additional seeds, allowing
    disconnected or incompletely traversed regions to be discovered. Interface
    networks use this reseeding behavior as well, which allows the same interface
    to be approached from multiple starting regions.

    Doublet vertices are maintained internally as primary/secondary solutions
    during discovery and expanded into separate dataframe rows only after the
    native search is complete.

    Parameters
    ----------
    net : Network
        Network being solved. The object provides ball coordinates, radii, group
        definitions, interface groups, search settings, spatial boxes, and build
        metrics.

    Notes
    -----
    The function updates ``net.verts`` in place, stores interface-native vertex
    state when applicable, updates the vertex-build timing metric, and writes the
    vertex output through ``write_verts``.
    """
    # Create the group indices
    if net.group is None:
        net.group = net.balls['num'].tolist()
    # Track group balls that have not yet been reached by vertex traversal.
    sphere_check_list = net.group.copy()

    cached_state = _load_cached_vertex_state(net)
    if cached_state is None:
        vert_ndxs = vlocs = vrads = vloc2s = vrad2s = averts = None
    else:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts = cached_state

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
        max_ball_rad = max(net.balls['rad'])
        # Iterate through the sphere check list
        for sphere in sphere_check_list:
            # Get the radius and location of the sphere
            sphere_rad, sphere_loc = net.balls['rad'][sphere], net.balls['loc'][sphere]
            # Create the sphere box
            sphere_box = box_search(sphere_loc)
            # Get the balls within the sphere box
            close_spheres = get_balls(sphere_box, dist=max_ball_rad - sphere_rad)
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

    if net.iface_grps is not None:
        interface_balls = net.iface_grps[0] | net.iface_grps[1]
        sphere_check_list[:] = [ball for ball in sphere_check_list if ball in interface_balls]

    # Reseed from any balls not reached by the initial traversal.
    while len(sphere_check_list) > 0:

        # Remove the seed before calling find_verts; b0 explicitly restores it
        # as the first seed candidate inside the new traversal.
        seed_ball = sphere_check_list.pop()
        # Find the vertices
        my_guuy = find_verts(
            b0=seed_ball,
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
    doublets = [0] * len(vert_ndxs)
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