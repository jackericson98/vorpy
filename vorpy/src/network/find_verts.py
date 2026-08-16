import time
import numpy as np
from numpy import sqrt
from vorpy.src.network.find_v0 import find_v0
from vorpy.src.network.fast import find_site_container
from vorpy.src.calculations import get_time
from vorpy.src.calculations import ndx_search
from vorpy.src.calculations import calc_vert


# Find network function. Keeps searching the network until all verts are found
def find_verts(locs, rads, max_vert, net_type, check_ndxs, b0=None, my_group=None,
               iface_grps=None, b_verts=None, vert_ndxs=None, vlocs=None, vrads=None, vloc2s=None, vrad2s=None,
               start_time=0, box=None, vert_box=None, group_box=None, tot_ball_num=None,
               printing=False, start_vert=0, split=False):
    """
    Traverse a network and discover vertices connected to an initial seed.

    The search begins from a verified vertex and follows each three-ball edge
    to locate neighboring four-ball vertices. Existing vertex state may be
    supplied so additional traversals can extend a partially solved network
    without duplicating previously discovered vertices.

    Normal group searches retain vertices containing at least one requested
    group ball. Interface searches require every retained vertex to contain
    at least one defining ball from each interface side. Interface seed
    candidates are prioritized by surface-to-surface proximity between the
    two groups.

    Parameters
    ----------
    locs : array-like
        Ball-center coordinates for the full system.
    rads : array-like
        Ball radii for the full system.
    max_vert : float
        Maximum permitted vertex radius/search extent.
    net_type : {'aw', 'pow', 'prm'}
        Network geometry being solved.
    check_ndxs : list of int
        Balls not yet reached by the current traversal.
    b0 : int, optional
        Preferred starting ball for seed discovery.
    my_group : collection of int, optional
        Ball indices defining the requested normal network or interface union.
    iface_grps : tuple, optional
        Two disjoint ball-index collections defining an interface.
    b_verts : list of list, optional
        Existing ball-to-vertex adjacency.
    vert_ndxs : list, optional
        Defining ball indices for previously discovered vertices.
    vlocs : list, optional
        Primary vertex locations.
    vrads : list, optional
        Primary vertex radii.
    vloc2s : list, optional
        Secondary vertex locations for doublets.
    vrad2s : list, optional
        Secondary vertex radii for doublets.
    start_time : float, optional
        Network build start time used for progress reporting.
    box : list, optional
        Global geometric bounds for accepted vertices.
    vert_box : list, optional
        Optional retained vertex bounds.
    group_box : list, optional
        Optional bounds used during seed discovery.
    tot_ball_num : int, optional
        Ball count used to estimate progress.
    printing : bool, optional
        Enable detailed diagnostic output.
    start_vert : int, optional
        Existing vertex count used in progress reporting.
    split : bool, optional
        Reserved for split-search workflows.

    Returns
    -------
    tuple or None
        ``(vert_ndxs, vlocs, vrads, vloc2s, vrad2s, check_ndxs, b_verts)``
        when a seed is found, otherwise ``None``.
    """

    # Calculate the maximum input ball radius
    max_ball_rad = max(rads)
    # Normalize the two interface selections.
    if iface_grps is not None:
        if len(iface_grps) != 2:
            raise ValueError(
                "iface_grps must contain exactly two ball-index collections."
            )

        iface_grps = tuple(
            set(group_indices)
            for group_indices in iface_grps
        )

        overlap = iface_grps[0].intersection(iface_grps[1])

        if overlap:
            raise ValueError(
                f"Interface groups overlap by {len(overlap)} balls."
            )
    # Get the group balls from which to check vertices against
    if my_group is None or len(my_group) == len(locs):
        my_group = list(range(len(locs)))
        # Calculate the rough number of vertices
        if tot_ball_num is None:
            tot_verts = int(6.6 * len(locs))
    # If a group was provided make sure to get its indices
    elif my_group is not None:
        # Calculate the number of vertices
        if tot_ball_num is None:
            tot_verts = int(6.6 * len(my_group) + int(60 * sqrt(len(my_group))))

    if tot_ball_num is not None:
        tot_verts = int(6.6 * tot_ball_num + int(60 * sqrt(tot_ball_num)))
    if b_verts is None:
        b_verts = [[] for _ in range(len(locs))]
    # Find the first verified vertex
    if len(my_group) == 1:
        v0 = find_v0(locs=locs, rads=rads, b_verts=b_verts, max_vert=max_vert, net_type=net_type, b0=my_group[0],
                     group_ndxs=my_group, vert_ndxs=vert_ndxs, group_box=group_box, iface_grps=iface_grps)

    elif len(my_group) == 4:
        v0_loc, v0_rad, v0_loc2, v0_rad2 = calc_vert(locs=[locs[_] for _ in my_group],
                                                     rads=[rads[_] for _ in my_group])
        v0 = {'balls': my_group, 'loc': v0_loc, 'rad': v0_rad, 'loc2': v0_loc2, 'rad2': v0_rad2}
    else:
        # In interface mode, seed from first-side balls ordered by their
        # minimum surface-to-surface separation from the opposite side.
        if iface_grps is not None:
            interface_side_1 = iface_grps[0]
            interface_side_2 = iface_grps[1]

            side_2_indices = np.array(sorted(interface_side_2), dtype=int)

            side_2_locs = np.asarray([locs[ball] for ball in side_2_indices], dtype=float)

            side_2_rads = np.asarray([rads[ball] for ball in side_2_indices], dtype=float)

            seed_cutoff = max_vert / 10

            seed_distances = []

            for ball in check_ndxs:
                if ball not in interface_side_1:
                    continue

                center_distances = np.linalg.norm(side_2_locs - np.asarray(locs[ball], dtype=float), axis=1)

                surface_distances = (center_distances - float(rads[ball]) - side_2_rads)

                minimum_surface_distance = float(np.min(surface_distances))

                if minimum_surface_distance <= seed_cutoff:
                    seed_distances.append((minimum_surface_distance, ball))

            # Try first-side balls closest to the opposite interface side first.
            # All eligible seeds remain available if the nearest candidates fail.
            seed_distances.sort(key=lambda item: item[0])

            seed_ndxs = [ball for _, ball in seed_distances]
            # For interface builds, estimate the total number of vertices from the
            # number of atoms that actually participate in the interface rather than
            # from the full union of both groups.
            if iface_grps is not None:
                tot_verts = max(50, 10 * len(seed_ndxs))

        else:
            seed_distances = None
            seed_ndxs = list(check_ndxs)

        # Force an explicitly supplied reseed ball to the front of the seed list.
        if b0 is not None:
            if iface_grps is None:
                valid_b0 = b0 in my_group
            else:
                valid_b0 = b0 in iface_grps[0] or b0 in iface_grps[1]

            if valid_b0:
                if b0 in seed_ndxs:
                    seed_ndxs.remove(b0)
                seed_ndxs.insert(0, b0)

        v0 = None

        for seed_ball in seed_ndxs:
            v0 = find_v0(locs=locs, rads=rads, b_verts=b_verts, max_vert=max_vert, net_type=net_type, b0=seed_ball,
                         group_ndxs=my_group, iface_grps=iface_grps, vert_ndxs=vert_ndxs,
                         group_box=group_box, box=box)

            if v0 is not None:
                break

    # Defensive verification: an interface seed must contain at least
    # one defining ball from each interface side.
    if v0 is not None and iface_grps is not None:
        v0_balls = set(v0["balls"])

        belongs_to_interface = all(bool(v0_balls.intersection(group_indices)) for group_indices in iface_grps)

        if not belongs_to_interface:
            v0 = None
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
    # Throttle progress output.
    last_print = 0
    # Traverse neighboring vertices depth-first.
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
            current_time = time.perf_counter()
            if current_time - last_print > 0.25:
                percentage = min((len(vlocs) / tot_verts) * 100, 100)
                my_time = current_time - start_time
                h, m, s = get_time(my_time)
                print("\rRun Time = {}:{:02d}:{:2.2f} - Process: finding vertices: {} verts - {:.2f} %"
                      .format(int(h), int(m), round(s, 2), len(vert_ndxs) + start_vert, percentage), end="")
                last_print = current_time
            # Get the edge from the top of the stack
            edge_balls, vert = e_stack.pop()
            # Find the next site in the network
            search_group = (
                tuple(iface_grps)
                if iface_grps is not None
                else my_group
            )

            vert_ndx_pr = find_site_container(edge_balls=edge_balls, locs=locs, rads=rads, b_verts=b_verts,
                                              vert_ndxs=vert_ndxs, max_vert=max_vert, net_type=net_type,
                                              vn_1=vert['balls'], box=box, vn_1_loc=vert['loc'],
                                              group_ndxs=search_group, printing=printing, max_ball_rad=max_ball_rad)


            # If the vertex is none continue
            if vert_ndx_pr is None:
                continue
            # Set the vertex and its index
            my_vert, metrics = vert_ndx_pr
            if iface_grps is not None:
                candidate_balls = set(my_vert["balls"])

                belongs_to_interface = all(
                    bool(candidate_balls.intersection(group_indices))
                    for group_indices in iface_grps
                )

                if not belongs_to_interface:
                    if printing:
                        print("\n[REJECTED NON-INTERFACE VERTEX]")
                        print(f"  balls = {my_vert['balls']}")

                    continue

            if my_vert['loc'] is None:
                continue
            if box is not None and any([box[0][k] > my_vert['loc'][k] or my_vert['loc'][k] > box[1][k] for k in range(3)]):
                continue
            if box is not None and 'loc2' in my_vert and my_vert['loc2'] is not None and any([box[0][k] > my_vert['loc2'][k] or my_vert['loc2'][k] > box[1][k] for k in range(3)]):
                my_vert['loc2'], my_vert['rad2'] = None, None
            # Queue the new vertex for traversal and append it to native state.
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
            # Update ball-to-vertex adjacency and mark reached balls as visited.
            for ball in my_vert['balls']:
                # noinspection PyTypeChecker
                b_vert_ndxs = [vert_ndxs[_] for _ in b_verts[ball]]
                # noinspection PyTypeChecker
                b_verts[ball].insert(ndx_search(b_vert_ndxs, my_vert['balls']), len(vert_ndxs) - 1)
                if ball in check_ndxs:
                    check_ndxs.remove(ball)
    # Return the values of the vertices
    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, check_ndxs, b_verts
