import time
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import calc_circ
from vorpy.src.network.fast import find_site_container
from vorpy.src.network.slow import find_site_container_slow
from vorpy.src.network.slow import find_site
from vorpy.src.calculations import verify_site


def is_interface_vertex(vertex_balls, iface_grps):
    """
    Check whether a vertex belongs to the requested interface.

    An interface vertex must contain at least one defining ball from
    each interface group.
    """
    if iface_grps is None:
        return True

    vertex_balls = set(vertex_balls)

    return all(
        bool(vertex_balls.intersection(group_indices))
        for group_indices in iface_grps
    )


def find_v0(locs, rads, b_verts, max_vert, net_type, b0=None, group_ndxs=None, iface_grps=None, metrics=None,
            vert_ndxs=None, group_box=None, box=None, timeout=None):
    """
    Find the initial verified vertex for a network.

    For a normal group network, this behaves as before.

    For an interface network, the returned seed vertex must contain at
    least one defining ball from each collection in iface_grps.
    """
    start_time = time.perf_counter()

    def timed_out():
        return timeout is not None and time.perf_counter() - start_time >= timeout
    # Make sure we have an existing-vertex list for the lower-level
    # site-finding functions.
    if vert_ndxs is None:
        vert_ndxs = []

    # Normalize the two interface selections so membership and
    # intersection checks are efficient.
    if iface_grps is not None:
        if len(iface_grps) != 2:
            raise ValueError(
                "find_v0 requires exactly two interface groups."
            )

        iface_grps = tuple(
            set(group_indices)
            for group_indices in iface_grps
        )

        # All balls that belong to either side of the interface.
        interface_balls = iface_grps[0] | iface_grps[1]

    else:
        interface_balls = None

    # If a specific starting ball was supplied for an interface search,
    # it must belong to one of the two interface groups.
    #
    # This does not by itself guarantee that the final vertex belongs
    # to the interface. The completed four-ball candidate is checked
    # later before it is returned.
    if (
        b0 is not None
        and interface_balls is not None
        and b0 not in interface_balls
    ):
        return None

    # Determine the initial spatial box from which to search.
    if b0 is not None:
        # Search around the explicitly supplied starting ball.
        my_box = box_search(locs[b0])

    elif group_box is not None:
        # Search around the center of the supplied group bounding box.
        group_center = [
            0.5 * abs(group_box[1][i] - group_box[0][i])
            + group_box[0][i]
            for i in range(3)
        ]

        my_box = box_search(group_center)

    elif group_ndxs is not None:
        # Search around a ball near the middle of the group-index list.
        middle_group_ball = group_ndxs[
            int(len(group_ndxs) / 2)
        ]

        my_box = box_search(locs[middle_group_ball])

    else:
        # Search around the middle ball in the complete system.
        my_box = box_search(
            locs[int(len(locs) / 2)]
        )

    # If no starting ball was supplied, find one in or near my_box.
    if b0 is None:
        b0s = []
        inc = 0

        # Expand the spatial search until at least one valid starting
        # ball is found.
        while len(b0s) < 1:
            b0s = get_balls([my_box], inc)

            # For a normal group search, the starting ball must belong
            # to the group and must not already have a known vertex.
            if group_ndxs is not None:
                b0s = [
                    ball
                    for ball in b0s
                    if ball in group_ndxs
                    and len(b_verts[ball]) == 0
                ]

            # For an interface search, the starting ball must belong to
            # one of the two interface sides.
            if interface_balls is not None:
                b0s = [
                    ball
                    for ball in b0s
                    if ball in interface_balls
                ]

            inc += 1

        # Use the first acceptable ball as the initial search ball.
        b0 = b0s[0]

    # Find nearby b1 candidates around b0.
    b1s = []
    inc = 0

    # Retain the existing behavior of gathering at least five nearby
    # candidates before sorting them by surface-to-surface distance.
    while len(b1s) < 5:
        b1s = get_balls([my_box], inc)
        inc += 1

    # Calculate the separation between b0 and every b1 candidate.
    b1_dists = [
        calc_dist(locs[b1], locs[b0])
        - (rads[b0] + rads[b1])
        for b1 in b1s
    ]

    # Sort b1 candidates from nearest to farthest.
    _, b1s_sorted = zip(
        *sorted(
            zip(b1_dists, b1s),
            key=lambda x: x[0],
        )
    )

    # Remove b0 itself.
    b1s_sorted = [
        ball
        for ball in b1s_sorted
        if ball != b0
    ]

    # For an interface search, try balls from the opposite interface
    # side first. Same-side candidates are retained as fallbacks because
    # b2 or the fourth defining ball may still complete the interface.
    if iface_grps is not None:
        if b0 in iface_grps[0]:
            opposite_group = iface_grps[1]

        elif b0 in iface_grps[1]:
            opposite_group = iface_grps[0]

        else:
            # This should already have been prevented above.
            return None

        opposite_b1s = [
            ball
            for ball in b1s_sorted
            if ball in opposite_group
        ]

        remaining_b1s = [
            ball
            for ball in b1s_sorted
            if ball not in opposite_group
        ]

        b1s_sorted = opposite_b1s + remaining_b1s

    # Preserve the original limit of five b1 candidates.
    b1s_sorted = b1s_sorted[:5]

    # Check each b1 candidate until an acceptable seed vertex is found.
    while len(b1s_sorted) > 0:
        if timed_out():
            return None
        b1 = b1s_sorted.pop(0)

        # Find the center between b0 and b1. Nearby b2 candidates are
        # searched around this location.
        b0_b1_com = calc_com(
            [locs[b0], locs[b1]]
        )

        b2s = []
        inc = 0
        b2_box = box_search(b0_b1_com)
        # Gather at least five b2 candidates near b0 and b1.
        while len(b2s) < 5:
            b2s = get_balls(b2_box, inc)
            inc += 1

        # b2 must be distinct from b0 and b1.
        b2s = [
            ball
            for ball in b2s
            if ball not in {b0, b1}
        ]

        # For interface mode, prioritize b2 candidates from a missing
        # interface side.
        #
        # This is only an ordering optimization. We do not remove the
        # other candidates because the fourth defining ball may be the
        # ball that completes the interface.
        if iface_grps is not None:
            current_balls = {b0, b1}

            has_group1 = bool(
                current_balls.intersection(iface_grps[0])
            )

            has_group2 = bool(
                current_balls.intersection(iface_grps[1])
            )

            if has_group1 and not has_group2:
                preferred_b2_group = iface_grps[1]

            elif has_group2 and not has_group1:
                preferred_b2_group = iface_grps[0]

            else:
                preferred_b2_group = None

            if preferred_b2_group is not None:
                preferred_b2s = [
                    ball
                    for ball in b2s
                    if ball in preferred_b2_group
                ]

                remaining_b2s = [
                    ball
                    for ball in b2s
                    if ball not in preferred_b2_group
                ]

                b2s = preferred_b2s + remaining_b2s

        # Construct valid three-ball circles from b0, b1, and each b2.
        my_circs = []

        for b2 in b2s:
            if timed_out():
                return None
            circle = [b0, b1, b2]

            circle_data = calc_circ(
                *[locs[ball] for ball in circle],
                *[rads[ball] for ball in circle],
            )

            if circle_data is not None:
                my_circs.append(
                    (circle, circle_data)
                )

        # Try the circles with the smallest absolute radius first.
        my_circs.sort(
            key=lambda x: abs(x[1][1])
        )

        for circ in my_circs:
            if timed_out():
                return None

            if net_type in {'prm', 'pow'}:
                my_vert = find_site_container(
                    circ[0],
                    locs=locs,
                    rads=rads,
                    b_verts=b_verts,
                    vert_ndxs=vert_ndxs,
                    max_vert=max_vert,
                    net_type=net_type,
                    group_ndxs=group_ndxs,
                    metrics=metrics,
                    box=box,
                )

            else:
                my_vert = find_site_container_slow(
                    circ[0],
                    locs=locs,
                    rads=rads,
                    b_verts=b_verts,
                    vert_ndxs=vert_ndxs,
                    max_vert=max_vert / 10,
                    net_type=net_type,
                    group_ndxs=group_ndxs,
                    iface_grps=iface_grps,
                    metrics=metrics,
                    box=box,
                )

            # No valid geometric site was found from this circle.
            if my_vert is None:
                continue

            # The site-container functions return the vertex and metrics.
            candidate = my_vert[0]

            # A seed vertex must have a valid primary location.
            if candidate['loc'] is None:
                continue

            # For an interface network, inspect the completed four-ball
            # vertex, not only the original three-ball circle.
            #
            # The candidate is accepted only when its defining balls
            # contain at least one ball from each interface group.
            if iface_grps is not None:
                candidate_balls = set(
                    candidate['balls']
                )

                belongs_to_interface = all(
                    bool(
                        candidate_balls.intersection(
                            group_indices
                        )
                    )
                    for group_indices in iface_grps
                )

                if not belongs_to_interface:
                    continue

            # Preserve the existing AW behavior: the initial seed should
            # not be a doublet.
            if (
                net_type == 'aw'
                and candidate.get('loc2') is not None
            ):
                continue

            # Reject a candidate whose primary location falls outside
            # the requested network box.
            if box is not None and any(
                box[0][axis] > candidate['loc'][axis]
                or candidate['loc'][axis] > box[1][axis]
                for axis in range(3)
            ):
                continue

            # This is a valid seed vertex. In interface mode, it has
            # already been confirmed to include both interface sides.
            return candidate

    # No acceptable seed vertex was found from this starting ball.
    return None


# still needs work


def find_v0_old(net, b_locs, b_rads, a0=None, group_balls=None):
    """
    Finds v0 using the ball finding functions to find a real verified site
    :param net: Network object to check from
    :param a0: The ball to seed from
    :param group_balls: List of balls for the building group based networks
    :return: V0 vertex
    """
    # Check to see if we need a group ball's box
    if a0 is not None:
        my_box = box_search(b_locs[a0])
    elif group_balls is not None:
        my_box = box_search(b_locs[group_balls[0]])
    else:
        # Find the middle sub_box of the set of boxes and
        mid = len(net.sub_boxes) // 2
        my_box = [mid, mid, mid]
    if a0 is None:
        a0s = []
        inc = 0
        # Keep grabbing balls until we have enough to get the current a0 increment
        while len(a0s) < 5:
            a0s = get_balls([my_box], inc)
            inc += 1
        # Pull an ball from the balls list
        a0 = a0s[0]
    a1s = []
    inc = 0
    # Get the 5 closest balls to a0
    while len(a1s) < 5:
        a1s = get_balls([my_box], inc)
        inc += 1
    # Set up the a2s lists
    a2s, j = [], 0
    # Check the a1s for verifiable
    while len(a1s) > 0:
        # Get the a1
        a1 = a1s.pop()
        # Add the circle check
        a2s.append([])
        inc = 0
        # Get the 20 closest balls to a0 and the current a1
        if len(b_locs) < 20:
            a2s[j] = [i for i in range(len(b_locs))]
        else:
            while len(a2s[j]) < 20:
                a2s[j] = get_balls([my_box], inc)
                inc += 1
        # Set up verified circles list for this a1
        verified_circles = []
        # Check each of the combinations for this a1
        for a2 in a2s[j]:
            # Use an edge object as a vehicle for calculating and verifying the inscribed circle
            circ = calc_circ(*[_.loc for _ in [a0, a1, a2]], *[_.rad for _ in [a0, a1, a2]])
            eloc, erad = None, None
            if circ is not None:
                eloc, erad = circ
            # If a circle can be made and the site does not overlap with any other balls, add it to the list
            if eloc is not None and erad < net.settings['max_vert'] and verify_site(eloc, erad, [a0, a1, a2], net,
                                                                                    net.settings['net_type']):
                verified_circles.append([a0, a1, a2])
        # Try to make a verified v0 site with the verified circles
        for circle in verified_circles:
            # Try to create a vertex
            my_vert = find_site(net, circle, group_ndxs=group_balls)
            # Check for a real site
            if my_vert is not None and my_vert[0].loc is not None:
                return my_vert[0]
        j += 1

