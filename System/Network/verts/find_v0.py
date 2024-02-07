from System.sys_funcs.calcs.calcs import calc_com, calc_dist
from System.sys_funcs.calcs.sorting import box_search, get_atoms
from System.sys_funcs.calcs.circle import calc_circ
from System.Network.verts.find_site.slow import find_site
from System.Network.verts.find_site.fast import find_site_container
from System.Network.verts.verify_site import verify_site
import numpy as np


def find_v0(alocs, arads, averts, max_vert, net_type, a0=None, group_atoms=None, metrics=None, vert_ndxs=None,
            group_box=None):
    """
    Finds v0 using the atom finding functions to find a real verified site
    """
    if vert_ndxs is None:
        vert_ndxs = []
    # Check to see if we need a group atom's box
    if a0 is not None:
        my_box = box_search(alocs[a0])
    elif group_box is not None:
        # Get the center of the group_box
        my_box = box_search([0.5 * abs(group_box[1][i] - group_box[0][i]) + group_box[0][i] for i in range(3)])
    elif group_atoms is not None:
        # Get the box for the
        my_box = box_search(alocs[group_atoms[int(len(group_atoms) / 2)]])
    else:
        # Find the middle sub_box of the set of boxes and
        my_box = box_search(alocs[int(len(alocs) / 2)])

    # If we still haven't found an a0
    if a0 is None:
        # Reset the a0 variables
        a0s = []
        inc = 0
        # Keep searching boxes until we find an atom
        while len(a0s) < 1:
            a0s = get_atoms([my_box], inc)
            if group_atoms is not None:
                a0s = [_ for _ in a0s if _ in group_atoms and len(averts[_]) == 0]
            inc += 1
        # Pull an atom from the atoms list
        a0 = a0s[0]
    # Reset the a1 variables
    a1s = []
    inc = 0
    # Get the 5 closest atoms to a0
    while len(a1s) < 5:
        a1s = get_atoms([my_box], inc)
        inc += 1
    # Sort the a1s
    a1_dists = [calc_dist(alocs[a1], alocs[a0]) - (arads[a0] + arads[a1]) for a1 in a1s]
    _, a1s_sorted = zip(*sorted(zip(a1_dists, a1s), key=lambda x: x[0]))
    a1s_sorted = [_ for _ in a1s_sorted if _ != a0][:5]
    # Set up the a2s lists
    a2s, j = [], 0
    # Check the a1s for verifiable
    while len(a1s_sorted) > 0:
        # Get the a1
        a1 = a1s_sorted.pop(0)
        # Find the center of mass for a0 and a1 locations
        a0_a1_com = calc_com([alocs[a0], alocs[a1]])

        inc = 0
        # Find a2s near a0 and a1
        while len(a2s) < 5:
            a2s = get_atoms(box_search(a0_a1_com), inc)
            inc += 1
        a2s = [_ for _ in a2s if _ not in {a0, a1}]

        my_circs = []
        # Check each of the combinations for this a1
        for a2 in a2s:
            # Set up the circle
            circle = [a0, a1, a2]
            my_circs.append((circle, calc_circ(*[alocs[_] for _ in circle], *[arads[_] for _ in circle])))
        my_circs.sort(key=lambda x: abs(x[1][1]))
        for circ in my_circs:
            circ[0].sort()
            # Try to create a vertex
            if net_type in ['del', 'pow']:
                my_vert = find_site_container(circ[0], alocs=alocs, arads=arads, averts=averts, vert_ndxs=vert_ndxs,
                                                 max_vert=max_vert, net_type=net_type, group_atoms=group_atoms,
                                                 metrics=metrics)
            else:

                my_vert, invalid_atoms = find_site(circ[0], alocs=alocs, arads=arads, averts=averts, vert_ndxs=vert_ndxs,
                                                   max_vert=max_vert, mv_inc=max_vert, net_type=net_type,
                                                   group_atoms=group_atoms, metrics=metrics)
            # Check for a real site that is not a doublet
            if my_vert is not None:
                if net_type == 'vor':
                    if my_vert[0]['loc'] is not None and my_vert[0]['loc2'] is None:
                        return my_vert[0]
                return my_vert[0]
        j += 1
