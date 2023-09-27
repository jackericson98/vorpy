from System.sys_funcs.calcs.calcs import calc_circ, box_search, get_atoms
from System.Network.verts.find_site.slow import find_site
from System.Network.verts.verify_site import verify_site
import numpy as np


def find_v0(alocs, arads, averts, max_vert, net_type, a0=None, group_atoms=None, metrics=None):
    """
    Finds v0 using the atom finding functions to find a real verified site
    """
    # Check to see if we need a group atom's box
    if a0 is not None:
        my_box = box_search(alocs[a0])
    elif group_atoms is not None:
        # Get the box for the
        my_box = box_search(alocs[group_atoms[0]])
    else:
        # Find the middle sub_box of the set of boxes and
        my_box = box_search(alocs[0])
    # Get the atoms in the given box and see if it's empty
    atoms = get_atoms([my_box])
    if len(atoms) > 0:
        a0 = atoms[0]
    # If we still haven't found an a0
    if a0 is None:
        # Reset the a0 variables
        a0s = []
        inc = 0
        # Keep searching boxes until we find an atom
        while len(a0s) < 1:
            a0s = get_atoms([my_box], inc)
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
    # Set up the a2s lists
    a2s, j = [], 0
    # Check the a1s for verifiable
    while len(a1s) > 0:
        # Get the a1
        a1 = a1s.pop()
        # Add the circle check
        a2s.append([])
        inc = 0
        # Get the 20 closest atoms to a0 and the current a1
        if len(alocs) < 20:
            a2s[j] = [_ for _ in range(len(alocs))]
        else:
            while len(a2s[j]) < 20:
                a2s[j] = get_atoms([my_box], inc)
                inc += 1
        # Set up verified circles list for this a1
        verified_circles = []
        # Check each of the combinations for this a1
        for a2 in a2s[j]:
            # Use an edge object as a vehicle for calculating and verifying the inscribed circle
            circ = calc_circ(alocs[a0], alocs[a1], alocs[a2], arads[a0], arads[a1], arads[a2])
            eloc, erad = None, None
            if circ is not None:
                eloc, erad = circ
            # Set up a list of atoms to test our edge atoms with
            # max_inc = int(max_vert / min(sub_box_size) - max_atom_rad) + 1
            # Grab the atoms we want to test against
            my_boxes = [box_search(loc=alocs[a2])]
            test_atoms = get_atoms(cells=my_boxes, dist=max_vert)
            # If a circle can be made and the site does not overlap with any other atoms, add it to the list
            filtered_test_atoms = [_ for _ in test_atoms if _ not in [a0, a1, a2]]
            test_locs = [alocs[_] for _ in filtered_test_atoms]
            test_rads = [arads[_] for _ in filtered_test_atoms]
            if eloc is not None and erad < max_vert and verify_site(loc=np.array(eloc), rad=erad, test_locs=np.array(test_locs), test_rads=np.array(test_rads), net_type=net_type):
                verified_circles.append([a0, a1, a2])
        # Try to make a verified v0 site with the verified circles
        for circle in verified_circles:
            # Try to create a vertex
            my_vert = find_site(circle, alocs=alocs, arads=arads, averts=averts, vert_ndxs=[], max_vert=max_vert, net_type=net_type, group_atoms=group_atoms, metrics=metrics)
            # Check for a real site
            if my_vert is not None and my_vert[0]['loc'] is not None:
                return my_vert[0]
        j += 1
