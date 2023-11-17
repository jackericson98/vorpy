from System.sys_funcs.calcs.calcs import box_search, get_atoms, calc_dist
from System.Network.verts.verify_site import verify_site
from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
import bisect
import numpy as np
import time


def find_site_pd_container(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, net_type, vn_1=None, vn_1_loc=None,
                           group_atoms=None, metrics=None):
    """
    Cycles through larger and larger areas searching for
    """
    # Set up the vert and invalid atoms parameters
    invalid_atoms, vert = [], None
    # If the previous vertex has been provided, add the other atom to the not allowed atoms
    if vn_1 is None:
        vn_1 = edge_atoms

    # Check if the edge contains a group atom, to see if the next atom needs to be checked or not
    # Start with check atoms as false if no group is defined
    check_atoms = False
    if group_atoms is not None:
        # If a group exists default to checking each atom
        check_atoms = True
        # Go through the edge atoms checking if they are in the group --> any vert found from another atom is included
        for atom in edge_atoms:
            # Take the potential index of the atom in group
            my_index = bisect.bisect_left(group_atoms, atom)
            # If the index is in the list check if the atom matches the index's element
            if my_index != len(group_atoms) and group_atoms[my_index] == atom:
                # If the element is found no need to check the atoms and break the for loop
                check_atoms = False
                break

    # Find the 3 boxes the edge atoms are in
    my_boxes = [box_search(loc=alocs[edge_atoms[_]]) for _ in range(3)]
    # Gather the surrounding atoms or the entire list of atoms we could be comparing to
    surr_atoms = get_atoms(cells=my_boxes, dist=max_vert)
    # Set the initial vert size
    mv_inc = 0.05
    # Look for the vertex and keep increasing box size until the vert is found
    while mv_inc < max_vert:
        # Search for the vertx in the current range
        vert, invalid_atoms = find_site_pd(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type,
                                           check_atoms, surr_atoms, my_boxes, invalid_atoms, vn_1,
                                           group_atoms=group_atoms, metrics=metrics, vn_1_loc=vn_1_loc)
        # If a vertex is found exit the loop
        if vert is not None:
            break
        # Increment the range for the search
        mv_inc *= 10
    # Return the vertex if found
    return vert


def find_site_pd(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type, check_atoms, surr_atoms,
                 my_boxes, invalid_atoms, vn_1, vn_1_loc=None, group_atoms=None, metrics=None):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = edge_atoms[:]

    # Time printing metrics <-- Delete later
    start = time.perf_counter()

    # Time printing metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()

    # Get the atoms not in the invalid atoms that are within the range specified
    test_atoms = [_ for _ in get_atoms(cells=my_boxes, dist=mv_inc) if _ not in invalid_atoms]

    # Sort the test atoms to be in order by distance from the previous vert location
    if vn_1_loc is not None:
        dists = [calc_dist(np.array(alocs[_]), np.array(vn_1_loc)) for _ in test_atoms]
        test_atoms = [_ for x, _ in sorted(zip(dists, test_atoms))]

    # Gather atoms metrics <-- Delete later
    if metrics is not None:
        metrics['gather_atoms'] += time.perf_counter() - start
        start = time.perf_counter()

    # Instantiate the list for test vertices to be calculated later. This saves us from sorting the vertices atoms twice
    test_verts = []
    # Go through the surrounding atoms to look for vertices that have been found before and filter out edge atoms
    for atom in test_atoms:
        # If the atom is in the previous vertex move on
        if atom in vn_1:
            continue
        # Check if we need to check and if so check for the atom in the list
        if check_atoms and atom not in group_atoms:
            continue
        # If we have found the vertex before it is not the previous vertex return
        atom_ndxs = edge_ndxs + [atom]
        atom_ndxs.sort()
        # Get the vertex's index/insert index
        check_verts = [vert_ndxs[_] for _ in averts[atom_ndxs[0]]]
        # Take the potential index of the atom in group
        my_vert_ndx = bisect.bisect_left(check_verts, atom_ndxs)
        # If the index returned is larger than the list or the vertex at the index is not equal to the atom_ndxs were ok
        if my_vert_ndx < len(check_verts) and atom_ndxs == check_verts[my_vert_ndx]:
            return None, invalid_atoms
        # Add the vertex indices to the test_vertices for calculation
        test_verts.append((atom_ndxs, atom))

    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start

    # Go through each atom in the given test atoms. Extremely optimized
    for i, vert in enumerate(test_verts):

        # Add the vertex atom to the
        vert_atoms, atom = vert
        # Calculate the correct vertex values
        start = time.perf_counter()
        if net_type == 'pow':
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms], power=True)
        elif net_type == 'del':
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms], power=False)
        else:
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms])

        # Record the calculate vertex metrics
        if metrics is not None:
            metrics['calc_vert'] += time.perf_counter() - start

        # Catch the none location case
        if vert_loc is None:
            invalid_atoms.append(atom)
            continue

        # Restart ste start time to only record verify site time to the verify site metrics
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_atoms = [_ for _ in surr_atoms if _ not in vert_atoms]
        # Get the locations from the test atoms
        test_locs = np.array([alocs[_] for _ in filtered_test_atoms])
        test_rads = np.array([arads[_] for _ in filtered_test_atoms])
        # Compare the vertex to the maximum allowed vertex and verify it
        if vert_rad < max_vert and verify_site(loc=np.array(vert_loc), rad=vert_rad, test_locs=test_locs,
                                               test_rads=test_rads, net_type=net_type):
            # Add the time for verification to the verify_site metrics
            if metrics is not None:
                metrics['verify_site'] += time.perf_counter() - start
            # Return the validated atom and the invalidated ist
            return [{'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad}, metrics], invalid_atoms
        else:
            # Add the atom to the invalid atoms list if it isn't verified
            invalid_atoms.append(atom)
    # Return the non-vertex and invalid atoms
    return None, invalid_atoms
