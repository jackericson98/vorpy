from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
from System.Network.verts.verify_site import verify_site
from System.sys_funcs.calcs.calcs import calc_dist, calc_com
from System.sys_funcs.calcs.sorting import box_search, get_atoms
from System.sys_funcs.calcs.circle import calc_circ
import bisect
import numpy as np
import time


def find_site_container(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, net_type, vn_1=None, vn_1_loc=None,
                        group_atoms=None, metrics=None, vn_1_rad=None, printing=False):
    """
    Cycles through larger and larger areas searching for
    """
    # Set up the vert and invalid atoms parameters
    invalid_atoms, vert = [], None

    # If no vn_1 is provided set it to the edge_atoms
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
    # Se the initial vert size
    mv_inc = 0.45
    # Look for the vert and keep increasing box size until the vert is found
    while vert is None and mv_inc < max_vert:
        # Search for the vertx in the current range
        if net_type == 'vor':
            vert, invalid_atoms = find_site_vor(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type,
                                                check_atoms, surr_atoms, my_boxes, invalid_atoms, vn_1,
                                                vn_1_loc, group_atoms=group_atoms, metrics=metrics, printing=printing)
        else:
            vert, invalid_atoms = find_site_pd(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc,
                                               net_type, check_atoms, surr_atoms, my_boxes, invalid_atoms, vn_1,
                                               vn_1_loc, group_atoms=group_atoms, metrics=metrics)
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
    if vn_1_loc is None:
        vn_1_loc = calc_com([alocs[_] for _ in edge_ndxs])

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
        # Calculate the 181L vertex values
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


def find_site_vor(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type, check_atoms, surr_atoms,
                  my_boxes, invalid_atoms, vn_1, vn_1_loc, group_atoms=None, metrics=None, printing=False):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = edge_atoms[:]
    # Time printing metrics <-- Delete later
    start = time.perf_counter()

    # Box search metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()

    # Get the atoms not in the invalid atoms that are within the range specified
    test_atoms = [_ for _ in get_atoms(cells=my_boxes, dist=mv_inc) if _ not in invalid_atoms]

    # Sort the test atoms to be in order by distance from the previous vert location
    if net_type != 'vor' and vn_1_loc is not None:
        dists = [calc_dist(np.array(alocs[_]), np.array(vn_1_loc)) for _ in test_atoms]
        test_atoms = [_ for x, _ in sorted(zip(dists, test_atoms))]

    # Gather atoms metrics <-- Delete later
    if metrics is not None:
        metrics['gather_atoms'] += time.perf_counter() - start
        start - time.perf_counter()

    # Instantiate the list for test vertices to be calculated later. This saves us from sorting the vertices atoms twice
    new_test_atoms = []
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
        # Get the vertices for the first atom. All atoms will contain the vertex so only one atom needs to be checked
        check_verts = [vert_ndxs[_] for _ in averts[atom_ndxs[0]]]
        # Use the ndx_search function to quickly search the list of sorted vertices
        my_vert_ndx = bisect.bisect_left(check_verts, atom_ndxs)
        # If the index returned is larger than the list or the vertex at the index is not equal to the atom_ndxs were ok
        if my_vert_ndx < len(check_verts) and atom_ndxs == check_verts[my_vert_ndx]:
            return None, invalid_atoms
        # Add the vertex indices to the test_vertices for calculation
        new_test_atoms.append(atom)

    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
        start = time.perf_counter()

    # Instantiate the calculated vertices list
    calculated_verts = []
    # Go through each atom in the given test atoms. Extremely optimized
    for i, atom in enumerate(new_test_atoms):

        # Combine the new atom with the edge atoms and sort
        vert_atoms = edge_atoms + [atom]
        vert_atoms.sort()
        # Make sure the vertex values are defined
        vert_loc, vert_rad, vert_loc2, vert_rad2 = None, None, None, None
        # If the network type is power
        if net_type == 'pow':
            # Calculate the power vertex values
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms], power=True)
        # If the network type is Delaunay
        elif net_type == 'del':
            # Calculate the Delaunay vertex values
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms], power=False)
        # If the network type is Voronoi
        elif net_type == 'vor':
            # Calculate the Voronoi vertex values
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms])
        # Catch the none location and the too large vertex cases
        if vert_loc is None or vert_rad > max_vert:
            continue

        # Delete the second location for the vertex if it is too large
        if vert_rad2 is not None and vert_rad2 > max_vert:
            vert_loc2, vert_rad2 = None, None

        # Add the vertex to the list of calculated vertices
        calculated_verts.append({'atoms': vert_atoms, 'loc': np.array(vert_loc), 'rad': vert_rad, 'loc2': vert_loc2,
                                 'rad2': vert_rad2})

    # Calculate vertices metrics <-- Delete later
    if metrics is not None:
        metrics['calc_vert'] += time.perf_counter() - start
        start = time.perf_counter()

    # If no vertices survived return
    if len(calculated_verts) == 0:
        return None, invalid_atoms
    # If there is only one vertex left, no need to sort. Just verify it
    elif len(calculated_verts) == 1:
        my_doopy = choose_vert(calculated_verts[0], edge_ndxs, surr_atoms, alocs, arads, metrics, start, net_type)[0], invalid_atoms
        return my_doopy

    # Instantiate the left and right vertex lists
    filtered_verts_left, filtered_verts_right = [], []
    # Get the centers of the edge atoms
    c0, c1, c2 = [alocs[_] for _ in edge_ndxs]
    # Get the center of the inscribed circle
    edge_center, edge_radius = calc_circ(alocs[edge_ndxs[0]], alocs[edge_ndxs[1]], alocs[edge_ndxs[2]],
                                         arads[edge_ndxs[0]], arads[edge_ndxs[1]], arads[edge_ndxs[2]])

    # Calculate the edge normal  direction - take cross product of vector centers of edge atoms - a0 a1 X a1, a2
    edge_direction = np.cross(c0 - c1, c0 - c2)
    edge_normal = edge_direction / np.linalg.norm(edge_direction)

    # Calculate the projection of the previous vertex onto the edge normal (value) or edge_normal dot prev vert center
    pv_dist = np.dot(edge_normal, edge_center - vn_1_loc)
    # Go through the calculated vertices made by the edge atoms and the surrounding atoms - filtering process
    for vert in calculated_verts:
        # Get the vertex's projected distance
        vert_proj_dist = np.dot(edge_normal, edge_center - vert['loc'])
        # Calculate the distance to the previous vertex and assign it as a value in the vertex dictionary
        vert['d2pv'] = abs(pv_dist - vert_proj_dist)

        # If the other atoms projection (value1) is less than the previous vertex's projection (value)
        if pv_dist < vert_proj_dist:
            # Add the vertex to the list of filtered vertices
            filtered_verts_left.append(vert)
        else:
            # Add the vertex to the list of filtered vertices
            filtered_verts_right.append(vert)

        vert['d2pv2'] = None
        if vert['loc2'] is not None:
            vert_proj_dist = np.dot(edge_normal, edge_center - vert['loc2'])
            flipped_vert = {'atoms': vert['atoms'], 'loc': vert['loc2'], 'rad': vert['rad2'], 'd2pv': abs(pv_dist - vert_proj_dist), 'loc2': vert['loc'], 'rad2': vert['rad']}
            # If the other atoms projection (value1) is less than the previous vertex's projection (value)
            if pv_dist < vert_proj_dist:
                # Add the vertex to the list of filtered vertices
                filtered_verts_left.append(flipped_vert)
            else:
                # Add the vertex to the list of filtered vertices
                filtered_verts_right.append(flipped_vert)

    # Sort the filtered vertices by distance to the previous vertex
    filtered_verts_left.sort(key=lambda my_vert: my_vert['d2pv'])
    filtered_verts_right.sort(key=lambda my_vert: my_vert['d2pv'])

    # Set up the left neighbor and the right neighbor variables for assignment
    left_neighbor, right_neighbor = None, None
    # If all vertices lie on the left side of the previous vertex
    if len(filtered_verts_right) == 0:
        # Get the leftmost vertex and the rightmost vertex
        vl, vr = filtered_verts_left[-1]['loc'], vn_1_loc
        # Counter variable
        i = 0
        # Loop through the vertices looking for the left and right neighbor
        while (left_neighbor is None or right_neighbor is None) and i < len(filtered_verts_left) - 1:
            # Grab the current vertex in the loop
            vi = filtered_verts_left[i]
            # Calculate the determinant of the vertex and the left most and right most vertices
            my_det = np.linalg.det([vl, vr, vi['loc']])
            # If the edge is straight, verify/return the leftmost vertex on the right
            if my_det == 0:
                # Verification
                my_doopy =choose_vert(filtered_verts_left[0], edge_ndxs, surr_atoms, alocs, arads, metrics, start, net_type)[0], invalid_atoms
                return my_doopy
            # If the vertex falls in the lower hull it is the left neighbor
            elif my_det > 0 and left_neighbor is None:
                left_neighbor = vi
            # If the vertex falls in the upper hull it is the right neighbor
            elif my_det < 0 and right_neighbor is None:
                right_neighbor = vi
            # Increment the counter
            i += 1
        if left_neighbor is None:
            left_neighbor = filtered_verts_left[-1]
        elif right_neighbor is None:
            right_neighbor = filtered_verts_left[-1]
    # If all vertices lie on the right side of the previous vertex
    elif len(filtered_verts_left) == 0:
        # Get the leftmost vertex and the rightmost vertex
        vr, vl = filtered_verts_right[-1]['loc'], vn_1_loc
        # Counter variable
        i = 0
        # Loop through the vertices looking for the left and right neighbor
        while (left_neighbor is None or right_neighbor is None) and i < len(filtered_verts_right) - 1:
            # Grab the current vertex in the loop
            vi = filtered_verts_right[i]
            # Calculate the determinant of the vertex and the left most and right most vertices
            my_det = np.linalg.det([vl, vr, vi['loc']])
            # If the edge is straight, verify/return the leftmost vertex on the right
            if my_det == 0:
                # Verification
                return choose_vert(filtered_verts_right[0], edge_ndxs, surr_atoms, alocs, arads, metrics, start, net_type)[0], invalid_atoms
            # If the vertex falls in the upper hull it is the left neighbor
            elif my_det < 0 and left_neighbor is None:
                left_neighbor = vi
            # If the vertex falls in the lower hull it is the right neighbor
            elif my_det > 0 and right_neighbor is None:
                right_neighbor = vi
            # Increment the counter
            i += 1

        if left_neighbor is None:
            left_neighbor = filtered_verts_right[-1]
        elif right_neighbor is None:
            right_neighbor = filtered_verts_right[-1]
    # If there are vertices on either side
    else:
        # Find the left most and right most vertices
        vl, vr = filtered_verts_left[-1], filtered_verts_right[-1]
        vert_det = np.linalg.det([vl['loc'], vr['loc'], vn_1_loc])
        # Assign the left and right neighbor variables
        left_neighbor, right_neighbor = None, None
        # Counter variable
        i = 0
        # Go through the vertices on the left og the vertex
        while left_neighbor is None and i < len(filtered_verts_left):
            # Get the current vertex in the loop
            vi = filtered_verts_left[i]
            # Calculate the determinant of the left most, right most and current vertex
            my_det = np.linalg.det([vl['loc'], vr['loc'], vi['loc']])
            # If they share a sign, we have found the vertex
            if my_det <= 0 and vert_det <= 0 or my_det >= 0 and vert_det >= 0:
                left_neighbor = vi
            # Increment the counter
            i += 1
        # Reset the counter variable
        i = 0
        # Go through the vertices on the right of the previous vertex
        while right_neighbor is None and i < len(filtered_verts_right):
            # Get the current vertex in the loop
            vi = filtered_verts_right[i]
            # Calculate the determinant of the left most, right most and current vertex
            my_det = np.linalg.det([vl['loc'], vr['loc'], vi['loc']])
            # If they share a sign, we have found the vertex
            if my_det <= 0 and vert_det <= 0 or my_det >= 0 and vert_det >= 0:
                right_neighbor = vi
            # Increment the counter
            i += 1

    # Check the left neighbor vertex
    if left_neighbor is not None:
        my_vert, extra_atom = choose_vert(left_neighbor, edge_ndxs, surr_atoms, alocs, arads, metrics, start, net_type)

        if my_vert is not None:
            return my_vert, invalid_atoms
        invalid_atoms.append(extra_atom)
    # Check the right neighbor vertex
    if right_neighbor is not None:
        my_vert, extra_atom = choose_vert(right_neighbor, edge_ndxs, surr_atoms, alocs, arads, metrics, start, net_type)
        if my_vert is not None:
            return my_vert, invalid_atoms
        invalid_atoms.append(extra_atom)
    return None, invalid_atoms


def choose_vert(my_vert, edge_ndxs, test_atoms, alocs, arads, metrics, start, net_type):
    # Create the extra atom variable
    extra_atom = None
    # Get the atoms surrounding the vertex, not including the vertex atoms
    my_check_atoms = [_ for _ in test_atoms if _ not in my_vert['atoms']]
    # Gather the locations and radii of the atoms
    test_locs = np.array([alocs[_] for _ in my_check_atoms])
    test_rads = np.array([arads[_] for _ in my_check_atoms])
    # Check the first location for the vertex
    if verify_site(np.array(my_vert['loc']), my_vert['rad'], test_locs, test_rads):
        # Check the second location if it exists, if it is within the allowed size range and if it is verified
        if my_vert['rad2'] is None or not verify_site(np.array(my_vert['loc2']), my_vert['rad2'], test_locs, test_rads):
            my_vert['loc2'], my_vert['rad2'] = None, None

        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        # Return what is left of the left vertex
        return [my_vert, metrics], extra_atom

    # If the first site is unverified try the other vertex site
    elif my_vert['loc2'] is not None and verify_site(loc=np.array(my_vert['loc2']), rad=my_vert['rad2'],
                                                     test_locs=test_locs, test_rads=test_rads, net_type=net_type):
        # Reset the left_vert variable with the other location and return it
        my_vert = {'atoms': my_vert['atoms'], 'loc': my_vert['loc2'], 'rad': my_vert['rad2'], 'loc2': None,
                      'rad2': None}
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        return [my_vert, metrics], extra_atom
    # We still need to return invalid atoms if they are not included
    return None, [_ for _ in my_vert['atoms'] if _ not in edge_ndxs][0]

