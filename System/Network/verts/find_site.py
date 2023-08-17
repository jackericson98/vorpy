from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
from System.Network.verts.verify_site import verify_site
from System.sys_funcs.calcs.calcs import box_search, get_atoms, calc_circ, calc_dist, ndx_search, rotate_points
from Visualize.mpl_visualize import plot_verts, plot_atoms
import matplotlib.pyplot as plt
import numpy as np
import time


def find_site(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, net_type, vn_1=None, vn_1_loc=None, group_atoms=None, metrics=None):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = edge_atoms[:]
    # Check if the edge contains a group atom or not
    check_atoms = True
    for ndx in edge_ndxs:
        if group_atoms is not None and ndx in group_atoms:
            check_atoms = False
            break
    # If the previous vertex has been provided, add the other atom to the not allowed atoms
    vert_atom_ndxs = vn_1
    if vn_1 is None:
        vert_atom_ndxs = edge_ndxs

    # Time printing metrics <-- Delete later
    start = time.perf_counter()
    # Grab the atoms we want to test against
    my_boxes = [box_search(loc=alocs[edge_atoms[_]]) for _ in range(3)]
    # Time printing metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()

    test_atoms = get_atoms(cells=my_boxes, dist=max_vert)

    # Get the center of the inscribed circle
    edge_center, edge_radius = calc_circ(alocs[edge_ndxs[0]], alocs[edge_ndxs[1]], alocs[edge_ndxs[2]], arads[edge_ndxs[0]], arads[edge_ndxs[1]], arads[edge_ndxs[2]])
    test_atom_tuples = []
    for atom in test_atoms:
        test_atom_tuples.append((atom, calc_dist(alocs[atom], np.array(edge_center)) - arads[atom]))
    test_atom_tuples.sort(key=lambda a: a[1])
    sorted_atoms = [_[0] for _ in test_atom_tuples]

    if metrics is not None:
        metrics['gather_atoms'] += time.perf_counter() - start
    # First look for vertices that have been found before
    new_test_atoms = []
    start = time.perf_counter()
    for atom in sorted_atoms:
        # If the atom is in the previous vertex move on
        if atom in vert_atom_ndxs:
            continue
        # Check if we need to check and if so check for the atom in the list
        if check_atoms and atom not in group_atoms:
            continue
        # If we have found the vertex before it is not the previous vertex return
        atom_ndxs = edge_ndxs + [atom]
        atom_ndxs.sort()
        # Get the vertex's index/insert index
        check_verts = [vert_ndxs[_] for _ in averts[atom_ndxs[0]]]
        my_vert_ndx = ndx_search(check_verts, atom_ndxs)
        if my_vert_ndx < len(check_verts) and atom_ndxs == check_verts[my_vert_ndx]:
            return
        new_test_atoms.append(atom)
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
    # Instantiate the vertex list and the size limit for vertices found
    verts = []
    # Go through each atom in the given test atoms. Extremely optimized
    for i, atom in enumerate(new_test_atoms):
        # Create the vertex and calculate its value
        vert_atoms = edge_atoms + [atom]
        vert_atoms.sort()
        vert_loc2, vert_rad2 = None, None
        # Calculate the correct vertex values
        start = time.perf_counter()
        if net_type == 'pow':
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms], power=True)
        elif net_type == 'del':
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms], power=False)
        else:
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(locs=[alocs[_] for _ in vert_atoms], rads=[arads[_] for _ in vert_atoms])
        if metrics is not None:
            metrics['calc_vert'] += time.perf_counter() - start
        # Catch the none location case
        if vert_loc is None:
            continue
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_atoms = [_ for _ in sorted_atoms if _ not in vert_atoms]
        test_locs = np.array([alocs[_] for _ in filtered_test_atoms])
        test_rads = np.array([arads[_] for _ in filtered_test_atoms])
        if abs(vert_rad) < max_vert and verify_site(loc=np.array(vert_loc), rad=vert_rad, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            if len(verts) > 0 and verts[0]['rad'] < vert_rad:
                return verts[0], metrics
            verts.append({'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad})
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if vert_loc2 is not None and abs(vert_rad2) < max_vert and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
                verts[-1]['loc2'], verts[-1]['rad2'] = vert_loc2, vert_rad2
        # Check to see if the doublet's site is verified
        elif vert_loc2 is not None and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            verts.append({'atoms': vert_atoms, 'loc': vert_loc2, 'rad': vert_rad2})
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
    # If no verts have been found return
    if len(verts) == 0:
        return
    # If we find only 1 vertex, return it
    elif len(verts) == 1 or verts[0]['rad'] < verts[1]['rad']:
        return verts[0], metrics
    return verts[1], metrics


def find_site_fast(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, net_type, vn_1, vn_1_loc, group_atoms=None,
                   metrics=None, vn_1_rad=None):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """

    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = edge_atoms[:]

    # Check if the edge contains a group atom or not. If not each atom paired with the edge atoms needs to be checked
    check_atoms = True
    for ndx in edge_ndxs:
        if group_atoms is not None and ndx in group_atoms:
            check_atoms = False
            break
    #extra_verts = [[357, 373, 375, 1338], [417, 419, 444, 2349], [233, 253, 1305, 1306], [419, 438, 2350, 2351], [63, 819, 1013, 1031]]
    extra_verts = [[311, 315, 316, 1379], [158, 207, 212, 225], [636, 644, 646, 824], [649, 650, 659, 668], [255, 257, 278, 284], [634, 636, 646, 824]]
    # Time printing metrics <-- Delete later
    start = time.perf_counter()
    # Grab the atoms we want to test against
    my_boxes = [box_search(loc=alocs[edge_atoms[_]]) for _ in range(3)]
    # Box search metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()

    # Gather the atoms surrounding the edge atoms
    test_atoms = get_atoms(cells=my_boxes, dist=max_vert)

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
        # Get the vertices for the first atom. All atoms will contain the vertex so only one atom needs to be checked
        check_verts = [vert_ndxs[_] for _ in averts[atom_ndxs[0]]]
        # Use the ndx_search function to quickly search the list of sorted vertices
        my_vert_ndx = ndx_search(check_verts, atom_ndxs)
        # If the index returned is larger than the list or the vertex at the index is not equal to the atom_ndxs were ok
        if my_vert_ndx < len(check_verts) and atom_ndxs == check_verts[my_vert_ndx]:
            return
        # Add the vertex indices to the test_vertices for calculation
        test_verts.append(atom_ndxs)

    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
        start = time.perf_counter()

    calculated_verts = []
    # Go through each atom in the given test atoms. Extremely optimized
    for i, vert in enumerate(test_verts):
        # Make sure the vertex values are defined
        vert_loc, vert_rad, vert_loc2, vert_rad2 = None, None, None, None

        # If the network type is power
        if net_type == 'pow':
            # Calculate the power vertex values
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert], rads=[arads[_] for _ in vert], power=True)
        # If the network type is Delaunay
        elif net_type == 'del':
            # Calculate the Delaunay vertex values
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert], rads=[arads[_] for _ in vert], power=False)
        # If the network type is Voronoi
        elif net_type == 'vor':
            # Calculate the Voronoi vertex values
            vert_loc, vert_rad, vert_loc2, vert_rad2 = calc_vert(locs=[alocs[_] for _ in vert], rads=[arads[_] for _ in vert])

        # Catch the none location and the too large vertex cases
        if vert_loc is None or vert_rad > max_vert:
            continue

        # Delete the second location for the vertex if it is too large
        if vert_rad2 is not None and vert_rad2 > max_vert:
            vert_loc2, vert_rad2 = None, None

        # # Check the vertex for overlap with the previous atom
        # if calc_dist(alocs[other_atom], np.array(vert_loc)) - arads[other_atom] - abs(vert_rad) < 0:
        #     continue

        # Add the vertex to the list of calculated vertices
        calculated_verts.append({'atoms': vert, 'loc': np.array(vert_loc), 'rad': vert_rad, 'loc2': vert_loc2,
                                 'rad2': vert_rad2})

    # Calculate vertices metrics <-- Delete later
    if metrics is not None:
        metrics['calc_vert'] += time.perf_counter() - start
        start = time.perf_counter()

    # If no vertices survived return
    if len(calculated_verts) == 0:
        return
    # If there is only one vertex left, no need to sort. Just verify it
    elif len(calculated_verts) == 1:
        return choose_vert(calculated_verts, test_atoms, alocs, arads, metrics, start, net_type)

    # Instantiate the left and right vertex lists
    filtered_verts_left, filtered_verts_right = [], []
    # Get the centers of the edge atoms
    c0, c1, c2 = [alocs[_] for _ in edge_ndxs]
    # Calculate the inscribed circle of the edge
    edge_center, edge_radius = calc_circ(c0, c1, c2, *[arads[_] for _ in edge_ndxs])

    # Calculate the edge normal  direction - take cross product of vector centers of edge atoms - a0 a1 X a1, a2
    edge_direction = np.cross(c0 - c1, c0 - c2)
    edge_normal = edge_direction / np.linalg.norm(edge_direction)

    # Calculate the projection of the previous vertex onto the edge normal (value) or edge_normal dot prev vert center
    pv_dist = np.dot(edge_normal, edge_center - vn_1_loc)
    plotting = False
    # Go through the calculated vertices made by the edge atoms and the surrounding atoms - filtering process
    printing = False
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
            flipped_vert = {'atoms': vert['atoms'], 'loc': vert['loc2'], 'rad': vert['rad'], 'd2pv': abs(pv_dist - vert_proj_dist), 'loc2': vert['loc'], 'rad2': vert['rad']}
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
    if printing:
        print([_['atoms'] for _ in filtered_verts_left])
        print([_['atoms'] for _ in filtered_verts_right])

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
                return choose_vert([filtered_verts_left[0]], test_atoms, alocs, arads, metrics, start, net_type)
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
                return choose_vert([filtered_verts_right[0]], test_atoms, alocs, arads, metrics, start, net_type)
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
        my_vert = choose_vert([left_neighbor], test_atoms, alocs, arads, metrics, start, net_type)

        if my_vert is not None:
            if my_vert[0]['atoms'] in extra_verts:
                print('\n\nleft neighbor: {}, coming from: {}, test_atoms used: {}'.format(my_vert[0]['atoms'], vn_1, test_atoms))
            return my_vert
    # Check the right neighbor vertex
    if right_neighbor is not None:
        my_vert = choose_vert([right_neighbor], test_atoms, alocs, arads, metrics, start, net_type)
        if my_vert is not None:
            if my_vert[0]['atoms'] in extra_verts:
                print('\n\nright neighbor: {}, coming from: {}, test_atoms used: {}'.format(my_vert[0]['atoms'], vn_1, test_atoms))
            return my_vert


def choose_vert(lr_verts, test_atoms, alocs, arads, metrics, start, net_type):
    # Get the first vertex in the sorted list
    my_vert = lr_verts.pop(0)
    # Get the atoms surrounding the vertex, not including the vertex atoms
    my_check_atoms = [_ for _ in test_atoms if _ not in my_vert['atoms']]
    # Gather the locations and radii of the atoms
    test_locs = np.array([alocs[_] for _ in my_check_atoms])
    test_rads = np.array([arads[_] for _ in my_check_atoms])
    # printing = False
    # if my_vert['atoms'] == [6, 7, 3143, 3145]:
    #     printing = True
    #     print(my_check_atoms[201])
    # Check the first location for the vertex
    if verify_site(np.array(my_vert['loc']), my_vert['rad'], test_locs, test_rads):
        # Check the second location if it exists, if it is within the allowed size range and if it is verified
        if my_vert['rad2'] is None or not verify_site(np.array(my_vert['loc2']), my_vert['rad2'], test_locs, test_rads):
            my_vert['loc2'], my_vert['rad2'] = None, None

        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        # Return what is left of the left vertex
        return my_vert, metrics

    # If the first site is unverified try the other vertex site
    elif my_vert['loc2'] is not None and verify_site(loc=np.array(my_vert['loc2']), rad=my_vert['rad2'],
                                                     test_locs=test_locs, test_rads=test_rads, net_type=net_type):
        # Reset the left_vert variable with the other location and return it
        my_vert = {'atoms': my_vert['atoms'], 'loc': my_vert['loc2'], 'rad': my_vert['rad2'], 'loc2': None,
                      'rad2': None}
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        return my_vert, metrics


def plot_vertex_2d(calc_verts, oa, pv_loc, edge_atom_locs, edge_atom_rads, edge_normal, real_vert):
    my_edge = calc_circ(*edge_atom_locs, *edge_atom_rads)
    chosen_vert = [_ for _ in calc_verts if _['atoms'] == [6, 7, 3143, 3145]][0]

    # Calculate the edge plane
    ep_norm = np.cross(pv_loc - my_edge[0], edge_normal)
    new_chosen_vert_point = rotate_points(ep_norm, np.array([chosen_vert['loc']]))
    new_cv_points = rotate_points(ep_norm, np.array([_['loc'] for _ in calc_verts]))
    new_pv_point = rotate_points(ep_norm, np.array([pv_loc]))
    new_oa_point = rotate_points(ep_norm, np.array([oa]))
    new_real_vert_loc = rotate_points(ep_norm, np.array([real_vert['loc']]))
    plt.scatter([_[0] for _ in new_cv_points], [_[1] for _ in new_cv_points])
    plt.scatter([new_pv_point[0][0]], [new_pv_point[0][1]], marker='o', s=20)
    plt.scatter([new_oa_point[0][0]], [new_oa_point[0][1]])
    plt.scatter([new_real_vert_loc[0][0]], [new_real_vert_loc[0][1]], marker='x', s=20)
    plt.scatter([new_chosen_vert_point[0][0]], [new_chosen_vert_point[0][1]], marker='.', s=40)

    plt.show()


def plot_vert_situation(edge_atoms, my_vert, vn_1_loc, vn_1_rad, alocs, arads, a0=None, a1=None, a2=None, v0=None):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        # actual_vert_other_atom = [_ for _ in actual_vert['atoms'] if _ not in edge_ndxs][0]
        # edge atoms
        plot_atoms([alocs[_] for _ in edge_atoms], [arads[_] for _ in edge_atoms], fig=fig, ax=ax, colors=['r', 'r', 'r'])
        # other atom
        if a0 is not None:
            plot_atoms([alocs[a0]], [arads[1779]], fig=fig, ax=ax, colors=['pink'])
        # interfering atom
        if a1 is not None:
            plot_atoms([alocs[3144]], [arads[3144]], fig=fig, ax=ax, colors=['orange'])
        # # actual vert other atom
        if a2 is not None:
            plot_atoms([alocs[7]], [arads[7]], fig=fig, ax=ax, colors=['purple'])
        # actual vert
        plot_verts([my_vert['loc2']], [my_vert['rad2']], fig=fig, ax=ax, spheres=True, colors=['b'])
        # closest vert
        if v0 is not None:
            plot_verts([v0['loc']], [v0['rad']], fig=fig, ax=ax, spheres=True, colors=['white'])
        # previous vert
        plot_verts([vn_1_loc], [vn_1_rad], fig=fig, ax=ax, spheres=True, colors=['green'], Show=True)


