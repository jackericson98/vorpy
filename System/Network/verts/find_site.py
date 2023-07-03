from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
from System.Network.verts.verify_site import verify_site
from System.sys_funcs.calcs.calcs import box_search, get_atoms, ndx_search, calc_circ, calc_dist
import numpy as np
import time


def find_site(edge_atoms, alocs, arads, vert_ndxs, max_vert, net_type, vn_1=None, group_atoms=None, metrics=None):
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
    vert_ndx_list_locs = []
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
        vert_ndx = ndx_search(vert_ndxs, atom_ndxs)
        # If the vertex has been found before return
        if vert_ndx < len(vert_ndxs) and vert_ndxs[vert_ndx] == atom_ndxs:
            return
        # else add the vert_ndx to the list
        vert_ndx_list_locs.append(vert_ndx)
        new_test_atoms.append(atom)
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
    # Instantiate the vertex list and the size limit for vertices found
    verts = []
    new_vert_ndx_list_locs = []
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
                return verts[0], vert_ndx_list_locs[0], metrics
            verts.append({'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad})
            new_vert_ndx_list_locs.append(vert_ndx_list_locs[i])
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if vert_loc2 is not None and abs(vert_rad2) < max_vert and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
                verts[-1]['loc2'], verts[-1]['rad2'] = vert_loc2, vert_rad2
                new_vert_ndx_list_locs.append(vert_ndx_list_locs[i] + 1)
        # Check to see if the doublet's site is verified
        elif vert_loc2 is not None and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            verts.append({'atoms': vert_atoms, 'loc': vert_loc2, 'rad': vert_rad2})
            new_vert_ndx_list_locs.append(vert_ndx_list_locs[i])
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
    # If no verts have been found return
    if len(verts) == 0:
        return
    # If we find only 1 vertex, return it
    elif len(verts) == 1 or verts[0]['rad'] < verts[1]['rad']:
        return verts[0], new_vert_ndx_list_locs[0], metrics
    return verts[1], new_vert_ndx_list_locs[1], metrics


def find_site_fast(edge_atoms, alocs, arads, vert_ndxs, max_vert, net_type, vn_1, group_atoms=None, metrics=None):
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

    # get the edge normal
    edge_normal = np.cross(alocs[edge_ndxs[0]] - alocs[edge_ndxs[1]], alocs[edge_ndxs[1]] - alocs[edge_ndxs[2]])
    edge_normal = edge_normal / np.linalg.norm(edge_normal)
    # Get the direction of the normal based on the previous atom
    other_atom = [_ for _ in vn_1 if _ not in edge_ndxs][0]
    atom_dir = np.array(edge_center) - alocs[other_atom]
    atom_dist = np.dot(atom_dir, edge_normal)
    if atom_dist < 0:
        edge_normal = - edge_normal
        atom_dist = - atom_dist

    if metrics is not None:
        metrics['gather_atoms'] += time.perf_counter() - start
    # First look for vertices that have been found before
    vert_ndx_list_locs = []
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
        vert_ndx = ndx_search(vert_ndxs, atom_ndxs)
        # If the vertex has been found before return
        if vert_ndx < len(vert_ndxs) and vert_ndxs[vert_ndx] == atom_ndxs:
            return
        # else add the vert_ndx to the list
        vert_ndx_list_locs.append(vert_ndx)
        new_test_atoms.append(atom)
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
    filtered_test_verts = []
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
        if vert_loc is None or vert_rad > max_vert:
            continue
        # print(edge_center, vert_loc, edge_normal, loc_dist)
        if vert_atoms == [1, 22, 135, 532]:
            overlap = calc_dist(alocs[other_atom], np.array(vert_loc)) - arads[other_atom] - abs(vert_rad)
            print(edge_ndxs, atom_dist, loc_dist, overlap)
        if calc_dist(alocs[other_atom], np.array(vert_loc)) - arads[other_atom] - abs(vert_rad) < 0:
            continue

        # Get the projected distance
        loc_dist = np.dot(np.array(edge_center) - np.array(vert_loc), edge_normal)

        filtered_test_verts.append({'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad, 'loc2': vert_loc2, 'rad2': vert_rad2, 'ndx_list_loc': vert_ndx_list_locs[i], 'loc_dist': loc_dist})
    # print(filtered_test_verts)
    filtered_test_verts.sort(key=lambda a: abs(atom_dist - a['loc_dist']), reverse=True)
    if len(filtered_test_verts) > 0:
        wrong_verts = []
        for vert in reversed(filtered_test_verts):
            if vert['atoms'] not in [[0, 1, 3, 4], [0, 1, 3, 6], [0, 1, 4, 8], [0, 1, 6, 8], [1, 3, 4, 14], [1, 3, 6, 24], [1, 3, 14, 28], [1, 3, 20, 22], [1, 3, 20, 28], [1, 3, 22, 24], [1, 4, 8, 18], [1, 4, 14, 20], [1, 4, 18, 30], [1, 4, 20, 30], [1, 6, 8, 21], [1, 6, 21, 24], [1, 8, 18, 31], [1, 8, 21, 31], [1, 14, 20, 28], [1, 18, 21, 31], [1, 18, 21, 533], [1, 18, 30, 533], [1, 20, 22, 567], [1, 20, 30, 533], [1, 20, 199, 381], [1, 20, 199, 567], [1, 20, 381, 383], [1, 20, 383, 533], [1, 21, 22, 24], [1, 21, 22, 135], [1, 21, 22, 383], [1, 21, 135, 532], [1, 21, 383, 533], [1, 22, 135, 532], [1, 22, 199, 381], [1, 22, 199, 567], [1, 22, 381, 383]]:

                # overlap = calc_dist(alocs[other_atom], np.array(vert['loc'])) - arads[other_atom] - abs(vert['rad'])
                # print((vert['atoms'], atom_dist, vert['loc_dist'], overlap, vert['rad']))
                continue
            return vert, vert['ndx_list_loc'], metrics
    return None
