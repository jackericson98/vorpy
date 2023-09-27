from System.Network.verts.calc_vert import calc_flat_vert
from System.sys_funcs.calcs.calcs import box_search, get_atoms, ndx_search, calc_dist
import time
import numpy as np


def find_site_pow_del(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, net_type, vn_1, vn_1_loc, group_atoms=None,
                      metrics=None):
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = edge_atoms[:]

    # Check if the edge contains a group atom or not. If not each atom paired with the edge atoms needs to be checked
    check_atoms = True
    for ndx in edge_ndxs:
        if group_atoms is not None and ndx in group_atoms:
            check_atoms = False
            break
    # Time printing metrics <-- Delete later
    start = time.perf_counter()
    # Grab the atoms we want to test against
    my_boxes = [box_search(loc=alocs[edge_atoms[_]]) for _ in range(3)]
    # Box search metrics <-- Delete later
    if metrics is not None:
        metrics['box_search'] += time.perf_counter() - start
        start = time.perf_counter()
    # If the previous vertex has been provided, add the other atom to the not allowed atoms
    vert_atom_ndxs = vn_1
    if vn_1 is None:
        vert_atom_ndxs = edge_ndxs
    # Gather the atoms surrounding the edge atoms
    test_atoms = get_atoms(cells=my_boxes, dist=max_vert)
    if metrics is not None:
        metrics['gather_atoms'] += time.perf_counter() - start
    # First look for vertices that have been found before
    test_verts = []
    start = time.perf_counter()
    for atom in test_atoms:
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
        test_verts.append(atom_ndxs)

    # Index search metrics <-- Delete later
    if metrics is not None:
        metrics['ndx_search'] += time.perf_counter() - start
        start = time.perf_counter()
    # Create the calculated vertices list
    calculated_verts = []
    # Go through each atom in the given test atoms. Extremely optimized
    for i, vert in enumerate(test_verts):
        # Make sure the vertex values are defined
        vert_loc, vert_rad, vert_loc2, vert_rad2 = None, None, None, None

        # If the network type is power
        if net_type == 'pow':
            # Calculate the power vertex values
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert], rads=[arads[_] for _ in vert],
                                                power=True)
        # If the network type is Delaunay
        elif net_type == 'del':
            # Calculate the Delaunay vertex values
            vert_loc, vert_rad = calc_flat_vert(locs=[alocs[_] for _ in vert], rads=[arads[_] for _ in vert],
                                                power=False)

        # Catch the none location and the too large vertex cases
        if vert_loc is None or vert_rad > max_vert:
            continue

        # Delete the second location for the vertex if it is too large
        if vert_rad2 is not None and vert_rad2 > max_vert:
            vert_loc2, vert_rad2 = None, None

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
        return calculated_verts[0], metrics

    # Find the closest vertex to the previous vertex
    dists = [calc_dist(np.array(_['loc']), np.array(vn_1_loc)) for _ in calculated_verts]
    dists, calculated_verts = zip(*sorted(zip(dists, calculated_verts)))

    # Return the closest vertex
    return calculated_verts[0], metrics

