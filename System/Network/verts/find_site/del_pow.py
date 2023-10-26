from System.sys_funcs.calcs.calcs import box_search, get_atoms, calc_circ, calc_dist, ndx_search
from System.Network.verts.verify_site import verify_site
from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
import numpy as np
import time


def find_site_pd_container(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, net_type, vn_1, vn_1_loc, group_atoms=None,
                             metrics=None, vn_1_rad=None):
    """
    Cycles through larger and larger areas searching for
    """
    # Define the global variables for the other functions to tap into
    global invalid_atoms, edge_anums
    edge_anums = edge_atoms
    invalid_atoms = []
    vert = None
    # Se the initial vert size
    mv_inc = 0.15
    # Look for the vert and keep increasing box size until the vert is found
    while vert is None and mv_inc < max_vert:
        vert = find_site_pd(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type, vn_1, vn_1_loc, group_atoms=group_atoms,
                              metrics=metrics)
        mv_inc *= 5
    # Las step find the vertex using the maximum size
    if vert is None:
        vert = find_site_pd(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, max_vert, net_type, vn_1, vn_1_loc,
                              group_atoms=group_atoms, metrics=metrics)
    # Return the vertex if found
    return vert


def find_site_pd(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type, vn_1=None, vn_1_loc=None, group_atoms=None, metrics=None):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """
    global invalid_atoms
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
    try:
        test_atoms = [_ for _ in get_atoms(cells=my_boxes, dist=mv_inc) if _ not in invalid_atoms]
    except NameError:
        invalid_atoms = []
        test_atoms = get_atoms(cells=my_boxes, dist=max_vert)
    surr_atoms = get_atoms(cells=my_boxes, dist=max_vert)

    if metrics is not None:
        metrics['gather_atoms'] += time.perf_counter() - start
    # First look for vertices that have been found before
    new_test_atoms = []
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
            invalid_atoms.append(atom)
            continue
        verts.append({'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad})
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_atoms = [_ for _ in surr_atoms if _ not in vert_atoms]
        test_locs = np.array([alocs[_] for _ in filtered_test_atoms])
        test_rads = np.array([arads[_] for _ in filtered_test_atoms])
        if vert_rad < max_vert and verify_site(loc=np.array(vert_loc), rad=vert_rad, test_locs=test_locs,
                                                    test_rads=test_rads, net_type=net_type):
            if metrics is not None:
                metrics['verify_site'] += time.perf_counter() - start
            return {'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad}, metrics
        else:
            invalid_atoms.append(atom)
