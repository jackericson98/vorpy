from System.sys_funcs.calcs.sorting import box_search, get_atoms, ndx_search
from System.Network.verts.verify_site import verify_site
from System.Network.verts.calc_vert import calc_flat_vert, calc_vert
import numpy as np
import time


def find_site(edge_atoms, alocs, arads, averts, vert_ndxs, max_vert, mv_inc, net_type, invalid_atoms=None,
              check_atoms=True, vn_1=None, vn_1_loc=None, group_atoms=None, metrics=None):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """
    if invalid_atoms is None:
        invalid_atoms = []
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = edge_atoms[:]

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

    test_atoms = [_ for _ in get_atoms(cells=my_boxes, dist=mv_inc) if _ not in invalid_atoms]

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
        if vert_ndxs is not None and len(vert_ndxs) > 0:
            check_verts = [vert_ndxs[_] for _ in averts[atom_ndxs[0]]]
            my_vert_ndx = ndx_search(check_verts, atom_ndxs)
            if my_vert_ndx < len(check_verts) and atom_ndxs == check_verts[my_vert_ndx]:
                return None, invalid_atoms
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
            invalid_atoms.append([_ for _ in vert_atoms if _ not in edge_ndxs])
            continue
        start = time.perf_counter()
        # Filter the vertex out if it is too large or not able to be made
        filtered_test_atoms = [_ for _ in surr_atoms if _ not in vert_atoms]
        test_locs = np.array([alocs[_] for _ in filtered_test_atoms])
        test_rads = np.array([arads[_] for _ in filtered_test_atoms])
        if abs(vert_rad) < max_vert and verify_site(loc=np.array(vert_loc), rad=vert_rad, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            if len(verts) > 0 and verts[0]['rad'] < vert_rad:
                return [verts[0], metrics], invalid_atoms
            verts.append({'atoms': vert_atoms, 'loc': vert_loc, 'rad': vert_rad, 'loc2': None, 'rad2': None})
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if vert_loc2 is not None and abs(vert_rad2) < max_vert and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
                verts[-1]['loc2'], verts[-1]['rad2'] = vert_loc2, vert_rad2
        # Check to see if the doublet's site is verified
        elif vert_loc2 is not None and verify_site(loc=np.array(vert_loc2), rad=vert_rad2, test_locs=test_locs, test_rads=test_rads, net_type=net_type):
            verts.append({'atoms': vert_atoms, 'loc': vert_loc2, 'rad': vert_rad2, 'loc2': None, 'rad2': None})
        if metrics is not None:
            metrics['verify_site'] += time.perf_counter() - start
        else:
            invalid_atoms.append([_ for _ in vert_atoms if _ not in edge_ndxs])
    # If no verts have been found return
    if len(verts) == 0:
        return None, invalid_atoms
    # If we find only 1 vertex, return it
    elif len(verts) == 1 or verts[0]['rad'] < verts[1]['rad']:
        return [verts[0], metrics], invalid_atoms
    return [verts[1], metrics], invalid_atoms
