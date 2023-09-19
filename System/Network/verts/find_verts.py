from System.Network.verts.calc_vert import calc_vert
from System.Network.verts.find_v0 import find_v0
from System.Network.verts.find_site import find_site
from System.Network.verts.find_site import find_site_fast
from System.sys_funcs.calcs.calcs import get_time, ndx_search
import time
from numpy import sqrt


# Find network function. Keeps searching the network until all verts are found
def find_verts(alocs, arads, max_vert, net_type, check_atoms, a0=None, my_group=None, averts=None, vert_ndxs=None,
               vlocs=None, vrads=None, vloc2s=None, vrad2s=None, start_time=0, print_metrics=False, vert_box=None):
    """
    Used a vertex and a combination of it's edge atoms to find the connecting vertex
    """
    # Metrics measuring < -- Deleting later
    metrics = None
    start = time.perf_counter()
    if print_metrics:
        metrics = {'ndx_search': 0, 'box_search': 0, 'gather_atoms': 0, 'verify_site': 0, 'calc_vert': 0, 'other': 0}
    # Get the group atoms from which to check vertices against
    if my_group is None or (my_group is not None and len(my_group.atoms) == len(alocs)):
        group_atoms = [i for i in range(len(alocs))]
        # Calculate the rough number of vertices
        tot_verts = 7 * len(alocs)
    # If a group was provided make sure to get its indices
    elif my_group is not None:
        group_atoms = my_group.atom_ndxs
        # Calculate the number of vertices
        tot_verts = 7 * len(group_atoms) + int(60 * sqrt(len(group_atoms)))
    else:
        return
    if averts is None:
        averts = [[] for _ in range(len(alocs))]
    # Find the first verified vertex
    if len(group_atoms) == 4:
        v0_loc, v0_rad, v0_loc2, v0_rad2 = calc_vert(locs=[alocs[_] for _ in group_atoms],
                                                     rads=[arads[_] for _ in group_atoms])
        v0 = {'atoms': group_atoms, 'loc': v0_loc, 'rad': v0_rad, 'loc2': v0_loc2, 'rad2': v0_rad2}
    else:
        v0 = find_v0(alocs=alocs, arads=arads, averts=averts, max_vert=max_vert, net_type=net_type, a0=a0,
                     group_atoms=group_atoms, metrics=metrics)
    # If no v0 is possible (e.g., a lone atom) return
    if v0 is None:
        return
    # Check if this is the first go around
    if vert_ndxs is None:
        for atom in v0['atoms']:
            # noinspection PyTypeChecker
            averts[atom].append(0)
        vert_ndxs = [v0['atoms']]
        vlocs = [v0['loc']]
        vrads = [v0['rad']]
        if 'loc2' in v0:
            vloc2s = [v0['loc2']]
            vrad2s = [v0['rad2']]
        else:
            vloc2s = [[None, None, None]]
            vrad2s = [None]
    else:
        for atom in v0['atoms']:
            avert_ndxs = [vert_ndxs[_] for _ in averts[atom]]
            # noinspection PyTypeChecker
            averts[atom].insert(ndx_search(avert_ndxs, v0['atoms']), len(vert_ndxs))
        vert_ndxs.append(v0['atoms'])
        vlocs.append(v0['loc'])
        vrads.append(v0['rad'])
        if 'loc2' in v0:
            vloc2s.append(v0['loc2'])
            vrad2s.append(v0['rad2'])
        else:
            vloc2s.append([None, None, None])
            vrad2s.append(None)
    # Set up the vertex stack
    vert_stack = [v0]
    # While the verts stack is not empty
    while vert_stack:
        # Get the vertex from the top of the stack
        vert = vert_stack.pop()
        # Set up the edge stack
        e_stack = [[[vert['atoms'][i], vert['atoms'][(i + 1) % 4], vert['atoms'][(i + 2) % 4]], vert] for i in range(4)]
        # While the edge stack is not empty
        while e_stack:
            # Get the percentage and print it
            percentage = min((len(vlocs) / tot_verts) * 100, 100)
            my_time = time.perf_counter() - start_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{:02d}:{:2.2f} - Process: finding vertices: {} verts - {:.2f} %"
                  .format(int(h), int(m), round(s, 2), len(vert_ndxs), percentage), end="")
            # Get the edge from the top of the stack
            edge_atoms, vert = e_stack.pop()
            # Find the next site in the network
            if net_type == 'pow':
                vert_ndx_pr = find_site(edge_atoms=edge_atoms, alocs=alocs, arads=arads, averts=averts,
                                             vert_ndxs=vert_ndxs, max_vert=max_vert, net_type=net_type,
                                             vn_1=vert['atoms'],
                                             vn_1_loc=vert['loc'], group_atoms=group_atoms, metrics=metrics)
            else:
                vert_ndx_pr = find_site_fast(edge_atoms=edge_atoms, alocs=alocs, arads=arads, averts=averts,
                                             vert_ndxs=vert_ndxs, max_vert=max_vert, net_type=net_type,
                                             vn_1=vert['atoms'],
                                             vn_1_loc=vert['loc'], group_atoms=group_atoms, metrics=metrics,
                                             vn_1_rad=vert['rad'])
            # If the vertex is none continue
            if vert_ndx_pr is None:
                continue

            # Set the vertex and its index
            my_vert, metrics = vert_ndx_pr
            # Check if there is a retaining box for the vertices and if the vertex is outside the box
            if vert_box is not None and any([vert_box[0][i] > my_vert[i] < vert_box[1][i] for i in range(3)]):
                continue
            # Add the vertex to the stack and the network
            vert_stack.append(my_vert)
            # Insert the vertices in order of increasing atom indices
            vert_ndxs.append(my_vert['atoms'])
            vlocs.append(my_vert['loc'])
            vrads.append(my_vert['rad'])
            if 'loc2' in my_vert:
                vloc2s.append(my_vert['loc2'])
                vrad2s.append(my_vert['rad2'])
            else:
                vloc2s.append([None, None, None])
                vrad2s.append(None)
            # Remove the atoms from the
            for atom in my_vert['atoms']:
                # noinspection PyTypeChecker
                avert_ndxs = [vert_ndxs[_] for _ in averts[atom]]
                # noinspection PyTypeChecker
                averts[atom].insert(ndx_search(avert_ndxs, my_vert['atoms']), len(vert_ndxs) - 1)
                if atom in check_atoms:
                    check_atoms.remove(atom)
    # Printing out metrics < --- Delete later
    if print_metrics:
        metrics['total'] = time.perf_counter() - start
        metrics['other'] = metrics['total'] - (metrics['ndx_search'] + metrics['box_search'] + metrics['gather_atoms'] +
                                               metrics['verify_site'] + metrics['calc_vert'])
        print('\n\nVertex Finding Time Metrics: \n'
              '  Index Search      = {:.3f} s\n'
              '  Box Search        = {:.3f} s\n'
              '  Get Atoms         = {:.3f} s\n'
              '  Verify Site       = {:.3f} s\n'
              '  Calculate Vertex  = {:.3f} s\n'
              '  Other             = {:.3f} s\n'
              '  Total             = {:.3f} s\n\n'
              .format(metrics['ndx_search'], metrics['box_search'], metrics['gather_atoms'], metrics['verify_site'],
                      metrics['calc_vert'], metrics['other'], metrics['total']))
    # Return the values of the vertices
    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, check_atoms
