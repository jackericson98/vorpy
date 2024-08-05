import time
import pandas as pd
from System.sys_funcs.calcs.sorting import global_vars, box_search, get_atoms
from System.sys_funcs.calcs.calcs import calc_dist
from System.Network.verts.find_verts import find_verts
from System.sys_funcs.output.net import write_verts


def find_net_verts(net):
    # Get the global variables
    global_vars(net.sub_boxes, net.box, net.settings['num_splits'], net.group.sys.max_atom_rad, net.sub_box_size)

    # Not sure what this does
    # vert_list_real = net.get_real_verts()
    # Create the group indices
    atom_nums = net.group.group_ndxs.copy()
    # Get the indices of the atoms in the network to keep track of the atoms that haven't been visited
    my_guuy = find_verts(alocs=net.spheres['loc'].to_numpy(), arads=net.spheres['rad'].to_numpy(),
                         max_vert=net.settings['max_vert'], net_type=net.settings['net_type'], check_atoms=atom_nums,
                         my_group=net.group.group_ndxs, start_time=net.start_time, print_metrics=net.settings['print_metrics'],
                         vert_box=net.group.sys.foam_box)
    if my_guuy is not None:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = my_guuy
    # Check to see if any of the atoms are encapsulated
    if len(atom_nums) > 0:
        skip_nums = []
        for atom in atom_nums:
            atom_rad, atom_loc = net.spheres['rad'][atom], net.spheres['loc'][atom]
            atom_box = box_search(atom_loc)
            near_atoms = get_atoms(atom_box, dist=net.group.sys.max_atom_rad - atom_rad)
            for atom2 in near_atoms:
                if calc_dist(atom_loc, net.spheres['loc'][atom2]) < abs(net.spheres['rad'][atom2] - atom_rad):
                    print("\nUh oh! Ball # {} is fully encapsulated by ball # {}! Skipping {}"
                          .format(atom, atom2, atom))
                    skip_nums.append(atom)
                    break
        for _ in skip_nums:
            atom_nums.pop(atom_nums.index(_))

    # Check for disconnects in the network
    threshold = 2
    if len(net.group.group_ndxs) <= 2:
        threshold = 0
    while len(atom_nums) > threshold:
        print("Atoms Disconnected: {}".format(atom_nums))
        a0 = atom_nums.pop()
        my_guuy = find_verts(a0=a0, alocs=net.spheres['loc'].to_numpy(), arads=net.spheres['rad'].to_numpy(),
                             max_vert=net.settings['max_vert'], net_type=net.settins['net_type'], check_atoms=atom_nums,
                             my_group=net.group.group_ndxs, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                             vloc2s=vloc2s, vrad2s=vrad2s, start_time=net.start_time, print_metrics=print_metrics,
                             vert_box=net.group.sys.foam_box, averts=averts)
        if my_guuy is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = my_guuy
        if net.group.sys.type == 'foam' and len(atom_nums) <= 0.25 * len(net.atoms['loc']):
            break
    # # Create the doublets list
    # if vert_list_real is not None and net.type == 'aw':
    #     missing_verts = [_ for _ in vert_list_real if _ not in vert_ndxs]
    #     print(missing_verts)
    #     extra_verts = [_ for _ in vert_ndxs if _ not in vert_list_real]
    #     print(extra_verts)
    doublets = [0 for _ in range(len(vert_ndxs))]
    # Incorporate the doublets into the vlocs, vatoms, vrads lists and lose the vloc2s and vrad2s
    i = 0
    while i < len(vlocs):
        # Check for doubletness
        if vrad2s[i] is not None:
            # Insert the relevant information into their respective lists
            vert_ndxs.insert(i + 1, vert_ndxs[i])
            vlocs.insert(i + 1, vloc2s[i])
            vrads.insert(i + 1, vrad2s[i])
            doublets.insert(i + 1, 1)
            # Preserve the relational aspects of vrad2s and vloc2s
            vrad2s.insert(i + 1, None)
            vloc2s.insert(i + 1, [None, None, None])
        i += 1

    # Make the dataframe
    net.verts = pd.DataFrame({"vatoms": vert_ndxs, 'vloc': vlocs, 'vrad': vrads, 'vdub': doublets})
    # Clear the print statement
    print("\r                                                                  ", end="")
    net.metrics['vert'] = time.perf_counter() - net.start_time
    write_verts(net)
