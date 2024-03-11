import time
import pandas as pd
from System.sys_funcs.calcs.sorting import global_vars, box_search, get_atoms
from System.sys_funcs.calcs.calcs import calc_dist
from System.Network.verts.find_verts import find_verts
from System.sys_funcs.output.net import write_verts


def find_net_verts(net, my_group=None, print_metrics=False):
    global_vars(net.sub_boxes, net.box, net.num_splits, net.sys.max_atom_rad, net.sub_box_size)
    # Check to see if a group has been provided
    if my_group is not None:
        atom_nums = my_group.atom_ndxs[:]
    else:
        atom_nums = [i for i in range(len(net.atoms))]
    vert_list_real = net.get_real_verts()
    # Get the indices of the atoms in the network to keep track of the atoms that haven't been visited
    net.atom_ndxs = [_ for _ in atom_nums]
    my_group_atom_ndxs = None
    if my_group is not None:
        my_group_atom_ndxs = my_group.atom_ndxs
    my_guuy = find_verts(alocs=net.atoms['loc'].to_numpy(), arads=net.atoms['rad'].to_numpy(),
                         max_vert=net.max_vert, net_type=net.type, check_atoms=atom_nums,
                         my_group=my_group_atom_ndxs, start_time=net.start_time, print_metrics=print_metrics,
                         vert_box=net.sys.foam_box)
    if my_guuy is not None:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = my_guuy
    # Check to see if any of the atoms are encapsulated
    if len(atom_nums) > 0:
        skip_nums = []
        for atom in atom_nums:
            atom_rad, atom_loc = net.atoms['rad'][atom], net.atoms['loc'][atom]
            atom_box = box_search(atom_loc)
            near_atoms = get_atoms(atom_box, dist=net.sys.max_atom_rad - atom_rad)
            for atom2 in near_atoms:
                if calc_dist(atom_loc, net.atoms['loc'][atom2]) < abs(net.atoms['rad'][atom2] - atom_rad):
                    print("\nUh oh! Ball # {} is fully encapsulated by ball # {}! Skipping {}"
                          .format(atom, atom2, atom))
                    skip_nums.append(atom)
                    break
        for _ in skip_nums:
            atom_nums.pop(atom_nums.index(_))

    # Check for disconnects in the network
    while len(atom_nums) > 0:
        a0 = atom_nums.pop()
        my_guuy = find_verts(a0=a0, alocs=net.atoms['loc'].to_numpy(), arads=net.atoms['rad'].to_numpy(),
                             max_vert=net.max_vert, net_type=net.type, check_atoms=atom_nums,
                             my_group=my_group.atom_ndxs, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                             vloc2s=vloc2s, vrad2s=vrad2s, start_time=net.start_time, print_metrics=print_metrics,
                             vert_box=net.sys.foam_box, averts=averts)
        if my_guuy is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = my_guuy
        if net.sys.type == 'foam' and len(atom_nums) <= 0.25 * len(net.atoms['loc']):
            break
    # Create the doublets list
    if vert_list_real is not None and net.type == 'vor':
        missing_verts = [_ for _ in vert_list_real if _ not in vert_ndxs]
        print(missing_verts)
        extra_verts = [_ for _ in vert_ndxs if _ not in vert_list_real]
        print(extra_verts)
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
