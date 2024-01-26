import os
import csv
import numpy as np
import pandas as pd
import time
from System.sys_funcs.calcs.sorting import divide_box, global_vars, ndx_search
from System.sys_funcs.calcs.calcs import calc_com, calc_dist
from System.Network.verts.find_verts import find_verts
from System.sys_funcs.output.net import add_metrics


def split_net(sys, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None, my_group=None,
              print_actions=None, num_atoms_sub_net=100, add_net_metrics=True, min_atom_split=30):
    # Sort the atoms in the main network
    sys.net.sort_atoms()
    # Calculate the group box
    group_box = sys.net.calc_box([sys.atoms['loc'][_] for _ in my_group.atoms],
                                 [sys.atoms['rad'][_] for _ in my_group.atoms], return_val=True, box_size=1.1)
    # Get the sub boxes
    sub_boxes = divide_box(group_box, round(len(my_group.atoms) / num_atoms_sub_net), c=0.1)
    print('num splits', len(sub_boxes))
    # Check for a max_vert that isn't defined
    if max_vert is None:
        max_vert = sys.net.max_vert

    # Sort the atoms into their sub_boxes
    atoms_lists = [[] for _ in range(len(sub_boxes))]
    atom_locs = sys.atoms['loc']
    # Loop through the atom locations and sort the atoms
    for atom in my_group.atoms:
        loc = atom_locs[atom]
        # Loop through the sub boxes to find the placement of the atom
        for j, sub_box in enumerate(sub_boxes):
            if [sub_box[0][k] <= loc[k] <= sub_box[1][k] for k in range(3)] == [True, True, True]:
                atoms_lists[j].append(atom)
    # If a list of atoms is too small add it to another
    skip_boxes = []
    for i, atoms_list in enumerate(atoms_lists):
        # If no atoms exist nothing to deal with
        if len(atoms_list) == 0:
            skip_boxes.append(i)
            continue
        if len(atoms_list) < min_atom_split:
            # Get the com of the atoms to find the closes sub_box to add to
            atoms_com = calc_com([atom_locs[_] for _ in atoms_list])
            min_dist = np.inf
            closest_sub_box = None
            for j, sub_box in enumerate(sub_boxes):
                # Make sure we aren't adding to a sub_box scheduled for deletion
                if j in skip_boxes or j == i:
                    continue
                # Calculate the distance of the com of the sub_box from the atoms_com
                my_dist = calc_dist(calc_com(sub_box), atoms_com)
                # Replace the variables if they are closer
                if my_dist < min_dist:
                    closest_sub_box, min_dist = j, my_dist
            # Add the atoms to the new sub_box
            atoms_lists[closest_sub_box] += atoms_list
            skip_boxes.append(i)
    for i, atom_list in enumerate(atoms_lists):
        # Skip the boxes to be skipped
        if i in skip_boxes:
            continue
    # Instantiate the global variables
    global_vars(sys.net.sub_boxes, sys.net.box, sys.net.num_splits, sys.max_atom_rad, sys.net.sub_box_size)
    vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = None, None, None, None, None, None, None
    # Create the subnetworks
    for i, atom_list in enumerate(atoms_lists):
        # Skip the boxes to be skipped
        if i in skip_boxes:
            continue
        # Get the atoms we are tying to find vertices for
        check_atoms = [_ for _ in atom_list if _ in my_group.atoms]
        atom_nums = check_atoms[:]
        # Find the initial vertices for the vertex group
        init_verts = find_verts(alocs=sys.atoms['loc'].to_numpy(), arads=sys.atoms['rad'].to_numpy(),
                                max_vert=max_vert, net_type=net_type, check_atoms=check_atoms,
                                my_group=atom_nums, start_time=sys.net.start_time,
                                vert_box=sys.foam_box, group_box=sub_boxes[i], vert_ndxs=vert_ndxs, vlocs=vlocs,
                                vrads=vrads, vloc2s=vloc2s, vrad2s=vrad2s, averts=averts,
                                tot_atom_num=len(my_group.atoms))
        # Check to see if find_verts fails
        if init_verts is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = init_verts

        # Check for disconnects in the network
        while len(atom_nums) > 0:
            # Grab the initial atom for the next search
            a0 = atom_nums.pop()

            # Find verts again
            more_verts = find_verts(a0=a0, alocs=sys.atoms['loc'].to_numpy(), arads=sys.atoms['rad'].to_numpy(),
                                    max_vert=max_vert, net_type=net_type, check_atoms=atom_nums,
                                    my_group=check_atoms, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                                    vloc2s=vloc2s, vrad2s=vrad2s, start_time=sys.net.start_time,
                                    vert_box=sys.foam_box, averts=averts, group_box=sub_boxes[i],
                                    tot_atom_num=len(my_group.atoms))
            # Check to see if find_verts fails
            if more_verts is not None:
                vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = more_verts
            # Every sphere needs a vert
            if sys.foam_box is not None and len(atom_nums) <= 0.25 * len(sys.atoms['loc']):
                break
    # Create the doublets list
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
    sys.net.verts = pd.DataFrame({"vatoms": vert_ndxs, 'vloc': vlocs, 'vrad': vrads, 'vdub': doublets})
    # Clear the print statement
    if sys.print_actions:
        print("\r                                                                  ", end="")
    sys.net.metrics['vert'] = time.perf_counter() - sys.net.start_time
    sys.net.build(surf_res=surf_res, max_vert=max_vert, box_size=box_size, build_surfs=build_surfs,
                  calc_verts=False, net_type=net_type, my_group=my_group, print_actions=print_actions)
    sys.net.metrics['splits'] = len(sub_boxes)
    if add_net_metrics:
        add_metrics(sys.net)


def split_net_slow(sys, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None, my_group=None,
                   print_actions=None, num_atoms_sub_net=50, add_net_metrics=True, min_atom_split=30):
    # Check to see if the pdb directory is suitable
    if sys.dir is None:
        sys.set_output_directory()
        os.mkdir(sys.dir + '/verts')
    # Sort the atoms in the main network
    sys.net.sort_atoms()
    # Calculate the group box
    group_box = sys.net.calc_box([sys.atoms['loc'][_] for _ in my_group.atoms],
                                 [sys.atoms['rad'][_] for _ in my_group.atoms], return_val=True, box_size=1.1)
    # Get the sub boxes
    sub_boxes = divide_box(group_box, round(len(my_group.atoms) / num_atoms_sub_net), c=1)
    print('num splits', len(sub_boxes))
    # Check for a max_vert that isn't defined
    if max_vert is None:
        max_vert = sys.net.max_vert

    # Sort the atoms into their sub_boxes
    atoms_lists = [[] for _ in range(len(sub_boxes))]
    atom_locs = sys.atoms['loc']
    # Loop through the atom locations and sort the atoms
    for atom in my_group.atoms:
        loc = atom_locs[atom]
        # Loop through the sub boxes to find the placement of the atom
        for j, sub_box in enumerate(sub_boxes):
            if [sub_box[0][k] <= loc[k] <= sub_box[1][k] for k in range(3)] == [True, True, True]:
                atoms_lists[j].append(atom)
    # If a list of atoms is too small add it to another
    skip_boxes = []
    for i, atoms_list in enumerate(atoms_lists):
        # If no atoms exist nothing to deal with
        if len(atoms_list) == 0:
            skip_boxes.append(i)
            continue
        if len(atoms_list) < min_atom_split:
            # Get the com of the atoms to find the closes sub_box to add to
            atoms_com = calc_com([atom_locs[_] for _ in atoms_list])
            min_dist = np.inf
            closest_sub_box = None
            for j, sub_box in enumerate(sub_boxes):
                # Make sure we aren't adding to a sub_box scheduled for deletion
                if j in skip_boxes or j == i:
                    continue
                # Calculate the distance of the com of the sub_box from the atoms_com
                my_dist = calc_dist(calc_com(sub_box), atoms_com)
                # Replace the variables if they are closer
                if my_dist < min_dist:
                    closest_sub_box, min_dist = j, my_dist
            # Add the atoms to the new sub_box
            atoms_lists[closest_sub_box] += atoms_list
            skip_boxes.append(i)
    # Instantiate the global variables
    global_vars(sys.net.sub_boxes, sys.net.box, sys.net.num_splits, sys.max_atom_rad, sys.net.sub_box_size)
    # Create the vertices file
    vert_file_name = sys.dir + '/verts/verts.txt'
    with open(vert_file_name, 'w') as verts_file, open(sys.dir + '/verts/averts.txt', 'w') as averts_file:
        # Header
        verts_file.write(sys.name + " vertices")
        averts_file.write((sys.name + ' atom vertices by box'))
        # Vertices
        for i, atom_list in enumerate(atoms_lists):
            # Reset the variables
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = None, None, None, None, None, None, None
            # Write the box info
            verts_file.write('\nBox {} - {}'.format(i, sub_boxes[i]))
            averts_file.write('\nBox {}'.format(i))
            # Skip the boxes to be skipped
            if i in skip_boxes:
                continue
            # Get the atoms we are tying to find vertices for
            check_atoms = [_ for _ in atom_list if _ in my_group.atoms]
            atom_nums = check_atoms[:]
            # Find the initial vertices for the vertex group
            init_verts = find_verts(alocs=sys.atoms['loc'].to_numpy(), arads=sys.atoms['rad'].to_numpy(),
                                    max_vert=max_vert, net_type=net_type, check_atoms=check_atoms,
                                    my_group=atom_nums, start_time=sys.net.start_time,
                                    vert_box=sys.foam_box, group_box=sub_boxes[i], vert_ndxs=vert_ndxs, vlocs=vlocs,
                                    vrads=vrads, vloc2s=vloc2s, vrad2s=vrad2s, averts=averts,
                                    tot_atom_num=len(my_group.atoms))
            # Check to see if find_verts fails
            if init_verts is not None:
                vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = init_verts
            # Check for disconnects in the network
            while len(atom_nums) > 0:
                # Grab the initial atom for the next search
                a0 = atom_nums.pop()

                # Find verts again
                more_verts = find_verts(a0=a0, alocs=sys.atoms['loc'].to_numpy(), arads=sys.atoms['rad'].to_numpy(),
                                        max_vert=max_vert, net_type=net_type, check_atoms=atom_nums,
                                        my_group=check_atoms, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                                        vloc2s=vloc2s, vrad2s=vrad2s, start_time=sys.net.start_time,
                                        vert_box=sys.foam_box, averts=averts, group_box=sub_boxes[i],
                                        tot_atom_num=len(my_group.atoms), printing=True if i == 5 else False)
                # Check to see if find_verts fails
                if more_verts is not None:
                    vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = more_verts
                # Every sphere needs a vert
                if sys.foam_box is not None and len(atom_nums) <= 0.25 * len(sys.atoms['loc']):
                    break

            # Write the vertex information into the vert file
            for j, vert in enumerate(vert_ndxs):
                # Write the line
                loc2 = 'None'
                if vloc2s[j] is not None:
                    loc2 = str([_ for _ in vloc2s[j]])
                verts_file.write('\n' + str(j) + ';' + str(vert) + ';' + str([_ for _ in vlocs[j]]) + ';' + str(vrads[j]) + ';' +
                                 loc2 + ';' + str(vrad2s[j]))

    # Combine the subnetworks
    vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts = combine_nets(vert_file_name, len(sys.atoms['loc']))
    # Create the doublets list
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
    sys.net.verts = pd.DataFrame({"vatoms": vert_ndxs, 'vloc': vlocs, 'vrad': vrads, 'vdub': doublets})
    # Clear the print statement
    if sys.print_actions:
        print("\r                                                                  ", end="")
    sys.net.metrics['vert'] = time.perf_counter() - sys.net.start_time
    sys.net.build(surf_res=surf_res, max_vert=max_vert, box_size=box_size, build_surfs=build_surfs,
                  calc_verts=False, net_type=net_type, my_group=my_group, print_actions=print_actions)
    sys.net.metrics['splits'] = len(sub_boxes)
    if add_net_metrics:
        add_metrics(sys.net)


def combine_nets(verts, num_atoms):
    # Create the dictionary based on the sub_boxes
    file_verts = []
    # Go through the regular vertices
    with open(verts, 'r') as verts_file:
        # Go through the vertices
        for i, line in enumerate(verts_file):
            if i == 0 or line[:3] == 'Box':
                continue
            # Split the line
            line = line.split(';')
            atoms = line[1][1:-1].split(',')
            loc = line[2][1:-1].split(',')
            if line[4][0] == 'N':
                loc2 = None
            else:
                loc2 = line[4][1:-1].split(',')
                loc2 = [float(_) for _ in loc2]
            if line[5][0] == 'N':
                rad2 = None
            else:
                rad2 = float(line[5])
            my_vert = {'atoms': [int(_) for _ in atoms if _ != ''], 'loc': [float(_) for _ in loc if _ != ''], 'rad': float(line[3]),
                               'loc2': loc2, 'rad2': rad2}
            if my_vert['atoms'] == []:
                print(my_vert)
            file_verts.append(my_vert)
    # Go through the sub_boxes one by one
    averts = [[] for _ in range(num_atoms)]
    verts = []

    for vert in file_verts:
        # Create a tracking variable for whether the vertex has been found or not
        vert_found = False
        if vert['atoms'] == []:
            print(vert)
        # Go through the atom vertices for the first atom in the vertex
        for my_vert in averts[vert['atoms'][0]]: ## This is where ndx search comes in
            if my_vert < len(verts) and verts[my_vert]['atoms'] == vert['atoms']:
                vert_found = True
                break
        if vert_found:
            continue
        else:
            # Add the vertex to the list of atoms
            verts.append(vert)
            vert_ndx = len(verts) - 1
            # Add the vertex to each atom in the vertex:
            for ndx in vert['atoms']:
                # Get the list of vertex indices

                atom_vert_ndxs = [verts[_]['atoms'] for _ in averts[ndx]]
                # find the place where to insert the vertex
                my_ndx = ndx_search(atom_vert_ndxs, vert['atoms'])
                averts[ndx].insert(my_ndx, vert_ndx)
    vert_ndxs, vlocs, vrads, vloc2s, vrad2s = [], [], [], [], []
    for i in range(len(verts)):
        vert_ndxs.append(verts[i]['atoms'])
        vlocs.append(verts[i]['loc'])
        vrads.append(verts[i]['rad'])
        vloc2s.append(verts[i]['loc2'])
        vrad2s.append(verts[i]['rad2'])
    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts
