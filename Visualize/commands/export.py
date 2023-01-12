from Visualize.commands.interpret import *
from System.sys_objs.group import Group


def group(sys, usr_npt, for_export=False):
    """
    Takes input strings interprets them and returns a group
    :param sys:
    :param for_export:
    :param usr_npt:
    :return:
    """
    # Initial selection check
    selections, selection_names = [], []
    for word in usr_npt:
        if word in ands or word.lower() in group_cmds or (for_export and word.lower() in export_cmds):
            selections.append([])
            selection_names.append("")
        else:
            selections[-1].append(word)
            selection_names[-1] += word
    # Make sure that the selections are long enough
    for selection in selections:
        selection += [None] * abs((2 - len(selection)))
    groups = []
    # Create the groups
    for i in range(len(selections)):
        # pull the selection variable
        selection = selections[i]
        # Get the first selections
        my_obj = get_obj(sys=sys, obj=selection[0], return_ndx=True)
        my_ndx = get_ndx(sys=sys, obj=my_obj, ndx=selection[1])
        # Check to see if the index provided is out of range
        checking_ndx, my_atoms = True, None
        while checking_ndx:
            # Try to create the selection
            try:
                my_atoms = [sys.mols, sys.residues, sys.atoms, sys.ndxs][my_obj - 1][my_ndx]
                checking_ndx = False
            except IndexError:
                my_ndx = get_ndx(sys=sys, obj=my_obj)
                continue
        # Check real quick for single atoms
        list_atoms = []
        if type(my_atoms) is list:
            for atom in my_atoms:
                if type(atom) is list:
                    for sub_atom in atom:
                        list_atoms.append(sub_atom)
                else:
                    list_atoms.append(atom)
        else:
            list_atoms = [my_atoms]
        # Create the group
        my_group = Group(sys, atoms=list_atoms, name="_".join(selection_names[i]))
        my_group.get_info()
        # Naming loop
        naming = True
        while naming:
            # Ask the user if they want to rename the group
            rename = input("Group name: \"{}\"\nconfirm >>>   ".format(my_group.name))
            # If they want to rename the group, let them
            if rename in ns:
                new_name = input("name >>>   ")
                new_name.replace(" ", "_")
                my_group.name = new_name
            # If they don't we are done
            elif rename in ys:
                naming = False
            elif rename in quits:
                return
            else:
                change_to_input = input("Group name: \'{}\' \nconfirm >>>   ".format(rename))
                if change_to_input in ys:
                    group.name = rename
        # Create the group
        groups.append(my_group)
    # If the groups have been made for export, do not add them to the system, just return them
    if for_export:
        return groups
    else:
        sys.groups += groups
        return groups



def export(sys, usr_npt):
    """
    Takes in input strings and exports them based on their option choices
    :param sys:
    :param usr_npt:
    :return:
    """
    if len(usr_npt) > 1 and usr_npt[1].lower() == 'all':
        sys.exports(network=True, pdb=True, surfaces=True, full_network_object=True, no_sol_network_object=True,
                    alter_atoms_script=True)
        return
    # Get the groups based off of what was specified
    groups = group(sys, usr_npt, for_export=True)
    # Go through the groups in the list
    for my_group in groups:
        my_group.get_info()
        my_group.exports()
    # Check to see if there is a second group
    if len(groups) > 1:
        # Ask the user if they want to export the interface
        exp_iface = input("Exporting interface between {} and {}\nconfirm >>>   ".format(groups[0].name, groups[1].name))
        if exp_iface.lower() in ys:
            groups[0].bff, groups[1].bff = groups[1], groups[0]
            groups[0].get_info()
            groups[1].get_info()
            groups[0].export_iface(groups[1], True, True)
            print("Groups {} and {}, and {}-{} interface exported".format(groups[0].name, groups[1].name, groups[0].name, groups[1].name))
        else:
            print("Groups {} and {} exported".format(groups[0].name, groups[1].name))
    else:
        print("Group {} exported".format(groups[0].name))
