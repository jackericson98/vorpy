from System.Group.group import Group
from Visualize.cmnd.interpret import *
from Visualize.cmnd.commands import *
from System.chemistry_interpreter import *


def group(sys, usr_npt, settings=None):
    # Create the object and index variables
    my_obj, my_ndx = None, None
    # Go through the different info things
    if usr_npt[0] == 'ns':
        my_obj = 'm'
        my_ndx = [0, len(sys.chains) - 2]
    # User only input "export"
    elif len(usr_npt) == 0:
        # Tell the user to pick an object and an index
        my_obj = get_obj(sys=sys)
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # User entered "export obj" and needs an index
    elif len(usr_npt) == 1:
        # Add a check for system
        # Check the object provided by the user
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj)
        # Get the group information
        obj_ndx = ['c', 'r', 'a', 'n'].index(my_obj)
        obj_list = [sys.chains, sys.residues, sys.balls, sys.ndxs][obj_ndx]
        # set up the name
        if my_obj == 'c':
            name = 'Chain_' + obj_list[obj_ndx].name
        elif my_obj == 'r':
            name = obj_list[obj_ndx].name + '_' + obj_list[obj_ndx].seq
        elif my_obj == 'a':
            name = obj_list.iloc[obj_ndx].res.name + '_' + str(obj_list.iloc[obj_ndx].res.seq) + '_' + \
                   str(obj_list.iloc[obj_ndx].name)
    # If the user input an object and an index of their own
    elif len(usr_npt) == 2:
        # Check the object
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj, ndx_npt=usr_npt[1])
        # Get the group information
        obj_ndx = ['c', 'r', 'a', 'n'].index(my_obj)
        obj_list = [sys.chains, sys.residues, sys.balls, sys.ndxs][obj_ndx]
        # set up the name
        if my_obj == 'c':
            name = 'Chain_' + obj_list[obj_ndx].name
        elif my_obj == 'r':
            name = obj_list[obj_ndx].name + '_' + obj_list[obj_ndx].seq
        elif my_obj == 'a':
            name = obj_list.iloc[obj_ndx].res.name + '_' + str(obj_list.iloc[obj_ndx].res.seq) + '_' + \
                   str(obj_list.iloc[obj_ndx].name)
    elif len(usr_npt) > 2:
        # Check the object
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj, ndx_npt=usr_npt[1:])
        # Get the group information
        name = my_obj + '_' + str(my_ndx[0]) + '_' + str(my_ndx[1])
    # Get the group information
    obj_ndx = ['c', 'r', 'a', 'n'].index(my_obj)
    obj_list = [sys.chains, sys.residues, sys.balls, sys.ndxs][obj_ndx]

    my_list = None
    # Get the slice and name of the group
    if len(my_ndx) == 1 and my_ndx[0] < len(obj_list):
        if obj_ndx == 2:
            my_list = [my_ndx[0]]
        else:
            my_list = [obj_list[my_ndx[0]]]
    elif len(my_ndx) <= 2:
        if obj_ndx == 2:
            my_list = [_ for _ in range(max(0, my_ndx[0]), min(len(obj_list), my_ndx[1] + 1))]
        else:
            my_list = obj_list[max(0, my_ndx[0]):min(len(obj_list), my_ndx[1] + 1)]
    # Create the group
    npt_list = [None] * 4
    npt_list[obj_ndx] = my_list
    my_group = Group(sys=sys, chains=npt_list[0], residues=npt_list[1], atoms=npt_list[2], name=name, settings=settings)
    return my_group


def get_group_spheres(atoms, identifier):
    """
    Takes in the atoms dataframe and the identifiers and returns a list of atom objects
    """
    # First see if the identifier is an atom index
    atom_indices = atoms['num'].to_list()
    try:
        my_atoms = atom_indices[int(identifier[0])]
        return [my_atoms]
    except ValueError:
        pass
    # Next see if it is a range:
    if '-' in identifier[0]:
        # Split the indices
        index1, index2 = identifier[0].split('-')
        # Make sure the two indices are intable
        try:
            index1, index2 = int(index1), int(index2)
        except ValueError:
            pass
        # Get all the atoms in the list
        my_atoms = []
        for i in range(index1, index2 + 1):
            my_atoms.append(atoms.iloc[i])
        return my_atoms
    # Next check to see if the identifier is a residue
    if identifier[0].lower() in residue_names and len(identifier) == 3:
        # Return the atom in atoms that has both the residue name and the residue sequence number and the atom name
        try:
            res_name, res_seq, atom_name = residue_names[identifier[0].lower()], int(identifier[1]), identifier[2]
            return atoms.loc[(atoms['res_name'] == res_name) &
                              (atoms['res_seq'] == res_seq) &
                              (atoms['name'] == atom_name.upper()), 'num'].to_lost()
        except Exception as e:
            pass
    # Check if the identifier is in the atom names
    if identifier[0].upper() in atoms['name'].values:
        return atoms.loc[atoms['name'] == identifier[0].upper(), 'num'].to_list()
    # Check if the identifier is an element
    if identifier[0].lower() in element_names:

        # Return the atoms in the atoms dataframe with the same name as the atom we want
        my_atoms = atoms.loc[atoms['element'] == element_names[identifier[0].lower()], 'num'].to_list()

        return my_atoms


def get_group_resids(resids, identifier):
    # First see if the identifier is an atom index
    try:
        my_res = resids[int(identifier[0])]
        return [my_res]
    except ValueError:
        pass
    # Next see if it is a range:
    if '-' in identifier[0]:
        # Split the indices
        index1, index2 = identifier[0].split('-')
        # Make sure the two indices are intable
        try:
            index1, index2 = int(index1), int(index2)
        except ValueError:
            pass
        # Get all the atoms in the list
        my_resids = []
        for i in range(index1, index2 + 1):
            my_resids.append(resids[i])
        return my_resids
    # Next look for residue names and sequences
    if len(identifier) == 2:
        try:
            seq = int(identifier[1])
            for res in resids:
                if res.name == identifier[0].upper() and res.seq == seq:
                    return [res]
        except ValueError:
            pass


def get_group_chains(chains, identifier):
    # First try the index option
    try:
        chain_index = float(identifier[0])
        if chain_index < len(chains):
            return [chains[chain_index]]
    except ValueError:
        pass
    # Loop through the chains and see if any have the identifier as their name
    for chain in chains:
        if chain.name == identifier[0].upper():
            return [chain]


def interpret_group_commands(my_sys, group_dict, command):
    """
    Takes in a system and a command and adds to the group dictionary for later group creation
    """
    # Check if the identifier is in the mols, chains, residues or atoms list
    all_dicts = [{_: 'c' for _ in chn_objs}, {_: 'r' for _ in res_objs}, {_: 'a' for _ in atom_objs}]
    # Set the identifier on em
    type_dict = {k: v for d in all_dicts for k, v in d.items()}
    # Check if the command is an atom command
    if command[0].lower() in type_dict and type_dict[command[0].lower()] == 'a':
        # interpret the command
        my_atoms = get_group_spheres(my_sys.balls, command[1:])
        # Check that the get_group_atoms function actually returned something
        if my_atoms is not None:
            # Add the atoms to the group dict
            group_dict['atoms'] += my_atoms
            return group_dict
    elif command[0].lower() in type_dict and type_dict[command[0].lower()] == 'r':
        # interpret the command
        my_resids = get_group_resids(my_sys.residues, command[1:])
        # Check that the get_group_atoms function actually returned something
        if my_resids is not None:
            # Add the atoms to the group dict
            group_dict['residues'] += my_resids
            return group_dict
    elif command[0].lower() in type_dict and type_dict[command[0].lower()] == 'c':
        # interpret the command
        my_chains = get_group_chains(my_sys.chains, command[1:])
        # Check that the get_group_atoms function actually returned something
        if my_chains is not None:
            # Add the atoms to the group dict
            group_dict['chains'] += my_chains
            return group_dict
    elif command[0].lower() in residue_names and len(command) >= 2:
        # First try the possibility of a number
        try:
            res_seq = float(command[1])
            for res in my_sys.residues:
                if res.name == residue_names[command[0].lower()] and res.seq == res_seq:
                    # Check for an atom identifier
                    if len(command) == 3 and command[2].upper() in residue_atoms[residue_names[command[0].lower()]]:
                        for atom_ndx in res.atoms:
                            atom = my_sys.balls.iloc[atom_ndx]

                            if atom['name'].strip() == command[2].upper().strip():
                                group_dict['atoms'].append(atom_ndx)
                                return group_dict
                    group_dict['residues'].append(res)
                    return group_dict
        except ValueError:
            # If no residue identifier is given and the
            if command[1].upper() in residue_atoms[residue_names[command[0].lower()]]:
                # Go through the residues in the system
                for res in my_sys.residues:
                    if res.name == residue_names[command[0].lower()]:
                        # Check for the atom name within the residue
                        for atom in res.atoms:
                            if my_sys.balls['name'][atom].upper() == command[1].upper():
                                group_dict['atoms'].append(atom)
                return group_dict
    elif command[0].lower() in residue_names:
        # Loop through the residues
        for res in my_sys.residues:
            if res.name == residue_names[command[0].lower()]:
                group_dict['residues'].append(res)
        return group_dict
    return group_dict


def ggroup(my_sys, group_commands, settings=None):
    """
    Takes in the raw grouping commands and makes the groups. No networks solved, just makes them.
    """
    # First check if we are comparing then we will call this same function with the settings altered
    if settings is not None and type(settings['net_type']) == list and settings['net_type'][0] == 'com':
        # Make a copy and set the network types for the two groups
        settings1, settings2 = settings.copy(), settings.copy()
        settings1['net_type'], settings2['net_type'] = settings['net_type'][1], settings['net_type'][2]
        # Make the group
        ggroup(my_sys, group_commands, settings1)
        ggroup(my_sys, group_commands, settings2)
        # Return so we dont keep making groups
        return
    # First case: if no groups are entered, then make the standard group (no sol for 'mol' or coarse and all for 'foam')
    if len(group_commands) == 0:
        # Given a foam, the group is going to be the whole set of spheres
        if my_sys.type == 'foam':
            # Make the foam group
            my_sys.groups = [Group(my_sys, name=my_sys.name + '_Network', atoms=my_sys.balls['num'].to_list(), settings=settings)]
        # If the given system is not a foam add only the residues to hold out the sol atoms
        else:
            my_sys.groups = [Group(my_sys, name=my_sys.name + '_Network', residues=my_sys.residues.copy(), settings=settings)]
        return
    # First check if there are specific names without identifiers, no sol, full
    if group_commands[0][0] in noSOL_objs:
        my_sys.groups = [Group(my_sys, name=my_sys.name + '_Network', residues=my_sys.residues.copy(), settings=settings)]
        return
    if group_commands[0][0] in full_objs:
        my_sys.groups = [Group(my_sys, name=my_sys.name + '_all_atoms_network', atoms=my_sys.balls['num'].to_list(), settings=settings)]
        return
    my_sys.groups = []
    # Loop through the names and identifiers
    for _ in group_commands:
        # Get the list of lists of groups
        my_grp_cmnds = group_commands[_]
        # Make a dictionary of items to add to
        group_dict = {'atoms': [], 'chains': [], 'residues': []}
        # Go through the individual grouping commands
        for sub_group in my_grp_cmnds:
            # Interpret the sub_group
            group_dict = interpret_group_commands(my_sys, group_dict, sub_group)
        name = '_and_'.join(['_'.join(_) for _ in my_grp_cmnds])
        # Finally make the group
        if group_dict is not None and sum([len(group_dict[_]) for _ in group_dict]) > 0:
            my_sys.groups.append(Group(my_sys, name=name, settings=settings, residues=group_dict['residues'],
                                       atoms=group_dict['atoms'], chains=group_dict['chains']))
