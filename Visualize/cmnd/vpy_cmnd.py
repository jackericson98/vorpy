import System.sys_objs.residue
from Visualize.cmnd.export import *
from Visualize.cmnd.set import *
from System.Group.group import Group


def load_atom_file(my_sys):
    # Keep asking for an atom file till one is loaded
    while True:
        # Set up the prompt
        usr_file = input("atom file >>>   ")
        # Check if the user entered load first
        usr_file = usr_file.split()
        if len(usr_file) > 1:
            usr_file = usr_file[1]
        elif len(usr_file) == 1:
            usr_file = usr_file[0]
        else:
            my_num = np.random.randint(8)
            usr_file = ['Na5', 'EDTA_Mg', 'cambrin', 'hairpin', 'DB1976', 'Na7', 'protein_ligand_complex', 'Complex_frame1'][my_num]
        # Check if the full path was loaded
        if os.path.exists(usr_file) and usr_file[-3:] in {"pdb", "gro", "mol", "cif"}:
            file_path = usr_file
        # Check if the file was loaded without the directory
        elif os.path.exists("./Data/test_data/" + usr_file):
            file_path = os.getcwd() + "/Data/test_data/" + usr_file
        # Check if the path exists without the directory or the extension
        elif os.path.exists("./Data/test_data/" + usr_file + ".pdb"):
            file_path = os.getcwd() + "/Data/test_data/" + usr_file + ".pdb"
        # Else try again
        else:
            print("{} is not a valid input file".format(usr_file))
            continue
        # Create the system and return
        my_sys.load_sys(file=file_path)
        my_sys.net = Network(sys=my_sys, atoms=my_sys.atoms)
        return


def load_another_file(my_sys):
    # Ask the usr what type of file they want to load
    file_type = input(
        "Enter file type: 1. Network (.csv)  2. Voronota Balls (.txt)  3. Voronota Vertices (.txt)  "
        "4. GROMACS index (.ndx)\nfile type (1-4) >>>   ")
    # If the input is a network file load it into the system
    if file_type.lower() in {"1", "1.", "one", "un", "uno", "1 "}:
        # Print the message that the network file type has been chosen
        print(
            "Network file type selected. Please enter a vorpy generated network file address "
            "(extension .csv) for \'{}\'".format(my_sys.name))
        # Ask the user to add the input file for the system name
        my_net_file = input("file address (.csv) >>>   ")
        # Check that the file is correct
        if my_net_file[-3:] == 'csv' and os.path.exists(my_net_file):
            # Load the network file
            my_sys.load_net(my_net_file)
        else:
            print("Bad file")
            return False
    # If the input is to load a ball file
    elif file_type.lower() in {"2", "2.", "two", "to", "too", "dos", "du", "due"}:
        # Print the message that the voronota balls file type has been chosen
        print("Voronota balls file type selected. Please enter a Voronota generated balls full file address "
              "(extension .txt) for \'{}\'".format(my_sys.name))
        # Ask the user to add the input file for the system name
        my_ball_file = input("file address (.txt) >>>   ".format(my_sys.name))
        # Check that the file is correct
        if my_ball_file[-3:] == 'txt' and os.path.exists(my_ball_file):
            # Load the network file
            my_sys.ball_file = my_ball_file
        else:
            print("Bad file")
            return False
    # If the input is to load a ball file
    elif file_type.lower() in {"3", "3.", "three", "tre", "tres"}:
        # Print the message that the voronota balls file type has been chosen
        print("Voronota vertices file type selected. Please enter a Voronota generated vertices file address "
              "(extension .txt) for \'{}\'".format(my_sys.name))
        # Ask the user to add the input file for the system name
        my_vert_file = input("file address (.txt) >>>   ".format(my_sys.name))
        # Check that the file is correct
        if my_vert_file[-3:] == 'txt' and os.path.exists(my_vert_file):
            # Load the network file
            my_sys.vert_file = my_vert_file
        else:
            print("Bad file")
            return False
    # If the input is to load a ball file
    elif file_type.lower() in {"4", "4.", "four", "quattro", "for", "4 "}:
        # Print the message that the voronota balls file type has been chosen
        print("GROMACS index file type selected. Please enter a GROMACS generated index file address "
              "(extension .ndx) for \'{}\'".format(my_sys.name))
        # Ask the user to add the input file for the system name
        my_ndx_file = input("file address (.ndx) >>>   ".format(my_sys.name))
        # Check that the file is correct
        if my_ndx_file[-3:] == 'ndx' and os.path.exists(my_ndx_file):
            # Load the network file
            my_sys.load_ndx(file=my_ndx_file)
        else:
            print("Bad file")
            return False
    else:
        print("Bad Number")
        return False
    return True


def create_group(my_sys, usr_npt):
    # Check for basic inputs
    if usr_npt[0].lower() == 'f':
        return Group(sys=my_sys, atoms=my_sys.atoms, name="{}_full".format(my_sys.name))
    # Check for no sol
    elif usr_npt[0].lower() == 'ns':
        return Group(sys=my_sys, chains=my_sys.chains, name=my_sys.name + "_no_SOL")

    # Create the object and index variables
    my_obj, my_ndx = None, None
    # User only input "export"
    if len(usr_npt) == 0:
        # Tell the user to pick an object and an index
        my_obj = get_obj(sys=my_sys)
        my_ndx = get_ndx(sys=my_sys, obj=my_obj)
    # User entered "export obj" and needs an index
    elif len(usr_npt) == 1:
        # Add a check for system
        # Check the object provided by the user
        my_obj = get_obj(sys=my_sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=my_sys, obj=my_obj)
    # If the user input an object and an index of their own
    elif len(usr_npt) >= 2:
        # Check the object
        my_obj = get_obj(sys=my_sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=my_sys, obj=my_obj, ndx_npt=usr_npt[1])
    # Get the group information
    obj_ndx = ['m', 'r', 'a', 'n'].index(my_obj)
    obj_list = [my_sys.chains, my_sys.residues, my_sys.atoms, my_sys.ndxs][obj_ndx]
    name_prfx = ['mol', 'resid', 'atom', 'ndx'][obj_ndx]
    my_list, name = None, None
    # Get the slice and name of the group
    if my_ndx is None:
        return
    elif len(my_ndx) == 1:
        my_list = [obj_list[my_ndx[0]]]
        name = name_prfx + '_' + str(my_ndx[0])
    elif len(my_ndx) <= 2:
        my_list = obj_list[max(0, my_ndx[0]):min(len(obj_list), my_ndx[1] + 1)]
        name = name_prfx + 's_' + str(my_ndx[0]) + '_' + str(my_ndx[1])
    # Create the group
    npt_list = [None] * 4
    npt_list[obj_ndx] = my_list
    return Group(sys=my_sys, chains=npt_list[0], residues=npt_list[1], atoms=npt_list[2], indices=npt_list[3],
                 name=name)


def vorpy(my_sys):
    """
    Main function that is looped. Checks the inputs and runs the correct functions
    :return:
    """
    print("Welcome to vorpy. For assistance type \'h\'. To quit type \'q\'")
    # Load the initial input file
    load_atom_file(my_sys)
    # Allow the user to keep loading files
    while True:
        # Ask the user if they have another file to load
        load_another = input("add files >>>   ")
        # Check if load another is requested
        if load_another.lower() in ns + dones + ['']:
            break
        # Give the user the interface to load files
        good_file = load_another_file(my_sys=my_sys)
        if good_file:
            continue
    # Get the number of atoms in the default grouping (if sol dne all atoms)
    atom_len = len(my_sys.atoms) - len(my_sys.sol.atoms) if my_sys.sol is not None else len(my_sys.atoms)
    # Print the default grouping information for the system
    print("Default group: {} atoms, {} residue{}, {} chain{}".format(atom_len, len(my_sys.residues), 's' if len(my_sys.residues) > 1 else '', len(my_sys.chains), 's' if (len(my_sys.chains) > 1) else ''))
    # Start the grouping loop
    while True:
        # Get an initial grouping input
        usr_npt = input("new group >>>   ")
        # Split the user input
        usr_npt = usr_npt.split()
        # Check that the user's input is valid
        if len(usr_npt) == 0 or usr_npt[0].lower() in my_objects:
            break
        elif usr_npt[0].lower() in show_cmds:
            show(my_sys, usr_npt)
            continue
        # Tell the user they f'd up
        print("Bad input")
    if len(usr_npt) == 0 or usr_npt[0] in ns:
        my_group = Group(sys=my_sys, residues=my_sys.residues, name=my_sys.name)
    else:
        # Create the group
        my_group = create_group(my_sys=my_sys, usr_npt=usr_npt)
        if my_group is not None:
            print("{} group created - {} atoms, {} residues, {} chains".format(my_group.name, len(my_group.atoms),
                                                                           len(my_group.residues), len(my_group.chains)))
    # Check if the network has been loaded
    if my_sys.net_file is None and my_sys.vert_file is None:
        net = my_sys.net
        # Keep asking for a setting to change
        while True:
            # Print the default settings
            print(u"Default settings: net type = {}, surf res = {:.2f} \u208B,  max vert  = {:.2f} \u208B,  "
                  u"box multiplier = {:.2f} x".format(my_group.sys.net.type, net.surf_res, net.max_vert, net.box_size))
            # Print the build settings and see if the user wants to change anything
            change_settings = input("alter set >>>   ")
            change_settings = change_settings.split()
            # If the user wants to change the settings:
            if len(usr_npt) == 0:
                break
            elif change_settings[0].lower() in ys:
                sett(my_sys, ["set"], vorpy2_set=True)
            # If the user changes the settings here, insert the inp-ut into the sett function
            elif change_settings[0].lower() in my_settings:
                sett(my_sys, change_settings, vorpy2_set=True)
            # If the user input is not a good one let them go again
            elif change_settings[0].lower() in ns + [""] + dones:
                break
        # Build the group
        my_sys.net.build(my_group=my_group, print_actions=True)
    # Check if both voronota files have been loaded
    elif my_sys.ball_file is not None and my_sys.vert_file is not None:
        my_sys.load_verts(file=my_sys.vert_file, vta_ball_file=my_sys.ball_file)

    # Export
    export(my_sys, usr_npt="e", my_group=my_group)

    # Exporting process
    while True:
        # Check if the user wants to export files
        export_files = input("export more files for {}? (y/n) >>>   ".format(my_group.name))
        # If the user wants to export files for the given group start the export process for the current group
        if export_files.lower() in ys:
            # Export
            export(my_sys, usr_npt="e", my_group=my_group)
        elif export_files.lower() in ns + quits:
            # Ask if the user wants to export another group
            export_another_group = input("export another group for {}? (y/n) >>>   ".format(my_sys.name))
            export_another_group = export_another_group.split()
            my_new_group = None
            if export_another_group[0] in ys:
                my_new_group = group(my_sys, [])
            elif export_another_group[0] in my_objects:
                my_new_group = group(my_sys, usr_npt)
            elif export_another_group in ns + quits:
                return
            my_group = my_new_group
            my_sys.net.build(my_group=my_group, print_actions=True)
