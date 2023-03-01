from System.system import *
from Visualize.commands.export import *
from Visualize.commands.load import *
from Visualize.commands.set import *


def load_atom_file():
    global sys
    # Keep asking for an atom file till one is loaded
    while True:
        # Set up the prompt
        usr_file = input("atom file >>>   ")
        # Check if the user entered load first
        usr_file = usr_file.split()
        if len(usr_file) > 1:
            usr_file = usr_file[1]
        else:
            usr_file = usr_file[0]
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
        sys = System(file=file_path)
        sys.net = Network(sys=sys, atoms=sys.atoms)
        print(sys.name + " loaded - {} atoms, {} chains, {} residues"
              .format(len(sys.atoms), len(sys.mols), len(sys.residues)))
        return


def load_another_file():
    # Load the global system variable
    global sys
    # Ask the usr what type of file they want to load
    file_type = input(
        "Enter file type: 1. Network (.csv)  2. Voronota Balls (.txt)  3. Voronota Vertices (.txt)  4. GROMACS index (.ndx)\nfile type (1-4) >>>   ")
    # If the input is a network file load it into the system
    if file_type.lower() in {"1", "1.", "one", "un", "uno", "1 "}:
        # Print the message that the network file type has been chosen
        print(
            "Network file type selected. Please enter a vorpy generated network file address (extension .csv) for \'{}\'".format(
                sys.name))
        # Ask the user to add the input file for the system name
        my_net_file = input("file address (.csv) >>>   ")
        # Check that the file is correct
        if my_net_file[-3:] == 'csv' and os.path.exists(my_net_file):
            # Load the network file
            sys.load_net(my_net_file)
        else:
            print("Bad file")
            return False
    # If the input is to load a ball file
    elif file_type.lower() in {"2", "2.", "two", "to", "too", "dos", "du", "due"}:
        # Print the message that the voronota balls file type has been chosen
        print(
            "Voronota balls file type selected. Please enter a Voronota generated balls full file address (extension .txt) for \'{}\'".format(
                sys.name))
        # Ask the user to add the input file for the system name
        my_ball_file = input("file address (.txt) >>>   ".format(sys.name))
        # Check that the file is correct
        if my_ball_file[-3:] == 'txt' and os.path.exists(my_ball_file):
            # Load the network file
            sys.ball_file = my_ball_file
        else:
            print("Bad file")
            return False
    # If the input is to load a ball file
    elif file_type.lower() in {"3", "3.", "three", "tre", "tres"}:
        # Print the message that the voronota balls file type has been chosen
        print(
            "Voronota vertices file type selected. Please enter a Voronota generated vertices file address (extension .txt) for \'{}\'".format(
                sys.name))
        # Ask the user to add the input file for the system name
        my_vert_file = input("file address (.txt) >>>   ".format(sys.name))
        # Check that the file is correct
        if my_vert_file[-3:] == 'txt' and os.path.exists(my_vert_file):
            # Load the network file
            sys.vert_file = my_vert_file
        else:
            print("Bad file")
            return False
    # If the input is to load a ball file
    elif file_type.lower() in {"4", "4.", "four", "quattro", "for", "4 "}:
        # Print the message that the voronota balls file type has been chosen
        print(
            "GROMACS index file type selected. Please enter a GROMACS generated index file address (extension .ndx) for \'{}\'".format(
                sys.name))
        # Ask the user to add the input file for the system name
        my_ndx_file = input("file address (.ndx) >>>   ".format(sys.name))
        # Check that the file is correct
        if my_ndx_file[-3:] == 'ndx' and os.path.exists(my_ndx_file):
            # Load the network file
            sys.load_ndx(file=my_ndx_file)
        else:
            print("Bad file")
            return False
    else:
        print("Bad Number")
        return False
    return True


def create_group(usr_npt):
    # Get the system variable
    global sys
    # Check for basic inputs
    if usr_npt[0].lower() == 'f':
        return Group(sys=sys, atoms=sys.atoms, name="{}_full".format(sys.name))
    # Check for no sol
    elif usr_npt[0].lower() == 'ns':
        return Group(sys=sys, mols=sys.mols[:-1], name=sys.name + "_no_SOL")

    # Create the object and index variables
    my_obj, my_ndx = None, None
    # User only input "export"
    if len(usr_npt) == 0:
        # Tell the user to pick an object and an index
        my_obj = get_obj(sys=sys)
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # User entered "export obj" and needs an index
    elif len(usr_npt) == 1:
        # Add a check for system
        # Check the object provided by the user
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # If the user input an object and an index of their own
    elif len(usr_npt) >= 2:
        # Check the object
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj, ndx_npt=usr_npt[1])
    # Get the group information
    obj_ndx = ['m', 'r', 'a', 'n'].index(my_obj)
    obj_list = [sys.mols, sys.residues, sys.atoms, sys.ndxs][obj_ndx]
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
    return Group(sys=sys, mols=npt_list[0], residues=npt_list[1], atoms=npt_list[2], indices=npt_list[3], name=name)


def vorpy():
    """
    Main function that is looped. Checks the inputs and runs the correct functions
    :return:
    """
    print("Welcome to vorpy. For assistance type \'h\'. To quit type \'q\'")
    global sys
    # Load the initial input file
    load_atom_file()
    # Allow the user to keep loading files
    while True:
        # Ask the user if they have another file to load
        load_another = input("load another file? (y/n) >>>   ")
        # Check if load another is requested
        if load_another.lower() in ys:
            # Give the user the interface to load files
            good_file = load_another_file()
            if good_file:
                continue
        elif load_another.lower() in ns:
            break
    # Start the grouping loop
    while True:
        # Get an initial grouping input
        usr_npt = input(
            "Create a group. (Full = \'f\', No SOL = \'ns\', mol = \'m\', res = \'r\', atom = \'a\', index = \'i\')\ngroup >>>   ")
        # Split the user input
        usr_npt = usr_npt.split()
        # Check that the user's input is valid
        if usr_npt[0].lower() in my_objects:
            break
        elif usr_npt[0].lower() in show_cmds:
            show(sys, usr_npt)
            continue
        # Tell the user they f'd up
        print("Bad input")
    # Create the group
    while True:
        my_group = create_group(usr_npt)
        if my_group is not None:
            # print("{} group created - {} atoms, {} residues, {} chains".format(my_group.name, len(my_group.atoms),
            #                                                                    len(my_group.resids), len(my_group.mols)))
            break
    # Check if the network has been loaded
    if sys.net_file is None and sys.vert_file is None:
        # Keep asking for a setting to change
        while True:
            # Print the default settings
            print(
                u"{} Build Settings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_multi = {:.2f} x,  build_surfs = {}, "
                u"flat_surfs = {}".format(my_group.name, sys.net.surf_res, sys.net.max_vert, sys.net.box_size,
                                          sys.net.build_surfs,
                                          sys.net.flat_Del))
            # Print the build settings and see if the user wants to change anything
            change_settings = input("change settings? (y/n) >>>   ")
            change_settings = change_settings.split()
            # If the user wants to change the settings:
            if change_settings[0].lower() in ys:
                sett(sys, ["set"], vorpy2_set=True)
            # If the user changes the settings here, insert the inp-ut into the sett function
            elif change_settings[0].lower() in my_settings:
                sett(sys, change_settings, vorpy2_set=True)
            # If the user input is not a good one let them go again
            elif change_settings[0].lower() in ns:
                break
        # Build the group
        sys.net.build(my_group=my_group)
    # Check if both voronota files have been loaded
    elif sys.ball_file is not None and sys.vert_file is not None:
        sys.load_verts(file=sys.vert_file, vta_ball_file=sys.ball_file)

    # Exporting process
    while True:
        # Check if the user wants to export files
        export_files = input("export files for {}? (y/n) >>>   ".format(my_group.name))
        # If the user wants to export files for the given group start the export process for the current group
        if export_files.lower() in ys:
            # Export
            export(sys, usr_npt="e", my_group=my_group)
        elif export_files.lower() in ns + quits:
            # Ask if the user wants to export another group
            export_another_group = input("export another group for {}? (y/n) >>>   ".format(sys.name))
            export_another_group = export_another_group.split()
            my_new_group = None
            if export_another_group[0] in ys:
                my_new_group = group(sys, [])
            elif export_another_group[0] in my_objects:
                my_new_group = group(sys, usr_npt)
            elif export_another_group in ns + quits:
                return
            my_group = my_new_group
            sys.net.build(my_group=my_group)

if __name__ == '__main__':
    # Welcome introduction
    sys = System()
    # Run vorpy
    vorpy()
