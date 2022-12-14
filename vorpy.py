from System.system import *

# Set the running variable
sys = None
loading = True
base_file, net_file, vert_file, ndx_file = "", "", "", ""
#
# clear_frame = "\r" + (" " * 100 + "\n") * 10

# Set up the yes and no variables
ys = ['y', 'yes', 'ya', 'yeet', 'yur', 'yoint', 'uhu', 'yup', 'jess', 'affirmative', 'yass']
ns = ['n', 'no', 'naur', 'nope', 'nonya', 'nope', 'nien', 'nada']
dones = ['', 'done', 'd', 'finished', 'finito', 'complete', 'doneso', 'don', 'fin', 'keep movin bruh', 'mf move on',
         'you still here?', 'go on', 'goo\'oon']




# Let the user load their files
while loading:

    # Reset the file type loaded
    sys_loaded, net_loaded, verts_loaded, ndx_loaded = False, False, False, False

    # Create the new loading header
    print("\r\n My Files:                                          \n\n\r" \
          "  1. System File  :  {}                                \n" \
          "  2. Network File :  {}                                \n" \
          "  3. Vertex File  :  {}                                \n" \
          "  4. Index File   :  {}                                \n"
          .format(base_file, net_file, vert_file, ndx_file), end="")

    # Prompt the user to load files
    new_file = input("\n\nLoad Files (when finished type \'done\'): \n\n  File Address:   ")

    # Check to see what type of file it is
    if new_file[-3:] == 'pdb' or new_file[-3:] == 'mol' or new_file[-3:] == 'gro' or new_file[-3:] == 'cif':
        sys_loaded = True
        if sys is not None:
            reset_sys = input("{} System is already loaded. Would you like to replace this system with: {} (Y/N)".format(sys.name, new_file))
            if reset_sys.lower() in ys:
                sys = System(new_file)
        else:
            sys = System(new_file)
            base_file = new_file

    # If the loaded file is a vertex or network file load them accordingly
    elif new_file[-3:] == 'txt':
        # If the new file is a vertex file load it
        if new_file[-9:-4] == 'verts':
            verts_loaded = True
            # If a vertex file has already been loaded make sure the user wants to load it if not load it
            if vert_file != "":
                replace_vert_file = input("There is a vertex file already loaded: {}\n "
                                          "Would you like to replace it? (Y/N):   ".format(vert_file))
                if replace_vert_file in ys or replace_vert_file in dones:
                    vert_file = new_file
            else:
                vert_file = new_file

        # If the new file is a network file load it
        elif new_file[-11:-4] == 'network':
            net_loaded = True
            # If a vertex file has already been loaded make sure the user wants to load it if not load it
            if net_file != "":
                replace_net_file = input("There is a Network file already loaded: {}\n "
                                          "Would you like to replace it? (Y/N):   ".format(net_file))
                if replace_net_file in ys or replace_net_file in dones:
                    net_file = new_file
            else:
                net_file = new_file

    elif new_file[-3:] == 'ndx':
        ndx_file = new_file

    # Move on if the user is done with the loading
    elif new_file == '' or new_file.lower() in dones:
        # Check to see that there is at least a system loaded
        if sys is None or sys.base_file is None:
            print("\rPlease load a system file before continuing", end="")
            continue
        else:
            # Make sure the user wants to continue and didn't accidentally hit the enter key
            are_you_sure = input("\nAre you finished loading files? (Y/N):   ")
            if are_you_sure.lower() in ys or are_you_sure.lower() in dones:
                loading = False
    # In all other case print an error and give the user a chance to try again
    else:
        print("\r\n\n{} Is not a valid input. Please provide the full address for one of the following file types: \n .pdb, "
              ".mol, .cif, .gro, _verts.txt (vorpy created), _network.txt (vorpy created), .ndx (GROMACS created)"
              .format(new_file), end="")

    # Load the new file
    if sys is not None:
        if verts_loaded:
            sys.load_verts(vert_file)
        elif net_loaded:
            sys.load_net(net_file)
        elif ndx_loaded:
            sys.load_ndx(ndx_file)
        elif sys_loaded:
            if net_file != "":
                sys.load_net(net_file)
            if vert_file != "":
                sys.load_verts(vert_file)
            if ndx_file != "":
                sys.load_ndx(ndx_file)

# Check to see if the network needs to be built
prepping = True
if sys.net is not None:
    print(" ")
    prepping = False

# Keep running the program until the user says not to
while prepping:

    # Create a network if one has not been created
    if sys.net is None:
        sys.net = Network(atoms=sys.atoms, sys=sys)

    # Create the header string
    prepping_header = "\r\nMy Settings:                                             \n\n" \
                      "  1. Surface Resolution    : {}                                \n" \
                      "  2. Maximum Vertex Radius : {}                                \n" \
                      "  3. Box Size              : {}                                \n" \
                      "  4. Solute Vertices       : {}                                \n"\
        .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts)


    # Update the user with the current network build settings
    print("\r" + prepping_header, end="")

    # Ask the user if they want to change any settings
    change_settings = input("\nAre you happy with the current build settings? (Y/N):   ")

    # If the user is ready to move on let them settings are requested to be change
    if change_settings.lower() not in ns:
        ready_to_build = input("\nBuild the network? (Y/N):   ")
        if ready_to_build.lower() in ys or ready_to_build.lower() in dones:
            prepping = False
        continue

    # Set up the settings changing tracking variable
    changing_settings = True
    # Create the settings changing loop
    while changing_settings:

        # Reset the setting variables
        settings = [None, None, None, None]
        # Ask the user which of the settings they want to change
        which_setting = input("\nWhich setting would you like to change (enter a number 1-4 or \'done\' to exit):   ")
        # Get the setting number through a try statement
        try:
            my_setting = int(which_setting)
        except ValueError:
            my_setting = 0

        # Change the settings
        if which_setting in dones:
            changing_settings = False
            continue
        elif my_setting == 0 or my_setting > 5:
            print("\'{}\' is not a valid input. Please try again".format(which_setting))
            continue
        elif my_setting == 1:
            settings[0] = input(u"\nSet the Surface Resolution (0.01 - 3.00 \u212B):  ")
        elif my_setting == 2:
            settings[1] = input(u"\nSet the Maximum Vertex Radius (1.0 - 20.0 \u212B):   ")
        elif my_setting == 3:
            settings[2] = input(u"\nSet the Retaining Box Multiplier (1.0 - 10.0 x):   ")
        elif my_setting == 4:
            settings[3] = input("\nCalculate the solute's vertices?: (Y/N)   ")
        # Go through the settings to see if they have been changed
        for i in range(3):
            # Try to convert the input to a float variable
            try:
                settings[i] = float(settings[i])
            except TypeError:
                pass
        # Check the solute vertices setting
        if settings[3] is not None:
            if settings[3].lower() in ys:
                settings[3] = True
            elif settings[3].lower() in ns:
                settings[3] = False

        # Check to see if any settings made it through
        if type(settings[0]) == float:
            sys.net.surf_res = settings[0]
        elif type(settings[1]) ==  float:
            sys.net.max_vert = settings[1]
        elif type(settings[2]) == float:
            sys.net.box_size = settings[2]
        elif type(settings[3]) == bool:
            sys.net.sol_verts = settings[3]
            print("HERE", sys.net.sol_verts)
        else:
            print("\rNo settings changed")

        # Create the header string
        prepping_header = "\r\nMy Settings:                                             \n\n" \
                          "  1. Surface Resolution    : {}                                \n" \
                          "  2. Maximum Vertex Radius : {}                                \n" \
                          "  3. Box Size              : {}                                \n" \
                          "  4. Solute Vertices       : {}                                \n" \
            .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts)
        print(prepping_header, end=" ")

        # Ask the user if they want to change another setting
        change_more_settings = input("\nChange another setting? (Y/N):   ")
        # Check their response
        if change_more_settings.lower() in ns:
            # Are you ready to build the network
            build_the_network = input("\nBuild the network? (Y/N):   ")
            if build_the_network.lower() in ys or build_the_network.lower() in dones:
                changing_settings = False
                prepping = False
        elif change_more_settings.lower() not in ys:
            print("\nIll take \'{}\' as a yes\n".format(change_more_settings), end="")

print("")
if sys.net_file is None:
    # Build the network
    sys.build_network()


def create_group(my_sys, name=None, my_group=None):

    # Create the group
    if my_group is None:
        my_group = Group(net=my_sys, name=name)
        my_group.mol_names, my_group.res_names, my_group.atom_names, my_group.ndx_names = [], [], [], []

    # Create the header for the group
    sys_header = "\r\n" \
                 "System elements: \n\n" \
                 "  1. Molecules:  {}\n".format([my_sys.mol_names[k] for k in range(len(my_sys.mols)) if k < 10]) + \
                 "  2. Residues :  {}\n".format([my_sys.res_names[k] for k in range(len(my_sys.residues)) if k < 10]) + \
                 "  3. Atoms    :  {}\n".format([my_sys.atom_names[k] for k in range(len(my_sys.atoms)) if k < 10]) + \
                 "  4. Indices  :  {}\n\n".format([my_sys.ndx_names[k] for k in range(len(my_sys.ndxs)) if k < 10]) + \
                 "Group: {}  {}  {}  {}\n\n"\
                     .format(my_group.mol_names, my_group.res_names, my_group.atom_names, my_group.ndx_names)
    # Print the system headers
    print(sys_header, sys_header)

    # Create the add elements tracking variable
    adding_elements = True
    # Keep allowing the user to add elements
    while adding_elements:

        # Ask the user what they want to add to the group
        add_elem = input("Choose an element to add (Enter 1-4, \'all\' for the full system or \'done\' to quit):   ")
        # Check to see if the user provided an integer or not
        try:
            change_elem = int(add_elem)
        except ValueError:
            change_elem = None
        # The general skip case
        if add_elem in dones or change_elem is None:
            pass
        # If the user wants to add molecules
        elif change_elem == 1:
            # Create the string variable to print the systems molecules
            my_sys_mols = "\r\nMolecules: \n"
            for j in range(len(sys.mols)):
                my_sys_mols += str(j + 1) + ". " + my_sys.mol_names[j] + "\n"
            print(my_sys_mols)
            # Ask what molecules the user wants to add to their group
            what_mols = input("Which molecule would you like to add? (enter a number 1 - {}):   "
                              .format(str(len(my_sys.mols) + 1)))
            # Figure out what molecule the user wants
            try:
                my_mol = int(what_mols) - 1
            except ValueError:
                my_mol = None
            # If the chosen molecule is not None, add the molecule to the group
            if my_mol is not None and my_mol < len(my_sys.mols):
                my_group.mol_names.append(my_sys.mol_names[my_mol])
                my_group.add_sele(my_sys.mols[my_mol], my_sys.mol_names[my_mol])
        # If the user wants to add residues
        elif change_elem == 2:
            # Create the string variable to print the systems residues
            my_sys_resids = "\r\nResidues : \n"
            for j in range(len(sys.residues)):
                my_sys_resids += str(j + 1) + ". " + my_sys.res_names[j] + "\n"
            print(my_sys_resids)
            # Ask what residues the user wants to add to their group
            what_resids = input("Which residue would you like to add? (enter a number 1 - {}):   "
                              .format(str(len(my_sys.residues))))
            # Figure out what molecule the user wants
            try:
                my_res = int(what_resids) - 1
            except ValueError:
                my_res = None
            # If the chosen residue is not None, add the residue to the group
            if my_res is not None and my_res < len(my_sys.residues):
                my_group.res_names.append(my_sys.res_names[my_res])
                my_group.add_sele(my_sys.residues[my_res], my_sys.res_names[my_res])
        # If the user wants to add atoms to the group
        elif change_elem == 3:
            # Create the string variable to print the systems atoms
            my_sys_atoms = "\r\nAtoms   : \n"
            for j in range(len(sys.atoms)):
                my_sys_atoms += str(j + 1) + ". " + my_sys.atom_names[j] + " - " + my_sys.atoms[j].res + " " + \
                                my_sys.atoms[j].res_seq + " " + "\n"
            print(my_sys_atoms)
            # Ask what atoms the user wants to add to their group
            what_atoms = input("Which atom would you like to add? (enter a number 1 - {}):   "
                              .format(str(len(my_sys.atoms))))
            # Figure out what atom the user wants
            try:
                my_atom = int(what_atoms) - 1
            except ValueError:
                my_atom = None
            # If the chosen atom is not None, add the atom to the group
            if my_atom is not None and my_atom < len(my_sys.atoms):
                my_group.atom_names.append(my_sys.atom_names[my_atom])
                my_group.add_sele([my_sys.atoms[my_atom]], my_sys.atom_names[my_atom])
        # If the user wants to add an index
        elif change_elem == 4:
            # Create the string variable to print the systems molecules
            my_sys_ndxs = "\r\nIndices : \n"
            for j in range(len(sys.mols)):
                my_sys_ndxs += str(j + 1) + ". " + my_sys.ndx_names[j] + "\n"
            print(my_sys_ndxs)
            # Ask what molecules the user wants to add to their group
            what_ndx = input("Which index would you like to add? (enter a number 1 - {}):   "
                              .format(str(len(my_sys.ndxs) + 1)))
            # Figure out what molecule the user wants
            try:
                my_ndx = int(what_ndx)
            except ValueError:
                my_ndx = None
            # If the chosen molecule is not None, add the molecule to the group
            if my_ndx is not None and my_ndx < len(my_sys.ndx_names):
                my_group.ndx_names.append(my_sys.ndx_names[my_ndx - 1])
                my_group.add_sele(my_sys.ndxs[my_ndx - 1], my_sys.ndx_names[my_ndx - 1])
        else:
            print("{} is not a valid entry".format(add_elem))

        # Create the header for the group
        group_header = "\r\n" \
                       "My Group: \n\n" \
                       "  1. Molecules           :  {}\n" \
                       "  2. Residues            :  {}\n" \
                       "  3. Atoms               :  {}\n" \
                       "  4. Indices             :  {}\n\n" \
            .format(my_group.mol_names, my_group.res_names, my_group.atom_names, my_group.ndx_names)
        print(group_header, end="")

        # Check if the user wants to add another element
        add_another_element = input("Add another element? (Y/N):   ")
        # If they specify no leave the adding elements loop
        if add_another_element.lower() in ns:
            adding_elements = False
    # Name the group
    if my_group.name is None:
        my_group.name = "_".join(my_group.mol_names + my_group.res_names + my_group.atom_names + my_group.ndx_names)
        my_group.name.replace(" ", "_")

    # After adding elements ask the user if they want to change the name of the group
    rename_group = input("\nGroup Name:  {}. Would you like to rename it? (Y/N):   ".format(my_group.name))
    if rename_group not in ns:
        new_group_name = input("New group name  (up to 10 characters, no spaces):   ")
        if new_group_name not in dones and new_group_name not in ns:
            my_group.name = new_group_name
            print("Group name changed to {}".format(my_group.name))
    # Lastly ask the user if they are finished making the group
    are_you_done = input("\nFinished making the group? (Y/N):   ")
    if are_you_done not in ns:
        my_sys.groups.append(my_group)
        my_sys.group_names.append(my_group.name)
    else:
        create_group(my_sys, my_group=my_group)



def get_sys_selects(my_sys):
    # Create the export objects header
    export_objects_header = "Export Objects: \n\n" \
                            "Prefixes (and lengths): i: Indices  = {}, g: Groups = {}, m: Molecules = {}), r: Residues = {}, a: Atoms = {}\n\n" \
                            "Instructions: Type the prefix of the object type followed by it's relative location in \n" \
                            "the list (e.g. the 30th atom = \'a30\'). To see a list of the objects in each \n" \
                            "classification type \'show\' prefix (e.g. \'show a\' to see all atoms). "\
        .format(len(my_sys.ndxs), len(my_sys.groups), len(my_sys.mols), len(my_sys.residues), len(my_sys.atoms))
    print(export_objects_header)
    # Set up the selection and it's name's variables
    my_selection, my_selection_name = None, None

    selecting = True
    while selecting:
        # Ask the user which objects they want exported
        export_select = input("\nPlease enter a selection for export:   ")
        # First check to see that the user entered a valid input
        if export_select[0].lower() in ['i', 'g', 's', 'm', 'r', 'a']:
            try:
                my_select_ndx = int(export_select[1:])
            except TypeError:
                print("{} is not a valid response".format(export_select))
                continue
        else:
            print("{} is not a valid response".format(export_select))
            continue

        # Get their selections
        if export_select[0].lower() == 'i' and my_select_ndx <= len(my_sys.ndxs):
            my_selection = my_sys.ndxs[my_select_ndx - 1]
            my_selection_name = my_sys.ndx_names[my_select_ndx - 1]
            selecting = False
        # Export groups
        elif export_select[0].lower() == 'g' and my_select_ndx <= len(my_sys.groups):
            my_selection = my_sys.groups[my_select_ndx - 1]
            my_selection_name = my_sys.group_names[my_select_ndx - 1]
            selecting = False
        # Export pre-determined system objects
        elif export_select[0].lower() == 's' and my_select_ndx < 3:
            if my_select_ndx == 1:
                my_selection = my_sys.atoms
                my_selection_name = "Full_Network"
            elif my_select_ndx == 2:
                my_selection = my_sys.sol
                my_selection_name = "No_SOL"
            selecting = False
        elif export_select[0].lower() == 'm' and my_select_ndx <= len(my_sys.mols):
            my_selection = my_sys.mols[my_select_ndx - 1]
            my_selection_name = my_sys.mol_names[my_select_ndx - 1]
            selecting = False
        elif export_select[0].lower() == 'r' and my_select_ndx <= len(my_sys.residues):
            my_selection = my_sys.residues[my_select_ndx - 1]
            my_selection_name = my_sys.res_names[my_select_ndx - 1]
            selecting = False
        elif export_select[0].lower() == 'a' and my_select_ndx <= len(my_sys.atoms):
            my_selection = [my_sys.atoms[my_select_ndx - 1]]
            my_selection_name = my_sys.atom_names[my_select_ndx - 1]
            selecting = False
        else:
            print(export_select, "is not a valid response. Please try again")

    return my_selection, my_selection_name


def export_sys_selects(my_sys):
    # Get the user's selections
    my_selection, my_selection_name = get_sys_selects(my_sys)


    # Create the exporting tracking variable
    exporting = True
    while exporting:
        # Create the group
        if type(my_selection) is list:
            my_selection = Group(atoms=my_selection, net=my_sys.net)
        my_selection.get_info()
        # Print the selected group 
        print("\n{} selected".format(my_selection_name))
        output_option = input("\n  1. All surfaces  2. Full Body (filled)  3. Full Body (empty)  4. Interface\n\n"
                               "Which of the above would you like to output? (choose 1 - 4):  ")
        # Check to make sure the user entered the correct value
        try:
            my_output_options = int(output_option)
        except ValueError:
            my_output_options = None
        # If the user wants all surfs
        if my_output_options == 1:
            # Make the surfaces file
            os.mkdir(my_sys.dir + "/Surfaces")
            os.chdir(my_sys.dir + "/Surfaces")
            # Export the surfaces one by one
            for k in range(len(sys.net.surfs)):
                # Export the
                surf = my_sys.net.surfs[k]
                # Get a random color
                my_color = np.random.rand(3)
                # Write each of the surfaces
                write_surfs([surf], "surf_" + str(surf.ndx[0]) + "_" + str(surf.ndx[1]), my_color)
            os.chdir(my_sys.dir)
        # If the user wants the body for their group
        elif my_output_options == 2:
            # Export the body
            write_surfs(my_selection.surfs, my_selection_name[0] + "_full", directory=my_sys.dir)
            print("Number of surfaces", len(my_selection.surfs))
        # If the user wants an empty body
        elif my_output_options == 3:
            # Export the body
            write_surfs(my_selection.body_surfs, my_selection_name[0] + "_empty", directory=my_sys.dir)
            print("Number of body surfaces", len(my_selection.body_surfs))
        # If the user wants an interface
        elif my_output_options == 4:
            # Finding interfaces tracking variable
            finding_interfaces = True
            while finding_interfaces:
                # Prompt the user to choose another selection
                print("Choose another selection for comparison\n")
                my_selection1, my_selection_name1 = get_sys_selects(my_sys)
                # Track whether the user wants to export the atoms around the interface
                export_iface_atoms_bool = False
                export_iface_atoms = input("Export interface atoms? (Y/N):   ")
                if export_iface_atoms in ys:
                    export_iface_atoms_bool = True
                # Change the new selection to a group object
                if type(my_selection1) is list:
                    my_selection1 = Group(net=my_sys.net, atoms=my_selection1, bff=my_selection, name=my_selection_name1)
                # Process the new group
                my_selection1.get_info()
                # Export the interface
                export_iface(groups=[my_selection, my_selection1], interface_atoms=export_iface_atoms_bool)
                # Ask the user if they want another interface
                find_another_iface = input("Calculate another interface with {}? (Y/N):   ".format(my_selection_name))
                if find_another_iface in ns:
                    finding_interfaces = False


# Create the analyzing variable
analyzing = True
# Set up the loop
while analyzing:
    # Get the data variables
    data = [sys.atoms, sys.mols, sys.residues, sys.net.verts, sys.net.edges, sys.net.surfs]
    a_len, m_len, r_len, v_len, e_len, s_len = [str(len(_)) for _ in data]

    # Create the analyzing header
    analyzing_header = "\r\n" + \
                       "My System:        " + " " * 25 +                                       "My Network:    \n\n" + \
                       "  Atoms        :  {}".format(a_len)        + " " * (25 - len(a_len)) + "  Vertices   :  {} \n" \
                       "  Molecules    :  {}".format(v_len, m_len) + " " * (25 - len(m_len)) + "  Edges      :  {} \n" \
                       "  Residues     :  {}".format(e_len, r_len) + " " * (25 - len(r_len)) + "  Surfaces   :  {} \n" \
                       "\nIndices  :  {}\nGroups   :  {}\n".format(s_len, sys.ndx_names, sys.group_names)


    # Clear the frame and print the system and network information
    print(analyzing_header, end="")

    # Print the System and network information
    create_new_group = input("\nCreate a new group? (Y/N):   ")
    # Check to see if the user wants to create a group
    if create_new_group.lower() not in ns:
        create_group(sys)
    # If the user inputted something not in the nos or the yeses lists
    elif create_new_group.lower() not in ns:
        print(create_new_group, "is not a valid response")
        continue

    # Check if the user wants to export any groups
    export_my_net_selects = input("\nReady to export? (Y/N):   ")

    # Check if the user wants to export selections
    if export_my_net_selects.lower() not in ns:
        export_sys_selects(sys)
    elif export_my_net_selects.lower() not in ns:
        print(export_my_net_selects, "is not a valid response")
        continue

    # Ask the user if they are finished
    quit_program = input("Quit? (Y/N)")

    # Check to see if the user wants to quit
    if quit_program.lower() in ns:
        continue
    elif quit_program in ys:
        print("Goodbye!")
        analyzing = False
    else:
        print(quit_program, "is not a valid response")
