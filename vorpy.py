from os import path
from System.system import *

"""
I want this to be a constantly running interface where the user can load, build, and export Voronoi network components. 
Has a running header that updates as commands are preformed. User can get system information based off of requests
"""



# Responses
ys = ['y', 'yes', 'ya', 'yeet', 'yur', 'yoint', 'uhu', 'yup', 'jess', 'affirmative', 'yass', '', 'yuss', 'yess',
      'yesss', 'yessss', 'yar', 'yuh', 'mhm']
ns = ['n', 'no', 'naur', 'nope', 'nonya', 'nope', 'nien', 'nada']
dones = ['done', 'd', 'finished', 'finito', 'complete', 'doneso', 'don', 'fin', 'keep movin bruh',
         'you still here?', 'go on', 'goo\'oon']
ands = ['&', 'and', 'nd', 'also', '+', '&&']
splitters = ['/', '-']

# General commands
quits = ['quit', 'q']
helps = ['h', 'help']

# Main commands
show_cmds = ['s', 'show', 'shw', 'sho', 'sh']
load_cmds = ['l', 'load', 'lod', 'laod', 'ld', 'lad', 'old']
set_cmds = ['st', 'set', 'assert', 'assign', 'make']
build_cmds = ['b', 'build', 'bld', 'bild', 'buld', 'bd']
group_cmds = ['g', 'group', 'gruop', 'grp', 'grup', 'grop', 'gp', 'gr']
export_cmds = ['e', 'export', 'xport', 'xprt', 'xpt', 'xp', 'expt', 'ext']

my_commands = quits + helps + show_cmds + load_cmds + set_cmds + build_cmds + group_cmds + export_cmds

# Objects
mol_objs = ['m', 'ms', 'molecule', 'molecules', 'mol', 'mols', 'ml', 'mls']
atom_objs = ['a', 'as', 'atom', 'atoms', 'at', 'ats', 'am', 'ams']
res_objs = ['r', 'rs', 'residue', 'residues', 'resid', 'resids', 'res', 'ress', 'reses', 'rdue', 'rdues']
ndx_objs = ['i', 'is',  'index', 'indexs', 'indexes', 'indices', 'ndx', 'ndxs', 'ndex']

my_objects = mol_objs + res_objs + atom_objs + ndx_objs

# Settings
surf_reses = ['surf_res', 'sr', 'surface_resolution', 'surface_res', 'surf_resolution']
max_verts = ['max_vert', 'mv', 'maximum_vertex', 'max_vertex', 'maximum_vert']
box_sizes = ['box_size', 'bs']
sol_vertses = ['sol_verts', 'sv']

my_settings = surf_reses + max_verts + box_sizes + sol_vertses


def create_header(sys):
    """
    Creates a header string to print in the terminal holding the loaded information
    :param sys:
    :return:
    """

    # Get the printable information

    # File strings
    my_file_strings = ["  System File :   {}".format(sys.base_file), "  Network File:   {}".format(sys.net_file),
                       "  Vertex File :   {}".format(sys.vert_file), "  Index File  :   {}".format(sys.ndx_file)]
    print("Files:\n")

    # System Strings
    my_sys_strings = ["System: {}".format(sys.name),
                      "  Atoms : {} - Molecules : {} - Residues : {}".format(len(sys.atoms), len(sys.mols), len(sys.residues)),
                      "  Additional info: {} ".format(sys.data[:min(len(sys.data) - 1, 5)])]

    # Network Strings
    my_net_strings = ["Network: "]
    if sys.net is not None:
        my_net_strings = ["Network: {}, {}, {}, {}".format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts)]

    # Index strings
    my_ndx_strings = ["Indexes: {}".format(sys.group_names)]

    print(my_file_strings, "\n\n", my_sys_strings, "\n\n", my_net_strings, "\n\n", my_ndx_strings)


def are_you_sure():
    ays = input("Are you sure? (y/n):   ")
    if ays in ys:
        return True
    return False


def invalid_input(string):
    if type(string) is list:
        string = " ".join(string)
    print("\'{}\' is not a valid input. try again or type \'h\' for help".format(string))


def get_help():
    """
    Shows a list of commands that the user has access to
    :return:
    """

    help_header = "Welcome to vorpy Help: ('h')"

    instructions_header = "Usage: Use a command and an object and its number (\'export mol 1\'), a setting and a value (\'set surf_res 0.1)\') or a\n" \
                          "       file (\'load /test_data/Na5.pdb\'). Use \'and\' to do multiple tasks or export interfaces (\'export mol 1 and group 3\')"



    commands_header = ["Commands:                                                                                                                   ",
                       "  1. load  : Loads file addresses for System (.pdb, .gro, .mol, .cif), Network (.txt), vertices (.txt) or index (.ndx) files",
                       "  2. set   : Sets build Settings with values (float, float, float, True)                                                    ",
                       "  3. build : Builds a Network for the System with the current settings                                                      ",
                       "  4. group : Groups together System objects                                                                                 ",
                       "  5. export: Exports System objects. Use \'and\' to export the interface between two groups (e.g. \'export ndx 1 and atom 3\')  ",
                       "  6. show  : Shows System elements in a given object category                                                               ",
                       "  7. quit  : Quits from the current process                                                                                 "]

    splitting_line = "--------------------------------------------------------------------------------------------------------------------------------"

    objects_header = ["Objects:                                           ",
                      "  1. mol : Molecule object from the current System ",
                      "           (Use the number or name of the Molecule)",
                      "  2. res : Residue object from the current System  ",
                      "  3. atom: Atom object from the current System     ",
                      "  4. ndx : Index loaded into the current System or ",
                      "           created by the user                     "]

    settings_header = ["Settings:                                                                   ",
                       "  1. surf_res : Surface Resolution (From 0.01 to 1 A, recommended 0.1 A)    ",
                       "  2. max_vert : Maximum Vertex Radius (From 0.10 to 20 A, recommended 7 A)  ",
                       "  3. box_size : Retaining Box Multiplier (From 1 to 10 A, recommended 1.5 A)",
                       "  4. sol_verts: Find the Vertices of all atoms or all atoms but the solute  ",
                       "                atoms (True/False, recommended True)                        ",
                       "                                                                             "]

    # Print everything
    print(splitting_line)
    print(help_header)
    print(splitting_line)
    print(instructions_header)
    print(splitting_line)
    for i in range(len(commands_header)):
        print(commands_header[i])
    print(splitting_line)
    for i in range(len(settings_header)):
        print(objects_header[i], "|", settings_header[i])
    print(splitting_line, "\n")


def get_ndx(ndx=None, list_len=None, obj=None):
    """
    Asks the user for the index of the object they specified
    :return:
    """
    if obj is not None and obj in my_objects:
        names = ["molecule", "residue", "atom", "index"]
        for i in range(4):
            if obj in [mol_objs, res_objs, atom_objs, ndx_objs][i]:
                obj = names[i]
    else:
        obj = "object"
    asking = True
    while asking:
        if ndx is None:
            extra = ""
            if list_len is not None:
                extra = " less than {}".format(list_len)
            prompt_str = "Enter {} index{} \nindex >>>   ".format(obj, extra)
            ndx = input(prompt_str)
        if ndx in quits or ndx in dones:
            return
        try:
            ndx = int(ndx)
            asking = False
        except ValueError:
            invalid_input(ndx)
            ndx = None
        if not asking and list_len is not None and ndx > list_len:
            invalid_input(ndx)
            asking = True
    return ndx


def get_obj(obj=None, return_ndx=True):
    """
    Makes the user type a proper object
    :return: 1-4 based on if it is a 1. molecule 2. residue 3. atom or 4. index
    """
    my_input, choosing = obj, False
    # If obj not in my objects
    if obj is None or obj.lower() not in my_objects:
        choosing = True

    # Keep asking the user to choose an object to export
    while choosing:
        # Prompt the user
        my_input = input("Enter an object type. (\'mol\', \'res\', \'atom\', or \'ndx\')\nobject type >>>   ")
        # Check to see if the user gave a valid response or not
        if my_input.lower() not in my_objects:
            # Tell the user they suck and try again
            invalid_input(my_input)
            continue
        # If they quit, then quit
        elif my_input.lower() in quits:
            return
        # Otherwise, we have a success
        else:
            choosing = False
    if return_ndx:
        # Go through and find the type of object we are getting
        objs = [mol_objs, res_objs, atom_objs, ndx_objs]
        for i in range(4):
            if my_input.lower() in objs[i]:
                return i + 1
    # As a failsafe
    return my_input


def get_file(file=None):
    if file is None:
        print("Enter a file address. (Use \'./\' to load a file from the \'.../vorpy\' directory):")
    checking_file = True
    while checking_file:
        # Get the file if None was specified
        if file is None:
            file = input("file address >>>   ")
            if file in quits:
                return
            test_file = file.split()
            if test_file[0] in load_cmds:
                file = file[len(test_file[0]) + 1:]
        if path.exists(file):
            checking_file = False
        else:
            invalid_input(file[0])
            file = None
            continue

    return file


def get_set(setting=None, val=None):
    """
        Makes the user type a proper Value
        :return: 1-4 based on if it is a 1. molecule 2. residue 3. atom or 4. index
        """
    my_input, choosing = val, False
    # If obj not in my objects
    if val is None or val.lower() not in my_settings:
        choosing = True

    # Keep asking the user to choose an object to export
    while choosing:
        my_input = setting
        if setting is None:
            # Prompt the user
            my_input = input("Enter setting type. (\'surf_res\', \'max_vert\', \'box_size\', or \'sol_verts\')\nsetting type >>>   ")
        # Check to see if the user gave a valid response or not
        if my_input.lower() not in my_settings:
            # Tell the user they suck and try again
            invalid_input(my_input)
            continue
        # If they quit, then quit
        elif my_input.lower() in quits:
            return
        # Otherwise, we have a success
        else:
            choosing = False
    # As a failsafe
    return my_input


def get_val(setting, val=None):
    """
        Asks the user for the index of the object they specified
        :return:
        """
    sett_name = ""
    if setting.lower() in my_settings:
        names = ["surface resolution", "maximum vertex", "box size", "solute vertices"]
        for i in range(4):
            if sett in [surf_reses, max_verts, box_sizes, sol_vertses][i]:
                sett_name = names[i]
    else:
        sett_name = "setting"
    asking = True
    while asking:
        if val is None:
            prompt_str = "Enter {} value \nvalue >>>   ".format(sett_name)
            val = input(prompt_str)
        if val in quits or val in dones:
            return
        if setting in sol_vertses:
            try:
                val = bool(val)
            except ValueError:
                val = None
        else:
            try:
                val = float(val)
            except ValueError:
                val = None
        if val is not None:
            asking = False
    return val



def print_list(names, list_name=None, width=150, height=30, cutoff=15):
    """
    Prints a long list in columns with numbers and allows the user to scroll through the list
    :param names:
    :param list_name:
    :param width:
    :param height:
    :param cutoff:
    :return:
    """
    # Check to see if a list name was provided
    if list_name is None:
        list_name = "My Objects"
    # First find the longest input in the list
    max_len = 0
    for name in names:
        if len(name) > max_len:
            max_len = len(name)
    if max_len > cutoff:
        max_len = cutoff

    # Figure out the columns. num cols = width / # of digits in index of last element + 2 ('. ') + max_len + 2 spaces
    num_cols = int(width / (2 + len(str(len(names) - 1)) + max_len + 2))
    # Print the first set of numbers
    i, row = 0, 0
    # Go through the names row by row, also print the header
    print(list_name)
    while row < height:
        row_str = ""
        for col in range(num_cols):
            if i >= len(names):
                row = np.inf
            else:
                row_str += str(i) + ". " + " " * (len(str(len(names) - 1)) - len(str(i))) + names[i] + " " * (
                            max_len - len(names[i])) + "  "
                i += 1
        print(row_str)
    # If that is all the data we are done and able to quit
    if len(names) < num_cols * height:
        return
    # In the case where the user wants to see a really long list, allow them to scroll
    scrolling = True
    while scrolling:
        my_response = input("Enter an index or a range or type 'q' to quit. (\'356\' or \'400-600\'")
        nums = None
        for i in range(len(my_response)):
            if my_response[i] in splitters:
                try:
                    nums = [int(my_response[:i], int(my_response[i + 1:]))]
                except ValueError:
                    nums = None
        # Check to see if a single number has been entered
        if nums is None:
            try:
                nums = [int(my_response)]
            except ValueError:
                nums = None
        # Print the lists
        if nums is not None and len(nums) == 1 and nums[0] < len(names):
            print_list(names[nums[0]:nums[0] + num_cols * height],
                       list_name=list_name + ": Elements " + str(nums[0]) + "-" + str(nums[0] + num_cols * height),
                       height=height, width=width, cutoff=max_len)
        elif nums is not None and nums[0] < nums[1] < len(names):
            # Check to see if the height needs to be changed
            new_height = height
            if nums[1] - nums[0] > num_cols * height:
                new_height = (nums[1] - nums[0]) // num_cols + 1
            print_list(names[nums[0]:nums[1]], list_name=list_name + ": Elements " + str(nums[0]) + "-" + str(nums[1]),
                       height=new_height, width=width, cutoff=max_len)
        else:
            invalid_input(my_response)



def show(sys, usr_npt):
    """
    Shows the input group type
    :return:
    """

    # If the user types 'Show' have a catch for it
    if len(usr_npt) == 1:
        show_var = get_obj(return_ndx=False).lower()
    # Get the list that the user wants to be shown if none was provided
    elif len(usr_npt) == 2 and usr_npt[1] in my_objects:
        show_var = usr_npt[1].lower()
    else:
        invalid_input(usr_npt)
        return

    # Get the correct list to show the user
    if show_var in mol_objs:
        show_name = "{} Molecules".format(sys.name)
        show_list = sys.mol_names
    elif show_var in res_objs:
        show_name = "{}".format(sys.name)
        show_list = sys.res_names
    elif show_var in atom_objs:
        show_name = "{}".format(sys.name)
        show_list = sys.atom_names
    elif show_var in ndx_objs:
        show_name = "{}".format(sys.name)
        show_list = sys.ndx_names
    else:
        show_name = ""
        show_list = []

    # Show the list
    if len(show_list) == 0:
        print("no objects to show")
        return
    else:
        print_list(show_list, show_name)



def load(sys, usr_npt):
    """
    Once one of the load commands is used try to load the rest of the string
    :param sys: System object to add the file to
    :param usr_npt:
    :return:
    """
    my_files = []
    if len(usr_npt) == 1:
        my_files.append(get_file())
    else:
        for file in usr_npt[1::2]:
            my_file = get_file(file)
            if my_file is None:
                return
            my_files.append(my_file)
    for file in my_files:
        # Check to see what type of file it is
        if file[-3:] == 'pdb' or file[-3:] == 'mol' or file[-3:] == 'gro' or file[-3:] == 'cif':
            if sys.name is not None and \
                    (sys.atoms is not None or sys.vert_file is not None or sys.net_file is not None):
                reset_sys = input("\n{} System is already loaded. Would you like to replace this system with: {} (y/n)"
                                  .format(sys.name, file))
                if reset_sys.lower() in ys:
                    my_sys = System(file)
                    print(my_sys.name, "loaded")
                    return my_sys
            else:
                my_sys = System(file)
                print(my_sys.name, "loaded")
                return my_sys
        # If the loaded file is a vertex or network file load them accordingly
        elif file[-3:] == 'txt':
            # If the new file is a vertex file load it
            if file[-9:-4] == 'verts':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.vert_file is not None and sys.vert_file != "":
                    replace_vert_file = input("\nThere is a vertex file already loaded: {}\n "
                                              "Would you like to replace it? (y/n):   ".format(sys.vert_file))
                    if replace_vert_file in ys or replace_vert_file in dones:
                        sys.load_verts(file)
                        print(sys.vert_file, "loaded")
                else:
                    sys.load_verts(file)
                    print(sys.vert_file, "loaded")

            # If the new file is a network file load it
            elif file[-11:-4] == 'network':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.net_file is not None or sys.net_file != "":
                    replace_net_file = input("\n{} is already loaded\n "
                                              "replace file >>>   ".format(sys.net_file))
                    if replace_net_file in ys:
                        sys.load_net(file)
                        print(sys.net_file, "loaded")
                    else:
                        return
                else:
                    sys.load_net(file)
                    print(sys.net_file, "loaded")
        # If the file is an index file load it accordingly
        elif file[-3:] == 'ndx':
            sys.load_ndx(file)
            print(sys.ndx_file, "loaded")

        # In all other case print an error and give the user a chance to try again
        else:
            print("\'{}\' is not a valid input. allowed file types: .pdb, .mol, .cif, .gro, .txt, .ndx. type "
                  "\'h\' for help".format(file))
            return


def sett(sys, usr_npt):
    """
    Set the network parameters
    :param sys:
    :param usr_npt:
    :return:
    """
    if len(usr_npt) == 1:
        my_set = get_set()
        my_val = get_val(my_set)
    elif len(usr_npt) == 2:
        my_set = get_set(usr_npt[1])
        my_val = get_val(my_set)
    elif len(usr_npt) <= 3:
        my_set = get_set(usr_npt[1])
        my_val = get_val(my_set, usr_npt[2])
    else:
        invalid_input(usr_npt)
        return
    # Check to see if a network has been created yet
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Set the surfaces resolution
    if my_set in surf_reses:
        sys.net.surf_res = my_val
        print(u"surface resolution set to {} \u212B".format(my_val))
    # Set the maximum vertex radius
    elif my_set in max_verts:
        sys.net.max_vert = my_val
        print(u"maximum vertex radius set to {} \u212B".format(my_val))
    # Set the box multiplier
    elif my_set in box_sizes:
        sys.net.box_size = my_val
        print("box size multiplier set to {} x".format(my_val))
    # Set the solute vertices
    elif my_set in sol_vertses:
        sys.net.sol_verts = my_val
        print("calculate solute vertices set to {}".format(my_val))
    else:
        invalid_input(usr_npt)


def build(sys):
    """
    Prints a pre-built header and asks the user if they are ready to build. Once confirmed prints a building header and
    builds
    :param sys:
    :return:
    """
    # Check to see if a network has been added
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Once the build command is used, the user is greeted with the build settings and asked if they are ready to build
    prepping_header = "\rsettings     -    surf_res = {:.2f}   max_vert  = {:.2f}   box_size = {:.2f}   sol_verts = {}"\
        .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts)
    print(prepping_header, end=" ")
    # The user is prompted to start the build - This could say eta and other build qualities
    pre_build_confirmation = input("\nbuild network (y/n) >>>   ")
    # If the user is ready to build, build the system
    if pre_build_confirmation.lower() in ys:
        sys.build_network()
    elif pre_build_confirmation.lower() in ns:
        print("note: to change build settings use the command line using the \'set\' command, a setting and a value ")
        return


def group(sys, usr_npt):
    """
    Takes input strings interprets them and returns a group
    :param sys:
    :param usr_npt:
    :return:
    """
    selections = []
    selection_names = []
    # Go through the user's input's grabbing the atoms, molecules, residues and indices
    for i in range(len(usr_npt)):
        # If the current value is and, continue
        if usr_npt[i] in ands:
            continue
        # If the user specified the index of the object us it
        if i + 1 < len(usr_npt):
            print(usr_npt, i + 1)
            my_ndx = int(usr_npt[i + 1])
        # If the user did not we need to ask for it
        else:
            my_ndx = get_ndx()
        # Exit this iteration if no ndx can be found
        if my_ndx is None:
            continue
        # Check to see if the object specified was a molecule
        elif usr_npt[i] in mol_objs and my_ndx < len(sys.mols):
            # Add the selections from the molecules and molecule names list of the sytem
            selections += sys.mols[my_ndx]
            selection_names.append(sys.mol_names[my_ndx])
        elif usr_npt[i] in res_objs:
            # Add the selections from the residues and residue names list in the system
            selections += sys.residues[my_ndx]
            selection_names.append(sys.res_names[my_ndx])
        elif usr_npt[i] in atom_objs:
            # Add the selections from the residues and residue names list in the system
            selections.append(sys.atoms[my_ndx])
            selection_names.append(sys.atom_names[my_ndx])
        elif usr_npt[i] in ndx_objs:
            # Add the selections from the residues and residue names list in the system
            selections += sys.ndxs[my_ndx]
            selection_names.append(sys.ndx_names[my_ndx])
    # Create a group and add it to the system
    my_group = Group(sys.net, atoms=selections, name="_".join(selection_names))
    sys.groups.append(my_group)
    naming = True
    while naming:
        # Ask the user if they want to rename the group
        rename = input("\ngroup name = {}\nchange name >>>   ".format(my_group.name))
        # If they want to rename the group, let them
        if rename in ys:
            new_name = input("What would you like to rename it to? (up to 10 characters, no spaces):   ")
            new_name.replace(" ", "_")
            my_group.name = new_name
        # If they don't we are done
        elif rename in ns:
            naming = False
        else:
            change_to_input = input("new group name {} \nchange name >>>  ".format(rename))
            if change_to_input in ys:
                group.name = rename



def export(sys, usr_npt):
    """
    Takes in input strings and exports them based on their option choices
    :param sys:
    :param usr_npt:
    :return:
    """
    my_obj, my_ndx, my_obj2, my_ndx2 = None, None, None, None
    if len(usr_npt) <= 1:
        my_obj = get_obj()
        my_ndx  = get_ndx()
    elif len(usr_npt) == 2:
        my_obj = get_obj(obj=usr_npt[1])
        my_ndx = get_ndx(obj=my_obj)
    elif len(usr_npt) == 3:
        my_obj = get_obj(obj=usr_npt[1])
        my_ndx = get_ndx(ndx=usr_npt[2], obj=my_obj)
    elif len(usr_npt) > 3 and usr_npt[3].lower() in ands:
        my_obj = get_obj(obj=usr_npt[1])
        my_ndx = get_ndx(ndx=usr_npt[2], obj=my_obj)
        if len(usr_npt) == 4:
            my_obj2 = get_obj()
            my_ndx2 = get_ndx()
        elif len(usr_npt) == 5:
            my_obj2 = get_obj(usr_npt[4])
            my_ndx2 = get_ndx(obj=my_obj2)
        else:
            my_obj2 = get_obj(usr_npt[4])
            my_ndx2 = get_ndx(usr_npt[5], obj=my_obj2)
    else:
        invalid_input(usr_npt)
        return

    # Get the atoms
    try:
        my_atoms = [sys.mols, sys.residues, sys.atoms, sys.ndxs][my_obj - 1][my_ndx]
    except IndexError:
        invalid_input(my_obj)
        return
    # In the case where we get one atom, make it a list
    if type(my_atoms) is not list:
        my_atoms = [my_atoms]

    # Create the first group, get its data and export it
    group1 = Group(atoms=my_atoms, net=sys.net)
    group1.get_info()
    export_body(group1, True, True)


    # Check to see if there is a second group
    if my_obj2 is not None:
        my_atoms = [sys.mols, sys.residues, sys.atoms, sys.ndxs][my_obj2 - 1][my_ndx2]
        if type(my_atoms) is not list:
            my_atoms = [my_atoms]
        group2 = Group(atoms=my_atoms, net=sys.net)
        group2.get_info()
        export_body(group2, True, True)
        # Ask the user if they want to export the interface
        exp_iface = input("Export interface between {} and {}? (y/n):   \nvorpy >>>")
        if exp_iface.lower() in ys:
            group1.bff, group2.bff = group2, group1
            group1.get_info()
            group2.get_info()
            export_iface([group1, group2], True, True)
            print("{}, {} and {}-{} interface exported".format(group1.name, group2.name, group1.name, group2.name))
        else:
            print("{} and {} exported".format(group1.name, group2.name))
    else:
        print("{} exported".format(group1.name))



def check_input(sys):
    """
    Main function that is looped. Checks the inputs and runs the correct functions
    :param sys:
    :return:
    """
    # Set up the prompt
    usr_npt = input("vorpy >>>   ")
    # Split it up by the spaces
    usr_npt = usr_npt.split()
    # Check to see if the initial input is a command
    if usr_npt[0] not in my_commands:
        invalid_input(usr_npt)
        return True

    ########################## Commands  ################################################

    # Check if the user's input is in loads
    if usr_npt[0].lower() in load_cmds:
        my_sys = load(sys=sys, usr_npt=usr_npt)
        if my_sys is not None:
            sys = my_sys
            return sys
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in set_cmds:
        sett(sys, usr_npt)
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in build_cmds:
        build(sys)
    # Check if the user's input is in groups
    elif usr_npt[0].lower() in group_cmds:
        group(sys, usr_npt)
    # Check if the user's input is in exports
    elif usr_npt[0].lower() in export_cmds:
        export(sys, usr_npt)
    # Check if the user's input is in shows
    elif usr_npt[0].lower() in show_cmds:
        show(sys, usr_npt)
    # Check if the user wants help
    elif usr_npt[0].lower() in helps:
        get_help()
    # Check to see if the user's input includes quits
    elif usr_npt[0].lower() in quits:
        if are_you_sure():
            return False

    # Unless the user quits, keep running the program
    return True


if __name__ == '__main__':
    # Welcome introduction
    print("Welcome to vorpy. For assistance type \'h\'")
    # Create the system
    mySys = System()
    # Set up the running variable
    running = True
    # Run the program
    while running:
        # create_header(mySys)
        running = check_input(mySys)
        if type(running) is not bool:
            mySys = running
            running = True
