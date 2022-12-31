import os
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
ndx_objs = ['i', 'is',  'index', 'indexs', 'indexes', 'indices', 'ndx', 'ndxs', 'ndex', 'group', 'g', 'grp']

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
    ays = input("confirm >>>   ")
    if ays.lower() in ys:
        return True
    elif ays.lower() in helps:
        help_()
    elif ays.lower() in quits:
        return
    return False


def invalid_input(string):
    if type(string) is list:
        string = " ".join(string)
    print("\'{}\' is not a valid input. try again or type \'h\' for help".format(string))


def help_():
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
    obj_name = None
    names = ["molecule", "residue", "atom", "index"]
    # Get the object and make sure it's real
    if obj is not None and obj in my_objects:
        for i in range(4):
            if obj in [mol_objs, res_objs, atom_objs, ndx_objs][i]:
                obj =  i + 1
                obj_name = names[i]
    elif type(obj) is int and obj <= 4:
        obj_name = names[obj - 1]
    else:
        obj = get_obj(obj=obj)
        obj_name = names[obj - 1]
    # Start asking for the index of the object
    asking = True
    while asking:
        max_ndx = len([sys.mols, sys.residues, sys.atoms, sys.groups][obj - 1])
        if max_ndx == 0:
            print("No {} to choose from. Choose another object or type \'h\' for help")
        elif ndx is None:
            extra = ""
            if list_len is not None:
                extra = " less than {}".format(list_len)
            prompt_str = "Enter {} index{} between 0 and {}\nindex >>>   ".format(obj_name, extra, max_ndx)
            ndx = input(prompt_str)
        if ndx.lower() in quits or ndx.lower() in dones:
            return
        elif ndx.lower() in helps:
            help_()
            continue
        elif ndx.lower() in show_cmds:
            show(obj_name)
            continue
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
    if obj is None or (type(obj) is str and obj.lower() not in my_objects) or (type(obj) is int and obj > 4):
        choosing = True
    # Keep asking the user to choose an object to export
    while choosing:
        # Prompt the user
        my_input = input("Enter an object type. (\'mol\', \'res\', \'atom\', or \'ndx\')\nobject >>>   ")
        # Check to see if the user gave a valid response or not
        if my_input.lower() in quits:
            return
        elif my_input.lower() in helps:
            help_()
        elif my_input.lower() not in my_objects:
            # Tell the user they suck and try again
            invalid_input(my_input)
            continue
        # Otherwise, we have a success
        else:
            choosing = False
    if return_ndx:
        # If the input is already an integer return it
        if type(my_input) is int and my_input <= 4:
            return my_input
        # Go through and find the type of object we are getting
        objs = [mol_objs, res_objs, atom_objs, ndx_objs]
        for i in range(4):
            if my_input.lower() in objs[i]:
                if len([sys.mols, sys.residues, sys.atoms, sys.ndxs][i]) > 0:
                    return i + 1
                else:
                    print("No {} in the system. Try again or typ \'h\' for help"
                          .format(["molecules", "residues", "atoms", "groups"][i]))

    # As a failsafe
    return my_input


def get_file(file=None):
    # Check if there is a file provided
    if file is None:
        print("Enter a file address. (Use \'./\' to load a file from the \'.../vorpy\' directory):")
    # Check the file
    checking_file = True
    while checking_file:
        # Get the file if None was specified
        if file is None:
            file = input("file >>>   ")
            if file.lower() in quits:
                return
            elif len(file) == 0:
                checking_file = True
                continue
            elif file.lower() in helps:
                help_()
            test_file = file.split()
            if test_file[0] in load_cmds:
                file = file[len(test_file[0]) + 1:]
        # Check if the initial file works
        if path.exists(file) and len(file) > 0:
            checking_file = False
        # Check if the file is in the ./Data/test_data folder
        elif path.exists("./Data/test_data/" + file) and len(file) > 0:
            file = os.getcwd() + "/Data/test_data/" + file
            checking_file = False
        # Check if it is just the raw name
        elif path.exists("./Data/test_data/" + file + ".pdb") and len(file) > 0:
            file = os.getcwd() + "/Data/test_data/" + file + ".pdb"
            checking_file = False
        # Otherwise, tell the user to try again
        else:
            invalid_input(file)
            file = None
            continue
    # Return the file
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
            my_input = input("Enter setting type. (\'surf_res\', \'max_vert\', \'box_size\', or \'sol_verts\')\nsetting >>>   ")
        # Check to see if the user gave a valid response or not
        if my_input.lower() not in my_settings:
            # Tell the user they suck and try again
            invalid_input(my_input)
            continue
        # If they quit, then quit
        elif my_input.lower() in quits:
            return
        elif my_input.lower() in helps:
            help_()
        # Otherwise, we have a success
        else:
            choosing = False
    # As a failsafe
    return my_input


def get_val(setting=None, val=None):
    """
        Asks the user for the index of the object they specified
        :return:
        """
    # If the setting is not in the settingsd
    if setting is None or setting.lower() not in my_settings:
        setting = get_set()
    # Set up the list of possible names
    names = ["surface resolution", "maximum vertex", "box size", "solute vertices"]
    # Find the value for the setting
    asking = True
    while asking:
        # If no val has been provided
        if val is None:
            prompt_str = "Enter \'{}\' value \nvalue >>>   ".format(setting)
            val = input(prompt_str)

        if val.lower() in quits or val.lower() in dones:
            return
        elif val.lower() in helps:
            help_()
        if setting in sol_vertses:
            if val.lower() in ['t', 'true', 'tr'] + ys:
                val = True
            elif val.lower() in ['f', 'false', 'flse', 'fl', 'fa', 'fs', 'fls'] + ns:
                val = False
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
        my_response = input("Enter an index or a range or type 'q' to quit. (\'356\' or \'400-600\')\nindex >>>   ")
        if my_response.lower() in quits:
            return
        elif my_response.lower() in helps:
            help_()
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



def show(usr_npt=None):
    """
    Shows the input group type
    :return:
    """
    global sys
    # If the user types 'Show' have a catch for it
    if usr_npt is None or len(usr_npt) == 1:
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



def load(usr_npt):
    """
    Once one of the load commands is used try to load the rest of the string
    :param usr_npt:
    :return:
    """
    global sys
    my_files = []
    if len(usr_npt) == 1:
        my_files.append(get_file())
        if my_files[-1] is None or my_files[-1].lower() in quits:
            return

    else:
        for file in usr_npt[1::2]:
            my_file = get_file(file)
            if my_file is None or my_file.lower() in quits:
                return
            my_files.append(my_file)
    for file in my_files:
        # Check to see what type of file it is
        if file[-3:] == 'pdb' or file[-3:] == 'mol' or file[-3:] == 'gro' or file[-3:] == 'cif':
            if sys.name is not None and \
                    (sys.atoms is not None or sys.vert_file is not None or sys.net_file is not None):
                reset_sys = input("Replacing {} with {}\nconfirm >>>   "
                                  .format(sys.name, file))
                if reset_sys.lower() in ys:
                    sys = System(file)
                    print(sys.name + " loaded - {} atoms, {} molecules, solute: {}"
                          .format(len(sys.atoms), len(sys.mols), sys.sol_name))
                    return sys
                elif reset_sys.lower() in helps:
                    help_()
                elif reset_sys.lower() in quits:
                    return
            else:
                sys = System(file)
                print(sys.name + " loaded - {} atoms, {} molecules, solute: {}"
                      .format(len(sys.atoms), len(sys.mols), sys.sol_name))
                return sys
        # If the loaded file is a vertex or network file load them accordingly
        elif file[-3:] == 'txt':
            # If the new file is a vertex file load it
            if file[-9:-4] == 'verts':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.vert_file is not None and sys.vert_file != "":
                    replace_vert_file = input("Replacing {} with {}\n "
                                              "confirm >>>   ".format(sys.vert_file, file))
                    if replace_vert_file.lower() in ys or replace_vert_file.lower() in dones:
                        sys.load_verts(file)
                        print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
                    elif replace_vert_file.lower() in helps:
                        help_()
                    elif replace_vert_file.lower() in quits:
                        return
                else:
                    sys.load_verts(file)
                    print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                          .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))

            # If the new file is a network file load it
            elif file[-11:-4] == 'network':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.net_file is not None or sys.net_file != "":
                    replace_net_file = input("Replacing {} with {}\n "
                                              "confirm >>>   ".format(sys.net_file, file))
                    if replace_net_file in ys:
                        sys.load_net(file)
                        print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box"
                              " size: {} x\n".format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
                    elif replace_net_file in helps:
                        help_()
                    else:
                        return
                else:
                    # Load the file
                    sys.load_net(file)
                    if len(sys.net.surfs) > 0:
                        print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
                    else:
                        print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
        # Check to see if it is a new network file
        elif file[-3:] == 'csv':
            # Check to see that this is a network file
            if file[-7:-4].lower() == 'net':
                sys.load_net(file=file)
        # If the file is an index file load it accordingly
        elif file[-3:] == 'ndx':
            sys.load_ndx(file)
            print(sys.ndx_file + "loaded -  {}".format(sys.ndx_names[:min(len(sys.ndx_names) - 1, 10)]))
        # In all other case print an error and give the user a chance to try again
        else:
            print("\'{}\' is not a valid input. allowed file types: .pdb, .mol, .cif, .gro, .txt, .ndx. type "
                  "\'h\' for help".format(file))
            return


def sett(usr_npt):
    """
    Set the network parameters
    :param usr_npt:
    :return:
    """
    global sys
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
        print(u"Surface resolution set to {} \u212B".format(my_val))
    # Set the maximum vertex radius
    elif my_set in max_verts:
        sys.net.max_vert = my_val
        print(u"Maximum vertex radius set to {} \u212B".format(my_val))
    # Set the box multiplier
    elif my_set in box_sizes:
        sys.net.box_size = my_val
        print("Box size multiplier set to {} x".format(my_val))
    # Set the solute vertices
    elif my_set in sol_vertses:
        sys.net.sol_verts = my_val
        print("Calculate solute vertices set to {}".format(sys.net.sol_verts))
    else:
        invalid_input(usr_npt)


def build():
    """
    Prints a pre-built header and asks the user if they are ready to build. Once confirmed prints a building header and
    builds
    :return:
    """
    global sys
    # If no system has been loaded tell the user to fuck off
    if len(sys.atoms) == 0:
        print("No atoms in the system. Use the \'load\' command or type \'h\' for help")
        return
    # Check to see if a network has been added
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Once the build command is used, the user is greeted with the build settings and asked if they are ready to build
    print("\rSettings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_size = {:.2f} x,  sol_verts = {}"
        .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts), end="")
    # The user is prompted to start the build - This could say eta and other build qualities
    pre_build_confirmation = input("\nconfirm >>>   ")
    # If the user is ready to build, build the system
    if pre_build_confirmation.lower() in ys:
        sys.build_network()
    elif pre_build_confirmation.lower() in ns:
        print("Use the \'set\' command to change a setting and a value. Type \'h\' for help")
        return
    elif pre_build_confirmation.lower() in helps:
        help_()
    elif pre_build_confirmation.lower() in quits:
        return


def group(usr_npt, for_export=False):
    """
    Takes input strings interprets them and returns a group
    :param usr_npt:
    :return:
    """
    global sys
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
        my_obj = get_obj(obj=selection[0], return_ndx=True)
        my_ndx = get_ndx(obj=my_obj, ndx=selection[1])
        # Check to see if the index provided is out of range
        checking_ndx, my_atoms = True, None
        while checking_ndx:
            # Try to create the selection
            try:
                my_atoms = [sys.mols, sys.residues, sys.atoms, sys.ndxs][my_obj - 1][my_ndx]
                checking_ndx = False
            except IndexError:
                my_ndx = get_ndx(obj=my_obj)
                continue
        # Check real quick for single atoms
        if type(my_atoms) is not list:
            my_atoms = [my_atoms]
        # Create the group
        my_group = Group(sys.net, atoms=my_atoms, name="_".join(selection_names[i]))
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


def export(usr_npt):
    """
    Takes in input strings and exports them based on their option choices
    :param usr_npt:
    :return:
    """
    if len(usr_npt) > 1 and usr_npt.lower() == 'all':
        sys.exports(network=True, pdb=True, surfaces=True, full_network_object=True, no_sol_network_object=True,
                    alter_atoms_script=True)
        return
    # Get the groups based off of what was specified
    groups = group(usr_npt, for_export=True)
    # Go through the groups in the list
    for my_group in groups:
        export_body(my_group, True, True)
    # Check to see if there is a second group
    if len(groups) > 1:
        # Ask the user if they want to export the interface
        exp_iface = input("Exporting interface between {} and {}\nconfirm >>>   ".format(groups[0].name, groups[1].name))
        if exp_iface.lower() in ys:
            groups[0].bff, groups[1].bff = groups[1], groups[0]
            groups[0].get_info()
            groups[1].get_info()
            export_iface(groups, True, True)
            print("Groups {} and {}, and {}-{} interface exported".format(groups[0].name, groups[1].name, groups[0].name, groups[1].name))
        else:
            print("Groups {} and {} exported".format(groups[0].name, groups[1].name))
    else:
        print("Group {} exported".format(groups[0].name))



def check_input():
    """
    Main function that is looped. Checks the inputs and runs the correct functions
    :return:
    """
    global sys
    # Set up the prompt
    usr_npt = input("vorpy >>>   ")
    # Split it up by the spaces
    usr_npt = usr_npt.split()
    # Check to see if the initial input is a command
    if len(usr_npt) == 0 or usr_npt[0] not in my_commands:
        invalid_input(usr_npt)
        return True

    ########################## Commands  ################################################

    # Check if the user's input is in loads
    if usr_npt[0].lower() in load_cmds:
        my_sys = load(usr_npt=usr_npt)
        if my_sys is not None:
            sys = my_sys
            return sys
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in set_cmds:
        sett(usr_npt)
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in build_cmds:
        build()
    # Check if the user's input is in groups
    elif usr_npt[0].lower() in group_cmds:
        group(usr_npt)
    # Check if the user's input is in exports
    elif usr_npt[0].lower() in export_cmds:
        export(usr_npt)
    # Check if the user's input is in shows
    elif usr_npt[0].lower() in show_cmds:
        show(usr_npt)
    # Check if the user wants help
    elif usr_npt[0].lower() in helps:
        help_()
    # Check to see if the user's input includes quits
    elif usr_npt[0].lower() in quits:
        if are_you_sure():
            return False

    # Unless the user quits, keep running the program
    return True


if __name__ == '__main__':
    # Welcome introduction
    print("Welcome to vorpy. For assistance type \'h\'. To quit type \'q\'")
    # Create the system
    sys = System()
    # Set up the running variable
    running = True
    # Run the program
    while running:
        # create_header(mySys)
        running = check_input()
        if type(running) is not bool:
            sys = running
            running = True
