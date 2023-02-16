# Responses
ys = ['y', 'yes', 'ya', 'yeet', 'yur', 'yoint', 'uhu', 'yup', 'jess', 'affirmative', 'yass', '', 'yuss', 'yess',
      'yesss', 'yessss', 'yar', 'yuh', 'mhm']
ns = ['n', 'no', 'naur', 'nope', 'nonya', 'nope', 'nien', 'nada']
dones = ['done', 'd', 'finished', 'finito', 'complete', 'doneso', 'don', 'fin', 'keep movin bruh',
         'you still here?', 'go on', 'goo\'oon']
ands = ['&', 'and', 'nd', 'also', '+', '&&']
splitters = ['/', '-']

# General commands
quits = ['quit', 'q', 'qt', 'exit', 'ext']
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
ndx_objs = ['i', 'is',  'index', 'indexs', 'indexes', 'indices', 'ndx', 'ndxs', 'ndex', 'group', 'g', 'grp', 'n']

my_objects = mol_objs + res_objs + atom_objs + ndx_objs

# Settings
surf_reses = ['surf_res', 'sr', 'surface_resolution', 'surface_res', 'surf_resolution', 'surfs', 'surf', 'surfs_res', 'surfs_resolution', 'surfaces_resolution', 'surfaces_res']
max_verts = ['max_vert', 'mv', 'maximum_vertex', 'max_vertex', 'maximum_vert', 'verts', 'vs', 'vert_size', 'max_vert_size', 'mvs', 'vert_max', 'vertex_max', 'vertex_maximum']
box_sizes = ['box_size', 'bm', 'box', 'bx_sz', 'size_box', 'containing_box', 'containing_box_size', 'box_multi', 'box_multiplier']
build_surfses = ['build_surfs', 'build_surfaces', 'bs', 'bld_srfs', 'cs', 'calc_surfs', 'surfs_build', 'surfaces_build', 'build_surf', 'build_surf']
flat_surfses = ['flat_surfs', 'flat_surfaces', 'fs', 'flt_srfs', 'surfaces_flat', 'surfs_flat', 'flat_surf', 'flat_surface', 'ff']

my_settings = surf_reses + max_verts + box_sizes + build_surfses + flat_surfses


def are_you_sure():
    ays = input("confirm >>>   ")
    if ays.lower() in ys:
        return True
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



    commands_header = ["Commands:                                                                                                                      ",
                       "  1. load  : Loads file types by their extension - System: (.pdb/.gro/.mol/.cif), Network: (.csv), Surface (.off), Index (.ndx)",
                       "  2. set   : Sets network build settings with a setting - (see Settings) and a value - (float, float, float, T/F, T/F)         ",
                       "  3. build : Builds the network. Asks the user to confirm and shows the current settings before starting the build process.    ",
                       "  4. export: Exports network objects with a name - (see below) and (optionally) an index (integer or range separated with \'-\') ",
                       "  5. show  : Shows element information by name (see below) for reference in a command (load/set/build/export)                  ",
                       "                                             (for more type \'c\')                                                             "]

    splitting_line = "--------------------------------------------------------------------------------------------------------------------------------"

    objects_header = ["Objects:                                           ",
                      "  1. mol : Molecule object from the current System ",
                      "  2. res : Residue object from the current System  ",
                      "  3. atom: Atom object from the current System     ",
                      "  4. ndx : Index loaded into the current System or ",
                      "           created by the user                     ",
                      "                                                   "]

    settings_header = ["Settings:                                                                   ",
                       "  1. surf_res : Surface Resolution (From 0.01 to 1 A, recommended 0.1 A)    ",
                       "  2. max_vert : Maximum Vertex Radius (From 0.10 to 20 A, recommended 7 A)  ",
                       "  3. box_size : Retaining Box Multiplier (From 1 to 10 A, recommended 1.5 A)",
                       "  4. build_surfs: Calculate the network's surfaces (True/False)",
                       "  5. flat_surfs: Build the surfaces flat (True/False)     "]

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
                row = 100000000000
            else:
                row_str += "(" + str(i) + ") - " + " " * (len(str(len(names) - 1)) - len(str(i))) + names[i] + " " * (
                            max_len - len(names[i])) + ",  "
                i += 1
        print(row_str)
    # If that is all the data we are done and able to quit
    if len(names) < num_cols * height:
        return
    # In the case where the user wants to see a really long list, allow them to scroll
    scrolling = True
    while scrolling:
        my_response = input("enter an index or a range or type 'q' to quit. (\'356\' or \'400-600\')\nindex >>>   ")
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
                       list_name=list_name + ": elements " + str(nums[0]) + "-" + str(nums[0] + num_cols * height),
                       height=height, width=width, cutoff=max_len)
        elif nums is not None and nums[0] < nums[1] < len(names):
            # Check to see if the height needs to be changed
            new_height = height
            if nums[1] - nums[0] > num_cols * height:
                new_height = (nums[1] - nums[0]) // num_cols + 1
            print_list(names[nums[0]:nums[1]], list_name=list_name + ": elements " + str(nums[0]) + "-" + str(nums[1]),
                       height=new_height, width=width, cutoff=max_len)
        else:
            invalid_input(my_response)



def get_obj(sys, obj=None, return_ndx=True):
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
        my_input = input("enter an object type. (\'mol\', \'res\', \'atom\', or \'ndx\')\nobject >>>   ")
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
                    print("no {} in the system. try again or typ \'h\' for help"
                          .format(["molecules", "residues", "atoms", "groups"][i]))

    # As a failsafe
    return my_input


def show(sys, usr_npt=None):
    """
    Shows the input group type
    :return:
    """
    # If the user types 'Show' have a catch for it
    if usr_npt is None or len(usr_npt) == 1:
        show_var = get_obj(sys=sys, return_ndx=False).lower()
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
