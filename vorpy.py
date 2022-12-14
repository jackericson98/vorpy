import numpy as np

from System.system import *

"""
I want this to be a constantly running interface where the user can load, build, and export Voronoi network components. 
Has a running header that updates as commands are preformed. User can get system information based off of requests
"""



# Responses
ys = ['y', 'yes', 'ya', 'yeet', 'yur', 'yoint', 'uhu', 'yup', 'jess', 'affirmative', 'yass', '', 'yuss', 'yess',
      'yesss', 'yessss', 'yar', 'yuh', 'mhm']
ns = ['n', 'no', 'naur', 'nope', 'nonya', 'nope', 'nien', 'nada']
dones = ['done', 'd', 'finished', 'finito', 'complete', 'doneso', 'don', 'fin', 'keep movin bruh', 'mf move on',
         'you still here?', 'go on', 'goo\'oon']
ands = ['&', 'and', 'nd', 'also', '+', '&&']

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
    ays = input("Are you sure? (Y/N):   ")
    if ays in ys:
        return True
    return False


def invalid_input(string):
    if type(string) is list:
        string = " ".join(string)
    print("\'{}\' is not a valid input. Please try again or type \'h\' for help".format(string))


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


def get_ndx(ndx_item=""):
    """
    Asks the user for the index of the object they specified
    :return:
    """
    asking = True
    my_ndx = None
    while asking:
        my_ndx = input("Please enter a(n) {} index".format(ndx_item))
        if my_ndx in quits or my_ndx in dones:
            return
        try:
            my_ndx = int(my_ndx) - 1
            asking = False
        except ValueError:
            invalid_input(my_ndx)
    return my_ndx


def get_obj(obj=None):
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
        my_input = input("Enter an object type. (\'mol\', \'res\', \'atom\', or \'ndx\'):   ")
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
    # Go through and find the type of object we are getting
    objs = [mol_objs, res_objs, atom_objs, ndx_objs]
    for i in range(4):
        if my_input.lower() in objs[i]:
            return i + 1
    # As a failsafe
    return my_input


def show(sys, usr_npt):
    """
    Shows the input group type
    :return:
    """
    # Get the list that the user wants to be shown if none was provided
    if len(usr_npt) >= 1 and usr_npt[1] in my_objects:
        show_var = usr_npt[1]
    else:
        finding_show_var = True
        show_var = None
        while finding_show_var:
            show_var = input("Which object list would you like to see? (enter \'h\' for assistance):   ")
            if show_var.lower() in helps:
                get_help()
                return
            elif show_var.lower() not in my_objects:
                invalid_input(show_var)
            else:
                finding_show_var = False

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
    if len(show_list) <= 30:
        print(show_name)
        for i in range(len(show_list)):
            print(str(i + 1) + ". " + show_list[i])
    else:
        print(show_name)
        for i in range(len(show_list[:30])):
            print(str(i + 1) + ". " + show_list[i])
        print("The list is too long")
        showing = True
        while showing:
            usr_rng = input("Specify a range of numbers to show separated by a dash and spaces between 1 and {}".format(len(show_list)))
            # If the user inputs a help command
            if usr_rng.lower() in helps:
                get_help()
                continue
            # If the user quits
            elif usr_rng.lower() in quits:
                if are_you_sure():
                    return
            # Split their inputs by spaces and try to create integers out of the inputs
            usr_rng.split()
            ndx1, ndx2 = None, None
            try:
                ndx1 = int(usr_rng[0])
                ndx2 = int(usr_rng[2])
            except ValueError or IndexError:
                invalid_input(usr_rng)
            if ndx1 is not None and ndx2 is not None and 0 < ndx1 < ndx2:
                if ndx2 >= len(show_list):
                    ndx2 = len(show_list)
                for i in range(show_list[ndx1 - 1:ndx2 - 1]):
                    print(str(i + ndx1) + ". " + show_list[i])



def load(sys, usr_npt):
    """
    Once one of the load commands is used try to load the rest of the string
    :param sys: System object to add the file to
    :param usr_npt:
    :return:
    """
    for file in usr_npt[1::2]:
        # Check to see what type of file it is
        if file[-3:] == 'pdb' or file[-3:] == 'mol' or file[-3:] == 'gro' or file[-3:] == 'cif':
            if sys.name is not None and (sys.atoms is not None or sys.vert_file is not None or sys.net_file is not None):
                reset_sys = input("\n{} System is already loaded. Would you like to replace this system with: {} (Y/N)".format(sys.name, file))
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
                verts_loaded = True
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.vert_file is not None and sys.vert_file != "":
                    replace_vert_file = input("\nThere is a vertex file already loaded: {}\n "
                                              "Would you like to replace it? (Y/N):   ".format(sys.vert_file))
                    if replace_vert_file in ys or replace_vert_file in dones:
                        sys.load_verts(file)
                        print(sys.vert_file, "loaded")
                else:
                    sys.load_verts(file)
                    print(sys.vert_file, "loaded")

            # If the new file is a network file load it
            elif file[-11:-4] == 'network':
                net_loaded = True
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.net_file is not None or sys.net_file != "":
                    replace_net_file = input("\nThere is a Network file already loaded: {}\n "
                                              "Would you like to replace it? (Y/N):   ".format(sys.net_file))
                    if replace_net_file in ys or replace_net_file in dones:
                        sys.load_net(file)
                        print(sys.net_file, "loaded")
                else:
                    sys.load_net(file)
                    print(sys.net_file, "loaded")
        # If the file is an index file load it accordingly
        elif file[-3:] == 'ndx':
            sys.load_ndx(file)
            print(sys.ndx_file, "loaded")

        # In all other case print an error and give the user a chance to try again
        else:
            print("\r\n\n{} Is not a valid input. Please provide the full address for one of the following file types: \n .pdb, "
                  ".mol, .cif, .gro, _verts.txt (vorpy created), _network.txt (vorpy created), .ndx (GROMACS created)"
                  .format(file), end="")
            return


def sett(sys, usr_npt):
    """
    Set the network parameters
    :param sys:
    :param usr_npt:
    :return:
    """
    # Check to see if a network has been created yet
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Set the surfaces resolution
    if usr_npt[1] in surf_reses:
        # Check to see that the user entered a value:
        try:
            sys.net.surf_res = float(usr_npt[2])
        except ValueError or IndexError:
            pass
    # Set the maximum vertex radius
    elif usr_npt[1] in max_verts:
        # Check to see that the user entered a value:
        try:
            sys.net.max_vert = float(usr_npt[2])
        except ValueError or IndexError:
            pass
    # Set the box multiplier
    elif usr_npt[1] in box_sizes:
        # Check to see that the user entered a value:
        try:
            sys.net.box_size = float(usr_npt[2])
        except ValueError or IndexError:
            pass
    # Set the solute vertices
    elif usr_npt[1] in sol_vertses:
        # Check to see that the user entered a value:
        try:
            sys.net.sol_verts = bool(usr_npt[2])
        except ValueError or IndexError:
            pass


def build(sys, usr_npt=None):
    """
    Prints a pre-built header and asks the user if they are ready to build. Once confirmed prints a building header and
    builds
    :param sys:
    :param usr_npt:
    :return:
    """
    # Check to see if a network has been added
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Once the build command is used, the user is greeted with the build settings and asked if they are ready to build
    prepping_header = "\rSettings     -     surf_res = {:.2f}   max_vert  = {:.2f}\n" \
                        "                   box_size = {:.2f}   sol_verts = {}"\
        .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts)
    print(prepping_header, end=" ")
    # The user is prompted to start the build - This could say eta and other build qualities
    pre_build_confirmation = input("\nBuild network? (Y/N):   ")
    # If the user is ready to build, build the system
    if pre_build_confirmation.lower() in ys:
        sys.build_network()
    elif pre_build_confirmation.lower() in ns:
        print("Note To change build settings use the command line using the \'set\' command, a setting and a value ")
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
    # Ask the user if they want to rename the group
    rename = input("\n{}\n\nAre you ok with this name? (Y/N):   ".format(my_group.name))
    if rename in ns:
        new_name = input("What would you like to rename it to? (up to 10 characters, no spaces):   ")
        new_name.replace(" ", "_")
        my_group.name = new_name


def export(sys, usr_npt):
    """
    Takes in input strings and exports them based on their option choices
    :param sys:
    :param usr_npt:
    :return:
    """

    # Separate the ands
    my_selects = [[]]
    for npt in usr_npt[1:]:
        if npt in ands:
            my_selects.append([])
        else:
            my_selects[-1].append(npt)
    # Create the export groups variable
    my_groups = []
    # Go through the selections
    for select in my_selects:
        # Give the selection a place-holder
        if len(select) == 1:
            usr_npt.append("")
        # If the user gave an invalid selection make them choose one
        my_obj = get_obj(select[0])
        if my_obj is None:
            return
        # Get the index
        try:
            my_ndx = int(select[1])
        except ValueError:
            my_ndx = None
        if my_ndx is None:
            my_ndx = get_ndx()
        # Get the atoms
        my_atoms = [sys.mols, sys.residues, sys.atoms, sys.ndxs][my_obj - 1][my_ndx - 1]
        if type(my_atoms) is not list:
            my_atoms = [my_atoms]
        # Create the group
        my_group = Group(atoms=my_atoms, net=sys.net)
        my_group.get_info()
        # Add it to the groups list
        my_groups.append(my_group)

    # Export the groups
    for my_group in my_groups:
        export_body(my_group, info_file=True, outer_atoms=True)
    # Export the interfaces
    if len(my_groups) == 2:
        my_iface_groups = my_groups
    # If there are more than two groups specified ask the user which of the two they want
    elif len(my_groups) > 2:
        # Prompt the user
        g1_ndx = get_ndx("Group 1")
        g2_ndx = get_ndx("Group 2")
        my_iface_groups = [sys.groups[g1_ndx], sys.groups[g2_ndx]]
    else:
        return
    # Prompt the user to export the 2 groups as an interface or not
    export_iface_prompt = input("\nExport interface between {} and {}? (Y/N):   "
                                .format(my_iface_groups[0].name, my_iface_groups[1].name))
    # If the user wants to export the interface between the two selected groups, do it
    if export_iface_prompt in ys:
        export_iface(groups=my_groups, info_file=True, interface_atoms=True)


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
        build(sys, usr_npt)
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
