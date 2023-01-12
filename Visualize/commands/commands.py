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
calc_surfses = ['calc_surfs', 'cs']

my_settings = surf_reses + max_verts + box_sizes + calc_surfses


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



    commands_header = ["Commands:                                                                                                                   ",
                       "  1. load  : Loads file addresses for System (.pdb, .gro, .mol, .cif), Network (.txt), vertices (.txt) or index (.ndx) files",
                       "  2. set   : Sets build Settings with values (float, float, float, True)                                                    ",
                       "  3. build : Builds a Network for the System with the current settings                                                      ",
                       "  4. export: Exports System objects. Use \'and\' to export the interface between two groups (e.g. \'export ndx 1 and atom 3\')  ",
                       "  5. show  : Shows System elements in a given object category                                                               ",
                       "  6. quit  : Quits from the current process                                                                                 "]

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
                       "  4. calc_surfs: Calculate the network's surfaces (True/False/Group,        ",
                       "                 recommended False)                                         "]

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
