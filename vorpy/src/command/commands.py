"""
This module defines command-related constants and command lists used throughout the Voronoi tessellation program.

The module contains several categories of command lists:

1. Response Lists:
   - Affirmative responses (ys)
   - Negative responses (ns)
   - None responses (nones)
   - Boolean responses (trues, falses)
   - Completion indicators (dones)
   - Conjunctions (ands)
   - Path separators (splitters)
   - File browser commands (browses)

2. General Commands:
   - Quit commands (quits)
   - Help commands (helps)

3. Main Commands:
   - Show commands (show_cmds)
   - Load commands (load_cmds)
   - Set commands (set_cmds)
   - Build commands (build_cmds)
   - Group commands (group_cmds)
   - Export commands (export_cmds)

4. Object Types:
   - Full system objects (full_objs)
   - No-solvent objects (noSOL_objs)
   - Chain/molecule objects (chn_objs)
   - Atom objects (atom_objs)
   - Residue objects (res_objs)
   - Index/group objects (ndx_objs)

5. Settings:
   - Surface resolution settings (surf_reses)
   - Maximum vertex settings (max_verts)
   - Box size settings (box_sizes)
   - Network type settings (net_types)
   - Surface color settings (surf_colors)
   - Surface scheme settings (surf_schemes)
   - Atom radius settings (atom_radii)
   - Surface factor settings (surf_factors)

6. Setting Values:
   - Power values (power_vals)
   - Voronoi values (voronoi_vals)

These command lists are used for command parsing and validation throughout the program, providing a consistent set of recognized commands and responses.
"""


# Responses
ys = ['y', 'yes', 'ya', 'yeet', 'yur', 'yoint', 'uhu', 'yup', 'jess', 'affirmative', 'yuss', 'yess', 'jeth',
      'yesss', 'yessss', 'yar', 'yuh', 'mhm', 'crabsolutely', 'dolphinitely', 'shell ya', 'whale of course']
ns = ['n', 'no', 'naur', 'nope', 'nonya', 'nope', 'nien', 'nada']
nones = ['none', 'noneya']
trues = ['t', 'true', 'tr', 'tru', 'truth', 'tu']
falses = ['f', 'false', 'flse', 'fl', 'fa', 'fs', 'fls']
dones = ['done', 'd', 'finished', 'finito', 'complete', 'doneso', 'don', 'fin']
ands = ['&', 'and', 'nd', 'also', '+', '&&']
splitters = ['/', '-']
browses = ['choose', 'select', 'browse', 'chs', 'brwse', 'brows', 'chose', 'get', 'file', 'folder']

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
full_objs = ['f', 'full', 'fl', 'ful', 'fs']
noSOL_objs = ['ns', 'nosol', 'no_sol', 'nos', 'nsol']
chn_objs = ['m', 'ms', 'molecule', 'molecules', 'mols', 'ml', 'mls', 'c', 'cs', "chain"]
atom_objs = ['a', 'as', 'atom', 'atoms', 'at', 'ats', 'am', 'ams']
res_objs = ['r', 'rs', 'residue', 'residues', 'resid', 'resids', 'res', 'ress', 'reses', 'rdue', 'rdues']
ndx_objs = ['i', 'is',  'index', 'indexs', 'indexes', 'indices', 'ndx', 'ndxs', 'ndex', 'group', 'g', 'grp', 'n']

my_objects = full_objs + noSOL_objs + chn_objs + res_objs + atom_objs + ndx_objs

# Settings
surf_reses = ['surf_res', 'sr', 'surface_resolution', 'surface_res', 'surf_resolution', 'surfs', 'surf', 'surfs_res', 'surfs_resolution', 'surfaces_resolution', 'surfaces_res']
max_verts = ['max_vert', 'mv', 'maximum_vertex', 'max_vertex', 'maximum_vert', 'verts', 'vs', 'vert_size', 'max_vert_size', 'mvs', 'vert_max', 'vertex_max', 'vertex_maximum']
box_sizes = ['box_size', 'bm', 'box', 'bx_sz', 'size_box', 'containing_box', 'containing_box_size', 'box_multi', 'box_multiplier', 'bs']
net_types = ['st', 'surf_type', 'net_type', 'nt']
surf_colors = ['sc', 'surf_colors', 'surf_color_map', 'surf_col', "scm"]
surf_schemes = ['ss', 'sirf_scheme', 'surf_scm']
atom_radii = ['ar', 'atom_radius', 'element_radius']
surf_factors = ['sf', 'surf_factor', 'surf_fac', 'surff']
build_surfses = ['bs', 'build_surfs', 'build_surf', 'build_surfses', 'build_surfs_es']
conc_cols = ['cc', 'conc_col', 'concave_col', 'concave_color', 'concave_surfaces', 'concave_colors']
vert_cols = ['vc', 'vert_col', 'vertex_color', 'vertex_colors']
edge_cols = ['ec', 'edge_col', 'edge_color', 'edge_colors']
round_tos = ["rt", "round", "round_to", "decimals", "decimal_places"]

# Settings vals
power_vals = ['pow', 'power', 'p', 'pwr', 'pwizzle']
voronoi_vals = ['vor', 'voronoi', 'vnoi', 'voron', 'vn', 'aw', 'additively_weighted', 'add_weight', 'awizzle']
delaunay_vals = ['del', 'dl', 'delaunay', 'dlny', 'dny', 'prm', 'primitive' 'prim', 'ptiv', 'ptizzle']
compare_vals = ['c', 'com', 'compare', 'cpr', 'compar', 'cum']

surf_scheme_curv_vals = ['curv', 'c', 'curvature']
surf_scheme_mean_vals = ['mean', 'mc', 'mean_curvature', 'm']
surf_scheme_gaus_vals = ['gaussian', 'gc', 'gauss', 'g']
surf_scheme_dist_vals = ['distance', 'dist', 'd']
surf_scheme_nout_vals = ['in_out', 'nout', 'no', 'ins_out']

surf_factor_vals = {
    **{_: 'lin' for _ in {'linear', 'lin', 'line'}},
    **{_: 'log' for _ in {'log', 'logarithmic', 'log'}},
    **{_: 'sqr' for _ in {'sqr', 'square', 'sq', 'squared'}},
    **{_: 'cub' for _ in {'cub', 'cube', 'cubed', 'cubaroonski'}}
}


file_types = ['net', 'vert', 'ball', 'ndx']
browse_names = {'browse', 'choose', 'brwse', 'chse', 'get', 'find'}


net_type_dict = {'pow': "Power", 'del': "Primitive", 'vor': "Additively Weighted"}

settings_dict = {'sr': 'Surface Resolution', 'mv': 'Maximum Vertex', 'bm': 'Box Multiplier', 'bs': 'Build Surfaces?',
                 'nt': 'Network Type', 'sc': 'Surface Color Map', 'ss': 'Surface Coloring Scheme'}


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


def print_help():
    print("""
Vorpy Help
---------------------------------

Usage:
  python vorpy.py <file> [options]

File:
  The first argument after 'vorpy.py' should be the file address of the ball or atom file.
  If the file is located in the 'Data/test_data' folder, specify the file name without the path or extension.
  Accepted file extensions include .pdb, .mol, .gro, .cif.

Options:
  -l <file>
    Load additional files like vertex files from previous runs, log files, Voronota vertex files, or GROMACS index files.

  -s <setting value>
    Adjust various simulation parameters:
      sr - Surface Resolution: Default = 0.2
      nt - Network Type: Default = Additively Weighted 'aw', Power 'pow', Primitive 'prm', or Compare 'com 'type1' 'type2''
      mv - Maximum Vertex: Default = 40
      bm - Box Multiplier: Default = 1.25
      sc - Surface Color: Default = 'viridis', 'plasma', 'rainbow', or any other matplotlib colormap
      ss - Surface Scheme: Default = curvature 'curv', inside vs outside spheres 'nout', distance from center 'dist'
      sf - Surface Coloring Scale: Default = linear 'lin', log 'log', squared 'square', cube 'cube'
      ar - Adjust Radii: 'element' 'value' or 'atom name' 'value' or 'residue' 'atom name' 'value'

  -g <identifier>
    Select specific balls or molecular elements using identifiers such as:
      b - Ball Identifier with index or range 'index1-index2'
      a - Atom Identifier with element, name, index, or range 'index1-index2'
      r - Residue Identifier with name, sequence number, index, or range 'index1-index2'
      c - Chain Identifier with name, index, or range 'index1-index2'
    Note: Use 'and' to combine components within the same group. Use multiple -g flags for multiple groups.

  -c <calculation_type>
    Identify additional calculations like interfaces between groups (iface, ifc, i) or calculating vertices up to specified layers (layers)

  -e <export_type>
    Specify the intensity and type of exports:
      Options include: small, medium, large, all
      Export choices: pdb, mol, cif, gro, set_atoms, info, logs, surfs, sep_surfs, edges, sep_edges, verts, sep_verts, shell, shell_edges, shell_verts, atoms, surr_atoms

Examples:
  Solve the network for tyrosine 2 and methionine 1 in the cambrin molecule, calculate their interface, and export large type results:
    python vorpy.py cambrin -s sr 0.05 and mv 80 -g tyr 2 -g met 1 -c iface -e large

  Calculate the primitive and power networks for the mg atom in the EDTA_Mg molecule and compare the difference:
    python vorpy.py EDTA_Mg -s nt compare prm pow -g mg

  Solve the network for hairpin and export the shell with inside and outside parts of the surfaces highlighted at high resolution:
    python vorpy.py hairpin -s ss nout and sr 0.01 -e shell and pdb

Note:
  Each option flag and its arguments must be separated by spaces.
  To use multiple commands for a single option, use 'and' or repeat the flag (except for groups to avoid creating multiple groups).
""")


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
            return 'q'
        elif my_response.lower() in helps:
            print_help()
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
            quit_check = print_list(names[nums[0]:nums[0] + num_cols * height],
                                    list_name=list_name + ": elements " + str(nums[0]) + "-" + str(nums[0] + num_cols * height),
                                    height=height, width=width, cutoff=max_len)
            if quit_check in quits:
                return
        elif nums is not None and nums[0] < nums[1] < len(names):
            # Check to see if the height needs to be changed
            new_height = height
            if nums[1] - nums[0] > num_cols * height:
                new_height = (nums[1] - nums[0]) // num_cols + 1
            quit_check = print_list(names[nums[0]:nums[1]], list_name=list_name + ": elements " + str(nums[0]) + "-" + str(nums[1]),
                                    height=new_height, width=width, cutoff=max_len)
            if quit_check in quits:
                return
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
            print_help()
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
        objs = [chn_objs, res_objs, atom_objs, ndx_objs]
        for i in range(4):
            if my_input.lower() in objs[i]:
                if len([sys.chains, sys.residues, sys.atoms, sys.ndxs][i]) > 0:
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

    # Get the 181L list to show the user
    if show_var in chn_objs:
        show_name = "{} Chains".format(sys.name)
        show_list = sys.chn_names
    elif show_var in res_objs:
        show_name = "{} Residues".format(sys.name)
        show_list = sys.res_names
    elif show_var in atom_objs:
        show_name = "{} Atoms".format(sys.name)
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
