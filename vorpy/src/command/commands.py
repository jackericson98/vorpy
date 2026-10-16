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
   - Edge width settings (edge_widths)
   - Vertex size settings (vertex_sizes)
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
max_verts = ['max_vert', 'mv', 'maximum_vertex', 'max_vertex', 'maximum_vert', 'verts', 'vert_size', 'max_vert_size', 'mvs', 'vert_max', 'vertex_max', 'vertex_maximum']
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
mesh_formats = ["ft", "file_type", "file_format", "format", "mesh_format"]
edge_widths = ["ew", "edge_width", "edge_radius"]
vertex_sizes = ["vs", "vertex_size", "vertex_radius", "vert_radius"]

mesh_format_vals = {
    "off": "off",
    "ply": "ply",
    "vtp": "vtp",
    "vtk": "vtp",
}

# Settings vals
power_vals = ['pow', 'power', 'p', 'pwr', 'pwizzle']
voronoi_vals = ['vor', 'voronoi', 'vnoi', 'voron', 'vn', 'aw', 'additively_weighted', 'add_weight', 'awizzle']
delaunay_vals = ['del', 'dl', 'delaunay', 'dlny', 'dny', 'prm', 'primitive', 'prim', 'ptiv', 'ptizzle']
compare_vals = ['c', 'com', 'compare', 'cpr', 'compar', 'cum']

surf_scheme_curv_vals = ['curv', 'c', 'curvature']
surf_scheme_mean_vals = ['mean', 'mc', 'mean_curvature', 'm']
surf_scheme_gaus_vals = ['gaussian', 'gc', 'gauss', 'g']
surf_scheme_dist_vals = ['distance', 'dist', 'd']
surf_scheme_nout_vals = ['in_out', 'nout', 'no', 'ins_out']

# Integrated curvature surface-coloring schemes
surf_scheme_int_mean_vals = [
    'int_mean_curv', 'int_mean', 'integrated_mean_curvature',
    'integrated_mean', 'imc'
]
surf_scheme_int_mean_sq_vals = [
    'int_mean_curv_sq', 'int_mean_sq', 'integrated_mean_curvature_squared',
    'integrated_mean_squared', 'imcs'
]
surf_scheme_int_gauss_vals = [
    'int_gauss_curv', 'int_gauss', 'integrated_gaussian_curvature',
    'integrated_gauss', 'igc'
]

surf_scheme_energy_vals = [
    'surf_energy',
    'surface_energy',
    'energy',
    'se'
]

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
                 'nt': 'Network Type', 'sc': 'Surface Color Map', 'ss': 'Surface Coloring Scheme',
                 'ew': 'Edge Width', 'vs': 'Vertex Size'}


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
                       "  5. flat_surfs: Build the surfaces flat (True/False)     ",
                       "  6. ew : Exported edge tube radius (default 0.05 A)",
                       "  7. vs : Exported vertex sphere radius (default 0.08 A)"]

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
VorPy command-line help
---------------------------------

Usage:
  python -m vorpy <input> [options]

Input:
  A structure accepted by the legacy loader (.pdb, .mol, .gro, or .cif), or a
  solved .vpy network archive.  Use --load-network before an archive path when
  loading an archive explicitly.

Workflow options:
  -l, --load <file> [and <file> ...]
    Load additional vertex, log, Voronota, or GROMACS index files.

  -s, --settings, --set <setting> <value>
    Set solver/export values such as sr (surface resolution), nt (network type),
    mv (maximum vertex), bm (box multiplier), sc, ss, sf, ft, ar, ew, and vs.
    The current CLI execution default for mv is 5 A; lower-level legacy helpers
    still default to 40 A.  The accepted mv range is 0.5-5000 A.  This split is
    intentionally unchanged pending scientific-owner confirmation.

  -g, --group <selector> [and <selector> ...]
    Select balls, atoms, residues, or chains.  Repeat -g for multiple groups.

  -b, --build
    Request the legacy build command.

  -i, --interface
    Request interface calculations from selected groups.

  -c, --calculate <target>
    Request system, group, interface, or chain calculations.
    For ordered typed A/B Results, combine -i with -c system; without a
    validated full-system source, the interface status remains UNRESOLVED.

  -e, --export <type> [and <type> ...]
    Export a preset or component, such as small, medium, large, all, shell,
    logs, surfs, edges, verts, or no_archive.

Analysis and execution:
  --profile geometry|interface|analysis|full
    Select the execution profile.

  --apollonius
    Request the AW Apollonius/dual interface representation.

  --alpha-complex <alpha>
    Request alpha-complex filtering at the given value in A.

  --power-interface cazals
    Run the Cazals-compatible Power interface workflow with two chain groups.
  --power-probe-radius <A>  --power-alpha <value>  --power-M <value>
  --power-vs-aw
    Configure or compare the Power interface workflow.

  --aw-alpha-experimental <alpha>
    Export experimental AW alpha-filtration data; repeat for multiple values.
  --aw-alpha-pairs-only-experimental
    Limit that experimental export to pair-overlap analysis.

Archive, frames, and diagnostics:
  --fail-if-exists
    Refuse existing output directories, archives, timing sidecars, and cache
    profiles; reserve final targets atomically before/at save finalization.
  --save-network <path.vpy>
    Save the solved network archive at the requested path.
  --load-network <path.vpy>
    Load a solved network archive; this option must precede the archive path.
  --profile-cache-export
    Profile archive-load/export work without rebuilding geometry (.vpy only).
  --timing-report <path.json>
    Write opt-in machine-readable load, network, interface, analysis, export,
    archive, and total wall-time measurements.
  --all-frames
    Process every frame in a PDB trajectory sequentially without prompting.
  --parallel-frames <workers>
    Process every PDB frame with up to this many worker processes.
    Implies --all-frames. Each frame has its own output folder and run.log.
  Input naming convention
    X.pdb is the solvated post-MD input, X_no_sol.pdb its dry pair, and
    X_frames/X_fN.pdb a named frame; the CLI reports this context without
    treating filenames as scientific certification.
  -v, --verbose
    Enable verbose build output.
  --diagnose-edges
    Audit constructed AW edges without changing exported geometry.

Boundary and solvent:
  --boundary none|shell
  --probe-radius <A>  --shell-spacing <A>  --shell-offset <A>
    Configure the optional virtual solvent shell boundary.
  --solvent-mode MOLECULE_ONLY|SOLVENT_CONTEXT|FULL_SOLVENT
  --solvent-context  --full-solvent
    Select the generator universe used for the solve.

Visualization:
  --export-dual-visualization
  --export-alpha-visualization <alpha>
    Export the requested AW dual or alpha visualization.

Use -h or --help to show this page.  Flags may be combined with the legacy
`and` separator where the selected command accepts multiple values.
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
