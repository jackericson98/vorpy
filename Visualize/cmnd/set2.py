import matplotlib as mpl
from Visualize.cmnd.interpret import *
from System.radii import element_radii, special_radii
from matplotlib._api.deprecation import MatplotlibDeprecationWarning as MPLDepWarn



def set_sr(surf_res, settings):
    # Quick catch if the max_vert value is in the form of a list
    if type(surf_res) is list:
        surf_res = surf_res[0]
    # try making the value a float value for use later
    try:
        # First set the value to a float value
        good_val = float(surf_res)
        # Check to see if it is within the range
        if not 0.01 <= good_val <= 10:
            print('surface resolution out of range (0.01 to 10 \u212B)')
            return settings['surf_res']
        # Print a confirmation that the setting has been changed
        print("surface resolution set to {} \u212B".format(good_val))
    except ValueError:
        # Tell the user they messed up and neet to get their life together
        print("\"{}\" is an invalid input for the surface resolution setting. Enter a float value "
              "(from 0.01 to 10 \u212B, recommended 0.1 \u212B)".format(surf_res))
        return settings['surf_res']
    return good_val


def set_mv(max_vert, settings):
    # Quick catch if the max_vert value is in the form of a list
    if type(max_vert) is list:
        max_vert = max_vert[0]
    # Try setting the maximum vertex value to a float for verification it works
    try:
        # First make it a float value
        good_val = float(max_vert)
        # Check to see if it is out of range
        if not 0.5 <= good_val <= 5000:
            print('maximum vertex out of range (0.5 to 5000 \u212B)')
            return settings['max_vert']
        print(u"maximum vertex radius set to {} \u212B".format(max_vert))
    except ValueError:
        print("\"{}\" is an invalid input for the maximum vertex radius setting. Enter a float value "
              "(From 0.10 to 20 A, recommended 7 A)".format(max_vert))
        return settings['max_vert']
    return good_val


def set_bs(box_size, settings):
    # Quick catch if the max_vert value is in the form of a list
    if type(box_size) is list:
        box_size = box_size[0]
    # Try setting the box size multiplier to a float value for verification it is the right user input
    try:
        # Make it a float value
        good_val = float(box_size)
        # Check that it is within range
        if not 1.0 < good_val < 10:
            print('box size multiplier out of range (1.0 to 10x)')
            return settings['box_size']
        print("box size multiplier set to {} x".format(good_val))
    except ValueError:
        print("\"{}\" is an invalid input for the box size multiplier setting. Enter a float value "
              "(From 1.0 to 10.0 X, recommended 1.5 X)".format(box_size))
        return settings['box_size']
    return good_val


def set_nt(net_type, settings):
    # Set up the list of different dictionaries
    all_dicts = [{_: 'aw' for _ in voronoi_vals}, {_: 'pow' for _ in power_vals}, {_: 'prm' for _ in delaunay_vals},
                 {_: 'com' for _ in compare_vals}]
    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}
    # If the net type is a list and the list contains the nets for comparison
    set_nets = []
    if type(net_type) is list:
        if len(net_type) > 1:
            set_nets = net_type[1:]
        net_type = net_type[0]
    # Make sure the net type is in the possible names
    if net_type not in interpreter:
        print('{} is not a valid network type. Please enter \'aw\', \'pow\', \'prm\', or \'com\''.format(net_type))
        return settings['net_type']
    # If we are comparing the network types
    if interpreter[net_type] == 'com':
        # Check to see if the set nets are available and at the very end add 'aw' and power so returned worst case
        set_nets = [interpreter[_] for _ in set_nets] + ['aw', 'pow']
        # Return the comparisons
        return [interpreter[net_type], set_nets[0], set_nets[1]]
    # Return the interpreted network type
    return interpreter[net_type]


def set_sc(surface_color, settings):
    try:
        my_cmap = mpl.colormaps.get_cmap(surface_color)
    except MPLDepWarn:
        my_cmap = mpl.cm.get_cmap(surface_color)
    except AttributeError:
        my_cmap = mpl.cm.get_cmap(surface_color)
    except Exception as e:
        print('{} is not a matplotlib colormap. Please choose a valid matplot lib color map (e.g. \"viridis\", '
              '\"plasma\", \"inferno\", \"cividis\", \"Greys\", \"Reds\", \"Greens\", \"Blues\", \"rainbow\"'
              .format(surface_color))
        return settings['surf_col']
    # Check that the surface color is in the possible lists of matplotlib colormaps
    print("surface color set to {}".format(surface_color))
    return surface_color


def set_ss(surf_scheme, settings):
    # Set up the list of different dictionaries
    all_dicts = [{_: 'curv' for _ in surf_scheme_curv_vals}, {_: 'dist' for _ in surf_scheme_dist_vals},
                 {_: 'ins_out' for _ in surf_scheme_nout_vals}, {_: 'none' for _ in nones}]
    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}
    # Check that the scheme entered is in the set of
    if surf_scheme not in interpreter:
        # Print a warning that the user has entered the wrong scheme
        print('{} is not a valid entry for surface coloring scheme. Please enter one of the following: \"curv\", '
              '\"dist\", \"ins_out\", or \"none\"'.format(surf_scheme))
        return settings['surf_scheme']
    return interpreter[surf_scheme]


def set_sf():
    pass


def set_ar(element_radius, settings):
    # Separate the element from the radius
    if len(element_radius) >= 3:
        residue = element



    elif len(element_radius) == 2:
        element = element_radius
    else:



    my_element, new_rad = my_val[0].strip().lower(), my_val[1]
    # Normalize 'element' column for comparison
    normalized_elements = sys.atoms['element'].str.strip().str.lower()
    matching_indices = sys.atoms.index[normalized_elements == my_element].tolist()
    count_changed = 0
    if matching_indices:
        # Record the old value of the first matching entry
        old_value = sys.atoms.loc[matching_indices[0], 'rad']

        # Replace 'rad' values where 'element' matches 'my_element'
        sys.atoms.loc[normalized_elements == my_element, 'rad'] = new_rad

        # Count and print the number of changed values
        count_changed = len(matching_indices)

        element_radii[my_val[0]] = my_val[1]
        print(u"{} atoms changed from {} to {}".format(count_changed, old_value, my_val[1]))
    else:
        print("No matching element found.")


def sett(setting, value, settings=None):
    """
    Take in all information after '-s' and return an updated dictionary
    """
    # Set the default settings
    if settings is None:
        settings = {'surf_res': 0.2, 'max_vert': 40, 'box_size': 1.25, 'net_type': 'aw', 'surf_col': 'plasma',
                    'surf_scheme': 'curv', 'scheme_factor': 'log', 'atom_rad': element_radii}
    # Set up the functions dictionary to return the value
    func_dict = {'surf_res': set_sr, 'max_vert': set_mv, 'box_size': set_bs, 'net_type': set_nt, 'surf_col': set_sc,
                 'surf_scheme': set_ss, 'scheme_factor': set_sf, 'atom_rad': set_ar}

    # Set up the interpretation dictionary
    all_dicts = [{_: 'surf_res' for _ in surf_reses}, {_: 'max_vert' for _ in max_verts},
                 {_: 'box_size' for _ in box_sizes}, {_: 'net_type' for _ in net_types},
                 {_: 'surf_col' for _ in surf_colors}, {_: 'surf_scheme' for _ in surf_schemes},
                 {_: 'scheme_factor' for _ in surf_schemes}, {_: 'atom_rad' for _ in atom_radii}]

    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}

    # Set the setting
    settings[interpreter[setting]] = func_dict[interpreter[setting]](value)

    # Return the settings
    return settings