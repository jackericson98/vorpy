import matplotlib as mpl
from Visualize.cmnd.interpret import *
from radii.Standard import element_radii, special_radii
from System.chemistry_interpreter import element_names, residue_names


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
    # First extract the value from the list if it is in fact a list
    if type(surface_color) is list:
        surface_color = surface_color[0]
    # Try each of the three possible options for surface coloring
    try:
        my_cmap = mpl.colormaps.get_cmap(surface_color)
        print("surface color set to {}".format(surface_color))
        return surface_color
    except Exception as e:
        pass
    # Try each of the three possible options for surface coloring
    try:
        my_cmap = mpl.cm.get_cmap(surface_color)
        print("surface color set to {}".format(surface_color))
        return surface_color
    except Exception as e:
        pass
    # Try each of the three possible options for surface coloring
    try:
        my_cmap = mpl.cm.get_cmap(surface_color)
        print("surface color set to {}".format(surface_color))
        return surface_color
    except Exception as e:
        pass
    # If none of the formatting options work print the error and return
    print('{} is not a matplotlib colormap. Please choose a valid matplotlib colormap (e.g. \"viridis\", '
          '\"plasma\", \"inferno\", \"cividis\", \"Greys\", \"Reds\", \"Greens\", \"Blues\", \"rainbow\"'
          .format(surface_color))
    return settings['surf_col']


def set_ss(surf_scheme, settings):
    # Make sure to extrac the surface scheme from the value
    if type(surf_scheme) is list:
        surf_scheme = surf_scheme[0]
    # Set up the list of different dictionaries
    all_dicts = [{_: 'gauss' for _ in surf_scheme_gaus_vals}, {_: 'dist' for _ in surf_scheme_dist_vals},
                 {_: 'mean' for _ in surf_scheme_mean_vals + surf_scheme_curv_vals},
                 {_: 'ins_out' for _ in surf_scheme_nout_vals}, {_: 'none' for _ in nones}]
    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}
    # Check that the scheme entered is in the set of
    if surf_scheme not in interpreter:
        # Print a warning that the user has entered the wrong scheme
        print('{} is not a valid entry for surface coloring scheme. Please enter one of the following: \"curv\", '
              '\"mean\", \"gaussian\", \"dist\", \"ins_out\", or \"none\"'.format(surf_scheme))
        return settings['surf_scheme']
    print('Surface Scheme set to {}'.format(interpreter[surf_scheme]))
    return interpreter[surf_scheme]


def set_sf(surface_factor, settings):
    if surface_factor[0].lower() in surf_factor_vals:
        return surf_factor_vals[surface_factor[0].lower()]
    print('{} is not a valid entry for surface coloring scale. Please enter one of the following \"lin\", \"log\", '
          '\"sqr\", \"cub\"'.format(surface_factor[0]))
    return settings['surf_factor']


def set_ar(element_radius, settings):

    # Create the changes list
    change_settings = {'element': {}, 'special': {}}
    if settings['atom_rad'] is not None:
        change_settings = settings['atom_rad']

    # Separate the element from the radius
    if len(element_radius) >= 3 and element_radius[0] not in atom_objs:
        # Get the residue
        residue, name, radius = element_radius[:3]
        # Check that this exists
        if residue.lower() in residue_names and residue_names[residue.lower()] in special_radii:
            if name.upper() in special_radii[residue_names[residue.lower()]]:
                try:
                    my_radius = float(radius)
                    print('All {} {} atoms radii changed from {} \u212B to {} \u212B'
                          .format(residue, name, special_radii[residue.upper()][name.upper()], radius))
                    # Add the radius
                    if residue_names[residue.lower()] in change_settings['special']:
                        change_settings['special'][residue_names[residue.lower()]][name.upper()] = my_radius
                        return change_settings
                    change_settings['special'][residue_names[residue.lower()]] = {name.upper(): my_radius}
                    return change_settings
                except ValueError:
                    print('{} is not a valid entry for radius. Please try a valid float entry')
                    return
            print('{} is not an atom in {}. Please try one of the following names: {}'
                  .format(name, residue, [_ for _ in special_radii[residue.upper()]]))
            return
        new_elem_rad = input('{} contains an invalid entry. Please re-enter your atom radius changing setting >>>   '.format(element_radius))
        new_elem_rad = new_elem_rad.split(' ')
        return set_ar(new_elem_rad, settings)

    # The case where the user wants to change just the element or all atoms with a certain name
    elif len(element_radius) == 2 or element_radius[0] in atom_objs:
        if element_radius[0] in atom_objs:
            element_radius = element_radius[1:]
        # If the changed name is in the regular elements use that
        if element_names[element_radius[0].lower()] in element_radii:
            # Check the value and that it is a float
            try:
                # Try creating a float value for the new radius
                my_radius = float(element_radius[1])
                print('All {} atoms radii changed from {} \u212B to {} \u212B.'
                      .format(element_names[element_radius[0].lower()], element_radii[element_names[element_radius[0].lower()]], element_radius[1]))
                change_settings['element'][element_radius[0].upper()] = my_radius
                return change_settings
            except ValueError:
                # Print the error saying the value is wrong
                print('{} is not a valid entry for radius'.format(element_radius[1]))
                return
        # Check special radii for specific changes (e.g. all alpha carbons)
        elif any([element_radius[0].upper() in special_radii[_] for _ in special_radii]):
            try:
                my_radius = float(element_radius[1])
            except ValueError:
                print('{} is not a valid entry for radius'.format(element_radius[1]))
                return
            # Loop through the special radii
            for residue in special_radii:
                if element_radius[0].upper() in residue:
                    if residue in change_settings['special']:
                        change_settings['special'][residue][element_radius[0].upper()] = my_radius
                    else:
                        change_settings['special'][residue] = {element_radius[0].upper(): my_radius}
            print('All {} radii changed to {}'.format(element_radius[0], element_radius[1]))
            return change_settings


def set_bt(build_type, settings):
    if build_type == 'logs':
        return 'logs'


def sett(setting, value, settings=None):
    """
    Take in all information after '-s' and return an updated dictionary
    """
    # Set the default settings
    if settings is None:
        settings = {'surf_res': 0.2, 'max_vert': 40, 'box_size': 1.25, 'net_type': 'aw', 'surf_col': 'plasma_r',
                    'surf_scheme': 'curv', 'scheme_factor': 'lin', 'atom_rad': None, 'bld_type': None}
    # Set up the functions dictionary to return the value
    func_dict = {'surf_res': set_sr, 'max_vert': set_mv, 'box_size': set_bs, 'net_type': set_nt, 'surf_col': set_sc,
                 'surf_scheme': set_ss, 'scheme_factor': set_sf, 'atom_rad': set_ar, 'bld_type': set_bt}

    # Set up the interpretation dictionary
    all_dicts = [{_: 'surf_res' for _ in surf_reses}, {_: 'max_vert' for _ in max_verts},
                 {_: 'box_size' for _ in box_sizes}, {_: 'net_type' for _ in net_types},
                 {_: 'surf_col' for _ in surf_colors}, {_: 'surf_scheme' for _ in surf_schemes},
                 {_: 'scheme_factor' for _ in surf_factors}, {_: 'atom_rad' for _ in atom_radii}]

    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}

    # Set the setting
    settings[interpreter[setting]] = func_dict[interpreter[setting]](value, settings)

    # Return the settings
    return settings
