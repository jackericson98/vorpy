import matplotlib as mpl
from vorpy.src.command.interpret import *
from vorpy.src.chemistry import element_radii
from vorpy.src.chemistry import special_radii
from vorpy.src.chemistry import element_names
from vorpy.src.chemistry import residue_names


def _vorpy_print_setting(setting_name, value, note=None):
    msg = f"Vorpy Setting - {setting_name} - set to {value}"
    if note:
        msg += f" {note}"
    print(msg)


def _vorpy_print_error(setting_name, value, valid_desc, example=None):
    msg = f"Vorpy Setting - {setting_name} - invalid input \"{value}\". Enter {valid_desc}"
    if example:
        msg += f" (e.g., {example})"
    print(msg)


def set_sr(surf_res, settings, print_change=False):
    """
    Sets the surface resolution parameter for the system.

    This function validates and sets the surface resolution value, which controls
    the level of detail in the generated surface mesh. The resolution determines
    how finely the surface is sampled, with higher values resulting in more detailed
    but computationally expensive meshes.

    Parameters:
    -----------
    surf_res : float or list
        The desired surface resolution value in angstroms. If a list is provided,
        the first element will be used.
    settings : dict
        Dictionary containing current system settings

    Returns:
    --------
    float
        The validated surface resolution value if successful, or the current
        surface resolution value if validation fails

    Notes:
    ------
    - Valid range: 0.01 to 10 angstroms
    - Recommended value: 0.1 angstroms
    - Values outside the valid range will result in an error message and
      return the current setting
    """
    # Quick catch if the max_vert value is in the form of a list
    if type(surf_res) is list:
        surf_res = surf_res[0]
    # try making the value a float value for use later
    try:
        # First set the value to a float value
        good_val = float(surf_res)
        # Check to see if it is within the range
        if not 0.001 <= good_val <= 10:
            if print_change:
                _vorpy_print_error(
                    "Surface Resolution",
                    good_val,
                    "a float value from 0.001 to 10 \u212B (recommended 0.1 \u212B)"
                )
            return settings['surf_res']
        # Print a confirmation that the setting has been changed
        if print_change:
            _vorpy_print_setting("Surface Resolution", f"{good_val} \u212B")
        return good_val
    except ValueError:
        _vorpy_print_error(
            "Surface Resolution",
            surf_res,
            "a float value from 0.001 to 10 \u212B (recommended 0.1 \u212B)"
        )
        return settings['surf_res']


def set_mv(max_vert, settings, print_change=False):
    """
    Sets the maximum vertex radius parameter for the system.

    This function validates and sets the maximum vertex radius value, which controls
    the size of the largest allowed vertex in the generated surface mesh. The maximum
    vertex radius determines the coarseness of the mesh, with larger values resulting
    in fewer but larger vertices.

    Parameters:
    -----------
    max_vert : float or list
        The desired maximum vertex radius value in angstroms. If a list is provided,
        the first element will be used.
    settings : dict
        Dictionary containing current system settings

    Returns:
    --------
    float
        The validated maximum vertex radius value if successful, or the current
        maximum vertex radius value if validation fails

    Notes:
    ------
    - Valid range: 0.5 to 5000 angstroms
    - Recommended value: 7 angstroms
    - Values outside the valid range will result in an error message and
      return the current setting
    """
    # Quick catch if the max_vert value is in the form of a list
    if type(max_vert) is list:
        max_vert = max_vert[0]
    # Try setting the maximum vertex value to a float for verification it works
    try:
        # First make it a float value
        good_val = float(max_vert)
        # Check to see if it is out of range
        if not 0.5 <= good_val <= 5000:
            if print_change:
                _vorpy_print_error(
                    "Maximum Vertex",
                    good_val,
                    "a float value from 0.5 to 5000 \u212B (recommended 7 \u212B)"
                )
            return settings['max_vert']
        if print_change:
            _vorpy_print_setting("Maximum Vertex", f"{good_val} \u212B")
        return good_val
    except ValueError:
        _vorpy_print_error(
            "Maximum Vertex",
            max_vert,
            "a float value from 0.5 to 5000 \u212B (recommended 7 \u212B)"
        )
        return settings['max_vert']


def set_bs(box_size, settings, print_change=False):
    """
    Sets the box size multiplier parameter for the system.

    This function validates and sets the box size multiplier value, which controls
    the size of the containing box relative to the molecular system. The box size
    multiplier determines how much empty space surrounds the system, with larger
    values resulting in more space for surface generation.

    Parameters:
    -----------
    box_size : float or list
        The desired box size multiplier value. If a list is provided,
        the first element will be used.
    settings : dict
        Dictionary containing current system settings

    Returns:
    --------
    float
        The validated box size multiplier value if successful, or the current
        box size multiplier value if validation fails

    Notes:
    ------
    - Valid range: 1.0 to 10.0
    - Recommended value: 1.5
    - Values outside the valid range will result in an error message and
      return the current setting
    """
    # Check if box_size is a list, and if so, use its first element
    if isinstance(box_size, list):
        box_size = box_size[0]
    try:
        # Attempt to convert box_size to a float
        good_val = float(box_size)
        # Validate that the value is within the allowed range [1.0, 10.0]
        if not (1.0 <= good_val <= 10.0):
            # If not valid and print_change is True, print an error message
            if print_change:
                _vorpy_print_error(
                    "Box Size Multiplier",
                    good_val,
                    "a float value from 1.0 to 10.0 (recommended 1.5)"
                )
            # Return the current setting if invalid
            return settings['box_size']
        # If print_change is True, print a confirmation message
        if print_change:
            _vorpy_print_setting("Box Size Multiplier", f"{good_val} x")
        # Return the validated value
        return good_val
    except (ValueError, TypeError):
        # If conversion fails, print an error message
        _vorpy_print_error(
            "Box Size Multiplier",
            box_size,
            "a float value from 1.0 to 10.0 (recommended 1.5)"
        )
        # Return the current setting if conversion fails
        return settings['box_size']


def set_nt(net_type, settings, print_change=False):
    """
    Sets the network type parameter for the system.

    This function handles the configuration of the network type used for surface generation.
    It supports multiple network types including:
    - Additively weighted Voronoi ('aw')
    - Power diagram ('pow')
    - Primitive/Delaunay ('prm')
    - Comparison mode ('com') for comparing different network types

    Parameters:
    -----------
    net_type : str or list
        The desired network type specification. If a list is provided,
        the first element specifies the type and subsequent elements
        specify networks to compare in comparison mode.
    settings : dict
        Dictionary containing current system settings

    Returns:
    --------
    str or list
        - For single network type: Returns the validated network type code
        - For comparison mode: Returns a list containing ['com', net1, net2]
          where net1 and net2 are the network types to compare

    Notes:
    ------
    - Valid network types: 'aw', 'pow', 'prm', 'com'
    - Comparison mode requires at least two network types to compare
    - If invalid network types are specified, defaults to ['aw', 'pow']
    """
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
        _vorpy_print_error(
            "Network Type",
            net_type,
            "one of: 'aw', 'pow', 'prm', or 'com'",
            example="aw"
        )
        return settings['net_type']
    # If we are comparing the network types
    if interpreter[net_type] == 'com':
        # Check to see if the set nets are available and at the very end add 'aw' and power so returned worst case
        set_nets = [interpreter[_] for _ in set_nets] + ['aw', 'pow']
        # Return the comparisons
        if print_change:
            _vorpy_print_setting("Network Type", "\"compare\"")
        return [interpreter[net_type], set_nets[0], set_nets[1]]
    # Return the interpreted network type
    if print_change:
        _vorpy_print_setting("Network Type", interpreter[net_type])
    return interpreter[net_type]


def set_sc(surface_color, settings, print_change=False):
    """
    Configures the surface color scheme for visualization.

    This function handles the setting of surface colors by:
    1. Attempting to validate the provided color map name against matplotlib's colormaps
    2. Supporting both new (mpl.colormaps) and legacy (mpl.cm) matplotlib colormap access methods
    3. Providing helpful error messages with examples of valid colormap names

    Parameters:
    -----------
    surface_color : str or list
        The desired colormap name. If a list is provided, uses the first element.
    settings : dict
        Dictionary containing current system settings including the default surface color

    Returns:
    --------
    str
        - If valid: Returns the validated colormap name
        - If invalid: Returns the current surface color from settings

    Notes:
    ------
    - Common valid colormaps include: 'viridis', 'plasma', 'inferno', 'cividis'
    - Also supports basic color maps like 'Greys', 'Reds', 'Greens', 'Blues'
    - Returns the current setting if an invalid colormap is specified
    """
    # First extract the value from the list if it is in fact a list
    if type(surface_color) is list:
        surface_color = surface_color[0]
    # Try each of the two possible options for surface coloring
    try:
        my_cmap = mpl.colormaps.get_cmap(surface_color)
        if print_change:
            _vorpy_print_setting("Surface Color", surface_color)
        return surface_color
    except Exception:
        pass
    try:
        my_cmap = mpl.cm.get_cmap(surface_color)
        if print_change:
            _vorpy_print_setting("Surface Color", surface_color)
        return surface_color
    except Exception:
        pass
    # If none of the formatting options work print the error and return
    _vorpy_print_error(
        "Surface Color",
        surface_color,
        "a valid matplotlib colormap name",
        example="viridis, plasma, inferno, cividis, Greys, Reds, Greens, Blues, rainbow"
    )
    return settings['surf_col']


def set_ss(surf_scheme, settings, print_change=False):
    """
    Configures the surface coloring scheme for the system.

    This function handles the setting of surface coloring schemes by:
    1. Validating the provided scheme against predefined options
    2. Supporting multiple aliases for each scheme type
    3. Providing helpful error messages with examples of valid schemes

    Parameters:
    -----------
    surf_scheme : str or list
        The desired surface coloring scheme. If a list is provided, uses the first element.
    settings : dict
        Dictionary containing current system settings including the default surface scheme

    Returns:
    --------
    str
        - If valid: Returns the validated scheme name
        - If invalid: Returns the current surface scheme from settings

    Notes:
    ------
    - Valid schemes include:
      - 'gauss' or 'gaussian' for Gaussian curvature
      - 'dist' for distance-based coloring
      - 'mean' or 'curv' for mean curvature
      - 'ins_out' for inside/outside coloring
      - 'none' for no special coloring
      - 'int_mean_curv' for integrated mean curvature
      - 'int_mean_curv_sq' for integrated squared mean curvature
      - 'int_gauss_curv' for integrated Gaussian curvature
    """
    # Make sure to extract the surface scheme from the value
    if type(surf_scheme) is list:
        surf_scheme = surf_scheme[0]
    # Set up the list of different dictionaries
    all_dicts = [
        {_: 'gauss' for _ in surf_scheme_gaus_vals},
        {_: 'dist' for _ in surf_scheme_dist_vals},
        {_: 'mean' for _ in surf_scheme_mean_vals + surf_scheme_curv_vals},
        {_: 'ins_out' for _ in surf_scheme_nout_vals},
        {_: 'none' for _ in nones},

        # Integrated curvature coloring
        {_: 'int_mean_curv' for _ in surf_scheme_int_mean_vals},
        {_: 'int_mean_curv_sq' for _ in surf_scheme_int_mean_sq_vals},
        {_: 'int_gauss_curv' for _ in surf_scheme_int_gauss_vals},
        {_: 'surf_energy' for _ in surf_scheme_energy_vals},
    ]
    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}
    # Check that the scheme entered is in the set of
    if surf_scheme not in interpreter:
        # Use _vorpy_print_error for invalid entry
        _vorpy_print_error(
            "Surface Scheme",
            surf_scheme,
            "a valid surface coloring scheme",
            example="curv, mean, gaussian, int_mean_curv, int_mean_curv_sq, int_gauss_curv, dist, ins_out, none"
        )
        return settings['surf_scheme']
    if print_change:
        _vorpy_print_setting("Surface Scheme", interpreter[surf_scheme])
    return interpreter[surf_scheme]


def set_sf(surface_factor, settings, print_change=False):
    """
    Configures the surface factor scaling for the system.

    This function handles the setting of surface factor scaling by:
    1. Validating the provided factor against predefined options
    2. Supporting multiple aliases for each factor type
    3. Providing helpful error messages with examples of valid factors

    Parameters:
    -----------
    surface_factor : list
        List containing the desired surface factor scaling type
    settings : dict
        Dictionary containing current system settings including the default surface factor

    Returns:
    --------
    str
        - If valid: Returns the validated factor name
        - If invalid: Returns the current surface factor from settings

    Notes:
    ------
    - Valid factors include:
      - 'lin' or 'linear' for linear scaling
      - 'log' or 'logarithmic' for logarithmic scaling
      - 'sqr' or 'square' for square root scaling
      - 'cub' or 'cube' for cubic root scaling
    """
    # Use _vorpy_print_setting and _vorpy_print_error for user feedback
    if surface_factor[0].lower() in surf_factor_vals:
        if print_change:
            _vorpy_print_setting("Surface Factor", surf_factor_vals[surface_factor[0].lower()])
        return surf_factor_vals[surface_factor[0].lower()]
    _vorpy_print_error(
        "Surface Factor",
        surface_factor[0],
        "a valid surface factor scaling",
        example="lin, linear, log, logarithmic, sqr, square, cub, cube"
    )
    return settings['surf_factor']


def set_rt(round_to, settings, print_change=False):
    if isinstance(round_to, list):
        round_to = round_to[0]

    round_to = int(round_to)

    if print_change:
        _vorpy_print_setting("Round To", round_to)

    return settings['round_to']


def set_ft(file_type, settings, print_change=False):
    """Validate and return the canonical geometry export format."""
    if isinstance(file_type, list):
        if len(file_type) == 0:
            return settings.get('file_type', 'off')
        file_type = file_type[0]

    value = str(file_type).strip().lower().lstrip('.')
    value = mesh_format_vals.get(value)
    if value is None:
        _vorpy_print_error(
            "Geometry File Format",
            file_type,
            "off, ply, or vtp",
            example="-s ft ply",
        )
        return settings.get('file_type', 'off')

    if print_change:
        _vorpy_print_setting("Geometry File Format", value.upper())
    return value


def set_ar(element_radius, settings, print_change=False):
    """
    Configures the atomic radii settings for the system.

    This function handles the setting of atomic radii by:
    1. Supporting two types of radius modifications:
       - Element-specific radii (e.g., 'C 1.7')
       - Residue-specific atom radii (e.g., 'ALA CA 1.7')
    2. Validating input against predefined element and residue names
    3. Supporting special cases for residue-specific atoms
    4. Providing helpful error messages for invalid inputs

    Parameters:
    -----------
    element_radius : list
        List containing the element/residue name and desired radius value
    settings : dict
        Dictionary containing current system settings including default atomic radii

    Returns:
    --------
    dict or None
        - If valid: Returns updated radius settings dictionary
        - If invalid: Returns None and prints error message

    Notes:
    ------
    - For element-specific changes, use format: 'element radius'
    - For residue-specific changes, use format: 'residue atom_name radius'
    - Radii must be valid float values
    - Element and residue names must match predefined values
    """
    print(element_radius)
    # Create the changes list
    change_settings = {'element': {}, 'special': {}}
    if settings['atom_rad'] is not None:
        change_settings = settings['atom_rad']

    # Separate the element from the radius
    if len(element_radius) >= 3 and element_radius[0] not in atom_objs:
        print(1)
        # Get the residue
        residue, name, radius = element_radius[:3]
        # Check that this exists
        if residue.lower() in residue_names and residue_names[residue.lower()] in special_radii:
            if name.upper() in special_radii[residue_names[residue.lower()]]:
                try:
                    my_radius = float(radius)
                    if print_change:
                        _vorpy_print_setting(
                            "Atomic Radius",
                            f"All {residue} {name} atoms radii changed from {special_radii[residue.upper()][name.upper()]} \u212B to {radius} \u212B"
                        )
                    # Add the radius
                    if residue_names[residue.lower()] in change_settings['special']:
                        change_settings['special'][residue_names[residue.lower()]][name.upper()] = my_radius
                        return change_settings
                    change_settings['special'][residue_names[residue.lower()]] = {name.upper(): my_radius}
                    return change_settings
                except ValueError:
                    _vorpy_print_error(
                        "Atomic Radius",
                        radius,
                        "a valid float value for radius",
                        example="1.7"
                    )
                    return
            _vorpy_print_error(
                "Atomic Radius",
                name,
                f"an atom in {residue}",
                example=", ".join([_ for _ in special_radii[residue.upper()]]),
            )
            return
        _vorpy_print_error(
            "Atomic Radius",
            " ".join(element_radius),
            "a valid residue/atom/radius entry",
            example="ALA CA 1.7"
        )
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
                if print_change:
                    _vorpy_print_setting(
                        "Atomic Radius",
                        f"All {element_names[element_radius[0].lower()]} atoms radii changed from {element_radii[element_names[element_radius[0].lower()]]} \u212B to {element_radius[1]} \u212B."
                    )
                change_settings['element'][element_radius[0].upper()] = my_radius
                return change_settings
            except ValueError:
                # Print the error saying the value is wrong
                _vorpy_print_error(
                    "Atomic Radius",
                    element_radius[1],
                    "a valid float value for radius",
                    example="1.7"
                )
                return
        # Check special radii for specific changes (e.g. all alpha carbons)
        elif any([element_radius[0].upper() in special_radii[_] for _ in special_radii]):
            try:
                my_radius = float(element_radius[1])
            except ValueError:
                _vorpy_print_error(
                    "Atomic Radius",
                    element_radius[1],
                    "a valid float value for radius",
                    example="1.7"
                )
                return
            # Loop through the special radii
            for residue in special_radii:
                if element_radius[0].upper() in residue:
                    if residue in change_settings['special']:
                        change_settings['special'][residue][element_radius[0].upper()] = my_radius
                    else:
                        change_settings['special'][residue] = {element_radius[0].upper(): my_radius}
            if print_change:
                _vorpy_print_setting(
                    "Atomic Radius",
                    f"All {element_radius[0]} radii changed to {element_radius[1]}"
                )
            return change_settings


def set_bt(build_type, settings, print_change=False):

    if isinstance(build_type, list):
        if len(build_type) == 0:
            return settings["bld_type"]

        build_type = build_type[0]

    build_type = str(build_type).strip().lower()

    if build_type in {"logs", "log"}:
        if print_change:
            _vorpy_print_setting("Build Type", "logs")

        return "logs"

    if build_type in {"normal", "build", "default", "fresh"}:
        if print_change:
            _vorpy_print_setting("Build Type", "normal")

        return "normal"

    _vorpy_print_error(
        "Build Type",
        build_type,
        "a valid build type",
        example="logs"
    )

    return settings["bld_type"]


def set_cc(conc_col, settings, print_change=False):
    """
    Updates the concave colors setting that will allow shells to show whether the surface comes inward or outward


    """
    if conc_col.lower in trues or conc_col:
        settings['conc_col'] = True
        if print_change:
            _vorpy_print_setting("Concave Colors", 'set to "True"')
        return True
    else:
        settings['conc_col'] = False
        if print_change:
            _vorpy_print_setting("Concave Colors", 'set to "False"')
        return False


def set_vc(vert_col, settings):
    """
    Configures the vertex color settings for the system.

    This function handles the setting of vertex colors by:
    1. Validating the provided color against predefined options
    2. Supporting multiple aliases for each color type
    3. Providing helpful error messages with examples of valid colors

    Parameters:
    -----------
    vert_col : str or list
        The desired vertex color. If a list is provided, uses the first element.
    settings : dict
        Dictionary containing current system settings including the default vertex color
    """

    return settings['vert_col']


def set_ec(edge_col, settings):
    """
    Configures the edge color settings for the system.

    This function handles the setting of edge colors by:
    1. Validating the provided color against predefined options
    2. Supporting multiple aliases for each color type
    3. Providing helpful error messages with examples of valid colors

    Parameters:
    -----------
    edge_col : str or list
        The desired edge color. If a list is provided, uses the first element.
    settings : dict
        Dictionary containing current system settings including the default edge color

    """
    return settings['edge_col']


def sett(setting, value, settings=None):
    """
    Updates system settings based on user input parameters.

    This function processes setting changes by:
    1. Initializing default settings if none provided
    2. Mapping user input to appropriate setting functions
    3. Interpreting various input formats for each setting type
    4. Applying the changes through specialized setting functions
    5. Returning the updated settings dictionary

    Parameters:
    -----------
    setting : str
        The setting name or alias to be modified
    value : str or list
        The new value(s) for the setting
    settings : dict, optional
        Current settings dictionary. If None, initializes with defaults.

    Returns:
    --------
    dict
        Updated settings dictionary with the new configuration
    """
    # Set the default settings
    if settings is None:
        settings = {'surf_res': 0.2, 'max_vert': 40, 'box_size': 1.25, 'net_type': 'aw', 'surf_col': 'plasma',
                    'surf_scheme': 'mean', 'scheme_factor': 'log', 'atom_rad': None, 'bld_type': None, 'conc_col': True,
                    'vert_col': 'red', 'edge_col': 'grey', 'round_to': 6, 'file_type': 'off'}
    # Set up the functions dictionary to return the value
    func_dict = {'surf_res': set_sr, 'max_vert': set_mv, 'box_size': set_bs, 'net_type': set_nt, 'surf_col': set_sc,
                 'surf_scheme': set_ss, 'scheme_factor': set_sf, 'atom_rad': set_ar, 'bld_type': set_bt,
                 'conc_col': set_cc, 'vert_col': set_vc, 'edge_col': set_ec, 'round_to': set_rt,
                 'file_type': set_ft}

    build_types = ["bt", "build_type", "bld_type", "build", "bld"]

    all_dicts = [{_: 'surf_res' for _ in surf_reses}, {_: 'max_vert' for _ in max_verts},
                 {_: 'box_size' for _ in box_sizes}, {_: 'net_type' for _ in net_types},
                 {_: 'surf_col' for _ in surf_colors}, {_: 'surf_scheme' for _ in surf_schemes},
                 {_: 'scheme_factor' for _ in surf_factors}, {_: 'atom_rad' for _ in atom_radii},
                 {_: 'bld_type' for _ in build_types},
                 {_: 'conc_col' for _ in conc_cols}, {_: 'vert_col' for _ in vert_cols},
                 {_: 'edge_col' for _ in edge_cols}, {_: 'round_to' for _ in round_tos},
                 {_: 'file_type' for _ in mesh_formats}]

    # Put all interpretations into one dictionary for convenience
    interpreter = {k: v for d in all_dicts for k, v in d.items()}

    # Set the setting
    settings[interpreter[setting]] = func_dict[interpreter[setting]](value, settings)
    # Return the settings
    return settings
