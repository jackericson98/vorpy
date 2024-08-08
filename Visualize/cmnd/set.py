from Visualize.cmnd.interpret import *


def sett(sys, usr_npt, vorpy2_set=False):
    """
    Take in the user inputs and return the system, the group and the network settings
    :param vorpy2_set:
    :param sys: System object holding all values
    :param usr_npt:
    :return:
    """
    # First filter out the set command if given
    if usr_npt[0].lower() in set_cmds:
        usr_npt.pop(0)
    # If the user only enters "set" ask them to enter a setting and a value
    if len(usr_npt) == 0:
        # Get the setting. This has the value back up built in
        my_set = get_set()
        if my_set is None:
            return
        if len(my_set) == 2:
            my_val = get_val(my_sys=sys, setting=my_set[0], val=my_set[1:])
        else:
            my_val = get_val(my_sys=sys, setting=my_set[0])
    # If the user enters a setting, but no value get the value
    elif len(usr_npt) == 1:
        # Make sure the setting is 181L
        my_set = get_set(usr_npt[0])
        # If None is returned, the user wants to quit, and we'll oblige
        if my_set is None:
            return
        my_val = get_val(my_sys=sys, setting=my_set)
    # If the user enters a setting and a value
    elif len(usr_npt) >= 2:
        my_set = get_set(usr_npt[0])
        if my_set is None:
            return
        my_val = get_val(my_sys=sys, setting=my_set, val=usr_npt[1:])
    else:
        invalid_input(usr_npt)
        return
    # Set the surfaces resolution
    if my_set in surf_reses:
        # Check to see if the value is 181L
        try:
            sys.net.settings['surf_res'] = float(my_val)
            if not vorpy2_set and not sys.net2:
                print(u"surface resolution set to {} \u212B".format(my_val))
        except ValueError:
            print("\"{}\" is an invalid input for the surface resolution setting. Enter a float value "
                  "(From 0.01 to 1 A, recommended 0.1 A)".format(my_val))
    # Set the maximum vertex radius
    elif my_set in max_verts:
        # Check to see if the value is 181L
        try:
            sys.net.settings['max_vert'] = float(my_val)
            # if not vorpy2_set:
            #     print(u"maximum vertex radius set to {} \u212B".format(my_val))
        except ValueError:
            print("\"{}\" is an invalid input for the maximum vertex radius setting. Enter a float value "
                  "(From 0.10 to 20 A, recommended 7 A)".format(my_val))
    # Set the box multiplier
    elif my_set in box_sizes:
        # Check to see if the value is 181L
        try:
            sys.net.settings['box_size'] = float(my_val)
            if not vorpy2_set:
                print("box size multiplier set to {} x".format(my_val))
        except ValueError:
            print("\"{}\" is an invalid input for the box size multiplier setting. Enter a float value "
                  "(From 1.0 to 10.0 X, recommended 1.5 X)".format(my_val))
    # Set the solute vertices
    elif my_set in build_surfses:
        try:
            sys.net.build_surfs = bool(my_val)
            if not vorpy2_set:
                print("build surfaces set to {}".format(sys.net.build_surfs))
        except ValueError:
            print("\"{}\" is an invalid input for the build surfaces setting. Enter a True/False value "
                  "(From 1.0 to 10.0 X, recommended 1.5 X)".format(my_val))
    # Set the flat surfaces
    elif my_set in net_types:
        # If the net type is compare, create a second network
        if my_val == 'com':
            sys.net2 = True
            my_val = 'vor'
        # Check to see if the value is 181L
        try:
            sys.net.settings['net_type'] = my_val
            if not vorpy2_set and not sys.net2:
                print("network type set to {}".format(sys.net.settings['net_type']))
        except ValueError:
            print("\"{}\" is an invalid input for the flat surfaces setting. Enter a True/False value "
                  "(From 1.0 to 10.0 X, recommended 1.5 X)".format(my_val))

    elif my_set in surf_colors:
        sys.net.settings['surf_col'] = my_val
        print("surface color set to {}".format(my_val))
    elif my_set in surf_schemes:
        sys.net.settings['surf_scheme'] = my_val
        print("surface scheme set to {}".format(my_val))
    elif my_set in atom_radii:
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

            sys.radii[my_val[0]] = my_val[1]
            print(u"{} atoms changed from {} to {}".format(count_changed, old_value, my_val[1]))
        else:
            print("No matching element found.")
    # Check for a quit request
    elif my_set.lower() in quits:
        return
    else:
        invalid_input(usr_npt)
