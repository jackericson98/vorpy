from Visualize.commands.interpret import *
from System.Network.network import Network


def sett(sys, usr_npt):
    """
    Set the network parameters
    :param sys: System object holding all values
    :param usr_npt:
    :return:
    """
    # If the user only enters "set" ask them to enter a setting and a value
    if len(usr_npt) == 1:
        # Get the setting. This has the value back up built in
        my_set = get_set()
        if my_set is None:
            return
        if len(my_set) == 2:
            my_val = get_val(setting=my_set[0], val=my_set[1])
        else:
            my_val = get_val(setting=my_set[0])
    # If the user enters a setting, but no value get the value
    elif len(usr_npt) == 2:
        # Make sure the setting is correct
        my_set = get_set(usr_npt[1])
        # If None is returned, the user wants to quit, and we'll oblige
        if my_set is None:
            return
        # If the user provided the value as well as the setting check the value
        if len(my_set) == 2:
            my_val = get_val(setting=my_set[0], val=my_set[1])
        # Otherwise, make them get the value
        else:
            my_val = get_val(setting=my_set[0])
    # If the user enters a setting and a value
    elif len(usr_npt) <= 3:
        my_set = get_set(usr_npt[1])
        if my_set is None:
            return
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
    elif my_set in calc_surfses:
        sys.net.sol_verts = my_val
        print("Calculate surfaces vertices set to {}".format(sys.net.sol_verts))
    elif my_set.lower() in quits:
        return
    else:
        invalid_input(usr_npt)
