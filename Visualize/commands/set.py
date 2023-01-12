from Visualize.commands.interpret import *
from System.Network.network import Network


def sett(sys, usr_npt):
    """
    Set the network parameters
    :param sys:
    :param usr_npt:
    :return:
    """
    if len(usr_npt) == 1:
        my_set = get_set()
        if my_set is None:
            return
        my_val = get_val(my_set)
    elif len(usr_npt) == 2:
        my_set = get_set(usr_npt[1])
        if my_set is None:
            return
        my_val = get_val(my_set)
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
    elif my_set in sol_vertses:
        sys.net.sol_verts = my_val
        print("Calculate solute vertices set to {}".format(sys.net.sol_verts))
    elif my_set.lower() in quits:
        return
    else:
        invalid_input(usr_npt)
