from System.Network.network import Network


def argv_build(my_sys, usr_npt):
    # Check for a group network
    for my_group in my_sys.groups:
        if len(my_group.ball_ndxs) > 0:
            my_group.build_network()
