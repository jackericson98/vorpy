from System.Network.network import Network


def argv_build(my_sys, usr_npt):
    # Check for a network
    if my_sys.net is None:
        my_sys.net = Network(sys=my_sys, atoms=my_sys.atoms)

    for my_group in my_sys.groups:
        if len(my_group.atoms) > 0:
            my_sys.build_network(my_group=my_group)
