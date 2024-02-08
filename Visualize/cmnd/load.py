from Visualize.cmnd.interpret import get_file
from Visualize.cmnd.commands import *
from System.system import System


def load(sys, usr_npt):
    """
    Once one of the load commands is used try to load the rest of the string
    :param sys:
    :param usr_npt:
    :return:
    """

    my_files = []
    if len(usr_npt) == 1:
        my_files.append(get_file())
        if my_files[-1] is None or my_files[-1].lower() in quits:
            return
    else:
        for file in usr_npt[1::2]:
            my_file = get_file(file)
            if my_file is None or my_file.lower() in quits:
                return
            my_files.append(my_file)
    for file in my_files:
        # Check to see what type of file it is
        if file[-3:] == 'pdb' or file[-3:] == 'mol' or file[-3:] == 'gro' or file[-3:] == 'cif':
            if sys.name is not None and \
                    (sys.atoms is not None or sys.vert_file is not None or sys.net_file is not None):
                reset_sys = input("replacing {} with {}\nconfirm >>>   "
                                  .format(sys.name, file))
                if reset_sys.lower() in ys:
                    sys = System(file)
                    print(sys.name + " loaded - {} atoms, {} molecules, solute: {}"
                          .format(len(sys.atoms), len(sys.chains), sys.sol.name))
                    return sys
                elif reset_sys.lower() in helps:
                    help_()
                elif reset_sys.lower() in quits:
                    return
            else:
                sys.load_sys(file=file)
                # noinspection PyUnresolvedReferences
                sys.print_info()
                return sys
        # If the loaded file is a vertex or network file load them accordingly
        elif file[-3:] == 'txt':
            # If the new file is a vertex file load it
            if file[-9:-4].lower() == 'verts' or file[-12:-4].lower() == 'vertices':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.vert_file is not None and sys.vert_file != "":
                    replace_vert_file = input("replacing {} with {}\n "
                                              "confirm >>>   ".format(sys.vert_file, file))
                    if replace_vert_file.lower() in ys or replace_vert_file.lower() in dones:
                        sys.load_verts(file, vta_ball_file=sys.ball_file)
                        print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
                    elif replace_vert_file.lower() in helps:
                        help_()
                    elif replace_vert_file.lower() in quits:
                        return
                else:
                    sys.load_verts(file, vta_ball_file=sys.ball_file)
                    print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                          .format(sys.name, len(sys.net.vta_verts), sys.net.max_vert, sys.net.box_size))
            elif file[-9:-4].lower() == 'balls':
                sys.ball_file = file
            # If the new file is a network file load it
            elif file[-11:-4].lower() in 'network':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.net_file is not None or sys.net_file != "":
                    replace_net_file = input("replacing {} with {}\n "
                                              "confirm >>>   ".format(sys.net_file, file))
                    if replace_net_file in ys:
                        sys.load_net(file)
                        print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box"
                              " size: {} x\n".format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
                    elif replace_net_file in helps:
                        help_()
                    else:
                        return
                else:
                    # Load the file
                    sys.load_net(file)
                    if len(sys.net.surfs) > 0:
                        print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
                    else:
                        print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.max_vert, sys.net.box_size))
        # Check to see if it is a new network file
        elif file[-3:] == 'csv':
            # Check to see that this is a network file
            if file[-7:-4].lower() == 'net':

                sys.load_net(file=file)

        # If the file is an index file load it accordingly
        elif file[-3:] == 'ndx':
            sys.load_ndx(file)
            print(sys.ndx_file + "loaded -  {}".format(sys.ndx_names[:min(len(sys.ndx_names) - 1, 10)]))
        # In all other case print an error and give the user a chance to try again
        else:
            print("\'{}\' is not a valid input. allowed file types: .pdb, .mol, .cif, .gro, .txt, .ndx. type "
                  "\'h\' for help".format(file))
            return
