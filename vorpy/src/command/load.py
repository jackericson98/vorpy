from vorpy.src.system import System
from vorpy.src.command.commands import *
from vorpy.src.command.interpret import get_file


def load(sys, usr_npt, balls_file=False):
    """
    Loads molecular structure files and associated data into the system.

    This function handles the loading of various file types:
    - Molecular structure files (.pdb, .mol, .gro, .cif)
    - Vertex files (.txt with 'verts' or 'vertices' in name)
    - Network files (.txt with 'net' in name)

    The function provides interactive confirmation prompts when:
    - Replacing an existing system
    - Replacing existing vertex files
    - Replacing existing network files

    Parameters:
    -----------
    sys : System
        The system object to load data into
    usr_npt : list
        List of file specifications to load
    balls_file : bool, optional
        Flag indicating if the file should be treated as a molecular structure file
        regardless of extension. Default is False.

    Returns:
    --------
    System or None
        Returns the updated system object if successful, None if loading is cancelled
    """
    # Initialize a list to store the files to load
    my_files = []
    # Process each file specification in the user input
    for npt in usr_npt:
        # If there's only one file specification, get it from the user
        if len(npt) == 1:
            my_files.append(get_file())
            # If the user quits or cancels the file selection, return None
            if my_files[-1] is None or my_files[-1].lower() in quits:
                return
        # If there are multiple file specifications, process each one
        else:
            for file in npt[1::2]:
                # Get the file from the user
                my_file = get_file(file)
                # If the user quits or cancels the file selection, return None
                if my_file is None or my_file.lower() in quits:
                    return
                # Add the file to the list
                my_files.append(my_file)
    # Process each file in the list
    for file in my_files:
        # Check to see what type of file it is
        if file[-3:] == 'pdb' or file[-3:] == 'mol' or file[-3:] == 'gro' or file[-3:] == 'cif' or balls_file:
            # If the system already exists, prompt the user to confirm replacement
            if sys.name is not None and \
                    (sys.atoms is not None or sys.files['verts_file'] is not None or sys.files['net_file'] is not None):
                reset_sys = input("replacing {} with {}\nconfirm >>>   "
                                  .format(sys.name, file))
                # If the user confirms the replacement, create a new system
                if reset_sys.lower() in ys:
                    sys = System(file)
                    print(sys.name + " loaded - {} atoms, {} molecules, solute: {}"
                          .format(len(sys.atoms), len(sys.chains), sys.sol.name))
                    return sys
                # If the user requests help, print the help message
                elif reset_sys.lower() in helps:
                    print_help()
                # If the user quits, return None
                elif reset_sys.lower() in quits:
                    return
            # If the system does not exist, load the new system
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
                if sys.files['verts_file'] is not None and sys.vert_file != "":
                    replace_vert_file = input("replacing {} with {}\n "
                                              "confirm >>>   ".format(sys.files['verts_file'], file))
                    # If the user confirms the replacement, load the vertices
                    if replace_vert_file.lower() in ys or replace_vert_file.lower() in dones:
                        sys.load_verts(file, vta_ball_file=sys.ball_file)
                        print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.settings['max_vert'], sys.net.settings['box_size']))
                    # If the user requests help, print the help message
                    elif replace_vert_file.lower() in helps:
                        print_help()
                    # If the user quits, return None
                    elif replace_vert_file.lower() in quits:
                        return
                # If the vertex file has not been loaded, load it
                else:
                    sys.load_verts(file, vta_ball_file=sys.files['ball_file'])
                    # print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                    #       .format(sys.name, len(sys.net.vta_verts), sys.net.max_vert, sys.net.box_size))
            elif file[-9:-4].lower() == 'balls':
                sys.ball_file = file
            # If the new file is a network file load it
            elif file[-11:-4].lower() in 'network':
                # If a vertex file has already been loaded make sure the user wants to load it if not load it
                if sys.net_file is not None or sys.net_file != "":
                    replace_net_file = input("replacing {} with {}\n "
                                              "confirm >>>   ".format(sys.net_file, file))
                    # If the user confirms the replacement, load the network
                    if replace_net_file in ys:
                        sys.load_net(file)
                        print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box"
                              " size: {} x\n".format(sys.name, len(sys.net.verts), sys.net.settings['max_vert'],
                                                     sys.net.settings['box_size']))
                    # If the user requests help, print the help message
                    elif replace_net_file in helps:
                        print_help()
                    # If the user quits, return None
                    else:
                        return
                else:
                    # Load the file
                    sys.load_net(file)
                    if len(sys.net.surfs) > 0:
                        print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.settings['max_vert'], sys.net.settings['box_size']))
                    else:
                        print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                              .format(sys.name, len(sys.net.verts), sys.net.settings['max_vert'], sys.net.settings['box_size']))
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
