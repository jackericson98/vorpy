from Visualize.commands.commands import *
from System.Network.network import Network


def build(sys):
    """
    Prints a pre-built header and asks the user if they are ready to build. Once confirmed prints a building header and
    builds
    :return:
    """
    # If no system has been loaded tell the user to screw off
    if len(sys.atoms) == 0:
        print("No atoms in the system. Use the \'load\' command or type \'h\' for help")
        return
    # Check to see if a network has been added
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Once the build command is used, the user is greeted with the build settings and asked if they are ready to build
    print(u"Settings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_size = {:.2f} x,  sol_verts = {}"
          .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.sol_verts))
    # The user is prompted to start the build - This could say eta and other build qualities
    pre_build_confirmation = input("confirm >>>   ")
    # If the user is ready to build, build the system
    if pre_build_confirmation.lower() in ys:
        sys.build_network(flat_faces=sys.net.flat_faces)
    elif pre_build_confirmation.lower() in ns:
        print("Use the \'set\' command to change a setting and a value. Type \'h\' for help")
        return
    elif pre_build_confirmation.lower() in helps:
        help_()
    elif pre_build_confirmation.lower() in quits:
        return
