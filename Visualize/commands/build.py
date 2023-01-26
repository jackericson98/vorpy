from Visualize.commands.commands import *
from Visualize.commands.set import sett
from System.Network.network import Network


def build(sys):
    """
    Prints a pre-built header and asks the user if they are ready to build. Once confirmed prints a building header and
    builds
    :return:
    """
    # If no system has been loaded tell the user to screw off
    if len(sys.atoms) == 0:
        print("no atoms in the system. use the \'load\' command or type \'h\' for help")
        return
    # Check to see if a network has been added
    if sys.net is None:
        sys.net = Network(sys=sys, atoms=sys.atoms)
    # Check to see if the network has been built before
    if sys.net.verts is not None:
        rebuild_npt = input("{} network already constructed. would you like to rebuild?\nconfirm >>>   ".format(sys.name))
        if rebuild_npt.lower() in ys:
            sys.net = Network(sys=sys, atoms=sys.atom, surf_res=sys.net.surf_res, max_vert=sys.net.max_vert,
                              box_size=sys.net.box_size, build_surfs=sys.net.build_surfs, flat_surfs=sys.net.flat_surfs)
        else:
            return
    # Once the build command is used, the user is greeted with the build settings and asked if they are ready to build
    print(u"settings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_size = {:.2f} x,  build_surfs = {}, flat_surfs = {}"
          .format(sys.net.surf_res, sys.net.max_vert, sys.net.box_size, sys.net.build_surfs, sys.net.flat_surfs))
    # The user is prompted to start the build - This could say eta and other build qualities
    pre_build_confirmation = input("confirm >>>   ")
    # If the user is ready to build, build the system
    if pre_build_confirmation.lower() in ys:
        sys.net.build()
    elif pre_build_confirmation.lower() in ns:
        # Ask the user if they would like to change the settings
        chng_stngs_npt = input("change settings?\nconfirm >>>   ")
        chng_stngs_npt_lst = chng_stngs_npt.split()
        if chng_stngs_npt_lst[0].lower() in ys + set_cmds:
            sett(sys=sys, usr_npt=chng_stngs_npt)
        else:
            print("use the \'set\' command to change a setting and a value or type \'h\' for help")
            return
    elif pre_build_confirmation.lower() in helps:
        help_()
    elif pre_build_confirmation.lower() in quits:
        return
