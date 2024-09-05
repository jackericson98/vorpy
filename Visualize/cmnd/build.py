from Visualize.cmnd.commands import *
from Visualize.cmnd.set import sett
from Visualize.cmnd.group import group


def build(sys, usr_npt=None):
    """
    Prints a pre-built header and asks the user if they are ready to build. Once confirmed prints a building header and
    builds
    :return:
    """
    # If no system has been loaded tell the user to screw off
    if len(sys.atoms) == 0:
        print("no atoms in the system. use the \'load\' command or type \'h\' for help")
        return
    # Check to see if a group was specified
    my_group = None
    if usr_npt is not None:
        my_group = group(sys, usr_npt)
    # Once the build command is used, the user is greeted with the build settings and asked if they are ready to build
    print(u"settings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_size = {:.2f} x,  build_surfs = {}, net_type = {}"
          .format(sys.net.settings['surf_res'], sys.net.settings['max_vert'], sys.net.settings['box_size'], sys.net.build_surfs, sys.net.settings['net_type']))
    # The user is prompted to start the build - This could say eta and other build qualities
    pre_build_confirmation = input("confirm >>>   ")
    # If the user is ready to build, build the system
    if pre_build_confirmation.lower() in ys:
        sys.net.build(calc_verts=not build_vta, my_group=my_group)
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
        print_help()
    elif pre_build_confirmation.lower() in quits:
        return
