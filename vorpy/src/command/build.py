from vorpy.src.command.commands import *
from vorpy.src.command.set import sett
from vorpy.src.command.group import group


def build(sys, usr_npt=None, build_vts=False):
    """
    Initiates the surface building process for a molecular system.

    This function handles the build command workflow:
    1. Validates that atoms are present in the system
    2. Optionally processes group-specific build parameters
    3. Displays current build settings including:
       - Surface resolution
       - Maximum vertex distance
       - Box size
       - Surface building flags
       - Network type
    4. Prompts for user confirmation before proceeding
    5. Either:
       - Proceeds with building if confirmed
       - Allows settings modification if declined
       - Provides help information if requested
       - Exits if quit command received

    Parameters:
    -----------
    sys : System
        The system object containing atoms and build settings
    usr_npt : list, optional
        Optional group-specific build parameters

    Returns:
    --------
    None
    """
    # If no system has been loaded tell the user to screw off
    if len(sys.atoms) == 0:
        print("no atoms in the system. use the \'load\' command or type \'h\' for help")
        return
    # Check to see if a group was specified
    my_group = None
    # If a group was specified, create a new group
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
