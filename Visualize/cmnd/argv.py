import sys
import os
from Visualize.cmnd.load import load
from Visualize.cmnd.set2 import sett
from Visualize.cmnd.group import ggroup
from System.system import System
from System.Group.group import Group
from System.sys_funcs.output.output import export_min1, export_min2, export_med, export_large, export_all, other_exports, set_sys_dir
from Visualize.cmnd.commands import ands, helps, print_help
from copy import deepcopy


"""
Argv rules: 
1. Space delimited
2. flags (-l: load, -s: set, -g: group, -e: export)
3. For multiple inputs use &&
4. Defaults to no sol, default settings, export all
"""


def argv_export(my_sys, usr_npt, add_on=None):
    """
    Exports the specified elements from the system.

    Specifications:

        'default' - (optional) exports network, information, no Sol shell, pdb, set atoms from the system and atoms,
                    surrounding atoms, shell, info from each group
        'all' - exports everything from the full system and each group, needs a warning first
        'info' - exports only the information files for the system and each group
        'surfs' - exports the built surfaces individually from the system and all group surfaces + verts and edges

    :param my_sys:
    :param usr_npt:
    :return:
    """
    # Go through each of the inputs in the exports
    if len(usr_npt) == 0:
        usr_npt.append(['default'])
    # Export the specified exports
    for npt in usr_npt:
        if npt[0].lower() in {'dir', 'directory'} and len(npt) == 2 and (os.path.isdir(npt[1]) or npt[1] == 'gsu_logs'):
            if npt[1] == 'gsu_logs':
                set_sys_dir(my_sys, "C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/{}/logs".format(my_sys.name))
            else:
                if add_on is None:
                    my_sys.dir = npt[1]
                else:
                    my_sys.dir = npt[1] + add_on
    for npt in usr_npt:
        if npt[0].lower() in {'dir', 'directory'}:
            continue
        export_npt(my_sys, npt[0])


def export_npt(my_sys, usr_npt=None):
    """
    Exports the selection from the usr_npts
    :param my_sys:
    :param usr_npt:
    :return:
    """
    # If nothing is specified export the defaults
    if usr_npt is None or usr_npt.lower() in {'default', '2', 'medium', '', 'med'}:
        export_large(sys=my_sys)

    # Small export
    elif usr_npt.lower() in {"tiny", "i", "info", "0", "smallest"}:
        export_min1(my_sys)

    # Medium small export
    elif usr_npt.lower() in {"small", "s", "1"}:
        export_min2(my_sys)

    # Large Export
    elif usr_npt.lower() in {"large", "l", "3"}:
        export_large(my_sys)

    # Export all
    elif usr_npt.lower() in {'all', 'a', 'everything'}:
        export_all(my_sys)

    else:
        other_exports(my_sys, usr_npt)


def interpret_argvs():
    # Separate the rest of the argv args
    my_args = sys.argv[2:]
    # Set up the commands dictionary
    cmnds = {'npt': [], 'set': [], 'grp': {}, 'bld': [], 'xpt': [], 'ifc': []}
    # Set the arg to load as a default
    arg = '-l'
    group_counter = -1
    # Go through the arguments
    while my_args:
        # Remove the first argument flag
        if my_args[0] in ands:
            my_args.pop(0)
        else:
            arg = my_args.pop(0)
            if arg == '-g':
                group_counter += 1
                cmnds['grp'][group_counter] = []
        # Gather the cmnd and the flag
        arg_cmnds = []
        while True:
            if len(my_args) == 0 or my_args[0][0] == '-' or my_args[0] in ands:
                break
            else:
                # Keep gathering the commands for the flag
                arg_cmnds.append(my_args.pop(0))
        # Add the command to the 181L list
        if arg.lower() == '-l':
            cmnds['npt'].append(arg_cmnds)
        elif arg.lower() == '-s':
            cmnds['set'].append(arg_cmnds)
        elif arg.lower() == '-g':
            cmnds['grp'][group_counter].append(arg_cmnds)
        elif arg.lower() == '-b':
            cmnds['bld'].append(arg_cmnds)
        elif arg.lower() == '-e':
            if arg_cmnds == 'logs':
                cmnds['set'].append(['bt', 'logs'])
            cmnds['xpt'].append(arg_cmnds)
        elif arg.lower() == '-i':
            cmnds['ifc'].append(arg_cmnds)
    # Return the lists
    return cmnds


def argv(my_sys):
    # First check if the argv value is in helps
    if sys.argv[1].lower() in helps:
        print_help()
        return
    # Load the atom file
    load(my_sys, [["", sys.argv[1]]])

    # Interpret the commands
    cmnds = interpret_argvs()
    # Set the system commands
    my_sys.cmnds = cmnds
    # Go through each of the ls
    load(my_sys, cmnds['npt'])
    for commandaroonski in cmnds['xpt']:
        if commandaroonski[0] == 'dir' and os.path.isdir(commandaroonski[1]):
            my_sys.files['dir'] = commandaroonski[1]
            cmnds['xpt'].pop(cmnds['xpt'].index(commandaroonski))
    # Declare the settings variable
    settings = None
    # Go through the user inputs loading files
    for my_set in cmnds['set']:
        # Alter the settings
        settings = sett(my_set[0], my_set[1:], settings)
    # Update the sphere radii in the system
    if settings is not None and settings['atom_rad'] is not None:
        my_sys.set_radii(settings['atom_rad']['element'], settings['atom_rad']['special'])

    # compare the groups
    if my_sys.groups is None or len(my_sys.groups) == 0:
        ggroup(my_sys, cmnds['grp'], settings)
    else:
        verts = my_sys.groups[0].verts
        ggroup(my_sys, cmnds['grp'], settings)
        my_sys.groups[0].verts = verts
    if my_sys.groups is None or len(my_sys.groups) == 0:
        print('{} not a valid group command. Calculating whole molecule'.format(cmnds['grp']))
        ggroup(my_sys, [['ns']])
    # Build the groups
    if settings is not None and len(settings['net_type']) > 1 and settings['net_type'][0] == 'com':
        new_groups = []
        for grp in my_sys.groups:
            copy_group = deepcopy(grp)
            copy_group.name = copy_group.name + '_' + settings['net_type'][1]
            copy_group.settings['net_type'] = settings['net_type'][1]
            grp.settings['net_type'] = settings['net_type'][2]
            grp.name = grp.name + '_' + settings['net_type'][2]
            # Delete the vertices from the grp group because they are only aw
            grp.verts = None
            copy_group.build()
            grp.build()
            new_groups.append(copy_group)
            # Compare the two networks
            my_sys.compare_networks(group1=copy_group, group2=grp)
        my_sys.groups += new_groups

    else:
        for grp in my_sys.groups:
            grp.build()
    # Make the system's interfaces
    my_sys.make_interfaces()
    # Export everything
    argv_export(my_sys, cmnds['xpt'])


if __name__ == '__main__':
    mySys = System()
    argv(mySys)
