from Visualize.commands.load import load
from Visualize.commands.set import sett
from Visualize.commands.group import group
from Visualize.commands.commands import *
from Visualize.commands.export import export
from System.sys_objs.group import Group
from System.Network.network import Network
import sys

"""
Argv rules: 
1. Space delimited
2. flags (-l: load, -s: set, -g: group, -e: export)
3. For multiple inputs use &&
4. Defaults to no sol, default settings, export all
"""


def print_error():
    pass


def load_argv(my_sys, usr_npt):
    # Go through the user inputs loading files
    while usr_npt:
        # Pop the file descriptor
        descriptor = usr_npt.pop(0)
        # Check to see that it is a descriptor
        if descriptor.lower() not in file_types or len(usr_npt) == 0:
            return
        # Load the file
        load(my_sys, [descriptor, usr_npt.pop(0)])
        # If the next value is && go again
        if len(usr_npt) > 0 and usr_npt[0] == '&&':
            usr_npt.pop(0)


def set_argv(my_sys, usr_npt):
    # Go through the user inputs loading files
    while usr_npt:
        # Pop the file descriptor
        descriptor = usr_npt.pop(0)
        # Check to see that it is a descriptor
        if descriptor.lower() not in my_settings or len(usr_npt) == 0:
            return
        # Load the file
        sett(my_sys, [descriptor, usr_npt.pop(0)], vorpy2_set=False)
        # If the next value is && go again
        if len(usr_npt) > 0 and usr_npt[0] == '&&':
            usr_npt.pop(0)


def group_argv(my_sys, usr_npt):
    # Create a group variable
    my_group = None
    # Go through the user inputs loading files
    while usr_npt:
        # Pop the file descriptor
        descriptor = usr_npt.pop(0)
        # Check to see that it is a descriptor
        if descriptor.lower() not in my_objects:
            return
        elif descriptor.lower() == 'ns':
            my_group = Group(my_sys, mols=my_sys.mols[:-1], name=my_sys.name + "_no_SOL")
            continue
        elif descriptor.lower in full_objs:
            return Group(my_sys, atoms=my_sys.atoms, name=my_sys.name + "_full")
        # Load the file
        my_group = group(my_sys, [descriptor, usr_npt.pop(0)], my_group)
        # If the next value is && go again
        if len(usr_npt) > 0 and usr_npt[0] == '+':
            usr_npt.pop(0)
    # Return the group object
    return my_group


def interpret_argvs(my_sys):
    # Load the atom file
    load(my_sys, ["", sys.argv[1]])
    if my_sys.net is None:
        my_sys.net = Network(my_sys, my_sys.atoms)
    # Separate the rest of the argv args
    my_args = sys.argv[2:]
    my_group = None
    arg = my_args.pop(0)
    export_commands = []
    while my_args:
        # Gather the commands and the flag
        cmnd_args = [arg]
        while True:
            arg = my_args.pop(0)
            if arg[0] == '-' or len(my_args) == 0:
                break
            cmnd_args.append(arg)
        # Execute the commands
        if cmnd_args[0].lower() == '-s' and len(cmnd_args) >= 3:
            set_argv(my_sys, cmnd_args[1:])
        elif cmnd_args[0].lower() == '-g' and len(cmnd_args) >= 2:
            my_group = group_argv(my_sys, cmnd_args[1:])
        elif cmnd_args[0].lower() == '-e' and len(cmnd_args) >= 2:
            export_commands.append(cmnd_args)
    if my_group is None:
        my_group = Group(sys=my_sys, mols=my_sys.mols[:-1], name=my_sys.name + "_no_SOL")
    net = my_sys.net
    print(u"{} Build Settings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_multi = {:.2f} x,  "
          u"build_surfs = {}, surf_type = {}"
          .format(my_group.name, net.surf_res, net.max_vert, net.box_size, net.build_surfs,
                  {'vor': 'Voronoi', 'del': 'Delaunay', 'pow': 'Power'}[net.type]))

    my_sys.net.build(my_group=my_group, build_surfs=True, output=True, print_actions=True)
    export(my_sys, my_group=my_group, usr_npt=[])

