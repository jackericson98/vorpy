from Visualize.commands.load import load
from Visualize.commands.set import sett
from Visualize.commands.group import group
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

def load_argv(my_sys, input):
    pass

def set_argv(my_sys, input):
    # Check that the object is real
    pass


def group_argv(my_sys, input):
    pass

def export_argv(my_sys):
    pass



def interpret_argvs(my_sys):
    # Load the atom file
    load(my_sys, ["", sys.argv[1]])
    if my_sys.net is None:
        my_sys.net = Network(my_sys, my_sys.atoms)
    # Separate the rest of the argv args
    my_args = sys.argv[2:]
    my_group = None
    while my_args:
        arg = my_args.pop(0)
        if arg.lower() == '-s' and len(my_args) >= 2:
            setting = my_args.pop(0)
            value = my_args.pop(0)
            sett(my_sys, [setting, value], vorpy2_set=False)
        elif arg.lower() == '-g' and len(my_args) >= 1:
            grouping = my_args.pop(0)
            if len(my_args) >= 1 and my_args[0][0] != '-':
                my_ndx = my_args.pop(0)
                my_group = group(my_sys, ["", grouping, my_ndx])
    if my_group is None:
        my_group = Group(sys=my_sys, mols=my_sys.mols[:-1], name=my_sys.name + "_no_SOL")
    print(
        u"{} Build Settings - surf_res = {:.2f} \u208B,  max_vert  = {:.2f} \u208B,  box_multi = {:.2f} x,  build_surfs = {}, "
        u"surf_type = {}".format(my_group.name, my_sys.net.surf_res, my_sys.net.max_vert, my_sys.net.box_size,
                                 my_sys.net.build_surfs, 'voronoi' if my_sys.net.type == 'vor' else 'flat'))

    my_sys.net.build(my_group=my_group, build_surfs=True, output=True, print_actions=True)


