from Visualize.argv.argv_load import argv_load, argv_load_atoms
from Visualize.argv.argv_set import argv_sett
from Visualize.argv.argv_group import argv_group
from Visualize.argv.argv_build import argv_build
from Visualize.argv.argv_export import argv_export
from System.Network.network import Network
from System.system import System
import sys

"""
Argv rules: 
1. Space delimited
2. flags (-l: load, -s: set, -g: group, -e: export)
3. For multiple inputs use &&
4. Defaults to no sol, default settings, export all
"""


def interpret_argvs(my_sys):
    # Separate the rest of the argv args
    my_args = sys.argv[2:]
    npt_cmnds, set_cmnds, grp_cmnds, bld_cmnds, xpt_cmnds, ifc_cmnds = [], [], [], [], [], []
    # Go through the arguments
    while my_args:
        # Remove the first argument flag
        arg = my_args.pop(0)
        # Gather the cmnd and the flag
        arg_cmnds = []
        while True:
            if len(my_args) == 0 or my_args[0][0] == '-':
                break
            else:
                # Keep gathering the commands for the flag
                arg_cmnds.append(my_args.pop(0))
        # Add the command to the correct list
        if arg.lower() == '-l':
            npt_cmnds.append(arg_cmnds)
        elif arg.lower() == '-s':
            set_cmnds.append(arg_cmnds)
        elif arg.lower() == '-g':
            grp_cmnds.append(arg_cmnds)
        elif arg.lower() == '-b':
            bld_cmnds.append(arg_cmnds)
        elif arg.lower() == '-e':
            xpt_cmnds.append(arg_cmnds)
        elif arg.lower() == '-i':
            ifc_cmnds.append(arg_cmnds)
    # Return the lists
    return npt_cmnds, set_cmnds, grp_cmnds, bld_cmnds, xpt_cmnds, ifc_cmnds


def argv(my_sys):
    # Load the atom file
    argv_load_atoms(my_sys, ["", sys.argv[1]])
    # Interpret the commands
    files, settings, groups, builds, exports, ifaces = interpret_argvs(my_sys)
    # Go through each of the ls
    argv_load(my_sys, files)

    argv_sett(my_sys, settings)
    argv_group(my_sys, groups, bff=ifaces)
    argv_build(my_sys, builds)
    argv_export(my_sys, exports)


if __name__ == '__main__':
    mySys = System()
    argv(mySys)
