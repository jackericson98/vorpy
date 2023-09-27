from Visualize.argv.argv_load import argv_load, argv_load_atoms, argv_load_foam
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


def interpret_argvs():
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
    if sys.argv[1].lower() == 'foam':
        argv_load_foam(my_sys, sys.argv)
    else:
        argv_load_atoms(my_sys, ["", sys.argv[1]])
    # Interpret the commands
    files, settings, groups, builds, exports, ifaces = interpret_argvs()
    # Go through each of the ls
    argv_load(my_sys, files)
    argv_sett(my_sys, settings)
    max_vert = my_sys.net.max_vert
    argv_group(my_sys, groups, bff=ifaces)
    if my_sys.net2:
        my_sys.net = Network(sys=my_sys, atoms=my_sys.atoms, net_type='vor')
        atoms2 = my_sys.atoms.copy()
        for my_group in my_sys.groups:
            if len(my_group.atoms) > 0:
                my_sys.net.build(my_group=my_group, max_vert=max_vert)
        atom_vals_vor = []
        atom_nums_vor = []
        for i, atom in my_sys.net.atoms.iterrows():
            if atom['complete']:
                atom_vals_vor.append({'num': i, 'vol': atom['vol'], 'sa': atom['sa']})
                atom_nums_vor.append(i)
            else:
                atom_vals_vor.append({})
        my_sys.net2 = my_sys.net
        my_sys.net = Network(sys=my_sys, atoms=atoms2, net_type='pow')
        for my_group in my_sys.groups:
            if len(my_group.atoms) > 0:
                my_sys.net.build(my_group=my_group, max_vert=max_vert)
        atom_vals = []
        for i, atom in my_sys.net.atoms.iterrows():
            if atom['complete'] and i in atom_nums_vor:
                atom_vals.append({'num': i, 'vol_pow': atom['vol'], 'vol_vor': atom_vals_vor[i]['vol'],
                                  'sa': atom['sa'], 'sa_vor': atom_vals_vor[i]['sa'],
                                  'vol_diff': abs(atom['vol'] - atom_vals_vor[i]['vol']) / atom_vals_vor[i]['vol'],
                                  'sa_diff': abs(atom['sa'] - atom_vals_vor[i]['sa']) / atom_vals_vor[i]['sa']})
        print('vol avg diff', sum([_['vol_diff'] for _ in atom_vals]) / len(atom_vals), 'sa avg diff', sum([_['sa_diff'] for _ in atom_vals]) / len(atom_vals), 'num cells', len(atom_vals))
        return

    elif my_sys.net_file is None:
        argv_build(my_sys, builds)

    argv_export(my_sys, exports)


if __name__ == '__main__':
    mySys = System()
    argv(mySys)
