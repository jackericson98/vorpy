import time

from Visualize.argv.argv_load import argv_load, argv_load_atoms, argv_load_foam
from Visualize.argv.argv_set import argv_sett
from Visualize.argv.argv_group import argv_group
from Visualize.argv.argv_build import argv_build
from Visualize.argv.argv_export import argv_export
from System.Network.network import Network
from System.system import System
import sys
import os
import csv

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
    # Set up the commands dictionary
    cmnds = {'npt': [], 'set': [], 'grp': [], 'bld': [], 'xpt': [], 'ifc': []}
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
        # Add the command to the 181L list
        if arg.lower() == '-l':
            cmnds['npt'].append(arg_cmnds)
        elif arg.lower() == '-s':
            cmnds['set'].append(arg_cmnds)
        elif arg.lower() == '-g':
            cmnds['grp'].append(arg_cmnds)
        elif arg.lower() == '-b':
            cmnds['bld'].append(arg_cmnds)
        elif arg.lower() == '-e':
            cmnds['xpt'].append(arg_cmnds)
        elif arg.lower() == '-i':
            cmnds['ifc'].append(arg_cmnds)
    # Return the lists
    return cmnds


def argv(my_sys):
    my_sys.print_actions = False
    # Load the atom file
    if sys.argv[1].lower() == 'foam':
        argv_load_foam(my_sys, sys.argv)
    else:
        argv_load_atoms(my_sys, ["", sys.argv[1]])
    # Interpret the commands
    cmnds = interpret_argvs()
    my_sys.cmnds = cmnds
    # Go through each of the ls
    argv_load(my_sys, cmnds['npt'])
    argv_sett(my_sys, cmnds['set'])
    max_vert = my_sys.net.max_vert
    argv_group(my_sys, cmnds['grp'], bff=cmnds['ifc'])
    # If we are comparing two network types
    if my_sys.net2:
        start = time.perf_counter()
        my_sys.net = Network(sys=my_sys, atoms=my_sys.atoms, net_type='pow')
        my_sys.name = my_sys.name + '_pow'
        atoms2 = my_sys.atoms.copy()
        my_sys.set_output_directory(cmnds['xpt'][0][1] + '/pow')
        for my_group in my_sys.groups:
            if len(my_group.atoms) > 0:
                my_sys.net.build(my_group=my_group, max_vert=max_vert, print_vert_metrics=False, print_actions=False, surf_res=0.7)
        argv_export(my_sys, cmnds['xpt'], add_on='/pow')
        atom_vals_pow = []
        atom_nums_pow = []
        for i, atom in my_sys.net.atoms.iterrows():
            if atom['complete']:
                atom_vals_pow.append({'num': i, 'vol': atom['vol'], 'sa': atom['sa']})
                atom_nums_pow.append(i)
            else:
                atom_vals_pow.append({})
        my_sys.net2 = my_sys.net
        my_sys.net = Network(sys=my_sys, atoms=atoms2, net_type='vor')
        my_sys.set_output_directory(cmnds['xpt'][0][1] + '/vor')
        os.chdir(my_sys.dir)
        for my_group in my_sys.groups:
            if len(my_group.atoms) > 0:
                my_sys.net.build(my_group=my_group, max_vert=max_vert, print_vert_metrics=False, print_actions=False, surf_res=0.7)
        atom_vals = []
        for i, atom in my_sys.net.atoms.iterrows():
            if atom['complete'] and i in atom_nums_pow:

                vdp = abs(atom['vol'] - atom_vals_pow[i]['vol']) / atom_vals_pow[i]['vol']
                sdp = abs(atom['sa'] - atom_vals_pow[i]['sa']) / atom_vals_pow[i]['sa']
                vdv = abs(atom['vol'] - atom_vals_pow[i]['vol']) / atom['vol']
                sdv = abs(atom['sa'] - atom_vals_pow[i]['sa']) / atom['sa']

                if vdp > 10 or sdp > 10 or vdv > 10 or sdv > 10:
                    continue

                atom_vals.append({'num': i, 'vol_vor': atom['vol'], 'vol_pow': atom_vals_pow[i]['vol'],
                                  'sa': atom['sa'], 'sa_pow': atom_vals_pow[i]['sa'],
                                  'vol_diff_pow': vdp, 'sa_diff_pow': sdp, 'vol_diff_vor': vdv, 'sa_diff_vor': sdv})
        folder = os.path.dirname(my_sys.base_file)
        if my_sys.foam_data is None:
            my_sys.foam_data = [0, 0, 0, 0, 0]
        if len(atom_vals) > 0:
            my_line = "\r{}".format(folder), *my_sys.foam_data, sum([_['vol_diff_vor'] for _ in atom_vals]) / len(atom_vals), sum([_['sa_diff_vor'] for _ in atom_vals]) / len(atom_vals), sum([_['vol_diff_pow'] for _ in atom_vals]) / len(atom_vals), sum([_['sa_diff_pow'] for _ in atom_vals]) / len(atom_vals), len(atom_vals), round((time.perf_counter() - start), 3)
        else:
            my_line = "\r{}".format(folder), *my_sys.foam_data, 0, 0, 0, 0, len(atom_vals), round((time.perf_counter() - start), 3)
        print(*my_line, end="")

        current_dir = os.getcwd()
        os.chdir('../..')
        try:
            with open('foam_data.csv', 'a') as foam_file:
                foam_writer = csv.writer(foam_file)
                foam_writer.writerow(my_line)
        except PermissionError:
            with openn('foam_data1.csv', 'a') as foam_file:
                foam_writer = csv.writer(foam_file)
                foam_writer.writerow(my_line)
        os.chdir(current_dir)
        argv_export(my_sys, cmnds['xpt'], add_on='/vor')
        return
    elif my_sys.net_file is None:
        argv_build(my_sys, cmnds['bld'])

    argv_export(my_sys, cmnds['xpt'])


if __name__ == '__main__':
    mySys = System()
    argv(mySys)
