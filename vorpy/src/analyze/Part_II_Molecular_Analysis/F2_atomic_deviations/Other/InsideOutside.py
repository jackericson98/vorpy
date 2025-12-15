import os
import numpy as np
import tkinter as tk
from tkinter import filedialog

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)


from vorpy.src.system.system import System
from vorpy.src.group.group import Group
from vorpy.src.analyze.tools.compare.read_logs import read_logs
from vorpy.src.analyze.tools.plot_templates.bar import bar


# Function to determine if an atom is inside the given group or not
def inside(atom_neighbors, group_nums):
    for num in atom_neighbors:
        if num not in group_nums:
            return False
    else:
        return True


# Function for plotting the
def inside_out(my_systems, logs, val='volume', title=''):
    sys_data = {}
    for system in my_systems:
        sys_name = system.name
        vor_atoms = read_logs(logs[sys_name]['vor'], True, True)['atoms']
        pow_atoms = read_logs(logs[sys_name]['pow'], True, True)['atoms']
        del_atoms = read_logs(logs[sys_name]['del'], True, True)['atoms']
        in_vor, out_vor, in_pow, out_pow, in_del, out_del = [], [], [], [], [], []
        for i in range(len(vor_atoms)):
            # Get the group atoms
            grp_atms = system.groups[0].atoms
            # Grab the atom dicts
            va, pa, da = vor_atoms[i], pow_atoms[i], del_atoms[i]
            # First we need to know if the atom is inside or outside
            vor_in, pow_in, del_in = inside(va['neighbors'], grp_atms), inside(pa['neighbors'], grp_atms), inside(pa['neighbors'], grp_atms)
            # if the atom isn't similarly on the inside and outside skip it
            if not all([vor_in, pow_in, del_in]) and not all([not vor_in, not pow_in, not del_in]):
                continue
            if vor_in:
                in_vor.append(va[val])
            else:
                out_vor.append(va[val])
            if pow_in:
                in_pow.append(pa[val])
            else:
                out_pow.append(pa[val])
            if del_in:
                in_del.append(da[val])
            else:
                out_del.append(da[val])
        sys_data[system.name] = {'in_vor': in_vor, 'out_vor': out_vor, 'in_pow': in_pow, 'out_pow': out_pow,
                                 'in_del': in_del, 'out_del': out_del}

    sys_names = [_.name for _ in my_systems if all([len(sys_data[_.name][__]) > 0 for __ in sys_data[my_systems[0].name]])]

    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig'}[_] for _ in sys_names]

    # Get the average volume for an inside and outside bulk for each scheme
    bulk_avgs = [([sum(sys_data[_]['in_vor']) / len(sys_data[_]['in_vor']) for _ in sys_names],
                  [sum(sys_data[_]['out_vor']) / len(sys_data[_]['out_vor']) for _ in sys_names]),
                 ([sum(sys_data[_]['in_pow']) / len(sys_data[_]['in_pow']) for _ in sys_names],
                  [sum(sys_data[_]['out_pow']) / len(sys_data[_]['out_pow']) for _ in sys_names]),
                 ([sum(sys_data[_]['in_del']) / len(sys_data[_]['in_del']) for _ in sys_names],
                  [sum(sys_data[_]['out_del']) / len(sys_data[_]['out_del']) for _ in sys_names])]

    # Get the percent differences of the bulks
    pow_in_bulk_pds = [100 * abs(bulk_avgs[1][0][i] - bulk_avgs[0][0][i]) / bulk_avgs[0][0][i] for i in range(len(bulk_avgs[0][0]))]
    pow_out_bulk_pds = [100 * abs(bulk_avgs[1][1][i] - bulk_avgs[0][1][i]) / bulk_avgs[0][1][i] for i in range(len(bulk_avgs[0][0]))]
    del_in_bulk_pds = [100 * abs(bulk_avgs[2][0][i] - bulk_avgs[0][0][i]) / bulk_avgs[0][0][i] for i in range(len(bulk_avgs[0][0]))]
    del_out_bulk_pds = [100 * abs(bulk_avgs[2][1][i] - bulk_avgs[0][1][i]) / bulk_avgs[0][1][i] for i in range(len(bulk_avgs[0][0]))]

    # Get the percent difference by atom
    atom_diffs = []
    atom_ses = []
    for my_str in {'in', 'out'}:
        for comp_type in {'_pow', '_del'}:
            atom_diffs.append([])
            atom_ses.append([])
            for sys_name in sys_names:
                # Average it all out by volume
                vals = [abs(sys_data[sys_name]['{}_vor'.format(my_str)][i] - sys_data[sys_name][my_str + comp_type][i]) / sys_data[sys_name]['{}_vor'.format(my_str)][i] for i in range(len(sys_data[sys_name]['in_pow']))]
                atom_diffs[-1].append(100 * sum(vals)/len(vals))
                atom_ses[-1].append(100 * np.std(vals) / np.sqrt(len(vals)))
    # Plot the data
    # for i in range(3):
        # # Choose the data and title based on the scheme
        # scheme = ['Additively Weighted', 'Power', 'Primitive'][i]
        # val_name = {'volume': 'Volume', 'sa': 'Surface Area', 'max curv': 'Curvature'}[val]
        # Plot using the bar function
    bar(atom_diffs, x_names=graph_labels, errors=atom_ses,
        legend_names=['Inside Power', 'Inside Primitive', 'Outside Power', 'Outside Primitive'], legend_title='Scheme',
        title='Average Volume % Difference By Atom', x_axis_title='Model', y_axis_title='% Difference',
        Show=True)


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    systems = []
    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'pdb':
                my_sys = System(file=folder + '/' + file)
                my_sys.groups = [Group(sys=my_sys, residues=my_sys.residues)]
                systems.append(my_sys)

    # Sort atoms by number of atoms
    num_atoms = [len(_.atoms) for _ in systems]
    systems = [x for _, x in sorted(zip(num_atoms, systems))]

    # Create the logs dictionary
    my_sys_names = [__.name for __ in systems]
    my_logs = {_: {'vor': None, 'pow': None, 'del': None} for _ in my_sys_names}
    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'csv':
                my_logs[file[:-13]][file[-12:-9]] = folder + '/' + file

    inside_out(systems, my_logs, val='volume')

