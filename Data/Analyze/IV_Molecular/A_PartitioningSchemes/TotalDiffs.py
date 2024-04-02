import csv
import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.bar import bar


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    # Get the systems in the designated folder
    systems = []
    for root, directory, files in os.walk(folder):
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

    # Choose what we are plotting Vol or SA
    plotting = ''

    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'csv':
                with open(folder + '/' + file, 'r') as my_file:
                    my_reader = csv.reader(my_file)
                    print(file)
                    for i, line in enumerate(my_reader):
                        if i == 5:
                            if plotting == '':
                                print(line[2])
                                my_logs[file[:-13]][file[-12:-9]] = float(line[2])
                            else:
                                print(line[3])
                                my_logs[file[:-13]][file[-12:-9]] = float(line[3])
    vor_vals = [my_logs[_]['vor'] for _ in my_sys_names]
    pow_vals = [my_logs[_]['pow'] for _ in my_sys_names]
    del_vals = [my_logs[_]['del'] for _ in my_sys_names]
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP'}[_] for _ in my_sys_names]

    pow_diff = [100 * abs(vor_vals[i] - pow_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
    del_diff = [100 * abs(vor_vals[i] - del_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
    # Create the
    bar([pow_diff, del_diff], x_names=graph_labels, legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Total Volume Difference')

