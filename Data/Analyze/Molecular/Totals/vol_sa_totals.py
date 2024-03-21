import csv
import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.plot_templates.bar import bar


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    drop_box_folder = filedialog.askdirectory()
    folder = drop_box_folder + '/Jack/Vorpy/Data/Molecular/logs_and_pdbs/'
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
    plotting = 'Volume'

    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'csv':
                with open(folder + '/' + file, 'r') as my_file:
                    my_reader = csv.reader(my_file)
                    for i, line in enumerate(my_reader):
                        if i == 5:
                            if plotting == '':
                                my_logs[file[:-13]][file[-12:-9]] = float(line[2])
                            else:
                                my_logs[file[:-13]][file[-12:-9]] = float(line[3])
    vor_vals = [my_logs[_]['vor'] for _ in my_sys_names]
    pow_vals = [my_logs[_]['pow'] for _ in my_sys_names]
    del_vals = [my_logs[_]['del'] for _ in my_sys_names]

    pow_diff = [100 * abs(vor_vals[i] - pow_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
    del_diff = [100 * abs(vor_vals[i] - del_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
    bar([pow_diff, del_diff], x_names=my_sys_names, legend_names=['Power Difference', 'Primitive Difference'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Total Molecule Surface Area')

