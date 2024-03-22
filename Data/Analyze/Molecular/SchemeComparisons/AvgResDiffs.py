import csv
import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.plot_templates.bar import bar
from Data.Analyze.read_logs import read_logs
from Data.Analyze.residue import residue_data


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

    # Create the log file name dictionary
    my_log_files = {_: {'vor': folder + _.name + '_vor_logs.csv',
                   'pow': folder + _.name + '_pow_logs.csv',
                   'del': folder + _.name + '_del_logs.csv'} for _ in my_sys_names}

    # Create the log dictionary
    my_logs = {}

    # Get the log information
    vor_vals, pow_vals, del_vals = [], [], []
    for system in systems:
        vor_vals.append(residue_data(system, read_logs(my_log_files[system.name]['vor'], return_dict=True), vol=True,
                                     sa=True))

        pow_vals.append(residue_data(system, read_logs(my_log_files[system.name]['pow'], return_dict=True), vol=True,
                                     sa=True))

        del_vals.append(residue_data(system, read_logs(my_log_files[system.name]['del'], return_dict=True), vol=True,
                                     sa=True))



    # Choose what we are plotting Vol or SA
    plotting = 'Volume'

    # Reorganize
    vor_vals = [my_logs[_]['vor'] for _ in my_sys_names]
    pow_vals = [my_logs[_]['pow'] for _ in my_sys_names]
    del_vals = [my_logs[_]['del'] for _ in my_sys_names]

    pow_diff = [100 * abs(vor_vals[i] - pow_vals[i])/vor_vals[i] for i in range(len(vor_vals))]

    del_diff = [100 * abs(vor_vals[i] - del_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
    # Create the
    bar([pow_diff, del_diff], x_names=my_sys_names, legend_names=['Power', 'Primitive'],
        Show=False, y_axis_title='% Difference', x_axis_title='Model', title='Total Molecule Surface Area')

