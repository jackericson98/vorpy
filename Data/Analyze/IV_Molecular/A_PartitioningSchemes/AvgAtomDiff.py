import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from Data.Analyze.tools.compare.read_logs import read_logs
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.bar import bar


def atoms_per_diff(systems, logs, val='volume'):
    # Create the averages lists
    avg_pow_diffs, pow_ses, avg_del_diffs, del_ses = [], [], [], []
    # Go through the loaded systems
    for system in systems:
        sys_name = system.name
        vor_atoms = read_logs(logs[sys_name]['vor'], True, True)['atoms']
        pow_atoms = read_logs(logs[sys_name]['pow'], True, True)['atoms']
        del_atoms = read_logs(logs[sys_name]['del'], True, True)['atoms']
        print(system.name)
        print(vor_atoms)
        # print(sys_name, len(vor_atoms), len(pow_atoms), len(del_atoms))
        pow_diffs, del_diffs = [], []
        for i in range(len(vor_atoms)):
            pow_diffs.append(abs(vor_atoms[i][val] - pow_atoms[i][val])/vor_atoms[i][val])
            del_diffs.append(abs(vor_atoms[i][val] - del_atoms[i][val])/vor_atoms[i][val])
        # Calculate the averages
        avg_pow_diffs.append(100*sum(pow_diffs)/len(pow_diffs))
        avg_del_diffs.append(100*sum(del_diffs)/len(del_diffs))
        # Calculate the standard errors
        pow_ses.append(100*np.std(pow_diffs)/np.sqrt(len(pow_diffs)))
        del_ses.append(100*np.std(del_diffs)/np.sqrt(len(del_diffs)))
    # Return the values
    return avg_pow_diffs, pow_ses, avg_del_diffs, del_ses


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    # Create the systems
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
    # Get the average differences and the standard errors
    avg_pow_diffs, pow_ses, avg_del_diffs, del_ses = atoms_per_diff(systems, my_logs, val='sa')
    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig'}[_] for _ in my_sys_names]
    # Plot the Data
    bar(data=[avg_pow_diffs, avg_del_diffs], title='Average Atom Surface Area % Difference', legend_title='Scheme',
        y_axis_title='% Difference', x_names=graph_labels, legend_names=['Power', 'Primitive'], Show=True,
        x_axis_title='Model', errors=[pow_ses, del_ses])

