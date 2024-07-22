import os

import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from Data.Analyze.tools.compare.read_logs import read_logs
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.box_whisker import box_whisker


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
    vals = []
    for sys_name in my_sys_names:
        logs = read_logs(my_logs[sys_name]['vor'], return_dict=True)
        vals.append([_['curvature'] for _ in logs['surfs']])
    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig'}[_] for _ in my_sys_names]

    # Bin the values
    increments = np.linspace(0, max([max(_) for _ in vals]), 80)
    for i, _ in enumerate(vals):
        # Scale the number of values by the size of the group
        bins = [[__ for __ in _ if increments[j] <= __ < increments[j + 1]] for j in range(len(increments[:-1]))]
        # Count the number based on the number of total values
        scaled_bins = [len(__) / len(_) for __ in bins]
        plt.plot([(increments[j] + increments[j + 1]) / 2 for j in range(len(increments[:-1]))][2:], [100* __ for __ in scaled_bins[2:]])

    plt.legend(graph_labels, ncols=len(graph_labels))
    plt.title('Surface Curvature Distributions By Model', fontdict=dict(size=20))
    plt.ylabel('Percentage', fontdict=dict(size=15))
    plt.xticks()
    plt.xlabel('Curvature', fontdict=dict(size=15))
    plt.show()

