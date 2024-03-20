import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
import matplotlib.pyplot as plt
from Data.Analyze.read_logs import read_logs
from System.Group.group import Group


def atoms_per_diff(systems, logs, val='volume', title=''):
    avg_pow_diffs, avg_del_diffs = [], []
    for system in systems:
        sys_name = system.name
        vor_atoms = read_logs(logs[sys_name]['vor'], True, True)['atoms']
        pow_atoms = read_logs(logs[sys_name]['pow'], True, True)['atoms']
        del_atoms = read_logs(logs[sys_name]['del'], True, True)['atoms']
        # print(sys_name, len(vor_atoms), len(pow_atoms), len(del_atoms))
        pow_diffs, del_diffs = [], []
        for i in range(len(vor_atoms)):
            pow_diffs.append(abs(vor_atoms[i][val] - pow_atoms[i][val])/vor_atoms[i][val])
            del_diffs.append(abs(vor_atoms[i][val] - del_atoms[i][val])/vor_atoms[i][val])
        avg_pow_diffs.append((sum(pow_diffs)/len(pow_diffs)))
        avg_del_diffs.append((sum(del_diffs)/len(del_diffs)))
    ymax = max(avg_del_diffs + avg_pow_diffs) * 100
    bar_width = 0.35
    x = range(len(systems))
    plt.bar(x, [_*100 for _ in avg_pow_diffs], width=bar_width, label='Power', color='skyblue', edgecolor='black')
    plt.bar([i + bar_width for i in x], [_*100 for _ in avg_del_diffs], width=bar_width, label='Primitive', color='orange',
                   edgecolor='black')
    plt.ylabel('% Difference', fontdict=dict(size=15))

    plt.title(title, fontdict=dict(size=20))
    plt.xticks([i + bar_width / 2 for i in x], [_.name for _ in systems], rotation=45, ha='right', font=dict(size=10))
    plt.legend(title='Scheme')
    plt.ylim(0, 1.3 * ymax)
    plt.show()


if __name__ == '__main__':
    # root = tk.Tk()
    # root.withdraw()
    # root.wm_attributes('-topmost', 1)
    # folder = filedialog.askdirectory()
    folder = '/Users/jackericson/Documents/logs'
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

    atoms_per_diff(systems, my_logs, title='Average Surface Area Difference', val='sa')




