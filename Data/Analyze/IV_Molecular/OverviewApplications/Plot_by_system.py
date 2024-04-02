import os
from Data.Analyze.tools.compare.read_logs import read_logs
from System.system import System
import tkinter as tk
from System.Group.group import Group
from tkinter import filedialog
import matplotlib.pyplot as plt


"""
we want to be able to return a dictionary that is able to be analyzed that looks something like this 
{by1: {avg vol: ..., tot vol: ..., vol list: [...], avg sa: ..., tot sa: ..., sa list: [...], avg curv: ..., curv list: [...]}
"""


def get_vals(log_info, my_sys, by='element', fig=None, ax=None):

    if not my_sys.groups:
        my_group = Group(my_sys, residues=my_sys.residues)
        my_sys.groups = [my_group]
    my_log_atom = read_logs(log_info, True)
    my_dict = {}
    for sys_atom in my_sys.groups[0].atoms:
        # Get the log atom
        atom = my_sys.atoms.loc[my_sys.atoms['num'] == sys_atom]
        # add the values to the averages
        if atom[by][sys_atom] in my_dict:
            my_dict[atom[by][sys_atom]]['vol list'].append(my_log_atom['volume'][sys_atom])
            my_dict[atom[by][sys_atom]]['sa list'].append(my_log_atom['sa'][sys_atom])
            my_dict[atom[by][sys_atom]]['curv list'].append(my_log_atom['max curv'][sys_atom])
        else:
            print(my_log_atom)
            my_dict[atom[by][sys_atom]] = {
                'vol list': [my_log_atom['volume'][sys_atom]],
                'sa list': [my_log_atom['sa'][sys_atom]],
                'curv list': [my_log_atom['max curv'][sys_atom]]
            }
    for _ in my_dict:
        my_dict[_]['avg vol'] = sum(my_dict[_]['vol list'])/len(my_dict[_]['vol list'])
        my_dict[_]['tot vol'] = sum(my_dict[_]['vol list'])
        my_dict[_]['avg sa'] = sum(my_dict[_]['sa list']) / len(my_dict[_]['sa list'])
        my_dict[_]['tot sa'] = sum(my_dict[_]['sa list'])
        my_dict[_]['avg curv'] = sum(my_dict[_]['curv list']) / len(my_dict[_]['curv list'])
    return my_dict


def plot_distribution(data, title='', labels=None):
    return


def plot_vals(val_dict, plotting='tot vol', title='', show=False, by='Atom'):
    """
    We want a function that plots the values for the dictionary by element
    """
    labels = [_ for _ in val_dict]
    data = [val_dict[_][plotting] for _ in labels]
    y_label_dict = {'tot vol': 'Total Volume', 'avg vol': 'Average Volume', 'vol list': 'Volume',
                    'tot sa': 'Total Surface Area', 'avg sa': 'Average Surface Area', 'sa list': 'Surface Area'}

    # This is where we plot the distributions
    if 'list' in plotting:
        plot_distribution(data, title=title)
    else:
        ymax = max(data)
        # Bar width
        bar_width = 0.35

        # Index for the x-axis
        x = range(len(labels))

        # Create the bar graph
        plt.bar(x, data, width=bar_width, color='skyblue', edgecolor='black')
        # Add labels and title
        plt.ylabel(y_label_dict[plotting], fontdict=dict(size=15))
        plt.xlabel(by, fontdict=dict(size=15))
        plt.title(title, fontdict=dict(size=20))

        # Angle the labels and add values at the top of the bars
        plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right')
        # for i, v in enumerate(data1):
        #     plt.text(i, v / 2, str(v) + ' \u212B\u00B3', ha='center', va='center', rotation=90)
        # for i, v in enumerate(data2):
        #     plt.text(i + bar_width, v / 2, str(v) + ' \u212B\u00B3', ha='center', va='center', rotation=90)
        plt.ylim(0, 1.25 * ymax)
        # Add legend with appropriate layout
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2)
        plt.title(title, fontdict=dict(size=20))
        # If the user wants to show the graph
        if show:
            plt.show()


if __name__ == '__main__':
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
    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'csv':
                my_logs[file[:-13]][file[-12:-9]] = folder + '/' + file
    for system in systems:
        for name in {'vor', 'pow', 'del'}:
            for my_type in {'Residue', 'Atom'}:
                my_vals = get_vals(my_logs[system.name][name], system, by='res_name')
                for data in {'avg vol', 'tot vol', 'vol list', 'avg sa', 'tot sa', 'sa list'}:
                    plot_vals(my_vals, show=True, plotting=data, by=my_type, title=system.name + ' ' + my_type + ' ' + data)



