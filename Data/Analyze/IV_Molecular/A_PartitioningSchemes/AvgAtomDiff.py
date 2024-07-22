import os
import numpy as np
import scipy as sp
import tkinter as tk
from tkinter import filedialog
from System.system import System
from Data.Analyze.tools.compare.read_logs import read_logs
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.bar import bar
from Data.Analyze.tools.plot_templates.scatter import scatter
import matplotlib.pyplot as plt


def atoms_per_diff(systems, logs, val='volume'):
    # Create the averages lists
    avg_pow_diffs, pow_ses, avg_del_diffs, del_ses = [], [], [], []
    # Go through the loaded systems
    all_pow_diffs, all_del_diffs, all_atom_types, all_pow_atom_types, all_del_atom_types = [], [], [], [], []
    for system in systems:
        sys_name = system.name
        vor_atoms = read_logs(logs[sys_name]['vor'], True, True)['atoms']
        pow_atoms = read_logs(logs[sys_name]['pow'], True, True)['atoms']
        del_atoms = read_logs(logs[sys_name]['del'], True, True)['atoms']
        # print(sys_name, len(vor_atoms), len(pow_atoms), len(del_atoms))
        pow_diffs, del_diffs, atom_types, pow_atom_types, del_atom_types = [], [], [], {'C': [], 'H': [], 'P': [], 'S': [], 'Se': [], 'O': [], 'N': []}, {'C': [], 'H': [], 'P': [], 'S': [], 'Se': [], 'O': [], 'N': []}
        for i in range(len(vor_atoms)):
            if abs(100*(vor_atoms[i][val] - pow_atoms[i][val])/vor_atoms[i][val]) > 150:
                pow_diffs.append(100*(vor_atoms[i][val] - pow_atoms[i][val])/vor_atoms[i][val])
                del_diffs.append(100*(vor_atoms[i][val] - del_atoms[i][val])/vor_atoms[i][val])
                atom_types.append(system.atoms['element'].iloc[i].strip())
                pow_atom_types[system.atoms['element'].iloc[i].strip()].append(100*(vor_atoms[i][val] - pow_atoms[i][val])/vor_atoms[i][val])
                del_atom_types[system.atoms['element'].iloc[i].strip()].append(100*(vor_atoms[i][val] - del_atoms[i][val])/vor_atoms[i][val])
        # Add the elements to their lists
        all_atom_types.append(atom_types)
        # Add the pow diffs and the dell diffs to their lists
        all_pow_diffs.append(pow_diffs)
        all_del_diffs.append(del_diffs)
        # Add the pow and del atom_type diffs
        all_pow_atom_types.append(pow_atom_types)
        all_del_atom_types.append(del_atom_types)
        # Calculate the averages
        avg_pow_diffs.append(sum([abs(_) for _ in pow_diffs])/len(pow_diffs))
        avg_del_diffs.append(sum([abs(_) for _ in del_diffs])/len(del_diffs))
        # Calculate the standard errors
        pow_ses.append(sp.stats.sem([abs(_) for _ in pow_diffs]))
        del_ses.append(np.std([abs(_) for _ in del_diffs])/np.sqrt(len(del_diffs)))
    # Return the values
    return avg_pow_diffs, pow_ses, avg_del_diffs, del_ses, all_pow_diffs, all_del_diffs, all_atom_types


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
    avg_pow_diffs, pow_ses, avg_del_diffs, del_ses, all_pow_diffs, all_del_diffs, atom_types\
        = atoms_per_diff(systems, my_logs)
    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig',
                     '1BNA': '1BNA', 'DB1976': 'DB1976', 'BSA': 'BSA'}[_] for _ in my_sys_names]
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J', '1BNA': 'K', 'DB1976': 'L'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    def sort_4_lists(lista, listb, listc, listd, liste, listf, listg, listh):
        # Zipping lists together and sorting by the first list
        sorted_lists = sorted(zip(lista, listb, listc, listd, liste, listf, listg, listh), key=lambda x: x[0])

        # Unpacking the sorted lists
        lista, listb, listc, listd, liste, listf, listg, listh = zip(*sorted_lists)

        # Converting tuples back to lists if needed
        lista = list(lista)
        listb = list(listb)
        listc = list(listc)
        listd = list(listd)
        liste = list(liste)
        listf = list(listf)
        listg = list(listg)
        listh = list(listh)

        # Return the lists
        return lista, listb, listc, listd, liste, listf, listg, listh

    new_graph_labels, avg_pow_diffs, pow_ses, avg_del_diffs, del_ses, all_pow_diffs, all_del_diffs, atom_types = (
        sort_4_lists(new_graph_labels, avg_pow_diffs, pow_ses, avg_del_diffs, del_ses, all_pow_diffs, all_del_diffs, atom_types))

    # Plot the Data
    bar(data=[avg_pow_diffs, avg_del_diffs], legend_title='Scheme',
        y_axis_title='% Difference', x_names=new_graph_labels, legend_names=['Power', 'Primitive'], Show=True,
        x_axis_title='Model', errors=[pow_ses, del_ses], y_range=[0, 70], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)


    my_colors = []
    color_dict = {'C': 'grey', 'O': 'r', 'N': 'b', 'P': 'darkorange', 'H': 'pink', 'S': 'y', 'Se': 'sandybrown'}
    for system in atom_types:
        system_colors = []
        for atom in system:
            if atom in color_dict:
                system_colors.append(color_dict[atom])
            else:
                print(atom)
                system_colors.append('purple')
        my_colors.append(system_colors)

    # Plot the data
    scatter(xs=[[3 * i + 1 for _ in all_pow_diffs[i]] for i in range(len(all_pow_diffs))] +
               [[3 * i + 2 for _ in all_del_diffs[i]] for i in range(len(all_del_diffs))],
            ys=all_pow_diffs + all_del_diffs, Show=False,
            markers=[['${}$'.format(_) for _ in system] for system in atom_types] * 2, alpha=0.3,
            colors=my_colors + my_colors,
            x_tick_labels=new_graph_labels,
            x_tick_label_locs=[3 * i + 1.5 for i in range(len(all_pow_diffs))],
            xlabel_size=30,
            x_axis_title='Model',
            ylabel_size=30,
            y_axis_title='% Deviation',
            marker_size=100, xlabel_rotation=45, xtick_anchor='ha')

    plt.plot([1], [1], label='Power - Left')

    plt.plot([1], [1], label='Primitive - Right')

    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2,
                        prop={'size': 25})

    plt.show()
