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
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.plot_templates.bar import bar

from vorpy.src.analyze.tools.batch.get_files import get_all_files


def per_diffs_surfs(skip_keys=[]):
    """Plots the average percentage differences for the given systems"""
    # Get the files
    files = get_all_files()
    # Loop through the files and get the surface data
    for key, value in files.items():
        if key in skip_keys:
            continue
        # Get the voronoi surface info
        aw_logs = read_logs2(value['aw'], all_=False, surfs=True)
        # Get the power surface info
        pow_logs = read_logs2(value['pow'], all_=False, surfs=True)
        # Get the primitive surface info
        prm_logs = read_logs2(value['prm'], all_=False, surfs=True)
        # add the logs to the files dictionary
        files[key]['aw surfs'] = aw_logs['surfs']
        files[key]['pow surfs'] = pow_logs['surfs']
        files[key]['prm surfs'] = prm_logs['surfs']
        
    



def surfs_per_diff(systems, logs, val='sa'):
    # Create the averages lists
    avg_pow_diffs, pow_ses, avg_del_diffs, del_ses = [], [], [], []
    # Go through the loaded systems
    for system in systems:
        sys_name = system.name
        vor_surfs = read_logs2(logs[sys_name]['vor'], True, True)['surfs']
        pow_surfs = {str(_['atoms'][0]) + '_' + str(_['atoms'][1]): _
                     for _ in read_logs2(logs[sys_name]['pow'], True, True)['surfs']}
        del_surfs = {str(_['atoms'][0]) + '_' + str(_['atoms'][1]): _
                     for _ in read_logs2(logs[sys_name]['del'], True, True)['surfs']}
        # print(sys_name, len(vor_atoms), len(pow_atoms), len(del_atoms))
        pow_diffs, del_diffs = [], []
        for i in range(len(vor_surfs)):
            if vor_surfs[i][val] == 0:
                continue
            # Get the hash string
            hsh_str = str(vor_surfs[i]['atoms'][0]) + '_' + str(vor_surfs[i]['atoms'][1])
            if hsh_str in pow_surfs:
                pow_diffs.append(abs(vor_surfs[i][val] - pow_surfs[hsh_str][val])/vor_surfs[i][val])
            if hsh_str in del_surfs:
                if system.name == 'BSA' and abs(vor_surfs[i][val] - del_surfs[hsh_str][val])/vor_surfs[i][val] > 100:

                    print(hsh_str, round(abs(vor_surfs[i][val] - del_surfs[hsh_str][val])/vor_surfs[i][val], 3))
                    continue
                del_diffs.append(abs(vor_surfs[i][val] - del_surfs[hsh_str][val])/vor_surfs[i][val])
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
    avg_pow_diffs, pow_ses, avg_del_diffs, del_ses = surfs_per_diff(systems, my_logs, val='sa')
    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP',
                     'pl_complex': 'Prot-Lig', 'BSA': 'BSA'}[_] for _ in my_sys_names]
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    def sort_4_lists(lista, listb, listc, listd, liste):
        # Zipping lists together and sorting by the first list
        sorted_lists = sorted(zip(lista, listb, listc, listd, liste), key=lambda x: x[0])

        # Unpacking the sorted lists
        lista, listb, listc, listd, liste = zip(*sorted_lists)

        # Converting tuples back to lists if needed
        lista = list(lista)
        listb = list(listb)
        listc = list(listc)
        listd = list(listd)
        liste = list(liste)

        # Return the lists
        return lista, listb, listc, listd, liste

    new_graph_labels, avg_pow_diffs, pow_ses, avg_del_diffs, del_ses = (
        sort_4_lists(new_graph_labels, avg_pow_diffs, pow_ses, avg_del_diffs, del_ses))

    # Plot the Data
    bar(data=[avg_pow_diffs, avg_del_diffs], legend_title='Scheme',
        y_axis_title='% Difference', x_names=new_graph_labels, legend_names=['Power', 'Primitive'], Show=True,
        x_axis_title='Model', errors=[pow_ses, del_ses], xlabel_size=30, ylabel_size=30, ytick_label_size=25, xtick_label_size=25, tick_width=2, tick_length=12)

