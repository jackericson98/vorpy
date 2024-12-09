import csv
import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.bar import bar


if __name__ == '__main__':

    # Choices
    # Choose what we are plotting Vol or SA
    plotting = 'Vol'
    # Choose to make it abs difference or total difference
    diff = 'tot'

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

    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'csv':
                with open(folder + '/' + file, 'r') as my_file:
                    my_reader = csv.reader(my_file)
                    for i, line in enumerate(my_reader):
                        if i == 5:
                            if plotting == 'Vol':
                                my_logs[file[:-13]][file[-12:-9]] = float(line[2])
                            elif plotting == 'SA':
                                my_logs[file[:-13]][file[-12:-9]] = float(line[3])
    vor_vals = [my_logs[_]['vor'] for _ in my_sys_names]
    pow_vals = [my_logs[_]['pow'] for _ in my_sys_names]
    del_vals = [my_logs[_]['del'] for _ in my_sys_names]
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig',
                     'BSA': 'BSA', 'hammerhead': 'H-Head', '1BNA': '1BNA', 'DB1976': 'DB1976'}
                    [_] for _ in my_sys_names]
    rb_sa_vals_pymol = [{'EDTA_Mg': 290.72, 'cambrin': 3180.126, '1BNA': 4665.118, '181L': 10085.962, 'DB1976': 1732.180,
                         'p53tet': 9387.583, 'streptavidin': 19044.982, 'hammerhead': 13868.580, 'NCP': 86372.086,
                         'pl_complex': 9678.118, 'hairpin': 4410.125, 'BSA': 104827.703}[_] for _ in my_sys_names]

    #
    #
    # rb_sa_vals = [{'EDTA_Mg': 449.840784, 'cambrin': 3172.407028, '1BNA': 4633.077172, '181L': 9118.281616,
    #                'DB1976': 683.894436, 'p53tet': 8613.472816, 'streptavidin': 18483.418668,
    #                'hammerhead': 10496.329656, 'NCP': 80066.345331, 'pl_complex': 8770.004628, 'hairpin': 3742.351508,
    #                'BSA': 96167.953541}[_] for _ in my_sys_names]
    #
    # rb_vol_vals = [{'EDTA_Mg': 271.496, 'cambrin': 5814.376, '1BNA': 6716.032, '181L': 22995.384, 'DB1976': 369.488,
    #                 'p53tet': 18731.216, 'streptavidin': 63989.128, 'hammerhead': 18939.936, 'NCP': 228134.552,
    #                 'pl_complex': 23138.856, 'hairpin': 5132.768, 'BSA': 339343.08}[_] for _ in my_sys_names]
    rb_sa_vals = [{'EDTA_Mg': 433.839, 'cambrin': 3125.924744, '1BNA': 4616.918620, '181L': 9176.314456,
                   'DB1976': 662.105696, 'p53tet': 8743.906396, 'streptavidin': 18822.882868,
                   'hammerhead': 10531.117348, 'NCP': 82304.614883, 'pl_complex': 8809.698816, 'hairpin': 3709.187016,
                   'BSA': 99064.141949}[_] for _ in my_sys_names]

    rb_vol_vals = [{'EDTA_Mg': 270.376, 'cambrin': 5548.656, '1BNA': 6627.4, '181L': 22658.264, 'DB1976': 365.328,
                    'p53tet': 18364.912, 'streptavidin': 62884.944, 'hammerhead': 18625.552, 'NCP': 223831.0240,
                    'pl_complex': 22842.08, 'hairpin': 5059.176, 'BSA': 332509.296}[_] for _ in my_sys_names]
    if diff == 'abs':
        pow_diff = [100 * abs(vor_vals[i] - pow_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
        del_diff = [100 * abs(vor_vals[i] - del_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
        if plotting == 'SA':
            rb_diff = [100 * abs(vor_vals[i] - rb_sa_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
        elif plotting == 'Vol':
            rb_diff = [100 * abs(vor_vals[i] - rb_vol_vals[i])/vor_vals[i] for i in range(len(vor_vals))]
    else:
        pow_diff = [round(100 * (pow_vals[i] - vor_vals[i])/vor_vals[i], 3) for i in range(len(vor_vals))]
        del_diff = [round(100 * (del_vals[i] - vor_vals[i])/vor_vals[i], 3) for i in range(len(vor_vals))]
        if plotting == 'SA':
            rb_diff = [round(100 * (rb_sa_vals[i] - vor_vals[i])/vor_vals[i], 3) for i in range(len(vor_vals))]
        elif plotting == 'Vol':
            rb_diff = [round(100 * (rb_vol_vals[i] - vor_vals[i])/vor_vals[i], 3) for i in range(len(vor_vals))]

    # Set the label codes
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J', '1BNA': '1BNA', 'DB1976': 'DB1976'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    def sort_3_lists(lista, listb, listc, listd):
        # Zipping lists together and sorting by the first list
        sorted_lists = sorted(zip(lista, listb, listc, listd), key=lambda x: x[0])

        # Unpacking the sorted lists
        lista, listb, listc, listd = zip(*sorted_lists)

        # Converting tuples back to lists if needed
        lista = list(lista)
        listb = list(listb)
        listc = list(listc)
        listd = list(listd)

        # Return the lists
        return lista, listb, listc, listd

    new_graph_labels, pow_diff, del_diff, rb_diff = sort_3_lists(new_graph_labels, pow_diff, del_diff, rb_diff)

    # Create the
    bar([pow_diff, del_diff], x_names=new_graph_labels, legend_names=['Power', 'Primitive', 'Solvent Accessible'],
        Show=True, y_axis_title='Difference', x_axis_title='Model', print_vals_on_bars=False,
        legend_orientation='Vertical', xlabel_size=30, ylabel_size=30, tick_width=2, tick_length=12,
        xtick_label_size=25, ytick_label_size=25)

