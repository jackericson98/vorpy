import os

import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog
from vorpy.src.system.system import System
from Data.Analyze.tools.compare.read_logs import read_logs
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.box_whisker import box_whisker
from Data.Analyze.tools.plot_templates.line import line_plot

combos = {'HG22': 'HG', 'HG21': 'HG', 'HG23': 'HG', 'HD1': 'HD', 'HD2': 'HD', 'HB1': 'HB', 'HG11': 'HG', 'CG1': 'CG',
          'CG2': 'CG', 'HG12': 'HG', 'HG13': 'HG', 'HW1': 'HW', 'HW2': 'HW', 'HH12': 'HH1', 'HH11': 'HH1',
          'HH21': 'HH2', 'HH22': 'HH2', 'HB3': 'HB', 'HB2': 'HB', 'HD1': 'HD', 'HD2': 'HD', 'HD3': 'HD', 'HG2': 'HG',
          'HG1': 'HG', 'HA1': 'HA', 'HA2': 'HA', 'H5\'\'': 'H5\'', 'H2\'\'': 'H2\'', 'H22': 'H21', 'H62': 'H61',
          'H42': 'H41', 'O1P': 'OP', 'O2P': 'OP', 'HD12': 'HD11', 'HD13': 'HD11', 'HD22': 'HD21', 'HD23': 'HD21',
          'H5\'2': 'H5\'1', 'H2\'2': 'H2\'1', 'H72': 'H71', 'H73': 'H71'}


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
    surfs = []
    for sys in systems:
        sys_name = sys.name
        logs = read_logs(my_logs[sys_name]['vor'])
        for i, surf in logs['surfs'].iterrows():
            atom_indices = [int(_) for _ in list(surf['atoms'])]
            try:
                atom0 = logs['atoms'].loc[logs['atoms']['num'] == atom_indices[0]].iloc[0]
                atom1 = logs['atoms'].loc[logs['atoms']['num'] == atom_indices[1]].iloc[0]
            except IndexError:
                continue
            atom_names = [atom0['name'].strip(), atom1['name'].strip()]

            for i, atom in enumerate(atom_names):
                if atom in combos:
                    atom_names[i] = combos[atom]
            if 'H' in atom_names[0]:
                atom_names = [atom_names[1], atom_names[0]]
            elif 'H' not in atom_names[1]:
                atom_names.sort()
            combined_names = ' - '.join(atom_names)
            surfs.append({'name': combined_names, 'curvature': surf['curvature'], 'system': sys.name})

    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig',
                     'BSA': 'BSA', '1BNA': '1BNA'}[_] for _ in my_sys_names]

    # Bin the values
    increments = np.linspace(0, max([_['curvature'] for _ in surfs]), 300)

    # Set up the x values
    # Set the label codes
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J', '1BNA': 'K'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]
    scaled_bins = []
    all_sections = [0.01, 0.066, 0.131, 0.148, 0.195, 0.227, 0.342, 0.386, 0.427, 0.473, 0.967, 1.082, 1.145, 2]
    nucleic_sections = [0.01, 0.069, 0.112, 0.145, 0.171, 0.195, 0.226, 0.397, 0.455, 0.983, 1.18, 2.000]
    protein1_sections = [0.01, 0.020, 0.030, 0.045, 0.059, 0.068, 0.081, 0.091, 0.12, 0.2, 0.227, 0.342, 0.4, 0.5, 0.967, 1.09, 2]
    protein2_sections = [0.01, 0.020, 0.030, 0.040, 0.059, 0.069, 0.091, 0.099, 0.12, 0.2, 0.227, 0.342, 0.4, 0.5, 0.967, 1.09, 1.123, 2]
    sections = [(protein1_sections[i], protein1_sections[i+1]) for i in range(len(protein1_sections) - 1)]
    surf_name_list = [{} for _ in range(len(sections))]
    for surf in surfs:
        for i, section in enumerate(sections):
            if section[0] <= surf['curvature'] < section[1]:
                if surf['name'] in surf_name_list[i]:
                    surf_name_list[i][surf['name']]['count'] += 1
                    surf_name_list[i][surf['name']]['curvatures'].append(surf['curvature'])
                else:
                    surf_name_list[i][surf['name']] = {'count': 1, 'curvatures': [surf['curvature']]}
    for i, _ in enumerate(surf_name_list):
        print(sections[i], [(__, _[__]['count'], round(np.mean(_[__]['curvatures']), 3), round(np.std(_[__]['curvatures']), 4)) for __ in _ if _[__]['count'] >= 10])
    # for _ in vals[1:]:
    #
    #     # Scale the number of values by the size of the group
    #     bins = [[__ for __ in _ if increments[j] <= __ < increments[j + 1]] for j in range(len(increments[:-1]))]
    #     # Count the number based on the number of total values
    #     scaled_bins.append([len(__) / len(_) for __ in bins])

    # sorted_lists = sorted(zip(new_graph_labels[1:], scaled_bins), key=lambda x: x[0])
    #
    # # Unpacking the sorted lists
    # lista, listb = zip(*sorted_lists)
    # line_plot([[(increments[j] + increments[j + 1]) / 2 for j in range(len(increments[:-1]))][5:] for _ in range(len(scaled_bins))],
    #           [_[5:] for _ in listb], title='Surface Curvature Distributions By Model',
    #           x_label='Curvature', y_label='Distribution', legend_title='Model', labels=lista,
    #           title_size=25, x_label_size=20, y_label_size=20)
    #
    # plt.show()

