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
from Data.Analyze.tools.calcs import sort_lists

amino_bbs = ['CA', 'HA', 'HA1', 'HA2', 'N', 'HN', 'H', 'C', 'O', 'OC1', 'OC2', 'OT1', 'OT2', 'H1', 'H2', 'H3']
amino_scs = ['CB', 'HB', 'HB1', 'HB2', 'HB3',
             'SD', 'CD', 'CD1', 'CD2', 'ND1', 'ND2', 'OD1', 'OD2', 'HD1', 'HD2', 'HD3', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'
             , 'CE', 'CE1', 'CE2', 'CE3', 'OE1', 'OE2', 'NE', 'NE1', 'NE2', 'HE', 'HE1', 'HE2', 'HE3', 'HE21', 'HE22',
             'CG', 'CG1', 'CG2', 'OG', 'SG', 'OG1', 'HG', 'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23',
             'CH2', 'NH1', 'OH', 'HH', 'HH1', 'HH2', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
             'NZ', 'CZ', 'CZ1', 'CZ2', 'CZ3', 'NZ', 'HZ', 'HZ1', 'HZ2', 'HZ3']

nucleic_acids = {'DT', 'DA', 'DG', 'DC', 'DU', 'U', 'G', 'A', 'T', 'C', 'GDP', 'OMC'}

nucleic_nbase = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'O2', 'O4',
                 'O6', 'H8', 'H1', 'H2', 'H21', 'H22', 'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07', 'H08', 'H09',
                 'H10', 'H11', 'H12', 'H61', 'H62', 'H13', 'H5', 'H6', 'H3', 'H71', 'H72', 'H73', 'H41', 'H42']
nucleic_sugr = ['O3\'', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\'', 'O2\'', 'CM2', 'H5\'1', 'H5\'2', 'H4\'',
                'H3\'', 'H2\'', 'H2\'1', 'H1\'', 'H2\'1', 'H2\'2']
nucleic_pphte = ['P', 'O1P', 'O2P', 'OP1', 'OP2', 'PA', 'PB', 'O1A', 'O1B', 'O2A', 'O2B', 'O3A', 'O3B', 'H5T', 'H3T']


def atoms_per_diff(systems, logs, val='volume'):
    # Go through the loaded systems
    all_atom_types, all_pow_atom_types, all_del_atom_types = [], [], []
    # Set up the min and max values for diff
    miny, maxy = 0, 0
    atoms_not_in = []
    for system in systems:
        sys_name = system.name
        vor_atoms = read_logs(logs[sys_name]['vor'], True, True)['atoms']
        pow_atoms = read_logs(logs[sys_name]['pow'], True, True)['atoms']
        del_atoms = read_logs(logs[sys_name]['del'], True, True)['atoms']
        # print(sys_name, len(vor_atoms), len(pow_atoms), len(del_atoms))
        atom_types, pow_atom_types, del_atom_types = [], [], []
        for i, atom in enumerate(vor_atoms):

            sys_atom = system.atoms.loc[system.atoms['num'] == atom['num']].to_dict(orient='records')[0]

            # if sys_atom['name'].strip() not in amino_bbs + amino_scs + ['HW1', 'HW2', 'OW']:
            #     if sys_atom['name'].strip() not in atoms_not_in:
            #         atoms_not_in.append(sys_atom['name'].strip())
            # if sys_atom['name'].strip() not in amino_scs:
            #     continue
            if sys_atom['name'] in ['HW1', 'HW2', 'OW', 'Na']:
                continue
            # if sys_atom['element'].strip().lower() == 'p':
            #     print(sys_atom['name'].strip())

            pow_val = 100 * (pow_atoms[i][val] - vor_atoms[i][val]) / vor_atoms[i][val]
            del_val = 100 * (del_atoms[i][val] - vor_atoms[i][val]) / vor_atoms[i][val]
            maxy, miny = max(maxy, pow_val, del_val), min(miny, pow_val, del_val)
            pow_atom_types.append(pow_val)
            del_atom_types.append(del_val)

            atom_types.append(sys_atom['element'])

        # Add the atom types
        all_atom_types.append(atom_types)
        # Add the pow and del atom_type diffs
        all_pow_atom_types.append(pow_atom_types)
        all_del_atom_types.append(del_atom_types)
    print(system.name)
    print(atoms_not_in)
    # Return the values
    return all_pow_atom_types, all_del_atom_types, all_atom_types, miny, maxy


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
            if file[-3:] == 'csv' and 'res' not in file:
                my_logs[file[:-13]][file[-12:-9]] = folder + '/' + file
    # Get the average differences and the standard errors
    all_pow_diffs, all_del_diffs, all_atom_types, miny, maxy = atoms_per_diff(systems, my_logs)
    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig',
                     '1BNA': '1BNA', 'DB1976': 'DB1976', 'BSA': 'BSA'}[_] for _ in my_sys_names]
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J', '1BNA': 'K', 'DB1976': 'L'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    new_graph_labels, all_pow_diffs, all_del_diffs, all_atom_types = sort_lists(new_graph_labels, all_pow_diffs, all_del_diffs, all_atom_types)

    # Plot the Data
    # bar(data=[avg_pow_diffs, avg_del_diffs], legend_title='Scheme',
    #     y_axis_title='% Difference', x_names=new_graph_labels, legend_names=['Power', 'Primitive'], Show=True,
    #     x_axis_title='Model', errors=[pow_ses, del_ses], y_range=[0, 70], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)

    color_dict = {'C': 'grey', 'O': 'r', 'N': 'b', 'P': 'darkorange', 'H': 'pink', 'S': 'y', 'Se': 'sandybrown'}
    num_bins = 250
    y_diff = maxy - miny
    y_step = y_diff / num_bins
    # Set up the value, number and element lists for each scheme
    pow_vals, del_vals, pow_nums, del_nums, pow_elems, del_elems, pow_colors, del_colors = [], [], [], [], [], [], [], []
    # We are looking for each system to have a set of bins for power, primitive containing color, location & marker size
    for i in range(len(all_atom_types)):
        # Get the number of atoms so we know how to scale everything
        my_num_atoms = len(all_atom_types[i])
        # Set up the power and del bins
        pow_bins = [{'C': [], 'H': [], 'P': [], 'S': [], 'Se': [], 'O': [], 'N': []} for _ in range(num_bins + 1)]
        del_bins = [{'C': [], 'H': [], 'P': [], 'S': [], 'Se': [], 'O': [], 'N': []} for _ in range(num_bins + 1)]
        # Go through each bin
        for j in range(len(all_atom_types[i])):
            # Get the element, the power value and the del value
            element, pow_val, del_val = all_atom_types[i][j], all_pow_diffs[i][j], all_del_diffs[i][j]
            # Assign the power bins
            try:
                pow_bins[int((pow_val - miny) // y_step)][element].append(pow_val)
                del_bins[int((del_val - miny) // y_step)][element].append(del_val)
            except KeyError:
                print("No color for {}. Leaving out".format(element))
        # Set up the value, number and element lists for each scheme
        sys_pow_vals, sys_del_vals, sys_pow_nums, sys_del_nums, sys_pow_elems, sys_del_elems, sys_pow_colors, sys_del_colors = [], [], [], [], [], [], [], []
        # Now go through the bins taking the average for each element
        for j in range(num_bins):
            for elem in pow_bins[j]:
                # Check if it is empty

                if len(pow_bins[j][elem])/my_num_atoms > 0.001:
                    # Get the average for each bin
                    sys_pow_vals.append(np.mean(pow_bins[j][elem]))
                    sys_pow_nums.append(10000 * len(pow_bins[j][elem]) / my_num_atoms)
                    sys_pow_elems.append(elem)
                    sys_pow_colors.append(color_dict[elem])
                # Check if the del is empty
                if len(del_bins[j][elem])/my_num_atoms > 0.001:
                    # Get the average for each bin
                    sys_del_vals.append(np.mean(del_bins[j][elem]))
                    sys_del_nums.append(10000 * len(del_bins[j][elem]) / my_num_atoms)
                    sys_del_elems.append(elem)
                    sys_del_colors.append(color_dict[elem])
        pow_vals.append(sys_pow_vals)
        del_vals.append(sys_del_vals)
        pow_nums.append(sys_pow_nums)
        del_nums.append(sys_del_nums)
        pow_elems.append(sys_pow_elems)
        del_elems.append(sys_del_elems)
        pow_colors.append(sys_pow_colors)
        del_colors.append(sys_del_colors)
    # Plot the data
    # scatter(xs=[[3 * i + 1 for _ in pow_vals[i]] for i in range(len(pow_vals))] +
    #            [[3 * i + 2 for _ in del_vals[i]] for i in range(len(del_vals))],
    #         ys=pow_vals + del_vals, Show=False,
    #         markers='o',
    #         alpha=0.3,
    #         colors=pow_colors + del_colors,
    #         x_tick_labels=new_graph_labels,
    #         x_tick_label_locs=[3 * i + 1.5 for i in range(len(all_pow_diffs))],
    #         xlabel_size=30,
    #         x_axis_title='Model',
    #         ylabel_size=30,
    #         y_axis_title='% Deviation',
    #         marker_size=pow_nums + del_nums,
    #         y_range=[-150, 100])
    scatter(xs=[[i for _ in del_vals[i]] for i in range(len(del_vals))],
            ys=del_vals, Show=False,
            markers='o',
            alpha=0.3,
            colors=del_colors,
            x_tick_labels=new_graph_labels,
            x_tick_label_locs=[i for i in range(len(all_del_diffs))],
            xlabel_size=30,
            x_axis_title='Model',
            ylabel_size=30,
            y_axis_title='% Deviation',
            marker_size=del_nums,
            y_range=[-120, 100])
    # plt.plot([1], [1], label='Power - Left')
    #
    # plt.plot([1], [1], label='Primitive - Right')
    #
    # legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2,
    #                     prop={'size': 25})

    plt.show()
# [['${}$'.format(_) for _ in sys_elems] for sys_elems in pow_elems] +
#                     [['${}$'.format(_) for _ in sys_elems] for sys_elems in del_elems]
