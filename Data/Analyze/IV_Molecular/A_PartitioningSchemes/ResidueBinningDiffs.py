import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.scatter import scatter
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.get_res_data import residue_data
from System.sys_funcs.calcs.sorting import sort_lists
import matplotlib.pyplot as plt

amino_dict = {
    "ala": {'type': 'NP', 'letter':  "A"},  # Alanine
    "arg": {'type': '+', 'letter':  "R"},  # Arginine
    "asn": {'type': 'P', 'letter':  "N"},  # Asparagine
    "asp": {'type': '-', 'letter':  "D"},  # Aspartic acid
    "cys": {'type': 'NP', 'letter':  "C"},  # Cysteine
    "gln": {'type': 'P', 'letter':  "Q"},  # Glutamine
    "glu": {'type': '-', 'letter':  "E"},  # Glutamic acid
    "gly": {'type': 'NP', 'letter':  "G"},  # Glycine
    "his": {'type': '+', 'letter':  "H"},  # Histidine
    "ile": {'type': 'NP', 'letter':  "I"},  # Isoleucine
    "leu": {'type': 'NP', 'letter':  "L"},  # Leucine
    "lys": {'type': '+', 'letter':  "K"},  # Lysine
    "met": {'type': 'NP', 'letter':  "M"},  # Methionine
    "phe": {'type': 'NP', 'letter':  "F"},  # Phenylalanine
    "pro": {'type': 'NP', 'letter':  "P"},  # Proline
    "ser": {'type': 'P', 'letter':  "S"},  # Serine
    "thr": {'type': 'P', 'letter':  "T"},  # Threonine
    "trp": {'type': 'NP', 'letter':  "W"},  # Tryptophan
    "tyr": {'type': 'P', 'letter':  "Y"},  # Tyrosine
    "val": {'type': 'NP', 'letter':  "V"}  # Valine
}
nucs_dict = {
    'dc': 'C',  # Cytosine
    'dg': 'G',  # Guanine
    'dt': 'T',  # Tyrosine
    'da': 'A',  # Adenine
    'du': 'U',  # Uracil
}


def residue_per_diff(systems, logs, val='volume'):
    # Go through the loaded systems
    all_res_types, all_pow_res_types, all_del_res_types = [], [], []
    # Set up the min and max values for diff
    miny, maxy = 0, 0
    atoms_not_in = []
    for system in systems:
        sys_name = system.name

        vor_out, vor_in = folder + '/res_data/{}_vor_res.csv'.format(system.name), None
        if os.path.exists(folder + '/res_data/{}_vor_res.csv'.format(system.name)):
            vor_in, vor_out = folder + '/res_data/{}_vor_res.csv'.format(system.name), None
        pow_out, pow_in = folder + '/res_data/{}_pow_res.csv'.format(system.name), None
        if os.path.exists(folder + '/res_data/{}_pow_res.csv'.format(system.name)):
            pow_in, pow_out = folder + '/res_data/{}_pow_res.csv'.format(system.name), None
        del_out, del_in = folder + '/res_data/{}_del_res.csv'.format(system.name), None
        if os.path.exists(folder + '/res_data/{}_del_res.csv'.format(system.name)):
            del_in, del_out = folder + '/res_data/{}_del_res.csv'.format(system.name), None
        # print('45: system = {}'.format(system.name))
        pow_vols, del_vols, pow_sas, del_sas = [], [], [], []
        # Get the values from the residue function
        vor_reses = residue_data(system, read_logs(my_log_files[system.name]['vor']), get_all=True, read_file=vor_in,
                                 output_file=vor_out)
        # print('49: vor_reses = {}'.format(vor_reses))
        pow_reses = residue_data(system, read_logs(my_log_files[system.name]['pow']), get_all=True, read_file=pow_in,
                                 output_file=pow_out)
        # print('51: pow_reses = {}'.format(pow_reses))
        del_reses = residue_data(system, read_logs(my_log_files[system.name]['del']), get_all=True, read_file=del_in,
                                 output_file=del_out)

        # print(sys_name, len(vor_atoms), len(pow_atoms), len(del_atoms))
        res_types, pow_res_types, del_res_types = [], [], []

        for _ in vor_reses:
            # Sub class level
            for __ in vor_reses[_]:
                # Res_seq level
                for ___ in vor_reses[_][__]:
                    if vor_reses[_][__][___] == {}:
                        continue
                    if vor_reses[_][__][___]['vol'] == 0 or vor_reses[_][__][___]['sa'] == 0:
                        continue

                    pow_val = 100 * (pow_reses[_][__][___]['vol'] - vor_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol']
                    del_val = 100 * (del_reses[_][__][___]['vol'] - vor_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol']
                    maxy, miny = max(maxy, pow_val, del_val), min(miny, pow_val, del_val)
                    pow_res_types.append(pow_val)
                    del_res_types.append(del_val)

                    res_types.append(__)

        # Add the atom types
        all_res_types.append(res_types)
        # Add the pow and del atom_type diffs
        all_pow_res_types.append(pow_res_types)
        all_del_res_types.append(del_res_types)
    # Return the values
    return all_pow_res_types, all_del_res_types, all_res_types, miny, maxy



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
    my_log_files = {_: {__: folder + '/' + _ + '_{}_logs.csv'.format(__) for __ in {'vor', 'pow', 'del'}}
                    for _ in my_sys_names}

    all_pow_diffs, all_del_diffs, all_res_types, miny, maxy = residue_per_diff(systems, my_log_files)

    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig',
                     '1BNA': '1BNA', 'DB1976': 'DB1976', 'BSA': 'BSA'}[_] for _ in my_sys_names]
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J', '1BNA': 'K', 'DB1976': 'L'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    new_graph_labels, all_pow_diffs, all_del_diffs, all_res_types = sort_lists(new_graph_labels, all_pow_diffs, all_del_diffs, all_res_types)

    color_dict = {'NP': 'green', '+': 'red', 'P': 'cyan', '-': 'k'}
    
    num_bins = 250
    y_diff = maxy - miny
    y_step = y_diff / num_bins
    # Set up the value, number and element lists for each scheme
    pow_vals, del_vals, pow_nums, del_nums, pow_types, del_types, pow_colors, del_colors = [], [], [], [], [], [], [], []
    # We are looking for each system to have a set of bins for power, primitive containing color, location & marker size
    for i in range(len(all_res_types)):
        # Get the number of atoms so we know how to scale everything
        my_num_ress = len(all_res_types[i])
        # Set up the power and del bins
        pow_bins = [{'NP': [], '+': [], 'P': [], '-': []} for _ in range(num_bins + 1)]
        del_bins = [{'NP': [], '+': [], 'P': [], '-': []} for _ in range(num_bins + 1)]
        # Go through each bin
        for j in range(len(all_res_types[i])):
            # Get the element, the power value and the del value
            residue, pow_val, del_val = all_res_types[i][j], all_pow_diffs[i][j], all_del_diffs[i][j]
            # Assign the power bins
            try:
                pow_bins[int((pow_val - miny) // y_step)][amino_dict[residue]['type']].append(pow_val)
                del_bins[int((del_val - miny) // y_step)][amino_dict[residue]['type']].append(del_val)
            except KeyError:
                print("No color for {}. Leaving out".format(residue))
        # Set up the value, number and element lists for each scheme
        sys_pow_vals, sys_del_vals, sys_pow_nums, sys_del_nums, sys_pow_types, sys_del_types, sys_pow_colors, sys_del_colors = [], [], [], [], [], [], [], []
        # Now go through the bins taking the average for each element
        for j in range(num_bins):
            for type in pow_bins[j]:
                # Check if it is empty

                if len(pow_bins[j][type])/my_num_ress > 0.001:
                    # Get the average for each bin
                    sys_pow_vals.append(np.mean(pow_bins[j][type]))
                    sys_pow_nums.append(10000 * len(pow_bins[j][type]) / my_num_ress)
                    sys_pow_types.append(type)
                    sys_pow_colors.append(color_dict[type])
                # Check if the del is empty
                if len(del_bins[j][type])/my_num_ress > 0.001:
                    # Get the average for each bin
                    sys_del_vals.append(np.mean(del_bins[j][type]))
                    sys_del_nums.append(10000 * len(del_bins[j][type]) / my_num_ress)
                    sys_del_types.append(type)
                    sys_del_colors.append(color_dict[type])
        pow_vals.append(sys_pow_vals)
        del_vals.append(sys_del_vals)
        pow_nums.append(sys_pow_nums)
        del_nums.append(sys_del_nums)
        pow_types.append(sys_pow_types)
        del_types.append(sys_del_types)
        pow_colors.append(sys_pow_colors)
        del_colors.append(sys_del_colors)

    # Plot the data
    scatter(xs=[[3 * i + 1 for _ in pow_vals[i]] for i in range(len(pow_vals))] +
               [[3 * i + 2 for _ in del_vals[i]] for i in range(len(del_vals))],
            ys=pow_vals + del_vals, Show=False,
            markers='o',
            alpha=0.3,
            colors=pow_colors + del_colors,
            x_tick_labels=new_graph_labels,
            x_tick_label_locs=[3 * i + 1.5 for i in range(len(all_pow_diffs))],
            xlabel_size=30,
            x_axis_title='Model',
            ylabel_size=30,
            y_axis_title='% Deviation',
            marker_size=pow_nums + del_nums,
            y_range=[-10, 10])
    # scatter(xs=[[i for _ in del_vals[i]] for i in range(len(del_vals))],
    #         ys=del_vals, Show=False,
    #         markers='o',
    #         alpha=0.3,
    #         colors=del_colors,
    #         x_tick_labels=new_graph_labels,
    #         x_tick_label_locs=[i for i in range(len(all_del_diffs))],
    #         xlabel_size=30,
    #         x_axis_title='Model',
    #         ylabel_size=30,
    #         y_axis_title='% Deviation',
    #         marker_size=del_nums,
    #         y_range=[-120, 100])
    # plt.plot([1], [1], label='Power - Left')
    #
    # plt.plot([1], [1], label='Primitive - Right')
    #
    # legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2,
    #                     prop={'size': 25})

    plt.show()
# [['${}$'.format(_) for _ in sys_elems] for sys_elems in pow_elems] +
#                     [['${}$'.format(_) for _ in sys_elems] for sys_elems in del_elems]