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
    "ala": {'type': 'NP', 'letter':  "A", 'weight': 89.09},  # Alanine
    "arg": {'type': '+', 'letter':  "R", 'weight': 174.20},  # Arginine
    "asn": {'type': 'P', 'letter':  "N", 'weight': 132.12},  # Asparagine
    "asp": {'type': '-', 'letter':  "D", 'weight': 133.10},  # Aspartic acid
    "cys": {'type': 'NP', 'letter':  "C", 'weight': 121.16},  # Cysteine
    "gln": {'type': 'P', 'letter':  "Q", 'weight': 146.15},  # Glutamine
    "glu": {'type': '-', 'letter':  "E", 'weight': 147.13},  # Glutamic acid
    "gly": {'type': 'NP', 'letter':  "G", 'weight': 75.07},  # Glycine
    "his": {'type': '+', 'letter':  "H", 'weight': 155.16},  # Histidine
    "ile": {'type': 'NP', 'letter':  "I", 'weight': 131.17},  # Isoleucine
    "leu": {'type': 'NP', 'letter':  "L", 'weight': 131.17},  # Leucine
    "lys": {'type': '+', 'letter':  "K", 'weight': 146.19},  # Lysine
    "met": {'type': 'NP', 'letter':  "M", 'weight': 149.21},  # Methionine
    "phe": {'type': 'NP', 'letter':  "F", 'weight': 165.19},  # Phenylalanine
    "pro": {'type': 'NP', 'letter':  "P", 'weight': 115.13},  # Proline
    "ser": {'type': 'P', 'letter':  "S", 'weight': 105.09},  # Serine
    "thr": {'type': 'P', 'letter':  "T", 'weight': 119.12},  # Threonine
    "trp": {'type': 'NP', 'letter':  "W", 'weight': 204.22},  # Tryptophan
    "tyr": {'type': 'P', 'letter':  "Y", 'weight': 181.19},  # Tyrosine
    "val": {'type': 'NP', 'letter':  "V", 'weight': 117.15}  # Valine
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
    all_res_types, all_vor_res_types, all_pow_res_types, all_del_res_types = [], [], [], []
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
        vor_vols, pow_vols, del_vols, vor_sas, pow_sas, del_sas = [], [], [], [], [], []
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
        res_types, vor_res_types, pow_res_types, del_res_types = [], [], [], []

        for _ in vor_reses:
            # Sub class level
            for __ in vor_reses[_]:
                # Res_seq level
                for ___ in vor_reses[_][__]:
                    if vor_reses[_][__][___] == {}:
                        continue
                    if vor_reses[_][__][___]['vol'] == 0 or vor_reses[_][__][___]['sa'] == 0:
                        continue
                    vor_val = vor_reses[_][__][___]['vol']
                    pow_val = pow_reses[_][__][___]['vol']
                    del_val = del_reses[_][__][___]['vol']
                    maxy, miny = max(maxy, vor_val, pow_val, del_val), min(miny, vor_val, pow_val, del_val)
                    vor_res_types.append(vor_val)
                    pow_res_types.append(pow_val)
                    del_res_types.append(del_val)

                    res_types.append(__)

        # Add the atom types
        all_res_types.append(res_types)
        # Add the pow and del atom_type diffs
        all_vor_res_types.append(vor_res_types)
        all_pow_res_types.append(pow_res_types)
        all_del_res_types.append(del_res_types)
    # Return the values
    return all_vor_res_types, all_pow_res_types, all_del_res_types, all_res_types, miny, maxy


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
    print([_.name for _ in systems])
    # Sort atoms by number of atoms
    num_atoms = [len(_.atoms) for _ in systems]
    systems = [x for _, x in sorted(zip(num_atoms, systems))]
    # Create the logs dictionary
    my_sys_names = [__.name for __ in systems]

    # Create the log file name dictionary
    my_log_files = {_: {__: folder + '/' + _ + '_{}_logs.csv'.format(__) for __ in {'vor', 'pow', 'del'}}
                    for _ in my_sys_names}

    all_vor_diffs, all_pow_diffs, all_del_diffs, all_res_types, miny, maxy = residue_per_diff(systems, my_log_files)

    # Create the labels manually for the systems in question
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP', 'pl_complex': 'Prot-Lig',
                     '1BNA': '1BNA', 'DB1976': 'DB1976', 'BSA': 'BSA'}[_] for _ in my_sys_names]
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J', '1BNA': 'K', 'DB1976': 'L'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    new_graph_labels, all_vor_diffs, all_pow_diffs, all_del_diffs, all_res_types = sort_lists(new_graph_labels, all_vor_diffs, all_pow_diffs, all_del_diffs, all_res_types)

    color_dict = {'NP': 'green', '+': 'red', 'P': 'cyan', '-': 'k'}
    
    num_bins = 250
    y_diff = maxy - miny
    y_step = y_diff / num_bins
    # Set up the value, number and element lists for each scheme
    vor_vals, pow_vals, del_vals, vor_nums, pow_nums, del_nums, vor_ress, pow_ress, del_ress, vor_colors, pow_colors, del_colors = [], [], [], [], [], [], [], [], [], [], [], []
    # We are looking for each system to have a set of bins for power, primitive containing color, location & marker size
    for i in range(len(all_res_types)):
        # Get the number of atoms so we know how to scale everything
        my_num_ress = len(all_res_types[i])
        # Set up the power and del bins
        vor_bins = [{_: [] for _ in amino_dict} for _ in range(num_bins + 1)]
        pow_bins = [{_: [] for _ in amino_dict} for _ in range(num_bins + 1)]
        del_bins = [{_: [] for _ in amino_dict} for _ in range(num_bins + 1)]
        # Go through each bin
        for j in range(len(all_res_types[i])):
            # Get the element, the power value and the del value
            residue, vor_val, pow_val, del_val = all_res_types[i][j], all_vor_diffs[i][j], all_pow_diffs[i][j], all_del_diffs[i][j]
            # Assign the power bins
            try:
                vor_bins[int((vor_val - miny) // y_step)][residue].append(vor_val)
                pow_bins[int((pow_val - miny) // y_step)][residue].append(pow_val)
                del_bins[int((del_val - miny) // y_step)][residue].append(del_val)
            except KeyError:
                print("No color for {}. Leaving out".format(residue))
        # Set up the value, number and element lists for each scheme
        (sys_vor_vals, sys_pow_vals, sys_del_vals, sys_vor_nums, sys_pow_nums, sys_del_nums, sys_vor_ress,
         sys_pow_ress, sys_del_ress, sys_vor_colors, sys_pow_colors, sys_del_colors) = (
            [], [], [], [], [], [], [], [], [], [], [], [])
        # Now go through the bins taking the average for each element
        for j in range(num_bins):
            for residue in vor_bins[j]:
                # Check if it is empty
                if len(vor_bins[j][residue])/my_num_ress > 0.001:
                    # Get the average for each bin
                    sys_vor_vals.append(np.mean(vor_bins[j][residue]))
                    sys_vor_nums.append(10000 * len(vor_bins[j][residue]) / my_num_ress)
                    sys_vor_ress.append(residue)
                    print(residue)
                    sys_vor_colors.append(color_dict[amino_dict[residue]['type']])
                if len(pow_bins[j][residue])/my_num_ress > 0.001:
                    # Get the average for each bin
                    sys_pow_vals.append(np.mean(pow_bins[j][residue]))
                    sys_pow_nums.append(10000 * len(pow_bins[j][residue]) / my_num_ress)
                    sys_pow_ress.append(residue)
                    sys_pow_colors.append(color_dict[amino_dict[residue]['type']])
                # Check if the del is empty
                if len(del_bins[j][residue])/my_num_ress > 0.001:
                    # Get the average for each bin
                    sys_del_vals.append(np.mean(del_bins[j][residue]))
                    sys_del_nums.append(10000 * len(del_bins[j][residue]) / my_num_ress)
                    sys_del_ress.append(residue)
                    sys_del_colors.append(color_dict[amino_dict[residue]['type']])
        vor_vals.append(sys_vor_vals)
        pow_vals.append(sys_pow_vals)
        del_vals.append(sys_del_vals)
        vor_nums.append(sys_vor_nums)
        pow_nums.append(sys_pow_nums)
        del_nums.append(sys_del_nums)
        vor_ress.append(sys_vor_ress)
        pow_ress.append(sys_pow_ress)
        del_ress.append(sys_del_ress)
        vor_colors.append(sys_vor_colors)
        pow_colors.append(sys_pow_colors)
        del_colors.append(sys_del_colors)

    # Plot the data
    scatter(xs=[[4 * i + 1 for _ in vor_vals[i]] for i in range(len(vor_vals))] +
               [[4 * i + 2 for _ in pow_vals[i]] for i in range(len(pow_vals))] +
               [[4 * i + 3 for _ in del_vals[i]] for i in range(len(del_vals))],
            ys=vor_vals + pow_vals + del_vals, Show=False,
            alpha=0.3,
            colors=vor_colors + pow_colors + del_colors,
            x_tick_labels=new_graph_labels,
            x_tick_label_locs=[4 * i + 2 for i in range(len(all_pow_diffs))],
            xlabel_size=30,
            x_axis_title='Model',
            ylabel_size=30,
            y_axis_title='% Deviation',
            marker_size=vor_nums + pow_nums + del_nums,
            markers=[['${}$'.format(amino_dict[_]['letter']) for _ in __] for __ in vor_ress] +
                    [['${}$'.format(amino_dict[_]['letter']) for _ in __] for __ in pow_ress] +
                    [['${}$'.format(amino_dict[_]['letter']) for _ in __] for __ in del_ress])
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