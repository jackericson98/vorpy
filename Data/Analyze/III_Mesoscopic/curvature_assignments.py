"""
Outputs a list plot of the different atomic curvature assignments

"""

import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.compare.read_logs import read_logs
import matplotlib.pyplot as plt
import numpy as np
from System.system import System


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    my_logs_info = read_logs(filedialog.askopenfilename(title='Choose Logs'))
    # Assign the atoms in the logs to the specific aspects of the pdb
    my_coarse_sys = System(file=filedialog.askopenfilename(title='Choose Coarse PDB'))


    surf_type_dict = {}
    for i, surf in my_logs_info['surfs'].iterrows():
        atom_indices = [int(_) for _ in list(surf['atoms'])]
        try:
            atom0 = my_logs_info['atoms'].loc[my_logs_info['atoms']['num'] == atom_indices[0]].iloc[0]
            atom1 = my_logs_info['atoms'].loc[my_logs_info['atoms']['num'] == atom_indices[1]].iloc[0]
            sys_a0 = my_coarse_sys.atoms.loc[my_coarse_sys.atoms['num'] == atom_indices[0]].iloc[0]
            sys_a1 = my_coarse_sys.atoms.loc[my_coarse_sys.atoms['num'] == atom_indices[1]].iloc[0]
        except IndexError:
            continue
        if 'sc' in my_coarse_sys.name:
            seq_dic = {0: ' BB', 1: ' SC'}
            atom_names = [sys_a1['res'].name + ' SC' if atom1['name'].strip() == 'pb' else sys_a1['res'].name + ' BB',
                          sys_a0['res'].name + ' SC' if atom0['name'].strip() == 'pb' else sys_a0['res'].name + ' BB']
        else:
            atom_names = [sys_a1['res'].name, sys_a0['res'].name]
        if 'H' in atom_names[0] or 'SOL' in atom_names[0]:
            continue
            atom_names = [atom_names[1], atom_names[0]]
        elif 'H' not in atom_names[1]:
            atom_names.sort()
        combined_names = ' - '.join(atom_names)
        if combined_names in surf_type_dict:
            surf_type_dict[combined_names].append(surf['curvature'])
        else:
            surf_type_dict[combined_names] = [surf['curvature']]

    new_surf_dict = {}
    new_surf_dict1 = {}
    for _ in surf_type_dict:
        curv_avg = sum(surf_type_dict[_]) / len(surf_type_dict[_])
        if len(surf_type_dict[_]) > 10:
            # Sort the outliars: Get the mean and standard deviation
            my_mean, my_std = np.mean(surf_type_dict[_]), np.std(surf_type_dict[_])
            # Filter out the outliars (2 stds)

            new_surf_dict[_] = [__ for __ in surf_type_dict[_] if abs(my_mean - __) < 1 * my_std]
            # print(_, my_mean, my_std, new_surf_dict[_])
            if len(new_surf_dict[_]) > 5:
                new_surf_dict1[_] = new_surf_dict[_]

    surf_dict = dict(sorted(new_surf_dict1.items(), key=lambda item: np.mean(item[1]), reverse=True))

    # Prepare data for plotting
    labels, values = zip(*surf_dict.items())

    # Create the boxplot
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.boxplot(values[:min(35, len(values))], labels=labels[:min(35, len(values))], patch_artist=True)

    # Set plot title and labels
    ax.set_title('Distribution of Curvatures ({})'.format(my_logs_info['data']['name'].capitalize()), fontdict=dict(size=30))
    ax.set_xlabel('Surface Type', fontdict=dict(size=30))
    ax.set_ylabel('Curvature', fontdict=dict(size=30))

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, size=20)
    plt.yticks(size=20)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.tick_params(axis='both', labelsize=20, width=2, length=15)

    # Display the plot
    plt.show()
