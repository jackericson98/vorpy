import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs import read_logs

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_logs = filedialog.askopenfilename(title='Choose Logs')

my_logs_info = read_logs(my_logs)

combos = {'HG22': 'HG', 'HG21': 'HG', 'HG23': 'HG', 'HD1': 'HD', 'HD2': 'HD', 'HB1': 'HB', 'HG11': 'HG', 'CG1': 'CG',
          'CG2': 'CG', 'HG12': 'HG', 'HG13': 'HG', 'HW1': 'HW', 'HW2': 'HW', 'HH12': 'HH1', 'HH11': 'HH1',
          'HH21': 'HH2', 'HH22': 'HH2', 'HB3': 'HB', 'HB2': 'HB', 'HD1': 'HD', 'HD2': 'HD', 'HD3': 'HD', 'HG2': 'HG',
          'HG1': 'HG', 'HA1': 'HA', 'HA2': 'HA', 'H5\'\'': 'H5\'', 'H2\'\'': 'H2\'', 'H22': 'H21', 'H62': 'H61',
          'H42': 'H41', 'O1P': 'OP', 'O2P': 'OP', 'HD12': 'HD11', 'HD13': 'HD11', 'HD22': 'HD21', 'HD23': 'HD21',
          'H5\'2': 'H5\'1', 'H2\'2': 'H2\'1', 'H72': 'H71', 'H73': 'H71'}

surf_type_dict = {}
for i, surf in my_logs_info['surfs'].iterrows():
    atom_indices = [int(_) for _ in list(surf['atoms'])]
    try:
        atom0 = my_logs_info['atoms'].loc[my_logs_info['atoms']['num'] == atom_indices[0]].iloc[0]
        atom1 = my_logs_info['atoms'].loc[my_logs_info['atoms']['num'] == atom_indices[1]].iloc[0]
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

        new_surf_dict[_] = [__ for __ in surf_type_dict[_] if __ > my_mean]
        # print(_, my_mean, my_std, new_surf_dict[_])
        if len(new_surf_dict[_]) > 5:
            new_surf_dict1[_] = new_surf_dict[_]


surf_dict = dict(sorted(new_surf_dict1.items(), key=lambda item: np.mean(item[1]), reverse=True))


# Prepare data for plotting
labels, values = zip(*surf_dict.items())

# Create the boxplot
fig, ax = plt.subplots(figsize=(12, 8))
ax.boxplot(values[:min(15, len(values))], labels=labels[:min(15, len(values))], patch_artist=True)

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
