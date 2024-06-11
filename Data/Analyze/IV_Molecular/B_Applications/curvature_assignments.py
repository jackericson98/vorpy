"""
Outputs a list plot of the different atomic curvature assignments

"""

import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.compare.read_logs import read_logs
import matplotlib.pyplot as plt

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_logs = filedialog.askopenfilename()

my_logs_info = read_logs(my_logs)

surf_type_dict = {}
for i, surf in my_logs_info['surfs'].iterrows():
    atom_indices = [int(_) for _ in list(surf['atoms'])]
    atom0 = my_logs_info['atoms'].loc[my_logs_info['atoms']['num'] == atom_indices[0]].iloc[0]
    atom1 = my_logs_info['atoms'].loc[my_logs_info['atoms']['num'] == atom_indices[1]].iloc[0]
    atom_names = [atom0['name'].strip(), atom1['name'].strip()]
    atom_names.sort()
    combined_names = '_'.join(atom_names)
    if combined_names in surf_type_dict:
        surf_type_dict[combined_names].append(surf['curvature'])
    else:
        surf_type_dict[combined_names] = [surf['curvature']]

new_surf_dict = {}
for _ in surf_type_dict:
    if sum(surf_type_dict[_]) / len(surf_type_dict[_]) >= 0.01 and len(surf_type_dict[_]) > 10:
        new_surf_dict[_] = surf_type_dict[_]

# Prepare data for plotting
labels, values = zip(*new_surf_dict.items())

# Create the boxplot
fig, ax = plt.subplots(figsize=(12, 8))
ax.boxplot(values, labels=labels, patch_artist=True)

# Set plot title and labels
ax.set_title('Distribution of Curvature Types')
ax.set_xlabel('Curvature Type')
ax.set_ylabel('Curvature Value')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)
plt.tight_layout()

# Display the plot
plt.show()
