"""
Outputs a list plot of the different atomic curvature assignments

"""

import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.compare.read_logs import read_logs

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

for _ in surf_type_dict:
    if sum(surf_type_dict[_]) != 0 and len(surf_type_dict[_]) > 10:
        print(_, surf_type_dict[_])
