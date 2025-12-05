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
from vorpy.src.analyze.tools.plot_templates.table import table
from vorpy.src.calculations.sorting import get_sys_type


def get_system_data(systems):
    # Make the titles for the table
    columns = ['Name', 'Type', '# Atoms', '# Residues', '# Chains', '# SOL Atoms']
    # Instantiate the list for the system data
    sys_data = []
    # Go through the list of systems
    for sys in systems:
        # Create the system dictionary to add stuff to
        sys_data.append([sys.name, get_sys_type(sys), len(sys.groups[0].atoms), len(sys.residues), len(sys.chains),
                         len(sys.sol.atoms)])

    # Create the table
    table(column_names=columns, rows=sys_data, Show=True, transpose=True)


if __name__ == '__main__':
    # Get the Dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    # Go through the folder and collect the information
    my_systems = []
    for root, directory, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'pdb':
                my_sys = System(file=folder + '/' + file)
                my_sys.groups = [Group(sys=my_sys, residues=my_sys.residues)]
                my_systems.append(my_sys)
    # Sort the systems by the number of atoms in the main group
    num_atoms = [len(_.groups[0].atoms) for _ in my_systems]
    my_systems = [x for _, x in sorted(zip(num_atoms, my_systems))]
    # Create the table
    get_system_data(my_systems)

