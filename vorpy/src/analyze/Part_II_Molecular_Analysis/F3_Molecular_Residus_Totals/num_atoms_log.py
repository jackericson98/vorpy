import os
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
from vorpy.src.calculations.sorting import get_sys_type
from vorpy.src.analyze.tools.plot_templates.inline_plot import inline_plot


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
    num_atoms = [len(_.atoms) for _ in my_systems]
    systems = [x for _, x in sorted(zip(num_atoms, my_systems))]
    types = [get_sys_type(_) for _ in systems]
    # Plot the number of atoms inline
    inline_plot(num_atoms)



