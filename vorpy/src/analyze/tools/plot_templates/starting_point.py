import os
import tkinter as tk
from tkinter import filedialog
from vorpy.src.system.system import System
from vorpy.src.group.group import Group


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    drop_box_folder = filedialog.askdirectory()
    folder = drop_box_folder + '/Jack/Vorpy/Data/IV_Molecular/logs_and_pdbs/'
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