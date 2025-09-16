import os
from os import path
import tkinter as tk
from tkinter import filedialog


import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
sys.path.append(project_root)


from vorpy.src.analyze.tools.CleanData.check_mol_data import check_mol_data


def print_run_statements(file_directory=None, python_name='python'):
    # Try to open up the foam_gen user_data file
    if file_directory is None:
        file_directory = filedialog.askdirectory(title='Get frame folder')
    else:
        file_directory = file_directory

    for folder in os.listdir(file_directory):
        my_missing_atoms = check_mol_data(file_directory + '/' + folder, print_statement=False)
        for key in ('aw', 'pow', 'prm'):
            if len(my_missing_atoms[key]) > 0:
                atoms_to_fix = " and a ".join([str(atom) for atom in my_missing_atoms[key]])
                print(f"{python_name} vorpy {file_directory + '/' + folder + '/' + folder + '.pdb'} -g a {atoms_to_fix} -s nt {key} -e dir {file_directory + '/' + folder}")


if __name__ == "__main__":
    print_run_statements()
