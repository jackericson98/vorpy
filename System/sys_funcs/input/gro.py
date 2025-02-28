import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.sys_objs.atom import make_atom
from pandas import DataFrame


# Read gro method. Interprets the data from a .cif file type
def read_gro1(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'gro':
        file = sys.base_file
    atoms = []
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
        for i, line in enumerate(my_file):
            line.split()
            if 2 <= i < len(my_file):
                atoms.append(make_atom(location=[line[3], line[4], line[5]], system=sys, index=i, name=line[1]))
    sys.atoms = DataFrame(atoms)


def read_gro(sys, file=None):
    """
    Read GROMACS file

    Processes the GROMACS file format and returns a ball dataframe that matches the standard vorpy balls dataframe
    format.

    """
    # Get the file if the file is not specified
    if file is None:
        file = sys.files['base_file']
    # Create the dictionary that holds the balls and the additional information
    file_dict = {'balls': [], 'Additional Lines': []}
    # Line splits
    line_splits = [0, 5, 8, 15, 20, 28, 36, 44]
    val_types = [int, str, str, int, float, float, float]
    line_vals = ['res_seq', 'res_name', 'atom_name', 'index', 'x', 'y', 'z']
    # Open the file
    with open(file, 'r') as read_file:
        # Loop through the lines
        for line in read_file.readlines():
            try:
                # Split the line into its constituent parts
                ball = {line_vals[j]: val_types[j](line[line_splits[j]: line_splits[j + 1]].strip()) for j in range(7)}
                # Make an atom
                ball = make_atom(sys, location=np.array([ball['x'], ball['y'], ball['z']]), index=ball['index'],
                                 name=ball['atom_name'], res_name=ball['res_name'])
                # Add the atom to the list
                file_dict['balls'].append(ball)

            except ValueError:
                file_dict['Additional Lines'].append(line)
    # Add the information
    sys.balls, sys.data = DataFrame(file_dict['balls']), file_dict['Additional Lines']
    sys.chains, sys.residues = [], []


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    my_file = filedialog.askopenfilename()
    read_gro(sys=None, file=my_file)

