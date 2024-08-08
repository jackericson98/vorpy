from System.sys_funcs.calcs.sorting import get_radius
from System.sys_objs.atom import make_atom
from System.sys_objs.residue import Residue
from System.sys_objs.chain import Chain, Sol
import os.path as path
import csv
import numpy as np
from pandas import DataFrame


# Read cif function. Interprets the data in a cif file
def read_cif(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'cif':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Get the starting number for the line
    num = int(my_file[0].split()[0])
    # Go through each line of the file
    for i in range(len(my_file)):
        # Split the line
        line = my_file[i].split()
        # Add the atoms
        if line == int(num) and len(line) >= 7:
            sys.atoms.append(make_atom([line[9], line[10], line[11]], element=line[3], index=i))