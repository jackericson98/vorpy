from System.sys_funcs.calcs.sorting import get_radius
from System.sys_objs.atom import make_atom
from System.Network.network import Network
from System.sys_objs.residue import Residue
from System.sys_objs.chain import Chain, Sol
import os.path as path
import csv
import numpy as np
from pandas import DataFrame



# Read mol method. Interprets the data from a .mol file type
def read_mol(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'mol':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    atoms = []
    # Go through the lines in the file
    for i in range(len(my_file)):
        # Get the line
        line = my_file[i]
        # If the line is an atom line add the data
        if len(line) > 6:
            # Add the data
            atoms.append(make_atom([line[0], line[1], line[2]], element=line[3], index=i))
    sys.atoms = DataFrame(atoms)