import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.compare.read_logs2 import read_logs2


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)


def verify_pdb(file=None, print_metadata=False, simple=False):
    """
    Verifies the contents of a pdb file and returns the metadata of the file
    """
    if file is None:
        file = filedialog.askopenfilename(title="Open PDB File")

    # If the verification is simple we just need to check that the file is longer than 0
    if simple:
        with open(file, 'r') as my_file:
            if len(my_file.readlines()) > 0:
                return True
            return False

    # Go through the lines collecting data
    with open(file, 'r') as my_file:
        # Create the dictionary for storing the pdb information
        pdb_dict = {'': file, 'atoms': 0, 'remarks': 0, 'other': 0}
        for line in my_file.readlines():
            # Split the line information for harvesting
            line_info = line.split()
            # Check the line information
            if line_info[0] == 'ATOM':
                pdb_dict['atoms'] += 1
            elif line_info[0] == 'REMARK':
                pdb_dict['remarks'] += 1
            else:
                pdb_dict['other'] += 1
    if print_metadata:
        print(pdb_dict)
    return pdb_dict


def verify_logs(file=None, print_metadata=False, simple=False):
    """
    Reads the logs and verifies that they 1. Exist and 2. Are complete
    """
    if file is None:
        file = filedialog.askopenfilename(title="Get Logs File")

    # If the verification is simple we just need to check that the file is longer than 0
    if simple:
        with open(file, 'r') as my_file:
            if len(my_file.readlines()) > 0:
                return True
            return False

    # We need to read the logs
    my_file_logs = read_logs2()
