import tkinter as tk
from tkinter import filedialog
import numpy as np
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.chemistry_interpreter import amino_names, nucleo_names, ion_names


"""
Takes in the pow and aw logs and gathers the radii and percent difference as well as sphericity and outputs the values 
in a csv file
"""

def get_group_info(logs):
    """
    Classifies the atoms into their set group
    """
    # Read the logs
    logs_info = read_logs2(logs, all_=False, balls=True)['atoms']
    # Create the atoms dictionary
    classifs = {}
    group = []
    #  Go through the logs and get the name and the classifications
    for i, atom in logs_info.iterrows():
        # Get the name and the classification
        if atom['Residue'] in amino_names:
            # Classify the atom as an amino acid atom
            classifs[atom['Index']] = 'aa'
            # Add the index to the group list
            group.append(atom['Index'])

        elif atom['Residue'] in nucleo_names:
            # Classify the atom as a nucleo atom
            classifs['Index'] = 'na'
            # Add the index to the group list
            group.append(atom['Index'])
        else:
            # Classify the atom as sol
            classifs[atom['Index']] = 'ho'
    # Return the stuff
    return classifs, set(group), logs


def get_atoms_info(aw_logs, pow_logs, classifs):
    """
    Takes the location of the logs
    """
    # Get the power logs that are the atoms
    pow_logs_a = pow_logs['atoms']
    # Create the atoms dictionary
    atoms = {}
    # Now loop through and classify the atoms in the group as
    for i, aw_atom in aw_logs.iterrows():
        # Get the power atom
        pow_atom = pow_logs_a.loc[pow_logs_a['']]
        # Create the dictionary
        atoms[]




    return atoms


def get_rad_vals(aw_logs=None, pow_logs=None, output_folder=None):
    # Get the aw logs if none are specified
    if aw_logs is None:
        # Ask for the aw logs
        aw_logs = filedialog.askopenfilename(title="Get AW Logs")
    # Get the pow logs if none are specified
    if pow_logs is None:
        # Ask for the pow logs
        pow_logs = filedialog.askopenfilename(title="Get POW Logs")
    # Get the output folder if one is not specified
    if output_folder is None:
        # Ask for the output folder
        output_folder = filedialog.askdirectory(title="Get Output Folder")


    # Create the aw_atoms dictionary
    aw_dict = get_atoms_info(aw_logs)
