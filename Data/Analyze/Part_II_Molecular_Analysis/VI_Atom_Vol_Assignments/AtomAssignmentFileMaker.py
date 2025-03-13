import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import csv
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.chemistry_interpreter import amino_names, nucleo_names, ion_names
from matplotlib import colormaps


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
    resies = {}
    chains = {}
    group = []
    #  Go through the logs and get the name and the classifications
    for i, atom in logs_info.iterrows():
        # Add the atom residue to the resies
        resies[atom['Index']] = atom['Residue'], atom['Residue Sequence']
        # Add the chain to the atom chain assignments
        chains[atom['Index']] = atom['Chain']
        # Get the name and the classification
        if atom['Residue'] in amino_names:
            # Classify the atom as an amino acid atom
            classifs[atom['Index']] = 'aa'
            # Add the index to the group list
            group.append(atom['Index'])

        elif atom['Residue'] in nucleo_names:
            # Classify the atom as a nucleo atom
            classifs[atom['Index']] = 'na'
            # Add the index to the group list
            group.append(atom['Index'])
        else:
            # Classify the atom as sol
            classifs[atom['Index']] = 'ho'
    
    # Return the stuff
    return classifs, set(group), logs_info, chains, resies


def get_atoms_info(aw_logs, pow_logs, classifs, chains, resies):
    """
    Takes the logs from the aw and pow solve and returns for each atom:
    0. Index
    1. Ball radius
    2. aw volume
    3. pow volume
    4. aw sa
    5. pow sa
    6. association
    7. aw sphericity
    8. pow sphericity
    9. aw sol facing
    10. pow sol facing
    11. aw protein facing
    12. pow protein facing
    13. aw nucleic facing
    14. pow nucleic facing
    15. aw separate chain interface
    16. pow separate chain interface
    17. aw separate residue interface
    18. pow separate residue interface
    """
    # Get the power logs that are the atoms
    pow_logs_a = pow_logs['atoms']
    # Create the atoms dictionary
    atoms = {}
    # Now loop through and classify the atoms in the group as
    for i, aw_atom in aw_logs.iterrows():
        # Get the index so we dont have to keep referencing it
        ndx = aw_atom['Index']
        # Get the power atom
        pow_atom = pow_logs_a.loc[pow_logs_a['Index'] == ndx].iloc[0]
        # Get the aw neighbors and the pow neighbors
        aw_nbors_assns = []
        pow_nbors_assns = []
        for nbor in aw_atom['Neighbors']:
            try:    
                aw_nbors_assns.append(classifs[nbor])
            except KeyError:
                aw_nbors_assns.append('ho')
        for nbor in pow_atom['Neighbors']:
            try:
                pow_nbors_assns.append(classifs[nbor])
            except KeyError:
                pow_nbors_assns.append('ho')
        # Get the chain and residue information
        chain = chains[ndx]
        res = resies[ndx]
        # Create the dictionary
        atoms[aw_atom['Index']] = {
            'Index': aw_atom['Index'],
            'aw rad': aw_atom['Radius'],
            'aw vol': aw_atom['Volume'],
            'pow vol': pow_atom['Volume'],
            'aw sa': aw_atom['Surface Area'],
            'pow sa': pow_atom['Surface Area'],
            'association': classifs[aw_atom['Index']],
            'aw sphericity': aw_atom['Sphericity'],
            'pow sphericity': pow_atom['Sphericity'],
            'aw sol facing': any([_ == 'ho' for _ in aw_nbors_assns]),
            'pow sol facing': any([_ == 'ho' for _ in pow_nbors_assns]),
            'aw aa facing': any([_ == 'aa' for _ in aw_nbors_assns]),
            'pow aa facing': any([_ == 'aa' for _ in pow_nbors_assns]),
            'aw na facing': any([_ == 'na' for _ in aw_nbors_assns]),
            'pow na facing': any([_ == 'na' for _ in aw_nbors_assns]),
            'aw sep chain iface': any([chain != chains.get(_, chain) for _ in aw_atom['Neighbors']]),
            'pow sep chain iface': any([chain != chains.get(_, chain) for _ in pow_atom['Neighbors']]),
            'aw sep res iface': any([res != resies.get(_, res) for _ in aw_atom['Neighbors']]),
            'pow sep res iface': any([res != resies.get(_, res) for _ in pow_atom['Neighbors']])
        }
    return atoms


def get_rad_vals(aw_logs=None, pow_logs=None, output_folder=None, write_csv=False):
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
    classifications, group, aw_logs_info, chains, resies = get_group_info(aw_logs)
    # Get the power logs read for interpreting later
    pow_logs_info = read_logs2(pow_logs)
    # Get the atom dictionary
    my_dict = get_atoms_info(aw_logs_info, pow_logs_info, classifications, chains, resies)
    # Get information on the whole molecule

    # If write csv is chosen do it in the output folder
    if write_csv:
        with open(output_folder + '/atomic_comparisons.csv', 'w') as writing_file:
            wc = csv.writer(writing_file, lineterminator='\n')
            wc.writerow(['Index', 'aw rad', 'aw vol', 'pow vol', 'aw sa', 'pow sa', 'association', 'aw sphericity',
                         'pow sphericity', 'aw sol facing', 'pow sol facing', 'aw aa facing', 'pow aa facing',
                         'aw na facing', 'pow na facing', 'aw sep chain iface', 'pow sep chain iface',
                         'aw sep res iface', 'pow sep res iface'])
            for spleesh in my_dict:
                wc.writerow(my_dict[spleesh].values())

    # return the dictionary
    return my_dict


def bool_assign(val):
    return val.lower() == 'True'


def get_dict_from_file(file):
    # Create the dictionary for the values
    my_dict = {}
    # Create the list of dictionary terms
    my_vals = ['Index', 'aw rad', 'aw vol', 'pow vol', 'aw sa', 'pow sa', 'association', 'aw sphericity',
               'pow sphericity', 'aw sol facing', 'pow sol facing', 'aw aa facing', 'pow aa facing', 'aw na facing',
               'pow na facing', 'aw sep chain iface', 'pow sep chain iface', 'aw sep res iface',
               'pow sep res iface']
    # Create the assignments for the type of values that we are gonna get from the thing that we read
    my_ass = [int, float, float, float, float, float, str, float, float, bool_assign, bool_assign, bool_assign,
              bool_assign, bool_assign, bool_assign, bool_assign, bool_assign, bool_assign, bool_assign]
    # Open the file
    with open(file, 'r') as reading_file:
        rf = csv.reader(reading_file)
        counter = 0
        for line in rf:
            if counter == 0:
                counter += 1
                continue
            my_dict[int(line[0])] = {my_vals[i]: my_ass[i](line[i]) for i in range(len(line))}
    return my_dict


def plot_my_stuff(dict_file=None):
    """
    Plots my stuff

    """
    # Check if a file is given to us
    if dict_file is not None:
        my_dict = get_dict_from_file(dict_file)
    # If no file is specified
    else:
        my_dict = get_rad_vals(write_csv=True)
    # Get the color maps
    my_cmap = plt.cm.viridis
    my_max_guy = max([my_dict[_]['aw vol'] for _ in my_dict])
    colors = my_cmap([__ / my_max_guy for __ in [my_dict[_]['aw vol'] for _ in my_dict]])
    # Now that we have everything information wise, we need to plot the stuff in a way that is good and we like to do
    for i, entry in enumerate(my_dict):
        plt.scatter([my_dict[entry]['aw rad']], 
                    [(my_dict[entry]['pow sphericity'] - my_dict[entry]['aw sphericity']) / my_dict[entry]['aw sphericity']], 
                    marker='x' if my_dict[entry]['aw sol facing'] else 'o', c=colors[i])
    plt.show()
    return my_dict


def plot_atom_sphericity_diff_by_classification(dict_file=None, my_dict=None):
    """
    Plots the atom's sphericity difference with the aw volume and colors based on if it is outside facing or not
    """
    # Check if a file is given to us
    if dict_file is not None:
        my_dict = get_dict_from_file(dict_file)
    # If no file is specified
    else:
        my_dict = get_rad_vals(write_csv=True)
    xs, ys, classifs = [], [], []
    for atom in my_dict:
        xs.append(my_dict[atom]['aw vol'])
        ys.append((my_dict[atom]['pow sphericity'] - my_dict[atom]['aw sphericity']) / my_dict[atom]['aw sphericity'])
        classifs.append(my_dict[atom]['association'])
    plt.scatter(xs, ys, c=classifs)
    plt.show()


def plot_atom_vol_by_radius(dict_file=None):
    """
    Plots the atom's volume by the radius
    """
    # Check if a file is given to us
    if dict_file is not None:
        my_dict = get_dict_from_file(dict_file)
    # If no file is specified
    else:
        my_dict = get_rad_vals(write_csv=True)
    xs, ys = [], []
    for atom in my_dict:
        xs.append(my_dict[atom]['aw rad'])
        ys.append(my_dict[atom]['aw vol'])
    plt.scatter(xs, ys)
    plt.show()  


def plot_atom_by_ho_facing(dict_file=None):
    """
    Plots the atom's volume by the radius
    """
    # Check if a file is given to us
    if dict_file is not None:
        my_dict = get_dict_from_file(dict_file) 
    else:
        my_dict = get_rad_vals(write_csv=True)
    xs, ys = [], []
    for atom in my_dict:
        xs.append(my_dict[atom]['aw rad'])
        ys.append(my_dict[atom]['sol facing'])
    plt.scatter(xs, ys)
    plt.show()



if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    # my_file = filedialog.askopenfilename(title="Get CSV File")
    my_file = None
    plot_my_stuff(dict_file=my_file)
    # plot_atom_sphericity_diff_by_classification()
    # plot_atom_vol_by_radius()
    # plot_atom_by_ho_facing()
