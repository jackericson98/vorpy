import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.compare.residue import nucleics, proteins, ions, sols
from Data.Analyze.tools.compare.read_logs import read_logs

# We want a full list of all surfaces with SA, classification,
# Classes include: protein-protein, intra-protein, Protein-Ligand, Protein-Nucleic, intra-nucleic, Nucleic-Nucleic,
#                  Nucleic-Sol, Protein-Sol


def classify_surf(system, a1_res, a2_res):

    # Set up the surface info dictionary
    surf_info = {'res1': a1_res, 'res2': a2_res}
    # Classify the interaction
    if a1_res == a2_res:
        if a1_res.name in nucleics:
            surf_info['csf'] = 'Intra-Nucleic'
        elif a1_res.name in proteins:
            surf_info['csf'] = 'Intra-Protein'
        else:
            surf_info['csf'] = 'Intra-Other'
    elif a1_res.name in nucleics:
        if a2_res.name in proteins:
            surf_info['csf'] = 'Protein-Nucleic'
        elif a2_res.name in nucleics:
            surf_info['csf'] = 'Nucleic-Nucleic'
        elif a2_res.name in ions:
            surf_info['csf'] = 'Nucleic-Ion'
        elif a2_res.name in sols:
            surf_info['csf'] = 'Nucleic-Sol'
        else:
            surf_info['csf'] = 'Nucleic-Other'
    elif a1_res.name in proteins:
        if a2_res.name in proteins:
            surf_info['csf'] = 'Protein-Protein'
        elif a2_res.name in nucleics:
            surf_info['csf'] = 'Protein-Nucleic'
        elif a2_res.name in ions:
            surf_info['csf'] = 'Protein-Ion'
        elif a2_res.name in sols:
            surf_info['csf'] = 'Protein-Sol'
        else:
            surf_info['csf'] = 'Protein-Other'

    return surf_info


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    # Create the systems
    systems = []
    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'pdb':
                my_sys = System(file=folder + '/' + file)
                my_sys.groups = [Group(sys=my_sys, residues=my_sys.residues)]
                systems.append(my_sys)
    # Sort atoms by number of atoms
    num_atoms = [len(_.atoms) for _ in systems]
    systems = [x for _, x in sorted(zip(num_atoms, systems))]

    # Set the output folder
    output_folder = filedialog.askdirectory() + '/'
    # Create the outputs by system
    my_maxes = []
    for my_sys in systems:
        # Read the logs
        my_log_vals = read_logs(folder + '/' + my_sys.name + '_vor_logs.csv', return_dict=True)
