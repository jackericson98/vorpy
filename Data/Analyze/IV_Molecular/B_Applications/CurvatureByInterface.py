import os
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.compare.residue import nucleics, proteins, ions, sols
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.plot_templates.box_whisker import box_whisker

# We want a full list of all surfaces with SA, classification,
# Classes include: protein-protein, intra-protein, Protein-Ligand, Protein-Nucleic, intra-nucleic, Nucleic-Nucleic,
#                  Nucleic-Sol, Protein-Sol


def classify_surf(a1_res, a2_res):

    # Set up the surface info dictionary
    surf_info = {'res1': a1_res, 'res2': a2_res}
    # Classify the interaction
    if a1_res == a2_res:
        name = a1_res.name.strip().lower()
        if name in nucleics:
            surf_info['csf'] = 'Intra-Nucleic'
        elif name in proteins:
            surf_info['csf'] = 'Intra-Protein'
        else:
            surf_info['csf'] = 'Intra-Other'
        return surf_info
    names = []
    for name in {a1_res.name.strip().lower(), a2_res.name.strip().lower()}:
        if name in proteins:
            names.append('Protein')
        elif name in nucleics:
            names.append('Nucleic')
        elif name in ions:
            names.append('Ions')
        elif name in sols:
            names.append('SOL')
        else:
            names.append('Other')
    names.sort()
    surf_info['csf'] = '-'.join(names)

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

    # Create the outputs by system
    my_sys_csfs = []
    for my_sys in systems:
        # Read the logs
        my_log_vals = read_logs(folder + '/' + my_sys.name + '_vor_logs.csv', return_dict=True)
        # Create the surface dictionary for later reference
        sys_csf_dict = {}
        # Sort the surfaces
        for surf in my_log_vals['surfs']:
            # Classification
            a1_res, a2_res = my_sys.atoms['res'].iloc[surf['atoms'][0]], my_sys.atoms['res'].iloc[surf['atoms'][1]]
            sclass = classify_surf(a1_res, a2_res)
            # Sorting
            if sclass['csf'] in sys_csf_dict:
                sys_csf_dict[sclass['csf']].append(surf['curvature'])
            else:
                sys_csf_dict[sclass['csf']] = [surf['curvature']]
        my_sys_csfs.append(sys_csf_dict)
    # Create the labels manually for the systems in question
    sys_names = [_.name for _ in systems]
    graph_labels = []
    for name in sys_names:
        try:
            graph_labels.append({'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                                 'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP',
                                 'pl_complex': 'Prot-Lig', '1BNA': '1BNA', 'DB1976': 'DB1976', '181L': '181L'}[name])
        except KeyError:
            graph_labels.append(name)
    # Loop through the comparisons
    comparisons = ['Intra-Nucleic', 'Intra-Protein', 'Protein-Protein', 'Nucleic-Nucleic', 'Nucleic-Protein',
                   'Protein-SOL', 'Nucleic-SOL']
    for comparison in comparisons:
        # Create the list for system names and data
        new_sys_names, data = [], []
        # Get the data for each system
        for i, my_sys in enumerate(graph_labels):
            if comparison in my_sys_csfs[i]:
                new_sys_names.append(my_sys)
                data.append(my_sys_csfs[i][comparison])
        box_whisker(data, x_names=new_sys_names, Show=True, title='Average Curvature - ' + comparison,
                    y_axis_title='Curvature')
