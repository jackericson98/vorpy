import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.compare.get_res_data import nucleics, proteins, ions, sols
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.plot_templates.bar import bar

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
        vor_surfs = read_logs(folder + '/' + my_sys.name + '_vor_logs.csv', return_dict=True)['surfs']
        pow_surfs = {str(_['atoms'][0]) + '_' + str(_['atoms'][1]): _
                     for _ in read_logs(folder + '/' + my_sys.name + '_pow_logs.csv', True, True)['surfs']}
        del_surfs = {str(_['atoms'][0]) + '_' + str(_['atoms'][1]): _
                     for _ in read_logs(folder + '/' + my_sys.name + '_del_logs.csv', True, True)['surfs']}
        # Create the surface dictionary for later reference
        pow_diffs, del_diffs = {}, {}
        # Sort the surfaces
        for i, surf in enumerate(vor_surfs):
            if vor_surfs[i]['sa'] == 0:
                continue
            val = 'sa'
            # Get the hash string
            hsh_str = str(vor_surfs[i]['atoms'][0]) + '_' + str(vor_surfs[i]['atoms'][1])
            # Classification
            a1_res, a2_res = my_sys.atoms['res'].iloc[surf['atoms'][0]], my_sys.atoms['res'].iloc[surf['atoms'][1]]
            sclass = classify_surf(a1_res, a2_res)
            #
            if hsh_str in pow_surfs:
                if sclass['csf'] in pow_diffs:
                    pow_diffs[sclass['csf']].append(abs(vor_surfs[i][val] - pow_surfs[hsh_str][val])/vor_surfs[i][val])
                else:
                    pow_diffs[sclass['csf']] = [abs(vor_surfs[i][val] - pow_surfs[hsh_str][val]) / vor_surfs[i][val]]
            if hsh_str in del_surfs:
                if sclass['csf'] in del_diffs:
                    del_diffs[sclass['csf']].append(abs(vor_surfs[i][val] - del_surfs[hsh_str][val])/vor_surfs[i][val])
                else:
                    del_diffs[sclass['csf']] = [abs(vor_surfs[i][val] - del_surfs[hsh_str][val]) / vor_surfs[i][val]]

        my_sys_csfs.append((pow_diffs, del_diffs))
    # Create the labels manually for the systems in question
    sys_names = [_.name for _ in systems]
    graph_labels = []
    for name in sys_names:
        try:
            graph_labels.append({'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                                 'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP',
                                 'pl_complex': 'Prot-Lig', '1BNA': '1BNA', 'DB1976': 'DB1976', '181L': '181L'}[
                                    name])
        except KeyError:
            graph_labels.append(name)
    # Loop through the comparisons
    comparisons = ['Intra-Nucleic', 'Intra-Protein', 'Protein-Protein', 'Nucleic-Nucleic', 'Nucleic-Protein',
                   'Protein-SOL', 'Nucleic-SOL']
    for comparison in comparisons:
        # Create the list for system names and data
        new_sys_names, pow_data, del_data, pow_ses, del_ses = [], [], [], [], []
        # Get the data for each system
        for i, my_sys in enumerate(graph_labels):
            if comparison in my_sys_csfs[i][0]:
                new_sys_names.append(my_sys)
                pow_data.append(sum(my_sys_csfs[i][0][comparison]) / len(my_sys_csfs[i][0][comparison]))
                del_data.append(sum(my_sys_csfs[i][1][comparison]) / len(my_sys_csfs[i][1][comparison]))
                pow_ses.append(np.std(my_sys_csfs[i][0][comparison]) / np.sqrt(len(my_sys_csfs[i][0][comparison])))
                del_ses.append(np.std(my_sys_csfs[i][1][comparison]) / np.sqrt(len(my_sys_csfs[i][1][comparison])))
        bar([pow_data, del_data], errors=[pow_ses, del_ses], legend_names=['Power', 'Primitive'], x_names=new_sys_names, Show=True, title='Average Surface Area % Diff - ' + comparison,
            y_axis_title='% Difference', x_axis_title='Model')


