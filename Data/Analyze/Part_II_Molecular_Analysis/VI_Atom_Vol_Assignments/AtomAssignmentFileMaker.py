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
from System.sys_objs.atom import get_element
from matplotlib import colormaps


"""
Takes in the pow and aw logs and gathers the radii and percent difference as well as sphericity and outputs the values 
in a csv file
"""


residue_colors = {
    # Amino Acids
    "ALA": '#8CFF8C',
    "GLY": '#FFFFFF', 
    "LEU": '#455E45',
    "SER": '#FF7042',
    "VAL": '#FF8CFF',
    "THR": '#B84C00',
    "LYS": '#4747B8',
    "ASP": '#A00042',
    "ILE": '#004C00',
    "ASN": '#FF7C70',
    "GLU": '#660000',
    "PRO": '#525252',
    "ARG": '#00007C',
    "PHE": '#534C42',
    "GLN": '#FF4C4C',
    "TYR": '#8C704C',
    "HIS": '#7070FF',
    "CYS": '#FFFF70',
    "MET": '#B8A042',
    "TRP": '#4F4600',
    "ASX": '#FF00FF',  # Assuming ASX and GLX represent Asp/Asn and Glu/Gln ambiguous cases
    "GLX": '#FF00FF',
    "PCA": '#FF00FF',  # Rare in standard use, included for completeness
    "HYP": '#FF00FF',  # Rare in standard use, included for completeness
    "GDP": '#FFD700',  # Gold color for GDP
    "OMC": '#E6B800',  # Dark yellow color for OMC
    "JZ4": '#FFD700',

    # Nucleic Acids
    **{_: '#A0A0FF' for _ in {"DA", "A"}},
    **{_: '#FF8C4B' for _ in {"DC", "C"}},
    **{_: '#FF7070' for _ in {"DG", "G"}},
    **{_: '#A0FFA0' for _ in {"DT", "T"}},
    **{_: '#B8B8B8' for _ in {"DU", "U"}},

    # Other special cases
    "Backbone": '#B8B8B8',
    "Special": '#5E005E',
    "Default": '#FF00FF',
    'MOL': '#9370DB',  # Medium Purple
    # Sol colors
    "SOL": '#00FFFF',
    "HOH": '#00FFFF',
    # Single atom colors
    "MG": '#00FA6D',
    "C": '#C8C8C8',
    "O": '#F00000',
    "H": '#FFFFFF',
    "N": '#8F8FFF',
    "S": '#FFC832',
    "P": '#FFA500',
    "CL": '#00FF00',
    "BR": '#A52A2A',
    "ZN": '#A52A2A',
    "NA": '#0000FF',
    "FE": '#FFA500',
    "CA": '#808090'
}
amino_bbs = ['CA', 'HA', 'HA1', 'HA2', 'N', 'HN', 'H', 'C', 'O', 'OC1', 'OC2', 'OT1', 'OT2', 'H1', 'H2', 'H3']
amino_scs = ['CB', 'HB', 'HB1', 'HB2', 'HB3',
             'SD', 'CD', 'CD1', 'CD2', 'ND1', 'ND2', 'OD1', 'OD2', 'HD1', 'HD2', 'HD3', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'
             , 'CE', 'CE1', 'CE2', 'CE3', 'OE1', 'OE2', 'NE', 'NE1', 'NE2', 'HE', 'HE1', 'HE2', 'HE3', 'HE21', 'HE22',
             'CG', 'CG1', 'CG2', 'OG', 'SG', 'OG1', 'HG', 'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23',
             'CH2', 'NH1', 'OH', 'HH', 'HH1', 'HH2', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
             'NZ', 'CZ', 'CZ1', 'CZ2', 'CZ3', 'NZ', 'HZ', 'HZ1', 'HZ2', 'HZ3']

nucleic_nbase = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'O2', 'O4', 'O6']
nucleic_sugr = ['O3\'', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\'', 'O2\'', 'CM2']
nucleic_pphte = ['P', 'O1P', 'O2P', 'OP1', 'OP2', 'PA', 'PB', 'O1A', 'O1B', 'O2A', 'O2B', 'O3A', 'O3B']
bb_sc_colors = {**{_: 'r' for _ in amino_bbs}, **{_: 'y' for _ in amino_scs}, **{_: 'blue' for _ in nucleic_nbase}, **{_: 'purple' for _ in nucleic_sugr}, **{_: 'maroon' for _ in nucleic_pphte}}


def average_pairwise_distance(points):
    # Convert list of points to a numpy array for efficient computation
    points = np.array(points)

    # Calculate the pairwise distances between all points
    dist_matrix = np.sqrt(np.sum((points[:, np.newaxis, :] - points[np.newaxis, :, :]) ** 2, axis=-1))

    # Since the distance matrix is symmetric and the diagonal is 0, we can extract the upper triangle excluding the diagonal
    upper_triangle_indices = np.triu_indices_from(dist_matrix, k=1)
    pairwise_distances = dist_matrix[upper_triangle_indices]

    # Calculate the average of these distances
    average_distance = np.mean(pairwise_distances)

    return average_distance


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
        if atom['Residue'].lower() in amino_names:
            # Classify the atom as an amino acid atom
            classifs[atom['Index']] = 'aa'
            # Add the index to the group list
            group.append(atom['Index'])

        elif atom['Residue'].lower() in nucleo_names:
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
            'element': get_element(atom_name=aw_atom['Name']),
            'name': aw_atom['Name'],
            'residue': aw_atom['Residue'],
            'residue sequence': aw_atom['Residue Sequence'],
            'rad': aw_atom['Radius'],
            'vdw vol': aw_atom['Van Der Waals Volume'],
            'aw vol': aw_atom['Volume'],
            'pow vol': pow_atom['Volume'],
            'pow vol diff': (pow_atom['Volume'] - aw_atom['Volume']) / aw_atom['Volume'],
            'aw vol diff': (aw_atom['Volume'] - pow_atom['Volume']) / pow_atom['Volume'],
            'aw sa': aw_atom['Surface Area'],
            'pow sa': pow_atom['Surface Area'],
            'association': classifs[aw_atom['Index']],
            'aw sphericity': aw_atom['Sphericity'],
            'pow sphericity': pow_atom['Sphericity'],
            'pow sphereicity diff': (pow_atom['Sphericity'] - aw_atom['Sphericity']) / aw_atom['Sphericity'],
            'aw sphereicity diff': (aw_atom['Sphericity'] - pow_atom['Sphericity']) / pow_atom['Sphericity'],
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
            wc.writerow(['Index', 'element', 'name', 'residue', 'residue sequence', 'vdw vol', 'aw rad', 'aw vol',
                         'pow vol', 'pow vol diff', 'aw vol diff', 'aw sa', 'pow sa', 'association', 'aw sphericity',
                         'pow sphericity', 'pow sphereicity diff', 'aw sphereicity diff', 'aw sol facing',
                         'pow sol facing', 'aw aa facing', 'pow aa facing', 'aw na facing', 'pow na facing',
                         'aw sep chain iface', 'pow sep chain iface', 'aw sep res iface', 'pow sep res iface'])
            for spleesh in my_dict:
                wc.writerow(my_dict[spleesh].values())
        with open(output_folder + '/atomic_comparisons_simple.csv', 'w') as writing_file:
            wc = csv.writer(writing_file, lineterminator='\n')
            wc.writerow(['Index', 'name', 'residue', 'residue sequence', 'pow vol diff', 'aw vol diff',
                         'pow sphereicity diff', 'aw sphereicity diff', 'vdw vol'])
            for spleesh in my_dict:
                wc.writerow([my_dict[spleesh]['Index'], my_dict[spleesh]['name'], my_dict[spleesh]['residue'],
                             my_dict[spleesh]['residue sequence'], my_dict[spleesh]['pow vol diff'],
                             my_dict[spleesh]['aw vol diff'], my_dict[spleesh]['pow sphereicity diff'],
                             my_dict[spleesh]['aw sphereicity diff'], my_dict[spleesh]['vdw vol']])
        with open(output_folder + '/pow_model_training.csv', 'w') as writing_file:
            wc = csv.writer(writing_file, lineterminator='\n')
            wc.writerow(['pow vol diff', 'element', 'name', 'residue', 'rad', 'pow sa', 'association', 'pow sphericity',
                         'pow sol facing', 'pow aa facing', 'pow na facing', 'pow sep chain iface', 'pow sep res iface'])
            for spleesh in my_dict:
                wc.writerow([my_dict[spleesh]['pow vol diff'], my_dict[spleesh]['element'], my_dict[spleesh]['name'],
                             my_dict[spleesh]['residue'], my_dict[spleesh]['aw rad'], my_dict[spleesh]['pow sa'],
                             my_dict[spleesh]['association'], my_dict[spleesh]['pow sphereicity'],
                             my_dict[spleesh]['pow sol facing'], my_dict[spleesh]['pow aa facing'],
                             my_dict[spleesh]['pow na facing'], my_dict[spleesh]['pow sep chain iface'],
                             my_dict[spleesh]['pow sep res iface']])

    # return the dictionary
    return my_dict, output_folder


def bool_assign(val):
    return val.lower() == 'true'


def get_dict_from_file(file):
    # Create the dictionary for the values
    my_dict = {}
    # Create the list of dictionary terms
    my_vals = ['Index', 'element', 'name', 'residue', 'residue sequence', 'vdw vol', 'aw rad', 'aw vol',
               'pow vol', 'pow vol diff', 'aw vol diff', 'aw sa', 'pow sa', 'association', 'aw sphericity', 'pow sphericity', 'pow sphereicity diff', 'aw sphereicity diff',
               'aw sol facing', 'pow sol facing', 'aw aa facing', 'pow aa facing', 'aw na facing',
               'pow na facing', 'aw sep chain iface', 'pow sep chain iface', 'aw sep res iface',
               'pow sep res iface']
    # Create the assignments for the type of values that we are gonna get from the thing that we read
    my_ass = [int, str, str, str, str, float, float, float, float, float, float, float, float, str, float, float, float, float, bool_assign,
              bool_assign, bool_assign, bool_assign, bool_assign, bool_assign, bool_assign, bool_assign, bool_assign, bool_assign]
    # Open the file
    with open(file, 'r') as reading_file:
        rf = csv.reader(reading_file)
        counter = 0
        for line in rf:
            if counter == 0:
                counter += 1
                continue
            my_dict[int(line[0])] = {my_vals[i]: my_ass[i](line[i]) for i in range(len(line))}
    output_folder = os.path.dirname(file)
    return my_dict, output_folder


def plot_my_stuff(dict_file=None, file_name=''):
    """
    Plots my stuff

    """
    # Check if a file is given to us
    if dict_file is not None:
        my_dict, output_folder = get_dict_from_file(dict_file)
    # If no file is specified
    else:
        my_dict, output_folder = get_rad_vals(write_csv=True)
    # create the groupings
    carbon_groupings = {0.025: {}, 0.05: {}, 0.075: {}, 0.1: {}, 0.125: {}, 0.157: {}, 0.2: {}, 0.5: {}}
    nitrogen_groupings = {-0.075: {}, -0.0325: {}, 0.05: {}, 0.5: {}}
    oxygen_groupings = {-0.04: {}, -0.01: {}, 0.05: {}, 0.5: {}}

    def add_key_grouping(entry, dictionary, s_diff):
        # Record the values
        for grouping in dictionary:
            if s_diff < grouping:
                my_key = entry['name'] + '_' + entry['residue']
                if my_key in dictionary[grouping]:
                    dictionary[grouping][my_key] += 1
                else:
                    dictionary[grouping][my_key] = 1
                return dictionary
        return dictionary
    for color_by in ['element', 'residue', 'bb_sc']:
        if color_by == 'element':
            color_dict = {'C': 'grey', 'O': 'r', 'N': 'b', 'P': 'darkorange', 'H': 'pink', 'S': 'y', 'Se': 'sandybrown'}
            colors = [color_dict[my_dict[entry]['element']] for entry in my_dict]
        elif color_by == 'residue':
            colors = [residue_colors[my_dict[entry]['residue']] for entry in my_dict]
        elif color_by == 'bb_sc':
            colors = [bb_sc_colors[my_dict[entry]['name']] if my_dict[entry]['name'] in bb_sc_colors else 'black' for entry in my_dict]

        for x_val in ['sphericity', 'vdw_volume']:
            my_res_dict, my_res_dict_no_h = {}, {}
            for ass in ['na', 'aa']:
                # Now that we have everything information wise, we need to plot the stuff in a way that is good and we like to do
                for i, entry in enumerate(my_dict):
                    if my_dict[entry]['association'] != ass:
                        continue


                    if x_val == 'sphericity':
                        s_diff = (my_dict[entry]['pow sphericity'] - my_dict[entry]['aw sphericity']) / my_dict[entry]['aw sphericity']
                    elif x_val == 'vdw_volume':
                        s_diff = my_dict[entry]['vdw vol']
                    vol_diff = (my_dict[entry]['pow vol'] - my_dict[entry]['aw vol']) / my_dict[entry]['aw vol']
                    plt.scatter([s_diff], [vol_diff], marker='x' if my_dict[entry]['aw sol facing'] else 'o', c=colors[i],
                                alpha=0.2)

                    if my_dict[entry]['residue'] in my_res_dict:
                        my_res_dict[my_dict[entry]['residue']].append((s_diff, vol_diff))
                    else:
                        my_res_dict[my_dict[entry]['residue']] = [(s_diff, vol_diff)]
                    if my_dict[entry]['element'].lower() != 'h':
                        if my_dict[entry]['residue'] in my_res_dict_no_h:
                            my_res_dict_no_h[my_dict[entry]['residue']].append((s_diff, vol_diff))
                        else:
                            my_res_dict_no_h[my_dict[entry]['residue']] = [(s_diff, vol_diff)]

                    if my_dict[entry]['element'] == 'C':
                        carbon_groupings = add_key_grouping(my_dict[entry], carbon_groupings, s_diff)
                    if my_dict[entry]['element'] == 'N':
                        nitrogen_groupings = add_key_grouping(my_dict[entry], nitrogen_groupings, s_diff)
                    if my_dict[entry]['element'] == 'O':
                        oxygen_groupings = add_key_grouping(my_dict[entry], oxygen_groupings, s_diff)

                if x_val == 'sphericity':
                    plt.xlabel('Sphericity Difference')
                elif x_val == 'vdw_volume':
                    plt.xlabel('VDW Volume')
                plt.ylabel('Volume Difference')
                if color_by == 'element':
                    # Get unique elements that appear in the data
                    unique_elements = set(my_dict[entry]['element'] for entry in my_dict)
                    # Add legend entries only for elements that appear, using 'Other' for unknown elements
                    element_handles = []
                    for elem in unique_elements:
                        if elem in color_dict:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor=color_dict[elem], label=elem, markersize=8))
                        else:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='gray', label='Other', markersize=8))
                elif color_by == 'residue':
                    # Get unique residues that appear in the data
                    unique_residues = set(my_dict[entry]['residue'] for entry in my_dict)
                    # Add legend entries only for residues that appear, using 'Other' for unknown residues
                    element_handles = []
                    for res in unique_residues:
                        if res in residue_colors:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor=residue_colors[res], label=res, markersize=8))
                        else:
                            print('New residue type: ', res)
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='gray', label='Other', markersize=8))
                elif color_by == 'bb_sc':
                    # Get unique atom names that appear in the data
                    unique_names = set(my_dict[entry]['name'] for entry in my_dict)
                    # Track which groups we've added to avoid duplicates
                    added_groups = set()
                    element_handles = []

                    for name in unique_names:
                        if name in amino_bbs and 'Back Bone' not in added_groups:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='r', label='Back Bone', markersize=8))
                            added_groups.add('Back Bone')
                        elif name in amino_scs and 'Side Chain' not in added_groups:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='y', label='Side Chain', markersize=8))
                            added_groups.add('Side Chain')
                        elif name in nucleic_nbase and 'Nucleobase' not in added_groups:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='blue', label='Nucleobase', markersize=8))
                            added_groups.add('Nucleobase')
                        elif name in nucleic_pphte and 'Phosphate' not in added_groups:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='maroon', label='Phosphate', markersize=8))
                            added_groups.add('Phosphate')
                        elif name in nucleic_sugr and 'Sugar' not in added_groups:
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='purple', label='Sugar', markersize=8))
                            added_groups.add('Sugar')
                        elif name not in (amino_bbs + amino_scs + nucleic_nbase + nucleic_pphte + nucleic_sugr) and 'Other' not in added_groups:
                            print(f'Unidentified atom name {name}\n')
                            element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor='gray', label='Other', markersize=8))
                            added_groups.add('Other')
                        else:
                            print(f'Unidentified atom name {name}\n')

                # Get current figure and axes
                fig = plt.gcf()
                ax = plt.gca()

                # Adjust the plot area to make room for legend
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

                # Place legend to the right of the plot
                if len(element_handles) > 30:
                    ax.legend(handles=element_handles, loc='center left', bbox_to_anchor=(1, 0.5),
                             fontsize='xx-small', ncol=1)
                elif len(element_handles) > 20:
                    ax.legend(handles=element_handles, loc='center left', bbox_to_anchor=(1, 0.5),
                             fontsize='x-small', ncol=1)
                else:
                    ax.legend(handles=element_handles, loc='center left', bbox_to_anchor=(1, 0.5),
                             fontsize='small', ncol=1)
                # Add legend entries for markers
                marker_handles = [
                    plt.scatter([], [], marker='x', color='black', label='Sol Facing', s=75),
                    plt.scatter([], [], marker='o', color='black', label='Not Sol Facing', s=75)
                ]
                # Increase font size of axis labels, ticks and title
                plt.xlabel(plt.gca().get_xlabel(), fontsize=20)
                plt.ylabel(plt.gca().get_ylabel(), fontsize=20)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                # Create legend with two columns
                plt.legend(handles=element_handles + marker_handles, ncol=2)
                plt.title(f'{file_name} {ass}', fontsize=25)
                plt.tight_layout()
                plt.show()
            print("Average Pairwise Distance:\n")
            for res in my_res_dict:
                print(f"{res} = {average_pairwise_distance(my_res_dict[res])}")
            print("Average Pairwise Distance (no H):\n")
            for res in my_res_dict_no_h:
                print(f"{res} = {average_pairwise_distance(my_res_dict_no_h[res])}")
    return my_dict


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    # my_file = filedialog.askopenfilename(title="Get CSV File")
    my_file = None
    plot_my_stuff(dict_file=my_file, file_name='NCP')
