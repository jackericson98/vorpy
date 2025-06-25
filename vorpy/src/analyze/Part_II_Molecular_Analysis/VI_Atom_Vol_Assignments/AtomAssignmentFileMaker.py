import os
import sys
import csv
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from scipy.spatial import distance_matrix, ConvexHull
from sklearn.metrics import silhouette_score
import matplotlib.patches as patches

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.chemistry.chemistry_interpreter import amino_names, nucleo_names, ion_names, sol_names
from vorpy.src.objects.atom import get_element
from vorpy.src.calculations.calcs import calc_com


def clustering_measures(full_coords, subsets, names, print_info=False):
    results = {}

    if len(full_coords) == 0:
        return
    total_volume = estimate_volume(full_coords)

    for i, subset_coords in enumerate(subsets):
        name = names[i]
        subset_volume = estimate_volume(subset_coords)
        density_ratio = (len(subset_coords) / subset_volume) / (
                    len(full_coords) / total_volume) if subset_volume > 0 else np.inf

        results[name] = density_ratio

    if print_info:
        for _ in results:
            if _ == 'Silhouette Score':
                continue
            print(_, results[_])

    return results


def estimate_volume(coords):

    if len(coords[0]) == 2:  # 2D case
        hull = ConvexHull(coords)
        return hull.area
    elif len(coords[0]) == 3:  # 3D case
        hull = ConvexHull(coords)
        return hull.volume
    return np.inf


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

nucleic_nbase = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'O2', 'O4',
                 'O6', 'H2', 'H21', 'H22', 'H3', 'H4', 'H41', 'H42', 'H5', 'H6', 'H61', 'H62', 'H8', 'H71', 'H72', 'H73']
nucleic_sugr = ['O3\'', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\'', 'O2\'', 'CM2', 'H1\'', 'H2\'', 'H2\'\'',
                'H3\'', 'H4\'', 'H5\'', 'H5\'\'', 'H3T', 'H5T', 'H2\'1', 'H2\'2', 'H5\'1', 'H5\'2']
nucleic_pphte = ['P', 'O1P', 'O2P', 'OP1', 'OP2', 'PA', 'PB', 'O1A', 'O1B', 'O2A', 'O2B', 'O3A', 'O3B']
bb_sc_colors = {**{_: 'r' for _ in amino_bbs}, **{_: 'y' for _ in amino_scs}, **{_: 'blue' for _ in nucleic_nbase},
                **{_: 'purple' for _ in nucleic_sugr}, **{_: 'maroon' for _ in nucleic_pphte}}


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
        elif atom['Residue'].lower() in sol_names:
            # Classify the atom as sol
            classifs[atom['Index']] = 'ho'
        else:
            # Classify the atom as other
            classifs[atom['Index']] = 'other'
            # print(f"Other Atom - Res: {atom['Residue']}, Name: {atom['Name']}")
            # Add the index to the group list
            group.append(atom['Index'])
    
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
            'pow nbors': pow_atom['Neighbors'],
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
    return my_dict, output_folder


def write_files(my_dict, output_folder=None, write_csv=False, ):
    # If write csv is chosen do it in the output folder
    if write_csv:
        with open(output_folder + '/atomic_comparisons.csv', 'w') as writing_file:
            wc = csv.writer(writing_file, lineterminator='\n')
            labels = ['Index', 'vdw group', 'sphere group', 'element', 'name', 'residue', 'residue sequence', 'vdw vol',
                      'rad', 'aw vol', 'pow vol', 'pow vol diff', 'aw vol diff', 'aw sa', 'pow sa', 'association',
                      'aw sphericity', 'pow sphericity', 'pow sphereicity diff', 'aw sphereicity diff', 'aw sol facing',
                      'pow sol facing', 'aw aa facing', 'pow aa facing', 'aw na facing', 'pow na facing',
                      'aw sep chain iface', 'pow sep chain iface', 'aw sep res iface', 'pow sep res iface']
            wc.writerow(labels)
            for spleesh in my_dict:
                wc.writerow([my_dict[spleesh][_] for _ in labels])
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
            wc.writerow(['aw vol', 'cluster', 'pow vol', 'vdw vol', 'element', 'name', 'pow nbors', 'residue', 'rad',
                         'pow sa', 'association', 'pow sphericity', 'pow sol facing', 'pow aa facing', 'pow na facing',
                         'pow sep chain iface', 'pow sep res iface'])
            for spleesh in my_dict:
                wc.writerow([my_dict[spleesh]['aw vol'], my_dict[spleesh]['vdw group'], my_dict[spleesh]['pow vol'],
                             my_dict[spleesh]['vdw vol'], my_dict[spleesh]['element'], my_dict[spleesh]['name'],
                             len(my_dict[spleesh]['pow nbors']), my_dict[spleesh]['residue'], my_dict[spleesh]['rad'],
                             my_dict[spleesh]['pow sa'], my_dict[spleesh]['association'],
                             my_dict[spleesh]['pow sphericity'], my_dict[spleesh]['pow sol facing'],
                             my_dict[spleesh]['pow aa facing'], my_dict[spleesh]['pow na facing'],
                             my_dict[spleesh]['pow sep chain iface'], my_dict[spleesh]['pow sep res iface']])

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


def get_rects(dictin, combine_per=95, cushion_per=2, print_info=True, ):
    # Looking to put the coordinates in their bounding boxes and combine boxes if they overlap. We want to return three lists:
    # names, coordinates, and bounding boxes.

    # Set up the initial bbox
    my_bbox = [[np.inf, np.inf], [-np.inf, -np.inf]]
    # Get the full bounding box
    for entry in dictin:
        # Loop through the coordinates in the entry
        for coord in dictin[entry]:
            # Set the bounding box
            my_bbox = [[min(my_bbox[0][0], coord[0]), min(my_bbox[0][1], coord[1])],
                       [max(my_bbox[1][0], coord[0]), max(my_bbox[1][1], coord[1])]]
    # Get the total width and height of the bbox
    fw, fh = my_bbox[1][0] - my_bbox[0][0], my_bbox[1][1] - my_bbox[0][1]

    # Get the cushion
    cush = (cushion_per / 100) * fw, (cushion_per / 100) * fh
    # Create a new dictionary
    new_dictin = {}
    # Loop through the preset groups
    for entry in dictin:
        # Get the coords
        coords = dictin[entry]
        # Get the bounding box and add the cushion
        bbox = [[min([_[0] for _ in coords]) - cush[0], min([_[1] for _ in coords]) - cush[0]],
                [max([_[0] for _ in coords]) + cush[1], max([_[1] for _ in coords]) + cush[1]]]

        # Add the box to the dictin
        new_dictin[entry] = {'bbox': bbox, 'coords': coords, 'olap_keys': []}

    # Create the is_inside function
    def is_inside(point, my_bbox):
        # Check if the point is inside
        return my_bbox[0][0] <= point[0] <= my_bbox[1][0] and my_bbox[0][1] <= point[1] <= my_bbox[1][1]

    # Loop through the entries again to find overlaps
    for entry in new_dictin:
        for entry2 in new_dictin:
            if entry == entry2:
                continue
            per = 100 * len([_ for _ in new_dictin[entry2]['coords'] if is_inside(_, new_dictin[entry]['bbox'])]) / len(
                new_dictin[entry2]['coords'])
            if per > combine_per:
                new_dictin[entry]['olap_keys'].append(entry2)

    # Sort the dictin by the number of overlaps
    nw_dictin = dict(sorted(new_dictin.items(), key=lambda item: len(item[1]['olap_keys']), reverse=True))

    # Merge overlapping bounding boxes
    merged_names, merged_coords, merged_bboxes, added_vals = [], [], [], []
    # Loop through the new sorted dictionary
    for entry in nw_dictin:
        # Check to see if the group has been merged yet
        if entry in added_vals:
            continue
        # Find the overlaps that arent in another group
        new_addys = [_ for _ in nw_dictin[entry]['olap_keys'] if _ not in added_vals]
        # Store the new values so we can know which groups to skip later
        added_vals += new_addys
        # Make sure the group stores all names in the new dictionary
        merged_names.append([entry] + new_addys)
        # Merge the set of coordinates
        merged_coords.append([nw_dictin[entry]['coords']] + [nw_dictin[_]['coords'] for _ in new_addys])
        # Add the bboxes together to find the new bbox
        all_bboxes = [nw_dictin[entry]['bbox']] + [nw_dictin[_]['bbox'] for _ in new_addys]
        # Find the maximum of all the bboxes to get the new bounding box
        new_bbox = [[min([_[0][0] for _ in all_bboxes]) + cush[0], min([_[0][1] for _ in all_bboxes]) + cush[1]],
                    [max([_[1][0] for _ in all_bboxes]) - cush[0], max([_[1][1] for _ in all_bboxes]) - cush[1]]]
        # Store the bbox
        merged_bboxes.append(new_bbox)
    if print_info:
        for i in range(len(merged_names)):
            print(f"Cluster=r: {len(merged_names)}, Count: {len(merged_coords[i])}, Names: {merged_names[i]}")
    # Return the bbox
    return merged_names, merged_coords, merged_bboxes


def get_convex_hulls(dictin, combine_per=95, cushion_per=2, print_info=True, outlier_thresh=3.0, make_csv=None, folder=None):
    my_bbox = [[np.inf, np.inf], [-np.inf, -np.inf]]
    my_coords = []

    for entry in dictin:
        for coord in dictin[entry]['coords']:
            my_bbox = [[min(my_bbox[0][0], coord[0]), min(my_bbox[0][1], coord[1])],
                       [max(my_bbox[1][0], coord[0]), max(my_bbox[1][1], coord[1])]]
            my_coords.append(coord)

    fw, fh = my_bbox[1][0] - my_bbox[0][0], my_bbox[1][1] - my_bbox[0][1]
    new_dictin = {}

    for entry in dictin:
        coords = np.array(dictin[entry]['coords'])

        if len(coords) >= 3:
            center = np.mean(coords, axis=0)
            dists = np.linalg.norm(coords - center, axis=1)
            threshold = np.mean(dists) + outlier_thresh * np.std(dists)
            coords = coords[dists <= threshold]

        if len(coords) < 3:
            hull_points = coords
        else:
            hull = ConvexHull(coords)
            hull_points = coords[hull.vertices]
            center = np.mean(hull_points, axis=0)
            hull_points = center + (hull_points - center) * (1 + cushion_per / 100)
            hull_points = np.vstack([hull_points, hull_points[0]])

        new_dictin[entry] = {'hull': hull_points, 'coords': coords, 'olap_keys': [], 'indexes': dictin[entry]['indexes']}

    def is_inside(point, hull_points):
        from matplotlib.path import Path
        return Path(hull_points).contains_point(point)

    for entry in new_dictin:
        for entry2 in new_dictin:
            if entry == entry2:
                continue
            per = 100 * len([_ for _ in new_dictin[entry2]['coords'] if is_inside(_, new_dictin[entry]['hull'])]) / len(
                new_dictin[entry2]['coords'])
            if per > combine_per:
                new_dictin[entry]['olap_keys'].append(entry2)

    nw_dictin = dict(sorted(new_dictin.items(), key=lambda item: len(item[1]['olap_keys']), reverse=True))
    merged_names, merged_coords, merged_hulls, merged_ndxs, added_vals = [], [], [], [], []

    for entry in nw_dictin:
        if entry in added_vals:
            continue
        new_addys = [_ for _ in nw_dictin[entry]['olap_keys'] if _ not in added_vals]
        added_vals += new_addys

        merged_names.append([entry] + new_addys)

        merged_coords.append(np.vstack([nw_dictin[entry]['coords']] + [nw_dictin[_]['coords'] for _ in new_addys]))
        all_points = np.vstack([nw_dictin[entry]['coords']] + [nw_dictin[_]['coords'] for _ in new_addys])
        indexes = nw_dictin[entry]['indexes'].copy()
        for addy in new_addys:
            indexes += nw_dictin[addy]['indexes']
        merged_ndxs.append(indexes)

        if len(all_points) >= 3:
            center = np.mean(all_points, axis=0)
            dists = np.linalg.norm(all_points - center, axis=1)
            threshold = np.mean(dists) + outlier_thresh * np.std(dists)
            all_points = all_points[dists <= threshold]

        if len(all_points) < 3:
            merged_hull = all_points
            area = 0
        else:
            hull = ConvexHull(all_points)
            merged_hull = all_points[hull.vertices]
            area = hull.volume
            center = np.mean(merged_hull, axis=0)
            merged_hull = center + (merged_hull - center) * (1 + cushion_per / 100)
            merged_hull = np.vstack([merged_hull, merged_hull[0]])

        merged_hulls.append(merged_hull)

    if print_info:
        for i in range(len(merged_names)):
            name_counts = {name: len(dictin[name]) for name in merged_names[i]}
            name_counts_str = ", ".join(f"{name} ({count})" for name, count in name_counts.items())
            print(f"Cluster: {i + 1}, Count: {len(merged_coords[i])}, "
                  f"COM: {[round(_, 4) for _ in calc_com(merged_coords[i])]}, Names: {name_counts_str}, "
                  f"Area: {ConvexHull(merged_coords[i]).volume if len(merged_coords[i]) >= 3 else 0}")

    return merged_names, merged_coords, merged_hulls, merged_ndxs


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

    for color_by in ['element', 'residue', 'bb_sc']:
        if color_by == 'element':
            color_dict = {'C': 'grey', 'O': 'r', 'N': 'b', 'P': 'darkorange', 'H': 'pink', 'S': 'y', 'Se': 'sandybrown'}
            colors = [color_dict[my_dict[entry]['element']] for entry in my_dict]
        elif color_by == 'residue':
            colors = [residue_colors[my_dict[entry]['residue']] for entry in my_dict]
        elif color_by == 'bb_sc':
            colors = [bb_sc_colors[my_dict[entry]['name']] if my_dict[entry]['name'] in bb_sc_colors else 'black' for entry in my_dict]
        # Create the dictionary for the sphericity and vdw vol groupings
        sphere_groupings, vdw_groupings = {}, {}
        # Go through the two different plotting x axes
        for x_val in ['sphericity', 'vdw_volume']:
            # Set up the res dictionaries and the coordinates lists
            my_res_dict, my_res_dict_no_h, coords, coords_no_h = {}, {}, [], []
            # Plot the different types of possible assignments for the atoms
            for ass in ['na', 'aa', 'other']:
                dictalonius = {}
                # Plot the rectangles
                fig, ax = plt.subplots()
                # Go through the entries in the dictionary to get the information to plot
                for i, entry in enumerate(my_dict):
                    # Check that the association makes sense
                    if my_dict[entry]['association'] != ass:
                        continue

                    # Get the x value for the set plotting x axis and choose the correct dictionary
                    if x_val == 'sphericity':
                        s_diff = (my_dict[entry]['pow sphericity'] - my_dict[entry]['aw sphericity']) / my_dict[entry]['aw sphericity']
                        dictalonius = sphere_groupings
                    else:
                        s_diff = my_dict[entry]['vdw vol']
                        dictalonius = vdw_groupings
                    # Get the y value or the volume difference
                    vol_diff = (my_dict[entry]['pow vol'] - my_dict[entry]['aw vol']) / my_dict[entry]['aw vol']
                    # Add the coordinates to the list for later reference
                    coords.append([s_diff, vol_diff])

                    # Create a label for later grouping
                    label = (my_dict[entry]['residue'], my_dict[entry]['name'])
                    # Add the value to the high specificity group (res, name)
                    if label in dictalonius:
                        dictalonius[label]['coords'].append([s_diff, vol_diff])
                        dictalonius[label]['indexes'].append(my_dict[entry]['Index'])
                    else:
                        dictalonius[label] = {'coords': [[s_diff, vol_diff]], 'indexes': [my_dict[entry]['Index']]}

                    # Plot the values
                    ax.scatter([s_diff], [vol_diff], marker='x' if my_dict[entry]['aw sol facing'] else 'o',
                                c=colors[i], alpha=0.2)

                    if my_dict[entry]['residue'] in my_res_dict:
                        my_res_dict[my_dict[entry]['residue']].append((s_diff, vol_diff))
                    else:
                        my_res_dict[my_dict[entry]['residue']] = [(s_diff, vol_diff)]
                    if my_dict[entry]['element'].lower() != 'h':
                        coords_no_h.append([s_diff, vol_diff])
                        if my_dict[entry]['residue'] in my_res_dict_no_h:
                            my_res_dict_no_h[my_dict[entry]['residue']].append((s_diff, vol_diff))
                        else:
                            my_res_dict_no_h[my_dict[entry]['residue']] = [(s_diff, vol_diff)]

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
                        elif name not in (amino_bbs + amino_scs + nucleic_nbase + nucleic_pphte + nucleic_sugr):
                            print(f'Unidentified atom name {name}')
                            if 'Other' not in added_groups:
                                element_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                                    markerfacecolor='gray', label='Other', markersize=8))
                                added_groups.add('Other')
                to_merges, all_coordss, bboxs, indexes = (
                    get_convex_hulls(dictalonius, cushion_per=10, combine_per=50, folder=output_folder, outlier_thresh=1.5))

                for i in range(len(indexes)):
                    for ndx in indexes[i]:
                        if x_val == 'sphericity':
                            my_dict[ndx]['sphere group'] = i
                        else:
                            my_dict[ndx]['vdw group'] = i
                if x_val == 'vdw_volume' and ass == 'other' and color_by == 'bb_sc':
                    for entry in my_dict:
                        if 'sphere group' not in my_dict[entry]:
                            my_dict[entry]['sphere group'] = -1
                        if 'vdw group' not in my_dict[entry]:
                            my_dict[entry]['vdw group'] = -1

                    write_files(my_dict, output_folder, True)

                for i, bbox in enumerate(bboxs):
                    #     # # plot the rectangle
                    #     # verts = [bbox[0], [bbox[0][0], bbox[1][1]], bbox[1], [bbox[1][0], bbox[0][1]], bbox[0]]
                    #     # ax.plot([_[0] for _ in verts], [_[1] for _ in verts])
                    ax.plot([_[0] for _ in bbox], [_[1] for _ in bbox], linewidth=1)

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
                plt.title(f'{file_name}', fontsize=25)
                plt.tight_layout()
                plt.show()

            if color_by == 'element':
                my_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
                               'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'A', 'DA', 'T', 'DT', 'G', 'DG',
                               'C', 'DC', 'U', 'DU', 'MOL']
                print(ass)

                for res in my_residues:
                    if res in my_res_dict:
                        print(res)
                print("Average Pairwise Distance:")
                for res in my_residues:
                    if res in my_res_dict:
                        print(f"{average_pairwise_distance(my_res_dict[res])}")
                print("Average Pairwise Distance (no H):")
                for res in my_residues:
                    if res in my_res_dict:
                        print(f"{average_pairwise_distance(my_res_dict_no_h[res])}")
                den_rat = clustering_measures(coords, my_res_dict.values(), [_ for _ in my_res_dict], False)
                print("Density Ratio")
                for res in my_residues:
                    if res in den_rat:
                        print(f"{den_rat[res]}")
                den_rat_no_h = clustering_measures(coords_no_h, my_res_dict_no_h.values(), [_ for _ in my_res_dict_no_h], False)
                print("Density Ratio No H")
                for res in my_residues:
                    if res in den_rat_no_h:
                        print(f"{den_rat_no_h[res]}")
    return my_dict


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    # my_file = filedialog.askopenfilename(title="Get CSV File")
    my_file = None
    plot_my_stuff(dict_file=my_file, file_name='BSA')
