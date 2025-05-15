import os
import sys
import csv
# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)
# from curses.ascii import isdigit
import numpy as np
import tkinter as tk
import datetime
from tkinter import filedialog
from vorpy.src.calculations import calc_com, round_func, calc_total_inertia_tensor, combine_inertia_tensors


def parse_string_lists(string_list):
    if 'np.float64' in string_list:
        new_string_list = string_list.split('np.float64')
        new_string = ''.join(new_string_list)
    else:
        new_string = string_list
    my_string = ''
    for letter in new_string:
        if letter in {'(', ')'}:
            continue
        else:
            my_string += letter
    if '], [' in my_string:
        new_vals = my_string.split('], [')
        for i, val in enumerate(new_vals):
            new_val = ''
            for letter in val:
                if letter in {'[', ']'}:
                    continue
                else:
                    new_val += letter
            new_vals[i] = [float(_) for _ in new_val.split(',')]
    elif '] [' in my_string:
        new_vals = my_string.split('] [')
        for i, val in enumerate(new_vals):
            new_val = ''
            for letter in val:
                if letter in {'[', ']'}:
                    continue
                else:
                    new_val += letter
            new_vals[i] = [float(_) for _ in new_val.split(',')]
    else:
        my_new_string = ''
        for letter in my_string:
            if letter in {'[', ']'}:
                continue
            else:
                my_new_string += letter
        new_vals = [float(_) for _ in my_new_string.split(',')]
    return new_vals


def combine_build_information(output_file, build_logs):
    # Create the lines list
    lines = [['build information']]
    # Write the first lines
    row_titles = ['Name', 'Location', 'Completion Date', 'Network Type', 'Surface Resolution', 'Box Size',
                  'Maximum Allowable Vertex', 'Total Time', 'Vertex Time', 'Connect Time',
                  'Surface Building Time', 'Analysis time', 'Maximum Found Vertex']

    lines.append(row_titles)
    # Create the dictionary
    build_dict = {row_titles[i]: [] for i in range(len(row_titles))}
    # Add the values for each of the build logs into the dictionary, so we can add them together
    for logaroony in build_logs:
        # Loop through the row titles adding stuff from each of the build logs
        for i in range(len(row_titles)):
            build_dict[row_titles[i]].append(build_logs[logaroony][0][i])

    # Write the info
    lines.append([build_dict[row_titles[0]][0], output_file, datetime.datetime.now(), build_dict[row_titles[3]][0],
                  build_dict[row_titles[4]][0], build_dict[row_titles[5]][0], build_dict[row_titles[6]][0],
                  sum([float(_) for _ in build_dict[row_titles[7]]]), sum([float(_) for _ in build_dict[row_titles[8]]]),
                  sum([float(_) for _ in build_dict[row_titles[9]]]), sum([float(_) for _ in build_dict[row_titles[10]]]),
                  sum([float(_) for _ in build_dict[row_titles[11]]]), max([float(_) for _ in build_dict[row_titles[12]]])])
    # Return the lines
    return lines


def combine_group_information(group_logs, sa, moi, spatial_moment, round_to=3):
    # Get the round function
    r = round_func(round_to)
    # Write the first line
    lines = [['group information']]
    row_titles = ['Name', 'Volume', 'Surface Area', 'Mass', 'Density', 'Center of Mass', 'VDW Volume',
                  'VDW Center of Mass', 'Moment of Inertia', 'Spatial Moment of Inertia']
    lines.append(row_titles)
    # Create the dictionary
    build_dict = {row_titles[i]: [] for i in range(len(row_titles))}
    # Add the values for each of the build logs into the dictionary, so we can add them together
    for logaroony in group_logs:
        # Loop through the row titles adding stuff from each of the build logs
        for i in range(len(row_titles)):
            build_dict[row_titles[i]].append(group_logs[logaroony][0][i])
    # Get the total volume
    vols = [float(_) for _ in build_dict[row_titles[1]]]
    # Get the masses
    masses = [float(_) for _ in build_dict[row_titles[3]]]
    # Get the van der waals volumes
    vdw_vols = [float(_) for _ in build_dict[row_titles[6]]]
    # Get the center of mass
    com = calc_com([parse_string_lists(_) for _ in build_dict[row_titles[5]]],
                   vols)
    # Get the vander waals com
    vdw_com = calc_com([parse_string_lists(_) for _ in build_dict[row_titles[7]]],
                       [float(_) for _ in build_dict[row_titles[3]]])
    # Write the info
    lines.append([build_dict[row_titles[0]][0], r(sum(vols)), r(sa), r(sum(masses)), r(sum(vdw_vols) / sum(vols)),
                  [r(_) for _ in com], r(sum(vdw_vols)), [r(_) for _ in vdw_com],
                  [[float(r(__)) for __ in _] for _ in moi], [[float(r(__)) for __ in _] for _ in spatial_moment]])
    # Return the lines
    return lines


def combine_atoms_lines(atom_logs):
    # Get the initial two lines
    lines = [['Atoms'],
             ['Index', 'Name', 'Residue', 'Residue Sequence', 'Chain', 'Mass', 'X', 'Y', 'Z', 'Radius', 'Volume',
              'Van Der Waals Volume', 'Surface Area', 'Complete Cell?', 'Maximum Mean Curvature',
              'Average Mean Surface Curvature', 'Maximum Gaussian Curvature', 'Average Gaussian Surface Curvature',
              'Sphericity', 'Isometric Quotient', 'Inner Ball?', 'Number of Neighbors', 'Closest Neighbor',
              'Closest Neighbor Distance', 'Layer Distance Average', 'Layer Distance RMSD', 'Minimum Point Distance',
              'Maximum Point Distance', 'Number of Overlaps', 'Contact Area', 'Non-Overlap Volume', 'Overlap Volume',
              'Center of Mass', 'Moment of Inertia Tensor', 'Bounding Box', 'neighbors']]
    atoms = {}
    atoms_data = []
    # Loop through the atom lines
    for atom_log_set in atom_logs:
        # Get the group's atoms
        for atom in atom_logs[atom_log_set]:
            # Get the values we want
            index = int(atom[0])
            location = np.array([float(_) for _ in atom[6:9]])
            rad = float(atom[9])
            mass = float(atom[5])
            vdw_vol = float(atom[11])
            # Check if the atom has been found before
            if index in atoms:
                print("Repeat Atom!!!!", atom_log_set, index)
                continue
            # Add the data
            atoms[index] = atom
            atoms_data.append({'mass': mass, 'loc': location, 'rad': rad, 'vdw_vol': vdw_vol, 'vol': float(10),
                               'com': np.array(parse_string_lists(atom[32])), 'moi': np.array(parse_string_lists(atom[33]))})
    group = set(atoms.keys())
    sorted_atoms = [atoms[key] for key in sorted(atoms)]
    lines = lines + sorted_atoms
    # Get the vander waals com
    vdw_com = calc_com([_['loc'] for _ in atoms_data], [_['vdw_vol'] for _ in atoms_data])
    com = [sum([_['com'][i] * _['vol'] for _ in atoms_data]) / sum([_['vol'] for _ in atoms_data]) for i in range(len(atoms_data[0]['com']))]
    moi = calc_total_inertia_tensor(atoms_data, vdw_com)
    spatial_moment = combine_inertia_tensors([_['moi'] for _ in atoms_data], [_['com'] for _ in atoms_data], com,
                                             [_['vol'] for _ in atoms_data])
    return lines, group, com, vdw_com, moi, spatial_moment


def combine_surface_lines(surface_logs, group):
    # Create the surfaces dictionary for later sorting
    surfaces = {}
    sa = 0
    ndx = -2
    # Loop through the surface logs adding the surfaces that aren't repeats
    for file in surface_logs:
        # Get the dictionary from the surfaces dictionaries
        my_dict = surface_logs[file]
        # Loop through each of the surfaces in the file's dictionary
        for surf in my_dict:
            indices = tuple([int(_) for _ in surf[1:3]])
            if indices in surfaces:
                # Check the values
                # if not all([surfaces[surf][j] == my_dict[surf][j] for j in range(len(my_dict[surf]))]):
                #     print(f"Bad surf match {surf}, {surfaces[surf]} != {my_dict[surf]}")
                continue
            else:
                # Check if both indices are in the group
                if indices[0] not in group or indices[1] not in group:
                    sa += float(3)
                surf = [ndx] + surf[1:]
                surfaces[indices] = surf
                ndx += 1
    # Sort the surfaces
    sorted_surfaces = [surfaces[key] for key in sorted(surfaces)]
    # Create the list of lines
    lines = [['Surfaces'], ['Index', 'Ball 1', 'Ball 2', 'Surface Area', 'Mean Curvature', 'Gaussian Curvature',
                            'Ball 1 Volume Contribution', 'Ball 2 Volume Contribution', 'Contact Area', 'Overlap']]

    # Write the rows
    lines += sorted_surfaces
    # Return the lines
    return lines, sa


def combine_edges_lines(output_file, edge_logs):
    # Create the edges dictionary for later sorting
    edges = {}
    ndx = -2
    # Loop through the edge logs adding the edges that aren't repeats
    for file in edge_logs:
        # Get the dictionary from the edges dictionaries
        my_dict = edge_logs[file]
        # Loop through each of the edges in the file's dictionary
        for edge in my_dict:
            indices = tuple([int(_) for _ in edge[1:4]])
            if indices in edges:
                # Check the values
                # if not all([edges[edge][j] == my_dict[edge][j] for j in range(len(my_dict[edge]))]):
                #     print(f"Bad edge match {edge}, {edges[edge]} != {my_dict[edge]}")
                continue
            else:
                if ndx >= 0:
                    edge = [ndx] + edge[1:]
                edges[indices] = edge
                ndx += 1
    # Sort the edges
    sorted_edges = [edges[key] for key in sorted(edges)]
    # Add to the output file
    with open(output_file, 'a') as of:
        # Create the csv writer
        of_csv = csv.writer(of, lineterminator='\n')
        # Write the header
        of_csv.writerow(['Edges'])
        of_csv.writerow(['Index', 'Ball 1', 'Ball 2', 'Ball 3', 'Length'])
        # Write the rows
        count = 0
        for row in sorted_edges:
            of_csv.writerow([count] + row[1:])
            count += 1
    # Close the file
    of.close()


def combine_vertex_lines(output_file, vertex_logs):
    # Create the vertices dictionary for later sorting
    vertices = {}
    ndx = -2
    # Loop through the vertex logs adding the vertices that aren't repeats
    for file in vertex_logs:
        # Get the dictionary from the vertices dictionaries
        my_dict = vertex_logs[file]
        # Loop through each of the vertices in the file's dictionary
        for vert in my_dict:
            indices = tuple([int(_) for _ in vert[1:5]])
            if indices in vertices:
                # Check the values
                # if not all([vertices[vert][j] == my_dict[vert][j] for j in range(len(my_dict[vert]))]):
                #     print(f"Bad vertex match {vert}, {vertices[vert]} != {my_dict[vert]}")
                continue
            else:
                if ndx >= 0:
                    vert = [ndx] + vert[1:]
                vertices[indices] = vert
                ndx += 1
    # Sort the vertices
    sorted_vertices = [vertices[key] for key in sorted(vertices)]
    # Add to the output file
    with open(output_file, 'a') as of:
        # Create the csv writer
        of_csv = csv.writer(of, lineterminator='\n')
        # Write the header
        of_csv.writerow(['Vertices'])
        of_csv.writerow(['Index', 'Ball 1', 'Ball 2', 'Ball 3', 'Ball 4', 'x', 'y', 'z', 'r'])
        count = 0
        # Write the rows
        for row in sorted_vertices:
            of_csv.writerow([count] + row[1:])
            count += 1
    # Close the file
    of.close()


def combine_logs(list_of_logs=None, output_dir=None):
    """
    Combines the logs of separate split files.
    1. Get the files
    2. Read the files and make a dictionary of each set -
        (build information, group information, Atoms, Edges, Surfaces, Vertices)
    3. Write the file
    """
    # If files aren't loaded up
    if list_of_logs is None:
        # Create the list of log file addresses to be combined
        list_of_logs = []
        # Keep looping through until no file is selected
        while True:
            # Get the logs file
            logs = filedialog.askopenfilename(title='Get new file')
            print(logs[len(os.path.dirname(os.path.dirname(logs))) + 1:])
            # Check if it exists
            if os.path.exists(logs) and logs != '':
                list_of_logs.append(logs)
            else:
                break

    # Go through the logs and
    # Get the lines for each group: Build, Group, Atoms, Surfaces, Edges, Vertices
    build, group, atoms, surfs, edges, verts = {}, {}, {}, {}, {}, {}
    for file in list_of_logs:
        build[file], group[file], atoms[file], surfs[file], edges[file], verts[file] = [], [], [], [], [], []

        with open(file, 'r') as read_file:
            rf = csv.reader(read_file)
            skip = False
            add_dict = None
            for line in rf:
                if skip:
                    skip = False
                    continue
                if line[0] == 'build informaiton' or line[0] == 'build information':
                    add_dict = build
                    skip = True
                    continue
                if line[0] == 'group information':
                    add_dict = group
                    skip = True
                    continue
                if line[0] == 'Atoms':
                    add_dict = atoms
                    skip = True
                    continue
                if line[0] == 'Surfaces':
                    add_dict = surfs
                    skip = True
                    continue
                if line[0] == 'Edges':
                    add_dict = edges
                    skip = True
                    continue
                if line[0] == 'Vertices':
                    add_dict = verts
                    skip = True
                    continue
                add_dict[file].append(line)
    # Check that the output dir is None
    if output_dir is None:
        # Choose the output directory for the final file
        output_dir = filedialog.askdirectory(title="Choose the output folder for \'Total_logs.csv\'")
    # Create the output file name
    bld_lines = combine_build_information(output_dir + '/Total_logs.csv', build)
    # Get the atoms
    atm_lines, group_nums, com, vdw_com, moi, spatial_moment = combine_atoms_lines(atoms)
    # Get the surfaces
    srf_lines, sa = combine_surface_lines(surfs, group_nums)
    # Group logs
    grp_lines = combine_group_information(group, sa, moi, spatial_moment)
    # Write
    with open(output_dir + '/Total_logs.csv', 'w') as outpt_file:
        of = csv.writer(outpt_file, lineterminator='\n')
        for line in bld_lines:
            of.writerow(line)
        for line in grp_lines:
            of.writerow(line)
        for line in atm_lines:
            of.writerow(line)
        count = -2
        for line in srf_lines:
            if count >= 0:
                of.writerow([count] + line[1:])
            else:
                of.writerow(line)
            count += 1
    combine_edges_lines(output_dir + '/Total_logs.csv', edges)
    combine_vertex_lines(output_dir + '/Total_logs.csv', verts)


if __name__ == '__main__':
    # combine_logs(list_of_logs=['/Users/jackericson/PycharmProjects/vorpy/Data/user_data/EDTA_Mg_2/a_1_aw/aw_logs.csv',
    #                            '/Users/jackericson/PycharmProjects/vorpy/Data/user_data/EDTA_Mg_3/a_2_aw/aw_logs.csv'],
    #              output_dir='/Users/jackericson/PycharmProjects/vorpy/Data/user_data')
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    combine_logs()
