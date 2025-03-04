import os
import csv
import tkinter as tk
import datetime
from tkinter import filedialog
from System.sys_funcs.calcs.calcs import calc_com, round_func

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)


def parse_string_lists(string_list, apply_type=float):
    # Test if it is just one single list
    if string_list[1] != '[':
        listy, current_number = [], ''
        for letter in string_list:
            if letter.isdigit() or letter == '.':
                current_number += letter
            elif letter == ',':
                listy.append(apply_type(current_number))
                current_number = ''
    else:
        listy, current_number = [[]], ''
        for letter in string_list:
            if letter.isdigit() or letter == '.':
                current_number += letter
            elif letter == ',' and len(current_number) > 0:
                listy[-1].append(apply_type(current_number))
                current_number = ''
            elif letter == ']':
                listy.append([])
                current_number = ''
    return listy


def get_sa():
    return None


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
            build_dict[row_titles[i]].append(build_logs[logaroony][i])
    # Write the info
    lines.append([build_dict[row_titles[0]][0], output_file, datetime.datetime.now(), build_dict[row_titles[3]][0],
                  build_dict[row_titles[4]][0], build_dict[row_titles[5]][0], build_dict[row_titles[6]][0],
                  sum(build_dict[row_titles[7]]), sum(build_dict[row_titles[8]]), sum(build_dict[row_titles[9]]),
                  sum(build_dict[row_titles[10]]), sum(build_dict[row_titles[11]]), max(build_dict[row_titles[12]])])
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
            build_dict[row_titles[i]].append(group_logs[logaroony][i])
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


def combine_atoms_lines(output_file, atom_logs):
    pass


def combine_surface_lines(output_file, surface_logs):
    # Create the surfaces dictionary for later sorting
    surfaces = {}
    # Loop through the surface logs adding the surfaces that aren't repeats
    for file in surface_logs:
        # Get the dictionary from the surfaces dictionaries
        my_dict = surface_logs[file]
        # Loop through each of the surfaces in the file's dictionary
        for surf in my_dict:
            if surf in surfaces:
                # Check the values
                if not all([surfaces[surf][j] == my_dict[surf][j] for j in range(len(my_dict[surf]))]):
                    print(f"Bad edge match {surf}, {surfaces[surf]} != {my_dict[surf]}")
                continue
            else:
                surfaces[surf] = my_dict[surf]
    # Sort the surfaces
    sorted_surfaces = [surfaces[key] for key in sorted(surfaces)]
    # Add to the output file
    with open(output_file, 'a') as of:
        # Create the csv writer
        of_csv = csv.writer(of)
        # Write the header
        of_csv.writerow(['Surfaces'])
        of_csv.writerow(['Index', 'Ball 1', 'Ball 2', 'Surface Area', 'Mean Curvature', 'Gaussian Curvature',
                         'Ball 1 Volume Contribution', 'Ball 2 Volume Contribution', 'Contact Area', 'Overlap'])
        # Write the rows
        for row in sorted_surfaces:
            of_csv.writerow(row)
    # Close the file
    of.close()


def combine_edges_lines(output_file, edge_logs):
    # Create the edges dictionary for later sorting
    edges = {}
    # Loop through the edge logs adding the edges that aren't repeats
    for file in edge_logs:
        # Get the dictionary from the edges dictionaries
        my_dict = edge_logs[file]
        # Loop through each of the edges in the file's dictionary
        for edge in my_dict:
            if edge in edges:
                # Check the values
                if not all([edges[edge][j] == my_dict[edge][j] for j in range(len(my_dict[edge]))]):
                    print(f"Bad edge match {edge}, {edges[edge]} != {my_dict[edge]}")
                continue
            else:
                edges[edge] = my_dict[edge]
    # Sort the edges
    sorted_edges = [edges[key] for key in sorted(edges)]
    # Add to the output file
    with open(output_file, 'a') as of:
        # Create the csv writer
        of_csv = csv.writer(of)
        # Write the header
        of_csv.writerow(['Edges'])
        of_csv.writerow(['Index', 'Ball 1', 'Ball 2', 'Ball 3', 'Length'])
        # Write the rows
        for row in sorted_edges:
            of_csv.writerow(row)
    # Close the file
    of.close()


def combine_vertex_lines(output_file, vertex_logs):
    # Create the vertices dictionary for later sorting
    vertices = {}
    # Loop through the vertex logs adding the vertices that aren't repeats
    for file in vertex_logs:
        # Get the dictionary from the vertices dictionaries
        my_dict = vertex_logs[file]
        # Loop through each of the vertices in the file's dictionary
        for vert in my_dict:
            if vert in vertices:
                # Check the values
                if not all([vertices[vert][j] == my_dict[vert][j] for j in range(len(my_dict[vert]))]):
                    print(f"Bad vertex match {vert}, {vertices[vert]} != {my_dict[vert]}")
                continue
            else:
                vertices[vert] = my_dict[vert]
    # Sort the vertices
    sorted_vertices = [vertices[key] for key in sorted(vertices)]
    # Add to the output file
    with open(output_file, 'a') as of:
        # Create the csv writer
        of_csv = csv.writer(of)
        # Write the header
        of_csv.writerow(['Vertices'])
        of_csv.writerow(['Index', 'Ball 1', 'Ball 2', 'Ball 3', 'Ball 4', 'x', 'y', 'z', 'r'])
        # Write the rows
        for row in sorted_vertices:
            of_csv.writerow(row)
    # Close the file
    of.close()




def combine_logs():
    """
    Combines the logs of separate split files.
    1. Get the files
    2. Read the files and make a dictionary of each set -
        (build information, group information, Atoms, Edges, Surfaces, Vertices)
    3. Write the file
    """
    # Create the list of log file addresses to be combined
    list_of_logs = []
    # Keep looping through until no file is selected
    while True:
        # Get the logs file
        logs = filedialog.askopenfilename(title='Get new file')
        # Check if it exists
        if os.path.exists(logs) and logs != '':
            list_of_logs.append(logs)
        else:
            break

    # Go through the logs and
