import os
import csv
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)


def combine_build_information(output_file, build_logs):
    # Open the file and write the headers
    with open(output_file, 'w') as of:
        # open the writer
        of_csv = csv.writer(of)
        # Write the first lines
        of_csv.writerow(['build information'])
        row_titles = ['Name', 'Location', 'Completion Date', 'Network Type', 'Surface Resolution', 'Box Size',
                      'Maximum Allowable Vertex', 'Total Time', 'Vertex Time', 'Connect Time',
                      'Surface Building Time', 'Analysis time', 'Maximum Found Vertex']
        of_csv.writerow(row_titles)
        # Create the dictionary
        build_dict = {}
        # Add the values for each of the build logs into the dictionary, so we can add them together
        for logaroony in build_logs:
            pass





def combine_group_information(output_file, group_logs):
    pass


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
