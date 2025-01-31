import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd


def get_information():
    # First open the folder with all the data
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    vpy_fl = filedialog.askopenfilename(title="Get Vorpy Vertices")
    vta_fl = filedialog.askopenfilename(title="Get Voronota Vertices")
    vvv_fl = filedialog.askopenfilename(title="Get V Vertices")

    # # Go through the vorpy vertices
    # with open(vpy_fl, 'r') as vpy_file:
    #     # Create the dictionary
    #     vpy_vrts = {}
    #     # Create a line counter
    #     line_counter = -1
    #     # Loop through the lines
    #     for vpy_line in vpy_file.readlines():
    #         # Increment the line counter
    #         line_counter += 1
    #         # Split the line
    #         line_list = [_ for _ in vpy_line.split(" ") if _ != ""]
    #         # If the line list is not the correct length we need to skip it so it doesn't crash the program
    #         if len(line_list) != 9:
    #             print(vpy_line)
    #             continue
    #         # Read the first line
    #         if line_counter == 0:
    #             # Get some vorpy information
    #             vpy_vrts['info'] = {
    #                 'Number': int(line_list[2]),
    #                 'Max Vert': float(line_list[9][:-1])
    #             }
    #             # Continue onto the next line
    #             continue
    #         # Get the vertex indices
    #         ndxs = tuple([int(_) for _ in line_list[:4]])
    #         # Get the vertex location
    #         loc, rad, dub = tuple([float(_) for _ in line_list[4:7]]), float(line_list[7]), int(line_list[8][0])
    #         # Check if the vertex is a doublet
    #         if dub == 1:
    #             vpy_vrts[ndxs]['loc2'], vpy_vrts['rad2'], vpy_vrts['dub'] = loc, rad, True
    #             continue
    #         # Add the vertex to the dictionary
    #         vpy_vrts[ndxs] = {'loc': loc, 'rad': rad, 'loc2': None, 'rad2': None, 'dub': False}

    # Go through the vorpy vertices
    with open(vta_fl, 'r') as vta_file:
        # Create the dictionary
        vta_vrts = {}
        # Create a line counter
        line_counter = -1
        # Loop through the lines
        for vta_line in vta_file.readlines():
            # Increment the line counter
            line_counter += 1
            # Split the line
            line_list = [_ for _ in vta_line.split(" ") if _ != ""]
            # If the line list is not the correct length we need to skip it so it doesn't crash the program
            if len(line_list) != 9:
                print(vta_line)
                continue
            # Get the vertex indices
            ndxs = tuple([int(_) for _ in line_list[:4]])
            # Get the vertex location
            loc, rad = tuple([float(_) for _ in line_list[4:7]]), float(line_list[7])
            # Check if the vertex is a doublet
            if ndxs in vta_vrts:
                if rad < vta_vrts[ndxs]:
                    vta_vrts =
                vta_vrts[ndxs]['loc2'], vta_vrts[ndxs]['rad2'], vta_vrts[ndxs]['dub'] = loc, rad, True
                continue
            # Add the vertex to the dictionary
            vta_vrts[ndxs] = {'loc': loc, 'rad': rad, 'loc2': None, 'rad2': None, 'dub': False}


folder = filedialog.askdirectory(title='Get the overall directory')
i = 0

results = {}

for roott, dirs, files in os.walk(folder):
    if i == 0:
        for dire in dirs:
            results[dire] = {}


    # Check for 'Voronota' and 'Vorpy' directories and specific files
    if 'Voronota' in dirs:
        key = None
        for dire in results:
            if dire in roott and dire != 'V':
                key = dire
        print(roott, key)
        voronota_path = os.path.join(roott, 'Voronota')
        # print
        if os.path.exists(voronota_path) and 'vertices.txt' in os.listdir(voronota_path):
            results[key]['Voronota'] = {'file': os.path.join(voronota_path, 'vertices.txt')}
        vorpy_path = os.path.join(roott, 'Vorpy')
        if os.path.exists(vorpy_path) and 'aw_verts.txt' in os.listdir(vorpy_path):
            results[key]['Vorpy'] = {'file': os.path.join(vorpy_path, 'aw_verts.txt')}


print(results)
# Output the dictionary
for key, value in results.items():

    # Load the data from both files
    try:
        file1 = pd.read_csv(value['Voronota']['file'], delimiter=' ', header=None, dtype=str)
        file2 = pd.read_csv(value['Vorpy']['file'], delimiter=' ', header=None, dtype=str)
    except KeyError:
        continue
    # Assume the first four columns are vertex identifiers for both files
    columns_of_interest = [0, 1, 2, 3]

    # Create sets of tuples from both files for comparison
    set1 = {tuple(row) for row in file1[columns_of_interest].values}
    set2 = {tuple(row) for row in file2[columns_of_interest].values}
    # print(set1, set2)
    # Find intersections (common elements)
    intersection = set1.intersection(set2)
    vpy_miss_verts = len(set2) - len(intersection)
    vta_miss_verts = len(set1) - len(intersection)

    # Calculate similarity metrics
    similarity_percentage = len(intersection) / min(len(set1), len(set2)) * 100

    print(f"Number of overlapping vertex configurations: {len(intersection)}")
    print(f"Percentage of smaller dataset represented in the overlap: {similarity_percentage:.2f}%")
    print(f"Number of VorPy Vertices = {set2} Number of Voronota Vertices = {set1} \n"
          f"Number of Voronota Vertices Missing from VorPy = {vta_miss_verts}, "
          f"Number of VorPy Vertices Missing from Voronota = {vpy_miss_verts}")

# Next get the voronota vertices

# Next filter the vertices that dont pertain to the main groups

# Get the Vorpy Vertices

