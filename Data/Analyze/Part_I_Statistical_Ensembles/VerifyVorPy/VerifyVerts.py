import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd

# First open the folder with all the data
root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
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

