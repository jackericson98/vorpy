import os
import sys
import pandas
import csv
import numpy
import sklearn
import tkinter as tk
from tkinter import filedialog

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools import read_logs2
from vorpy.src.inputs import read_pdb_simple




def gather_file_names():
    # Get the folder path
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askdirectory(title="Select a folder")
    # check for pbc
    pbc = False
    if 'non' in file_path:
        pbc = True

    file_dict = {}
    # Gether the files and make the dictionary with the file_names as the keys and the file_paths as the values
    file_names = os.listdir(file_path)
    # Loop through the files
    for sub_folder in file_names:
        # Get the cv, number of balls, density, overlap, pbc/non-pbc, and the folder number
        values = sub_folder.split("_")
        mean, cv, num_balls, density, overlap, folder_num = values
        # Make the key for the dictionary
        key = (mean, cv, num_balls, density, overlap, pbc, folder_num)
        # Add the file_path to the dictionary
        file_dict[key] = {'aw': os.path.join(file_path, sub_folder, 'aw', 'aw_logs.csv'), 'pow': os.path.join(file_path, sub_folder, 'pow', 'pow_logs.csv'), 'pbc': os.path.join(file_path, sub_folder, 'balls.txt')}
    # return the dictionary
    return file_dict


def get_atoms_info(file_dict, csv_file, start=False):
    # Check to see if we need a header with the names of the columns
    if start:
        # Open the csv file and write the header
        with open(csv_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Index', 'Name', 'Residue', 'Residue Sequence', 'Chain', 'Mass', 'X', 'Y', 'Z', 'Radius', 'AW Volume', 'POW Volume', 'AW Van Der Waals Volume', 
                             'POW Van Der Waals Volume', 'AW Surface Area', 'POW Surface Area', 'AW Complete Cell?', 'POW Complete Cell?', 'Maximum Mean Curvature', 
                             'Average Mean Surface Curvature', 'Maximum Gaussian Curvature', 'Average Gaussian Surface Curvature', 'AW Sphericity', 'POW Sphericity', 
                             'AW Isometric Quotient', 'POW Isometric Quotient', 'AW Inner Ball?', 'POW Inner Ball?', 'AW Number of Neighbors', 'AW Closest Neighbor', 
                             'POW Closest Neighbor', 'AW Closest Neighbor Distance', 'POW Closest Neighbor Distance', 'AW Layer Distance Average', 'POW Layer Distance Average', 
                             'AW Layer Distance RMSD', 'POW Layer Distance RMSD', 'AW Minimum Point Distance', 'POW Minimum Point Distance', 'AW Maximum Point Distance', 
                             'POW Maximum Point Distance', 'AW Number of Overlaps', 'POW Number of Overlaps', 'AW Contact Area', 'POW Contact Area', 'AW Non-Overlap Volume', 
                             'POW Non-Overlap Volume', 'AW Overlap Volume', 'POW Overlap Volume', 'AW Center of Mass', 'POW Center of Mass', 'AW Moment of Inertia Tensor', 
                             'POW Moment of Inertia Tensor', 'AW Bounding Box', 'POW Bounding Box', 'AW neighbors', 'POW neighbors'])
    # Loop through the dictionary
    for key, value in file_dict.items():
        # Open the aw_logs file
        aw_logs = read_logs2(value['aw'], all_=False, balls=True)
        # Open the pow_logs file
        pow_logs = read_logs2(value['pow'], all_=False, balls=True)
        # Create the list for the row to be written to the csv file
        row = [aw_logs['Index'], aw_logs['Name'], aw_logs['Residue'], aw_logs['Residue Sequence'], aw_logs['Chain'], aw_logs['Mass'], aw_logs['X'], aw_logs['Y'], aw_logs['Z'], 
               aw_logs['Radius'], aw_logs['Volume'], pow_logs['Volume'], aw_logs['Van Der Waals Volume'], pow_logs['Van Der Waals Volume'], aw_logs['Surface Area'], pow_logs['Surface Area'], 
               aw_logs['Complete Cell?'], pow_logs['Complete Cell?'], aw_logs['Maximum Mean Curvature'], aw_logs['Average Mean Surface Curvature'], aw_logs['Maximum Gaussian Curvature'], 
               aw_logs['Average Gaussian Surface Curvature'], aw_logs['Sphericity'], pow_logs['Sphericity'], aw_logs['Isometric Quotient'], pow_logs['Isometric Quotient'], 
               aw_logs['Inner Ball?'], pow_logs['Inner Ball?'], aw_logs['Number of Neighbors'], aw_logs['Closest Neighbor'], pow_logs['Closest Neighbor'], 
               aw_logs['Closest Neighbor Distance'], pow_logs['Closest Neighbor Distance'], aw_logs['Layer Distance Average'], pow_logs['Layer Distance Average'], 
               aw_logs['Layer Distance RMSD'], pow_logs['Layer Distance RMSD'], aw_logs['Minimum Point Distance'], pow_logs['Minimum Point Distance'], aw_logs['Maximum Point Distance'], 
               pow_logs['Maximum Point Distance'], aw_logs['Number of Overlaps'], pow_logs['Number of Overlaps'], aw_logs['Contact Area'], pow_logs['Contact Area'], 
               aw_logs['Non-Overlap Volume'], pow_logs['Non-Overlap Volume'], aw_logs['Overlap Volume'], pow_logs['Overlap Volume'], aw_logs['Center of Mass'], pow_logs['Center of Mass'], 
               aw_logs['Moment of Inertia Tensor'], pow_logs['Moment of Inertia Tensor'], aw_logs['Bounding Box'], pow_logs['Bounding Box'], aw_logs['neighbors'], pow_logs['neighbors']]

        # Write the row to the csv file
        writer.writerow(row)

def get_surfs_info(file_dict, csv_file, start=False):
    # Check to see if we need a header with the names of the columns
    if start:
        # Open the csv file and write the header
        with open(csv_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow()




