import os
import numpy as np
import tkinter as tk
from tkinter import filedialog


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)


def get_txt_from_pdb(file, out_file):
    # Open the files
    with open(file, 'r') as pdb, open(out_file, 'w') as txt:
        # Create the counter variable
        counter = 0
        # Loop through the pdb file
        for line in pdb.readlines():
            # We only need tha atom file
            if line[:4].lower() == 'atom':
                # get the location and radius values from the line
                x, y, z = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                rad = float(line[60:66])
                # Write the atom information
                txt.write(f"{x} {y} {z} {round(rad, 4)} # {counter} \n")
                counter += 1


def clean_folder(folder=None):
    """
            Structure :
            > aw
                > aw_verts.txt
                > aw_logs.csv
            > pow
                > pow_verts.txt
                > pow_logs.csv
            > balls.pdb
            > balls.txt
            > set_balls.pml
            > retaining_box.off
            """
    # Get the folder if none has been selected yet
    if folder is None:
        folder = filedialog.askdirectory()
    # C
    for subfolder in os.listdir(folder):
        # Create a joined directory name for referencing
        sub = os.path.join(folder, subfolder)
        # Check that the aw and pow folder exist
        if not os.path.exists(os.path.join(sub, 'aw')) or os.path.exists(os.path.join(sub, 'pow')):
            print("No aw or pow folder - ", subfolder)
            continue
        # Check to see if the balls txt file is in the main directory
        if 'balls.txt' not in os.listdir(sub):
            if 'balls.pdb' not in os.listdir(sub):
                print("No pdb file - ",  subfolder)
            # Copy the pdb into the txt
            get_txt_from_pdb(os.path.join(sub, 'balls.pdb'), os.path.join(sub, 'balls.txt'))


