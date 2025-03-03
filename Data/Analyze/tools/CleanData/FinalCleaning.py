import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.CleanData.CreateVertsFromLogs import write_log_to_vert


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
    num_subs = len([_ for _ in os.listdir(folder)])
    for i, subfolder in enumerate(os.listdir(folder)):
        # Make a print
        print(f"\r Folder {i + 1}/{num_subs}", end="")
        # Create a joined directory name for referencing
        sub = os.path.join(folder, subfolder)
        # Check that the aw and pow folder exist
        if not os.path.exists(os.path.join(sub, 'aw')) or not os.path.exists(os.path.join(sub, 'pow')):
            print("No aw or pow folder - ", subfolder)
            continue
        # Check to see if the balls txt file is in the main directory
        if 'balls.txt' not in os.listdir(sub):
            if 'balls.pdb' not in os.listdir(sub):
                print("No pdb file - ",  subfolder)
            # Copy the pdb into the txt
            get_txt_from_pdb(os.path.join(sub, 'balls.pdb'), os.path.join(sub, 'balls.txt'))

        # Check if the verts don't exist and if not make them exist
        if not os.path.exists(sub + '/aw/aw_verts.txt'):
            # Check if there are no logs
            if not os.path.exists(sub + '/aw/aw_logs.csv'):
                print("No aw logs - ", subfolder)
                continue
            # Get the verts
            write_log_to_vert(sub + '/aw/aw_logs.csv', sub + '/aw/aw_verts.txt')

        # Check if the pow verts don't exist
        if not os.path.exists(sub + '/aw/aw_verts.txt'):
            # Check if there are no logs
            if not os.path.exists(sub + '/aw/aw_logs.csv'):
                print("No pow logs - ", subfolder)
                continue
            # Get the verts
            write_log_to_vert(sub + '/aw/aw_logs.csv', sub + '/aw/aw_verts.txt')

        # Check for the set atoms file
        if not os.path.exists(sub + '/set_balls.pml'):
            # See if we can just rename the set atoms file
            if os.path.exists(sub + '/set_atoms.pml'):
                os.rename(sub + '/set_atoms.pml', sub + '/set_balls.pml')
            # Otherwise let us know
            else:
                print("No set balls or set atoms - ", subfolder)
                continue

        # Make a list of the files we want
        files_wanted = {'aw_verts.txt', 'aw_logs.csv', 'pow_verts.txt', 'pow_logs.csv', 'set_balls.pml',
                        'retaining_box.off', 'set_balls.pml', 'balls.txt', 'balls.pdb'}
        # Now delete all files that aren't correct
        for rooot, direc, filess in os.walk(sub):
            for my_file in filess:
                if my_file not in files_wanted:
                    os.remove(rooot + '/' + my_file)


if __name__ == '__main__':
    clean_folder()
