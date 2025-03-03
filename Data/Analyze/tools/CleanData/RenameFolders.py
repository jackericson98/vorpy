import datetime
import os
import shutil
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.CleanData.VerifyFiles import verify_logs, verify_pdb


def get_the_best_folders(dict_list, num_needed=20):
    good, medium, trash = [], [], []
    # Loop through the folder and place them into three groups
    for my_dict in dict_list:
        if my_dict['pdb'] is not None:
            if my_dict['aw'] is not None and my_dict['pow'] is not None:
                good.append(my_dict)
            else:
                medium.append(my_dict)
        else:
            trash.append(my_dict)
    # Check if we have enough of the others to go ahead
    if len(good) >= num_needed:
        medium += good[num_needed:]
        new_good = good[:num_needed]
    else:
        new_good = good + medium
        if len(new_good) > num_needed:
            medium = new_good[num_needed:]
            new_good = new_good[:num_needed]
    # Return the three lists
    return new_good, medium, trash


def move_contents(source_folder, target_folder):
    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # List all the entries in the source folder
    entries = os.listdir(source_folder)

    for entry in entries:
        # Construct full file path
        source_path = os.path.join(source_folder, entry)
        target_path = os.path.join(target_folder, entry)

        # Move each entry to the target folder
        shutil.move(source_path, target_path)


def copy_contents(source_folder, target_folder):
    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # List all the entries in the source folder
    entries = os.listdir(source_folder)

    for entry in entries:
        source_path = os.path.join(source_folder, entry)
        target_path = os.path.join(target_folder, entry)

        # Check if the entry is a file or a folder
        if os.path.isdir(source_path):
            # Recursively copy an entire directory tree rooted at source_path
            if os.path.exists(target_path):
                # If the target directory already exists, shutil.copytree would raise an error
                # so we delete the existing target directory
                shutil.rmtree(target_path)
            shutil.copytree(source_path, target_path)
        else:
            # Copy each file to the target folder
            shutil.copy2(source_path, target_path)


def rename_folders(folder=None, include_mean=True, include_cv=True, include_ball_num=True, include_den=True,
                   include_olap=False, include_pbc=False, include_sar=False, verify_files=False):
    """
    Renames all folders in the main folder to a similar naming convention
    """
    # In no folder is provided, get one
    if folder is None:
        folder = filedialog.askdirectory(title="Select Data Directory")
    # Get the PBC and Olap information
    if 'nonPBC' in folder:
        pbc = False
    else:
        pbc = True
    if 'Olap_0.0' in folder:
        olap = 0.0
    elif 'Olap_0.5' in folder:
        olap = 0.5
    else:
        olap = 1.0
    print("Sorting the folders")
    # Create a dictionary for storing the information about the folders
    subfolder_dict = {}
    # Loop through the subfolders
    for subfolder in os.listdir(folder):
        # Get the information from the subfolder
        subfolder_info = subfolder.split('_')
        # Get the variables from the subfolder information
        if subfolder[-3:] == 'csv':
            continue
        mean, cv, num_balls, den = subfolder_info[:4]
        mean, cv, num_balls, den = float(mean), float(cv), int(num_balls), float(den)
        # Check for more information
        try:
            num = int(subfolder_info[-1])
        except ValueError:
            num = 0
        # Verify the files
        # Get the files from the subfolder
        pdb_fl, aw_fl, pow_fl = get_files(os.path.join(folder, subfolder))

        # Quick verification
        pdb_fl = pdb_fl if pdb_fl is not None and verify_pdb(pdb_fl, simple=True) else None
        aw_fl = aw_fl if aw_fl is not None and verify_logs(aw_fl, simple=True) else None
        pow_fl = pow_fl if pow_fl is not None and verify_pdb(pow_fl, simple=True) else None

        # Create the subfolder information dictionary
        sub_info_dict = {'subfolder': subfolder, 'olap': olap, 'pbc': pbc, 'pdb': pdb_fl, 'aw': aw_fl, 'pow': pow_fl,
                         'mean': mean, 'num': num, 'ball_num': num_balls, 'cv': cv, 'den': den}
        # Add the information to the subfolder dictionary
        if (cv, den) in subfolder_dict:
            subfolder_dict[(cv, den)].append(sub_info_dict)
        else:
            subfolder_dict[(cv, den)] = [sub_info_dict]
    print("Determining the quality of the folders")
    # We need to filter through the data now to get the best of the best
    want, keep, bad, counter = [], [], [], {}
    for cv in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        for den in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]:
            try:
                my_good, my_med, my_trash = get_the_best_folders(subfolder_dict[(cv, den)])
            except KeyError:
                print("Missing completely: ", cv, den)
                continue
            want += my_good
            keep += my_med
            bad += my_trash
            # Keep a counter for the directories
            counter[(cv, den)] = 0

    # Get the outer folder
    outer_folder = os.path.dirname(folder)
    # Create the new folders
    os.mkdir(outer_folder + '/New_Data3')
    os.mkdir(outer_folder + '/Keep_Data')
    os.mkdir(outer_folder + '/Bad_Data')

    # Create a ledger to keep track of everything
    with open(outer_folder + '/file_move_ledger.txt', 'w') as fml:
        fml.write(f"Folder cleaning ledger {datetime.datetime.now()} files in folder {folder}\n")
    print("Placing the folders")
    # Now we need to change the names of the good directories
    for my_dir in want:
        # First we need to make the new directory
        new_dir_list = []
        # Check if the values want to be added
        if include_mean:
            new_dir_list.append(my_dir['mean'])
        if include_cv:
            new_dir_list.append(my_dir['cv'])
        if include_ball_num:
            new_dir_list.append(my_dir['ball_num'])
        if include_den:
            new_dir_list.append(my_dir['den'])
        if include_olap:
            new_dir_list.append(my_dir['olap'])
        if include_pbc:
            new_dir_list.append(['pbc'])
        new_dir_list.append(counter[(my_dir['cv'], my_dir['den'])])
        # Create the directory name
        new_dir = outer_folder + "/New_Data3/" + "_".join([str(_) for _ in new_dir_list])
        # Create the new directory
        os.mkdir(new_dir)
        # Move all the stuff
        move_contents(os.path.join(folder, my_dir['subfolder']), new_dir)
        # Increment the counter
        counter[(my_dir['cv'], my_dir['den'])] += 1
        # Tell the ledger what was just performed
        with open(outer_folder + '/file_move_ledger.txt', 'a') as fml:
            fml.write(f"Good Files moved from {folder + '/' + my_dir['subfolder']} to {new_dir}\n")

    # Now go through the keep data
    for my_dir in keep:
        # Create the keep move to sub directory
        new_dir = outer_folder + '/Keep_Data/' + my_dir['subfolder']
        os.mkdir(new_dir)
        # Move the data
        move_contents(os.path.join(folder, my_dir['subfolder']), new_dir)
        # Tell the ledger what was just performed
        with open(outer_folder + '/file_move_ledger.txt', 'a') as fml:
            fml.write(f"Keep Files moved from {folder + '/' + my_dir['subfolder']} to {new_dir}\n")

    # Now go through the keep data
    for my_dir in bad:
        # Create the keep move to sub directory
        new_dir = outer_folder + '/Bad_Data/' + my_dir['subfolder']
        os.mkdir(new_dir)
        # Move the data
        move_contents(os.path.join(folder, my_dir['subfolder']), new_dir)
        # Tell the ledger what was just performed
        with open(outer_folder + '/file_move_ledger.txt', 'a') as fml:
            fml.write(f"Bad Files moved from {folder + '/' + my_dir['subfolder']} to {new_dir}\n")


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    rename_folders()
