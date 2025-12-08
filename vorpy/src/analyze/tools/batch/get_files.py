import os
from os import path
import tkinter as tk
from tkinter import filedialog


def get_files(folder):
    """Traverse the folder and extract required PDB and log files.
    Logs are identified based on subfolder names containing 'aw' or 'pow'.
    """
    pdb_file, aw_logs, pow_logs = None, None, None

    for root_dir, sub_dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.pdb') and 'atoms' not in file and 'diff' not in file:
                pdb_file = path.join(root_dir, file)

        # Check sub-subfolders for logs
        for sub_dir in sub_dirs:
            if 'aw' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        aw_logs = path.join(root_dir, sub_dir, file)
            elif 'pow' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        pow_logs = path.join(root_dir, sub_dir, file)
    return pdb_file, aw_logs, pow_logs


def get_all_files(folder=None):
    """Goes through the folder finding all the systems and gethers their logs, but only if they are full"""
    # Get the folder from the user if not provided
    if folder is None:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory()
    # Create a dictionary to store the files
    files = {}
    # Loop through the folder finding the systems
    for subfolder in os.listdir(folder):
        # Make the new dictionary entry
        sys_letter = subfolder[0]
        sys_name = subfolder.split('_')[1]

        if all(os.path.exists(os.path.join(folder, subfolder, f"{log}_logs.csv")) for log in ['aw', 'pow', 'prm']):
            files[sys_letter] = {log: os.path.join(folder, subfolder, f"{log}_logs.csv") for log in ['aw', 'pow', 'prm']}
            files[sys_letter]['name'] = sys_name
        elif all(os.path.exists(os.path.join(folder, subfolder, f"{log}/{log}_logs.csv")) for log in ['aw', 'pow', 'prm']):
            files[sys_letter] = {log: os.path.join(folder, subfolder, f"{log}/{log}_logs.csv") for log in ['aw', 'pow', 'prm']}
            files[sys_letter]['name'] = sys_name
    print(files)
    return files


if __name__ == "__main__":
    files = get_all_files()
    for key, value in files.items():
        print(key, value)