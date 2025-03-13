import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import shutil
import tkinter as tk
from tkinter import filedialog


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


def integrate_folders(subfolder1, subfolder2):
    # First move the vertices
    shutil.move(subfolder2 + '/pow_verts.txt', subfolder1 + '/pow/pow_verts.txt')
    shutil.move(subfolder2 + '/aw_verts.txt', subfolder1 + '/aw/aw_verts.txt')
    # Next move the logs and rename them
    shutil.move(subfolder2 + '/balls_Network_aw/balls_Network_aw_logs.csv', subfolder1 + '/aw/aw_logs.csv')
    shutil.move(subfolder2 + '/balls_Network_pow/balls_Network_pow_logs.csv', subfolder1 + '/pow/pow_logs.csv')
    # Last move the folder out of the main area
    new_bad_trash_dir = os.path.dirname(os.path.dirname(subfolder2)) + '/Bad_Data/' + os.path.basename(subfolder2)
    os.mkdir(new_bad_trash_dir)
    move_contents(subfolder2, new_bad_trash_dir)
    os.rmdir(subfolder2)


def adjust_after_runs(folder=None):
    if folder is None:
        folder = filedialog.askdirectory()
    dir_dict = {}
    for subfolder in os.listdir(folder):
        subdir_info = subfolder.split("_")
        key = tuple(subdir_info[:5])
        if key in dir_dict:
            try:
                if len(subdir_info) == 6:
                    integrate_folders(dir_dict[key], folder + '/' + subfolder)
                else:
                    integrate_folders(folder + '/' + subfolder, dir_dict[key])
            except FileNotFoundError:
                continue
        else:
            dir_dict[key] = folder + '/' + subfolder


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    adjust_after_runs()

