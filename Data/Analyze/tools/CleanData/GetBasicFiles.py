import os
import shutil
import tkinter as tk
from tkinter import filedialog


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


def get_basic_files(folder=None):
    if folder is None:
        folder = filedialog.askdirectory()
    place_directory = os.path.dirname(folder) + '/Basic_Data'
    os.mkdir(os.path.dirname(folder) + '/Basic_Data')
    num_folders = len([_ for _ in os.listdir(folder)])
    i = 1
    for subfolder in os.listdir(folder):

        print(f"\rCopying folder {i}/{num_folders}", end="")
        i += 1
        os.mkdir(place_directory + '/' + subfolder)
        os.mkdir(place_directory + '/' + subfolder + '/aw')
        os.mkdir(place_directory + '/' + subfolder + '/pow')

        # if os.path.exists(folder + '/' + subfolder + '/aw/edges.off') or os.path.exists(folder + '/' + subfolder + '/aw/atoms'):
        #
        for filename in ['balls.txt', 'balls.pdb', '/aw/aw_verts.txt', '/aw/aw_logs.csv', '/pow/pow_verts.txt', '/pow/pow_logs.csv', 'info.txt', 'retaining_box.off', 'set_atoms.pml', 'set_balls.pml']:
            try:
                shutil.copy2(folder + '/' + subfolder + '/' + filename, place_directory + '/' + subfolder + '/' + filename)
            except FileNotFoundError:
                continue


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    get_basic_files()

