import os
import tkinter as tk
from tkinter import filedialog
import shutil


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


def rename_pdb(folder=None, rename_to=None):

    if folder is None:
        folder = filedialog.askdirectory(title="Get PDB File")

    for rooot, folders, files in os.walk(folder):
        for file in files:
            # Change the main PDB name
            if file[-3:] == 'pdb' and 'atoms' not in file and 'surr' not in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'balls.pdb')

            elif file[-3:] == 'pdb' and 'surr' in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'surrounding_balls.pdb')

            elif file[-3:] == 'pdb' and 'atoms' in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'group_balls.pdb')

            elif file[-3:] == 'csv' and 'logs' in file and 'pow' in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'pow_logs.csv')

            elif file[-3:] == 'csv' and 'logs' in file and 'aw' in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'aw_logs.csv')

            elif file[-3:] == 'txt' and 'info' in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'info.txt')

            elif file[-3:] == 'off' and 'verts' in file and 'atoms' not in rooot and 'shell' not in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'verts.off')

            elif file[-3:] == 'pml' and 'set_atoms' in file:
                os.rename(rooot + '/' + file, rooot + '/' + 'set_balls.pml')


def get_nets_in_order(folder=None):

    if folder is None:
        folder = filedialog.askdirectory(title="Get folder to sort")

    for rooot, folders, files in os.walk(folder):
        os.mkdir(rooot + '/aw')
        os.mkdir(rooot + '/pow')

        for subfolder in folders:
            if 'aw' in subfolder:
                try:
                    os.rename(rooot + '/' + subfolder + '/atoms', rooot + subfolder + '/balls')
                except FileNotFoundError:
                    pass
                move_contents(rooot + '/' + subfolder, rooot + '/aw')
                os.rmdir(rooot + '/' + subfolder)

            if 'pow' in subfolder:
                try:
                    os.rename(rooot + '/' + subfolder + '/atoms', rooot + subfolder + '/balls')
                except FileNotFoundError:
                    pass
                move_contents(rooot + '/' + subfolder, rooot + '/pow')
                os.rmdir(rooot + '/' + subfolder)

        try:
            shutil.move(rooot + '/aw_verts.txt', rooot + '/aw/aw_verts.txt')
            shutil.move(rooot + '/pow_verts.txt', rooot + '/pow/pow_verts.txt')
        except FileNotFoundError:
            continue

    for file in os.listdir(folder + '/aw'):
        if 'logs' in file:
            os.rename(folder + '/aw/' + file, folder + '/aw/aw_logs.csv')

    for file in os.listdir(folder + '/pow'):
        if 'logs' in file:
            os.rename(folder + '/pow/' + file, folder + '/pow/pow_logs.csv')


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    main_folder = filedialog.askdirectory()
    length = len([name for name in os.listdir(main_folder)])
    i = 1
    for less_folder in os.listdir(main_folder):

        print(f"folder {i}/{length}")
        try:
            rename_pdb(main_folder + '/' + less_folder)
        except FileExistsError:
            pass


        try:
            get_nets_in_order(main_folder + '/' + less_folder)
        except FileExistsError:
            continue

        i += 1

