import tkinter as tk
from datetime import datetime
from tkinter import filedialog
import os
import platform


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)


# Try to open up the foam_gen user_data file
try:
    file_directory = filedialog.askdirectory(initialdir='../foam_gen/Data/user_data')
except:
    file_directory = filedialog.askdirectory()


my_dirs_unfiltered = []
# Get the directories in the data directory
for my_dir in os.listdir(file_directory):
    my_dirs_unfiltered.append(my_dir)

strings = []

# Get the directory that this is in

thine_dir = os.getcwd()

# Detect OS
if platform.system() == "Windows":
    OS = "windows"
else:
    OS = "linux"

run_dirs, numbers = [], []
# We want to create a script to run all of these
num_done = 0
tot = 0
for my_dir in my_dirs_unfiltered:
    # Get the settings to find the pdb within the directory
    settings = my_dir.split('_')
    try:
        dinky_winky = int(settings[-1])
        number = dinky_winky
        if dinky_winky > 19:
            continue
        new_file = '_'.join(settings[:-1])
        export_type = 'logs'
    except ValueError:
        export_type = 'large'
        new_file = '_'.join(settings)
        number = 0
    tot += 1
    run_dir = file_directory + '/' + my_dir + '/' + new_file + '.pdb'
    export_dir = file_directory + '/' + my_dir

    if len(settings) < 4:
        print(settings)
        continue
    if float(settings[3]) == 0.05:
        max_vert = 150
    elif float(settings[3]) <= 0.25:
        max_vert = 100
    elif float(settings[3]) <= 0.35:
        max_vert = 60
    elif float(settings[3]) <= 0.45:
        max_vert = 30
    else:
        max_vert = 25
    # Check if the folder for AW exists, aka the network is Done
    if (not os.path.exists(export_dir + '/chain_a_aw') and '.csv' not in export_dir) and not os.path.exists(export_dir + '/' + new_file + '_Network_aw'):
        # Check if the vertices have been solved
        if os.path.exists(export_dir + '/verts.txt'):
            strings.append('\npy vorpy.py {} -s mv {} -s nt compare -e dir {} -e {} -l verts {}'
                           .format(run_dir, max_vert, export_dir, export_type, export_dir + '/verts.txt'))
        else:
            strings.append('\npy vorpy.py {} -s mv {} -s nt compare -e dir {} -e {} -g chain a'.format(run_dir, max_vert, export_dir, export_type))
        numbers.append(number)
    else:
        num_done += 1

# Sort the strings by the las number on their
strings = [x for _, x in sorted(zip(numbers, strings), key=lambda _: _)]

# Define chunk size for how many strings per file
chunk_size = 400
num_files = (len(strings) + chunk_size - 1) // chunk_size  # Calculate number of files

# Initialize file writers and create the files
file_handles = []
for j in range(num_files):
    file_name = f"{thine_dir}/foam_runs_{j}.{'sh' if OS == 'linux' else 'bat'}"
    mode = 'w'  # Write mode for initial creation
    file_handles.append(open(file_name, mode))

# Write the strings evenly into the files
for i, string in enumerate(strings):
    # Distribute first 'important' strings sequentially across all files
    file_index = i % num_files if i < num_files else i // chunk_size
    foam_write = file_handles[file_index]

    # For Linux files, add a header only once per file
    if foam_write.tell() == 0 and OS == 'linux':
        foam_write.write('#!/bin/sh\n')

    # Write the current string to the appropriate file
    foam_write.write(string)

# Close all file handles
for foam_write in file_handles:
    foam_write.close()


print(f"{num_done}/{tot} finished at {datetime.now().strftime('%Y-%m-%d %I:%M %p')}")
