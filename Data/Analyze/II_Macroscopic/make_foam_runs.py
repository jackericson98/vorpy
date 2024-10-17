import tkinter as tk
from tkinter import filedialog
import os


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
file_directory = filedialog.askdirectory()


my_dirs_unfiltered = []
# Get the directories in the data directory
for my_dir in os.listdir(file_directory):
    my_dirs_unfiltered.append(my_dir)

strings = []

# Get the directory that this is in

thine_dir = os.getcwd()

OS = 'windows'

if OS == 'linux':
    # We want to create a script to run all of these
    num_done = 0
    tot = 0
    for my_dir in my_dirs_unfiltered:
        # Get the settings to find the pdb within the directory
        settings = my_dir.split('_')
        try:
            dinky_winky = int(settings[-1])
            if dinky_winky > 19:
                continue
            new_file = '_'.join(settings[:-1])

        except ValueError:
            new_file = '_'.join(settings)
        tot += 1
        run_dir = file_directory + '/' + my_dir +'/' + new_file + '.pdb'
        export_dir = file_directory + '/' + my_dir
        if settings[3] == '0.05':
            print('I am counting these')
            mv = 200
        elif settings[3] == '0.1':
            mv = 100
        else:
            mv = 50
        if os.path.exists(export_dir) and not os.path.exists(export_dir + '/' + new_file + '_Network_aw'):

            strings.append('\npython3 vorpy.py {} -s mv {}} -s nt compare -e dir {} -e logs -g chain a'.format(run_dir, mv, export_dir))
        else:
            num_done += 1

    j = 0
    for i in range(len(strings)):
        pass
        if j < i // 400 or i == 0:
            j = i // 400
            with open(thine_dir + '/foam_runs_{}.sh'.format(j), 'w') as foam_write:
                foam_write.write('#!/bin/sh\n')
                foam_write.write(strings[i])
        else:
            with open(thine_dir + '/foam_runs_{}.sh'.format(j), 'a') as foam_write:
                foam_write.write(strings[i])

elif OS == 'windows':
    # We want to create a script to run all of these
    num_done = 0
    tot = 0
    strings = []
    for my_dir in my_dirs_unfiltered:
        # Get the settings to find the pdb within the directory
        settings = my_dir.split('_')
        try:
            dinky_winky = int(settings[-1])
            if dinky_winky > 19:
                continue
            new_file = '_'.join(settings[:-1])

        except ValueError:
            new_file = '_'.join(settings)
        tot += 1
        run_dir = file_directory + '/' + my_dir + '/' + new_file + '.pdb'
        export_dir = file_directory + '/' + my_dir
        if settings[3] == '0.05':
            print('I am counting these')
            mv = 200
        elif settings[3] == '0.1':
            mv = 100
        else:
            mv = 50
        if os.path.exists(export_dir) and not os.path.exists(export_dir + '/' + new_file + '_Network_aw'):
            strings.append('\npy vorpy.py {} -s nt compare -e dir {} -s mv {} -g chain a'.format(run_dir, export_dir, mv))
        else:
            num_done += 1

    j = 0
    for i in range(len(strings)):
        if j < i // 220 or i == 0:
            j = i // 220
            with open(thine_dir + '/foam_runs_{}.bat'.format(j), 'w') as foam_write:
                foam_write.write(strings[i])
        else:
            with open(thine_dir + '/foam_runs_{}.bat'.format(j), 'a') as foam_write:
                foam_write.write(strings[i])


print('{}/{} finished'.format(num_done, tot))
