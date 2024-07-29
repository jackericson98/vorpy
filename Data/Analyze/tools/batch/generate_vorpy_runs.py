"""
1. Go through each directory in the directory we want
"""
import os
from os import path
import tkinter as tk
from tkinter import filedialog
import datetime


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
folder = filedialog.askdirectory()

keywords = ('lognormal', 'True')

num_tot, num_needed = 0, 0
type_dir = {}
for directory in os.listdir(folder):
    # Split the directory name by underscores
    my_dir_split = directory.split('_')
    # Look for specified keywords to only generate runs we are interested in
    skip_dir = False
    for word in keywords:
        if word not in my_dir_split:
            skip_dir = True
    if skip_dir:
        continue
    num_tot += 1
    run_file_num = num_needed // 341
    # Check to see if the path for the solved exists
    if not path.exists(folder + '/' + directory + '/vor'):

        num_needed += 1

        # split_dir = directory.split('_')
        # if split_dir[1] in type_dir:
        #     if split_dir[3] in type_dir[split_dir[1]]:
        #         type_dir[split_dir[1]][split_dir[3]]['count'] += 1
        #         type_dir[split_dir[1]][split_dir[3]]['dirs'].append(directory)
        #     else:
        #         type_dir[split_dir[1]][split_dir[3]] = {'count': 1, 'dirs': [directory]}
        # else:
        #     type_dir[split_dir[1]] = {split_dir[3]: {'count': 1, 'dirs': [directory]}}
        # # # # Generate the pdb file by separating the number if there is one.
        # try:
        #     int(my_dir_split[-1])
        #     pdb_file = '_'.join(my_dir_split[:-1]) + '.pdb'
        # except ValueError:
        #     pdb_file = directory + '.pdb'
        # pdb_file = folder + '/' + directory + '/' + pdb_file
        # with open('foam_runs_' + str(run_file_num) + '.bat', 'a') as write_file:
        #     write_file.write('py vorpy.py {} -s nt compare -e dir {} -e logs -s mv 75\n'.format(pdb_file, folder + '/' + directory))
        # # # print('py vorpy.py {} -s nt compare -e dir {} -e large'.format(pdb_file, folder + '/' + directory))

print('Number of directories solved: {}/{} at {}'.format(num_tot - num_needed, num_tot, datetime.datetime.now()))

# run_strs = []
# for _ in type_dir:
#     for __ in type_dir[_]:
#         my_dirs = type_dir[_][__]['dirs']
#         my_count = type_dir[_][__]['count']
#         print(_, __, type_dir[_][__]['count'])
#         temp_run_strs = []
#         for directory in my_dirs:
#             my_dir_split = directory.split('_')
#             try:
#                 int(my_dir_split[-1])
#                 pdb_file = '_'.join(my_dir_split[:-1]) + '.pdb'
#             except ValueError:
#                 pdb_file = directory + '.pdb'
#             pdb_file = folder + '/' + directory + '/' + pdb_file
#             temp_run_strs.append('py vorpy.py {} -s nt compare -e dir {} -e logs -s mv {}\n'
#                                  .format(pdb_file, folder + '/' + directory, 50 + 10 / float(my_dir_split[3])))
#
#         if my_count >= 19:
#             run_strs.insert(0, temp_run_strs[0])
#             run_strs.insert(0, temp_run_strs[1])
#             run_strs += temp_run_strs[2:]
#         else:
#             run_strs += temp_run_strs
#
# num_files = 12
# for i, _ in enumerate(run_strs):
#     run_file_num = i % 12
#     with open('foam_runs_' + str(run_file_num) + '.bat', 'a') as write_file:
#         write_file.write(_)
