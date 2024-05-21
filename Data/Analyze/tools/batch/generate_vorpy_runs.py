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
    run_file_num = num_tot // 500
    # Check to see if the path for the solved exists
    if not path.exists(folder + '/' + directory + '/vor'):

        num_needed += 1
        # Generate the pdb file by separating the number if there is one.
        # try:
        #     int(my_dir_split[-1])
        #     pdb_file = '_'.join(my_dir_split[:-1]) + '.pdb'
        # except ValueError:
        #     pdb_file = directory + '.pdb'
        # pdb_file = folder + '/' + directory + '/' + pdb_file
        # with open('foam_runs_' + str(run_file_num) + '.bat', 'a') as write_file:
        #     write_file.write('py vorpy.py {} -s nt compare -e dir {} -e large\n'.format(pdb_file, folder + '/' + directory))
        # print('py vorpy.py {} -s nt compare -e dir {} -e large'.format(pdb_file, folder + '/' + directory))

print('Number of directories solved: {}/{} at {}'.format(num_tot - num_needed, num_tot, datetime.datetime.now()))
