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

with open('C:/Users/i7-8700/PycharmProjects/vorpy/foam_runs.bat', 'w') as foam_file:
    # We want to create a script to run all of these
    num_done = 0
    tot = 0
    for my_dir in my_dirs_unfiltered:
        # Get the settings to find the pdb within the directory
        settings = my_dir.split('_')
        try:
            dinky_winky = int(settings[-1])
            if dinky_winky >= 19:
                continue
            new_file = '_'.join(settings[:-1])

        except ValueError:
            new_file = '_'.join(settings)
        tot += 1
        run_dir = file_directory + '/' + my_dir +'/' + new_file + '.pdb'
        export_dir = file_directory + '/' + my_dir
        if os.path.exists(run_dir) and not os.path.exists(run_dir + '/vor'):
            foam_file.write('\npython3 vorpy.py {} -s nt compare -e dir {} -e large'.format(run_dir, export_dir))
        else:
            num_done += 1

print('{}/{} finished'.format(num_done, tot))