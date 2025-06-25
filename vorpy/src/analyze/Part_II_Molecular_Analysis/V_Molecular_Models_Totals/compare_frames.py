import os
from os import path
import tkinter as tk
from tkinter import filedialog
import numpy as np

import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs import read_logs

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
logs_folder = filedialog.askdirectory()

my_logs = []
for file, directory, x in os.walk(logs_folder):
    for my_file in x:
        my_logs.append(read_logs(file + '/' + my_file, return_dict=True, new_logs=True))

# Get the totals
vols, sas = [], []
for logs in my_logs:
    vols.append(logs['group data']['volume'])
    sas.append(logs['group data']['sa'])

print('Average Volume = ', sum(vols)/len(vols), ' +- ', np.std(vols)/np.sqrt(len(vols)))
print('Average Surface Area = ', sum(sas)/len(sas), ' +- ', np.std(sas)/np.sqrt(len(sas)))
