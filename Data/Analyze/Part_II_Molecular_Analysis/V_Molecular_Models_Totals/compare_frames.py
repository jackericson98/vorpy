import os
from os import path
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.compare.read_logs import read_logs
import numpy as np


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
