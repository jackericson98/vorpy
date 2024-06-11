import os
from os import path
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.compare.read_logs import read_logs


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
logs_folder = filedialog.askdirectory()

my_logs = []
for file, directory, x in os.walk(logs_folder):
    print(file, directory, x)
    my_logs.append(read_logs(file, return_dict=True))

# Get the totals
vols, sas = [], []
for logs in my_logs:
    vols.append(logs['group data']['volume'])
    sas.append(logs['group data']['sa'])

print(vols)
print(sas)
