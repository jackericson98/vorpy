import csv
import tkinter as tk
from tkinter import filedialog
import os
from Data.Analyze.tools.compare.read_logs import read_logs


def get_logs_and_pdbs(make_file=True, output_file_name=None):
    logs_pdbs = {}

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    user_data = filedialog.askdirectory(title='Get User Data')
    for rroot, directories, files in os.walk(user_data):
        for directory in directories:
            if 'aw' in directory or 'pow' in directory:
                continue
            logs_pdbs[directory] = {}
            for rrooot, dircs, filese in os.walk(rroot + '/' + directory):
                for file in filese:
                    if file[-3:] == 'pdb':
                        logs_pdbs[directory]['pdb'] = rrooot + '/' + file
                # print(filese)
                for dircy_dirc in dircs:
                    if dircy_dirc[-2:] == 'aw':
                        for rootytooty, dincretories, flies in os.walk(rrooot):
                            for file in flies:
                                if file[-3:] == 'csv' and rootytooty[-2:] == 'aw':
                                    logs_pdbs[directory]['aw'] = rootytooty + '/' + file
                                    # print(rootytooty + '/' + file)
                    elif dircy_dirc[-3:] == 'pow':
                        for rootytooty, dincretories, flies in os.walk(rrooot):
                            for file in flies:
                                if file[-3:] == 'csv' and rootytooty[-3:] == 'pow':
                                    logs_pdbs[directory]['pow'] = rootytooty + '/' + file
                # print(logs_pdbs[directory])
    if make_file:
        if output_file_name is None:
            output_file_name = 'logs_pdbs.txt'
        with open(output_file_name, 'w') as loggy_woggys:
            for _ in logs_pdbs:
                # print(logs_pdbs[_])
                if 'pdb' in logs_pdbs[_] and 'aw' in logs_pdbs[_] and 'pow' in logs_pdbs[_]:
                    loggy_woggys.write(logs_pdbs[_]['pdb'] + '\n')
                    loggy_woggys.write(logs_pdbs[_]['aw'] + '\n')
                    loggy_woggys.write(logs_pdbs[_]['pow'] + '\n')

    return logs_pdbs


data = {}

loggan_paul = get_logs_and_pdbs(False)
for directory in loggan_paul:
    if 'aw' not in loggan_paul[directory] or 'pow' not in loggan_paul[directory]:
        continue
    vals = directory.split('_')
    cv, den = float(vals[1]), float(vals[3])
    # Open the aw logs and get the maximum vertex
    with open(loggan_paul[directory]['aw'], 'r') as aw_logs:
        aw_read = csv.reader(aw_logs)
        reading = False
        vert_rads = []
        for line in aw_read:

            if reading:
                vert_rads.append(float(line[8]))
            if len(line) == 9 and line[8] == 'r':
                reading = True
            else:
                reading = False
    # Record the maximum vertex for aw
    max_aw = max(vert_rads)
    len_aw_rads = len(vert_rads)
    # Open the aw logs and get the maximum vertex
    with open(loggan_paul[directory]['pow'], 'r') as pow_logs:
        pow_read = csv.reader(pow_logs)
        reading = False
        vert_rads = []
        for line in pow_read:

            if reading:
                vert_rads.append(float(line[8]))
            if len(line) == 9 and line[8] == 'r':
                reading = True
            else:
                reading = False
    # Record the maximum vertex for aw
    max_pow = max(vert_rads)
    len_pow_rads = len(vert_rads)

    if (cv, den) in data:
        data[(cv, den)]['pow'].append(max_pow)
        data[(cv, den)]['aw'].append(max_aw)
    else:
        data[(cv, den)] = {'pow': [max_pow], 'aw': [max_aw], 'pow_verts': len_pow_rads}

for _ in data:
    print(_, 'pow', data[_]['pow'])
    print(_, 'aw', data[_]['aw'])
