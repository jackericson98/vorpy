import csv
import os
from os import path
import tkinter as tk
from tkinter import filedialog
from scipy import stats
from matplotlib import pyplot as plt
from System.system import System
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from Data.Analyze.tools.compare.read_logs import read_logs
from System.sys_funcs.input.pdb import read_pdb_line
from System.sys_funcs.output.atoms import make_pdb_line
import numpy as np
from scipy.optimize import curve_fit


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
sets = {}
fig = plt.figure(figsize=(7, 5))

colors = ['purple', 'green', 'orange']
Names = {'devries.pdb': 'A. DeVries', 'gal': 'B. Gal-Or', 'lemlich.pdb': 'R. Lemlich'}
# counter = 0
# folders = []

#
# while True:
#     folders.append(filedialog.askdirectory(title='Choose the folder mf'))
#     plot_another = input('Plot Another?    ')
#     if plot_another == 'n':
#         break

counter = 0
folders = ['C:/Users/Optiplex_7060/GSU Dropbox Dropbox/John Ericson/Data/II_Macroscopic/Real_foams/1.0_0.0_1000_0.5_False_devries',
           'C:/Users\Optiplex_7060\GSU Dropbox Dropbox\John Ericson\Data\II_Macroscopic\Real_foams/1.0_0.0_1000_0.5_False_gal_or',
           'C:/Users\Optiplex_7060\GSU Dropbox Dropbox\John Ericson\Data\II_Macroscopic\Real_foams/1.0_0.0_1000_0.5_False_lemlich']
for folder in folders:
    for rroot, directys, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'pdb' and 'atoms' not in file and 'diff' not in file:
                pdb = rroot + '/' + file
            if file[-3:] == 'csv' and 'vor' in file:
                aw_logs = rroot + '/' + file
            if file[-3:] == 'csv' and 'pow' in file:
                pow_logs = rroot + '/' + file

    # Take in the PDB
    # pdb = filedialog.askopenfilename(title='Get the pdb file mf')
    my_split_pdb = pdb.split('/')[-1]
    density = Names[my_split_pdb.split('_')[5]]
    cv = my_split_pdb.split('_')[1]
    # Get the power logs
    # pow_logs = filedialog.askopenfilename(title='Get the pow logs mf')
    # Get the Voronoit Logs
    # vor_logs = filedialog.askopenfilename(title='Get the vor logs mf')


    old_logs = True
    if old_logs:
        pow = read_logs(pow_logs)
        aw = read_logs(aw_logs)
    else:
        pow = read_logs2(pow_logs)
        aw = read_logs2(aw_logs)
    my_sys = System(pdb, simple=True)
    rads, diffs, names = [], [], []
    # Loop through the pdb file
    with open(pdb, 'r') as pdb_reader, open(pdb[:-4] + '_vor_diff_colored.pdb', 'w') as pdb_writer:
        mini, maxi = np.inf, -np.inf
        # Loop through the lines
        for i, line in enumerate(pdb_reader.readlines()):

            if i == 0:
                pdb_writer.write(line)
                continue
            npl = read_pdb_line(line)
            # Check that the number is actually in the dataframe

            # Get the pow atom and the vor atom
            try:
                if old_logs:
                    pow_atom = pow['atoms'].loc[pow['atoms']['num'] == i - 1].to_dict('records')[0]
                    vor_atom = aw['atoms'].loc[aw['atoms']['num'] == i - 1].to_dict('records')[0]
                else:
                    pow_atom = pow['atoms'].loc[pow['atoms']['Index'] == i - 1].to_dict('records')[0]
                    vor_atom = aw['atoms'].loc[aw['atoms']['Index'] == i - 1].to_dict('records')[0]
            except IndexError:
                continue
            if old_logs:
                vol_diff = ((pow_atom['volume'] - vor_atom['volume'])) / vor_atom['volume']
            else:
                vol_diff = ((pow_atom['Volume'] - vor_atom['Volume'])) / vor_atom['Volume']

            # print(density, i, pow_atom['volume'], vor_atom['volume'])
            # Check for crazy volume difference and trigger a volume difference
            if vol_diff >= 10:
                continue
            diffs.append(vol_diff)
            rads.append(npl['temperature_factor'])
            names.append(npl['atom_name'])
            if vol_diff < mini:
                mini = vol_diff
            if vol_diff > maxi:
                maxi = vol_diff
    print(len(diffs))
    if counter == 0:
        max_rad = max(rads)
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    try:
        color = colors[counter]
    except IndexError:
        color = None
    # Plot the radius to difference values
    if counter == 0:
        plt.plot([min(rads), max(rads)], [0, 0], c='k')
    mean = np.mean([100 * _ for _ in diffs])
    abs_mean = np.mean([100 * abs(_) for _ in diffs])
    print(f"plot: {density} --- Mean = {mean}, Absolute Mean = {abs_mean}")
    # plt.plot([min(rads), max(rads)], [mean, mean], linestyle=':', c=color, linewidth=3)
    plt.plot([min(rads), max_rad], [abs_mean, abs_mean], linestyle='--', c=color, linewidth=1)
    popt, pcov = curve_fit(func, np.array(rads), np.array(diffs))
    # plt.text(s='y = {:.2f} * exp(-{:.2f} * x) + {:.2f}'.format(*popt), x=1.5, y=1.5, font=dict(size=10))
    plt.plot([max_rad] + rads, [100 * _ for _ in func(np.array([max_rad] + rads), *popt)], c=color, linewidth=1, label=density)

    plt.scatter(rads, [100 * _ for _ in diffs], s=2, alpha=0.3, c=color)

    counter += 1

plt.xlabel('Radius', fontdict=dict(size=30))
plt.ylabel('% Difference', fontdict=dict(size=30))
# plt.ti
plt.xticks(font=dict(size=20))
plt.yticks(font=dict(size=20))
# plt.ylim([-17, 27])
plt.title('AW vs Power Volume', font=dict(size=30))
plt.tick_params(axis='both', width=2, length=12)
plt.legend(fontsize=25)

plt.tight_layout()
plt.show()
