import os
import csv
import numpy as np
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.plot_templates.line import line_plot


# We need the two directories
root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)

# As for the open directory and then the closed directory
open_logs = filedialog.askopenfilename(title='Choose Open Logs')
closed_logs = filedialog.askopenfilename(title='Choose Closed Logs')


# We need to go through each directory and open the pdb files
def get_box_vols(foam_data_file):
    # Set up the dictionary
    value_dict = {}
    # Go through each foam_data file
    with open(foam_data_file, 'r') as foam_data:
        fdf = csv.reader(foam_data)
        for line in fdf:
            if len(line) <= 1:
                continue
            try:
                volume = float(line[1]) ** 3
            except ValueError:
                continue
            file_name = line[0]
            file_name = file_name.split('/')
            file_info = file_name[-1].split('_')
            try:
                cv, density = file_info[1], file_info[3]
            except IndexError:
                continue

            if cv in value_dict:
                if density in value_dict[cv]:
                    value_dict[cv][density].append(volume)
                else:
                    value_dict[cv][density] = [volume]
            else:
                value_dict[cv] = {density: [volume]}
    return value_dict


# We need to get the box volume and sort it into the correct designation


# for each cv we need to plot the percentage
open_dick = get_box_vols(open_logs)
closed_dick = get_box_vols(closed_logs)

combined_dick = {}
densities = []
# Go through the two dictionaries and combine the data for comparison
for _ in open_dick:
    if _ not in closed_dick:
        continue
    for __ in open_dick[_]:
        if __ not in densities:
            densities.append(__)
        if __ not in closed_dick[_]:
            closed_dick[_][__] = [np.mean(open_dick[_][__])]
        open_avg = np.mean(open_dick[_][__])
        closed_avg = np.mean(closed_dick[_][__])
        open_std_err = np.std(open_dick[_][__]) / np.sqrt(len(open_dick[_][__]))
        closed_std_err = np.std(closed_dick[_][__]) / np.sqrt(len(closed_dick[_][__]))

        diff = closed_avg - open_avg
        tot_err = closed_std_err + open_std_err

        per_diff = (closed_avg - open_avg) / closed_avg

        new_std_err = 400 * np.sqrt(((closed_avg ** 2) * open_std_err ** 2 + (open_avg ** 2) * closed_std_err ** 2) /
                                    (closed_std_err + open_std_err) ** 4)



        add_dick = {'open': open_dick[_][__], 'closed': closed_dick[_][__], 'open_avg': open_avg,
                    'closed_avg': closed_avg, 'open_err': open_std_err, 'closed_err': closed_std_err, 'diff': diff,
                    'per_diff': per_diff, 'tot_std_err': new_std_err, 'tot_err': tot_err}
        if _ in combined_dick:
            combined_dick[_][__] = add_dick
        else:
            combined_dick[_] = {__: add_dick}

# Sort the dictionary
cvs = [_ for _ in closed_dick]
cvs.sort()
densities.sort()
x_vals, y_vals, err_vals = [], [], []
for _ in cvs:
    y_vals.append([])
    err_vals.append([])
    x_vals.append([])
    for __ in densities:
        try:
            y_vals[-1].append(combined_dick[_][__]['diff'] if combined_dick[_][__]['diff'] != 0 else np.nan)
            err_vals[-1].append(combined_dick[_][__]['tot_err'])
            x_vals[-1].append(__)
        except KeyError:
            continue

# x_vals = [densities for _ in cvs]
# y_vals = [[combined_dick[_][__]['per_diff'] for __ in densities] for _ in cvs]
# err_vals = [[combined_dick[_][__]['tot_err'] for __ in densities] for _ in cvs]

# Plot the
line_plot(x_vals, y_vals)
