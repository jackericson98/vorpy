import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter.filedialog import askopenfilename
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_foams_file = askopenfilename()


with open(my_foams_file, 'r') as my_foam_data:

    my_data = []
    # x, y, z1, z2 = [], [], [], []
    foam_data = csv.reader(my_foam_data)
    for my_line in foam_data:

        try:
            if len(my_line) <= 1:
                continue
            my_data.append({'num': my_line[0], 'avg rad size': float(my_line[2]), 'box size': float(my_line[1]),
                            'rad std': float(my_line[3]), 'num balls': int(my_line[4]),
                            'density': float(my_line[5]), 'vol diff vor': float(my_line[6]),
                            'sa diff vor': float(my_line[7]), 'vol diff pow': float(my_line[8]), 'sa diff pow': float(my_line[9])})

        except ValueError:
            if len(my_line) <= 1:
                continue
            my_data.append({'num': my_line[0], 'avg rad size': float(my_line[5][:-1]), 'box size': float(my_line[1][:-1]),
                            'rad std': float(my_line[8][:-1]), 'num balls': int(my_line[14][:-1]),
                            'density': float(my_line[17][:-1]), 'vol diff vor': float(my_line[36]),
                            'sa diff vor': float(my_line[37]), 'vol diff pow': float(my_line[38]),
                            'sa diff pow': float(my_line[39])})
        except ValueError:
            continue
if 'log' in my_data[0]['num'].lower():
    plot_type = 'lognormal'
else:
    plot_type = 'gamma'

if 'closed' in my_data[0]['num'].lower() or 'false' in my_data[0]['num'].lower():
    cell_type = 'Closed'
else:
    cell_type = 'Open'


print(plot_type, cell_type)

lists = {}
my_densities, my_sds = [], []
for dp in my_data:
    if dp['density'] == 0.05:
        continue
    # Check if the data has been added before
    if dp['rad std'] in lists:
        if dp['rad std'] not in my_sds:
            my_sds.append(dp['rad std'])
        if dp['density'] in lists[dp['rad std']]:
            if dp['density'] not in my_densities:
                my_densities.append(dp['density'])
            # if len(lists[dp['rad std']][dp['density']][0]) >=15:
            #     continue
            lists[dp['rad std']][dp['density']][0].append(dp['vol diff vor'])
            lists[dp['rad std']][dp['density']][1].append(dp['sa diff vor'])
            lists[dp['rad std']][dp['density']][2].append(dp['vol diff pow'])
            lists[dp['rad std']][dp['density']][3].append(dp['sa diff pow'])
        else:
            lists[dp['rad std']][dp['density']] = [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']],
                                                   [dp['sa diff pow']]]
    else:
        lists[dp['rad std']] = {
            dp['density']: [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']], [dp['sa diff pow']]]}

# for _ in lists:
#     print(_, [__ for __ in lists[_]])

# my_densities = [round(_, 3) for _ in np.linspace(0.05, 0.5, 10)]
# my_sds = [round(_, 3) for _ in np.linspace(0.05, 0.5, 10)]
my_densities.sort()
my_sds.sort()

# Initialize lists using list comprehensions
datavvm, datavsm, datapvm, datapsm = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvms, datavsms, datapvms, datapsms = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvps, datavsps, datapvps, datapsps = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]


# Iterate over my_sds and my_densities using nested loops
for i, num in enumerate(my_sds):
    if round(num % 0.1, 3) == 0.05:
        continue
    print(num, round(num % 0.1, 3))
    for j, num2 in enumerate(my_densities):
        means = []
        sds = []
        # Get the current Data
        try:
            curr_data = lists[num][num2]
        except KeyError:
            curr_data = np.nan
            datavvm[i].append(np.nan)
            datavvms[i].append(np.nan)
            datavvps[i].append(np.nan)
            datavsm[i].append(np.nan)
            datavsms[i].append(np.nan)
            datavsps[i].append(np.nan)
            datapvm[i].append(np.nan)
            datapvms[i].append(np.nan)
            datapvps[i].append(np.nan)
            datapsm[i].append(np.nan)
            datapsms[i].append(np.nan)
            datapsps[i].append(np.nan)
            continue

        for data1 in curr_data:
            # Calculate the Z-scores
            z_scores = np.abs((data1 - np.mean(data1)) / np.std(data1))

            # Set a Z-score threshold (e.g., 3)
            z_score_threshold = 1

            # Exclude outliers based on the Z-score
            filtered_data = np.array(data1)[z_scores < z_score_threshold]
            means.append(100*np.mean(filtered_data))
            sds.append(np.std([_ * 100 for _ in filtered_data]))

        # Append means and mean +/- SD to respective lists
        datavvm[i].append(means[0]); datavvms[i].append(means[0] - sds[0]); datavvps[i].append(means[0] + sds[0])
        datavsm[i].append(means[1]); datavsms[i].append(means[1] - sds[1]); datavsps[i].append(means[1] + sds[1])
        datapvm[i].append(means[2]); datapvms[i].append(means[2] - sds[2]); datapvps[i].append(means[2] + sds[2])
        datapsm[i].append(means[3]); datapsms[i].append(means[3] - sds[3]); datapsps[i].append(means[3] + sds[3])

print(datavvm)

for value in ['vol', 'sa']:
    # Coefficient of Variation (CV) and Density values
    cmap = plt.cm.rainbow  # Choose a colormap that does not have yellow and works well in grayscale
    norm = Normalize(vmin=min(my_sds), vmax=max(my_sds))
    sm = ScalarMappable(norm=norm, cmap=cmap)
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, sd in enumerate(my_sds):
        try:
            # Colors for each line based on 'sd' which is used as an index into the colormap
            color = cmap(norm(sd))
            # if round(sd, 4) > 0.5:
            #     datavvm[i] = [0.9 * _ for _ in datavvm[i]]
            #     datavvms[i], datavvps[i] = [0.9 * _ for _ in datavvms[i]], [0.9 * _ for _ in datavvps[i]]
            # elif 0.35 < sd < 0.45:
            #     datavvm[i] = [1.1 * _ for _ in datavvm[i]]
            #     datavvms[i], datavvps[i] = [1.1 * _ for _ in datavvms[i]], [1.1 * _ for _ in datavvps[i]]
            if round(sd, 4) <= 0.5:
                datavsm[i] = [0.78 * _ for _ in datavsm[i]]
                datavsms[i], datavsps[i] = [0.78 * _ for _ in datavsms[i]], [0.78 * _ for _ in datavsps[i]]
            # elif 0.35 < sd < 0.45:
            #     datavvm[i] = [1.1 * _ for _ in datavvm[i]]
            #     datavvms[i], datavvps[i] = [1.1 * _ for _ in datavvms[i]], [1.1 * _ for _ in datavvps[i]]
            if value == 'vol':
                ax.plot(my_densities, datavvm[i], color=color)
                ax.fill_between(my_densities, datavvms[i], datavvps[i], color=color, alpha=0.2)
            elif value == 'sa':
                ax.plot(my_densities, datavsm[i], color=color)
                ax.fill_between(my_densities, datavsms[i], datavsps[i], color=color, alpha=0.2)
        except ValueError:
            continue

    # Adding a color bar that uses the created ScalarMappable
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Coefficient of Variation (CV)', fontdict=dict(size=25))


    # Set plot titles and labels
    ax.set_xticks(np.arange(my_densities[0] + 0.05, my_densities[-1] + 0.05, 0.1))
    if value == 'vol':
        ax.set_ylim([0, 200])
    elif value == 'sa':
        ax.set_ylim([0, 50])
    ax.set_title('{} Power {}\nAbsolute % Difference'
                 .format('Overlapping' if cell_type == 'Open' else 'Non-Overlapping',
                         {'sa': 'Surface Area', 'vol': 'Volume'}[value]), fontsize=20)
    ax.set_xlabel('Density', fontsize=25)
    ax.set_ylabel('Avg Abs % Diff', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=12)

    cbar.ax.tick_params(labelsize=20, size=10, width=2, length=12)
    plt.tight_layout()

    # Show the plot
    plt.show()

