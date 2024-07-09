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

open_adjustments = {
    '1000': [-0.42779658559670897, 0.9522701611684636, 0.003273616493599463],
    '100': [-0.5121075296037707, 0.936176070116948, 0.004823092120831488],
    '100000': [-0.3563740799873354, 0.9518810242087892, 0.0033141283705400313]
}


def invers_square(val, abc):
    return abc[0] * val ** 2 + abc[1] * val + abc[2]


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
                            'sa diff vor': float(my_line[7]),
                            'vol diff pow': float(my_line[8]), 'sa diff pow': float(my_line[9]),
                            'num cells': int(my_line[10])})
        except ValueError:
            continue
print(my_data[0]['num'])

if 'log' in my_data[0]['num'].lower():
    plot_type = 'lognormal'
else:
    plot_type = 'gamma'

if 'closed' in my_data[0]['num'].lower() or 'false' in my_data[0]['num'].lower():
    cell_type = 'Closed'
else:
    cell_type = 'Open'

lists = {}

for dp in my_data:
    # Check if the data has been added before
    if dp['rad std'] in lists:
        if dp['density'] in lists[dp['rad std']]:
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


if plot_type == 'gamma':
    my_densities = np.arange(0.025, 0.525, 0.025)
    my_sds = np.arange(2.0, 10.5, 0.5)

if plot_type == 'lognormal':
    my_densities = np.arange(0.025, 0.525, 0.025)
    my_sds = np.arange(0.1, 0.5, 0.025)

my_densities = [round(_, 3) for _ in my_densities]
my_sds = [round(_, 3) for _ in my_sds]

# Initialize lists using list comprehensions
datavvm, datavsm, datapvm, datapsm = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvms, datavsms, datapvms, datapsms = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvps, datavsps, datapvps, datapsps = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]


# Iterate over my_sds and my_densities using nested loops
for i, num in enumerate(my_sds):
    for j, num2 in enumerate(my_densities):
        means = []
        sds = []
        # Get the current Data
        try:
            curr_data = lists[num][num2]
        except KeyError:
            continue

        for data1 in curr_data:
            # Calculate the Z-scores
            z_scores = np.abs((data1 - np.mean(data1)) / np.std(data1))

            # Set a Z-score threshold (e.g., 3)
            z_score_threshold = 1

            # Exclude outliers based on the Z-score
            filtered_data = np.array(data1)[z_scores < z_score_threshold]
            means.append(100*np.mean(filtered_data))
            sds.append(np.std([_ * 100 for _ in filtered_data]) / np.sqrt(len(filtered_data)))

        # Append means and mean +/- SD to respective lists
        datavvm[i].append(means[0]); datavvms[i].append(means[0] - sds[0]); datavvps[i].append(means[0] + sds[0])
        datavsm[i].append(means[1]); datavsms[i].append(means[1] - sds[1]); datavsps[i].append(means[1] + sds[1])
        datapvm[i].append(means[2]); datapvms[i].append(means[2] - sds[2]); datapvps[i].append(means[2] + sds[2])
        datapsm[i].append(means[3]); datapsms[i].append(means[3] - sds[3]); datapsps[i].append(means[3] + sds[3])

if cell_type == 'Open':
    my_densities = [invers_square(_, open_adjustments['1000']) for _ in my_densities]
for value in {'vol', 'sa'}:
    # Coefficient of Variation (CV) and Density values
    cmap = plt.cm.rainbow  # Choose a colormap that does not have yellow and works well in grayscale
    norm = Normalize(vmin=min(my_sds), vmax=max(my_sds))
    sm = ScalarMappable(norm=norm, cmap=cmap)
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, sd in enumerate(my_sds):
        # Colors for each line based on 'sd' which is used as an index into the colormap
        color = cmap(norm(sd))

        try:
            ax.plot(my_densities[2:], datavsm[i][2:], color=color)
            ax.fill_between(my_densities[2:], datavsms[i][2:], datavsps[i][2:], color=color, alpha=0.2)
        except ValueError:
            print(sd, my_densities[2:], datavsm[i][2:])

    # Adding a color bar that uses the created ScalarMappable
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    if plot_type == 'lognormal':
        cbar.set_label('Coefficient of Variation (CV)', fontdict=dict(size=25))
    elif plot_type == 'gamma':
        cbar.set_label('\u03b2 Value', fontdict=dict(size=20))


    # Set plot titles and labels

    ax.set_xticks(np.arange(0, 0.55, 0.1))
    ax.set_ylim([0, 60])
    ax.set_xlim([0.05, 0.55])
    # ax.set_title('{} {} Power {} % Diff'.format(plot_type.capitalize(), cell_type, {'sa': 'Surface Area', 'vol': 'Volume'}[value]), fontsize=20)
    ax.set_xlabel('Density', fontsize=25)
    ax.set_ylabel('% Difference', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=12)

    cbar.ax.tick_params(labelsize=20, size=10, width=2, length=12)
    plt.tight_layout()

    # Show the plot
    plt.show()

# Create a single plot
# fig, ax = plt.subplots(figsize=(8, 6))
# for i in range(len(datavvm)):
#     ax.plot(my_densities[2:], datavvm[i][2:], label=str(my_sds[i]))
#     ax.fill_between(my_densities[2:], datavvms[i][2:], datavvps[i][2:], alpha=0.2)
#
#
# # Set plot title and legend
# ax.set_xticks(np.arange(my_densities[1], my_densities[-1] + 0.05, 0.05))
# ax.set_title('Power Volume Deviation (Closed)', fontsize=30)
# ax.set_xlabel('Density', fontsize=20)
# ax.set_ylabel('% Difference', fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=15)
# legend = ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1))
# legend.set_title('CV')
#
# # Adjust the right margin to make room for the legend
# plt.subplots_adjust(right=0.8)
#
# # Show the plot
# plt.show()
