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


# Choose the column you want to plot
plotting_choice = 22

# Make the dictionary for the title name
title_dict = {0: 'File', 1: 'Average # of Overlaps', 2: 'STD # of Overlaps', 3: 'Minimum # of Overlaps',
              4: 'Maximum # of Overlaps', 5: 'Average Neighbors (AW)', 6: 'Average Neighbors (Pow)',
              7: 'Average Absolute Difference\n in Neighbor # (AW base)',
              8: 'Standard Deviation in Absolute\n Difference in Neighbor # (AW base)',
              9: 'Average Absolute Difference\n in Neighbor # (Pow base)',
              10: 'Standard Deviation in Absolute\n Difference in Neighbor # (Pow base)',
              11: 'Average Sphericity (AW)', 12: 'Average Sphericity (Pow)',
              13: 'Average Absolute Difference\n in Sphericity (AW base)',
              14: 'Standard Deviation Absolute\n Difference in Sphericity (AW base)',
              15: 'Average Absolute Difference\n in Sphericity (Pow base)',
              16: 'Standard Deviation Absolute\n Difference in Sphericity (Pow base)',
              17: 'Average Maximum Spike Distance (AW)', 18: 'Average Maximum Spike Distance (Pow)',
              19: 'Average Absolute Difference\n in Max Spike Distance (AW base)',
              20: 'Standard Deviation Absolute\n Difference Max Spike Distance (AW base)',
              21: 'Average Absolute Difference\n Max Spike Distance (Pow base)',
              22: 'Standard Deviation Absolute\n Difference Max Spike Distance (Pow base)'}
y_label_dict = {0: 'File', 1: '# of Overlaps', 2: 'STD # of Overlaps', 3: '# of Overlaps',
                4: '# of Overlaps', 5: 'Neighbors', 6: 'Neighbors',
                7: 'Absolute % Difference\n Neighbor #', 8: 'Standard Deviation in Absolute\n Difference in Neighbor #',
                9: 'Absolute % Difference\n Neighbor #',
                10: 'Standard Deviation in Absolute\n Difference in Neighbor #',
                11: 'Sphericity', 12: 'Sphericity', 13: 'Absolute % Difference Sphericity',
                14: 'Standard Deviation Absolute\n % Difference in Sphericity',
                15: 'Absolute % Difference Sphericity',
                16: 'Standard Deviation Absolute\n % Difference in Sphericity',
                17: 'Spike Distance', 18: 'Spike Distance', 19: 'Absolute % Difference \nSpike Distance',
                20: 'Standard Deviation Absolute\n % Difference Max Spike Distance',
                21: 'Absolute % Difference \nSpike Distance',
                22: 'Standard Deviation Absolute\n % Difference Max Spike Distance'}

with open(my_foams_file, 'r') as my_foam_data:

    my_data = []
    # x, y, z1, z2 = [], [], [], []
    foam_data = csv.reader(my_foam_data)
    for i, my_line in enumerate(foam_data):
        if i == 0:
            continue
        try:
            my_vals = my_line[0].split('_')
            if len(my_vals) == 1:
                continue
            my_data.append({'rad std': float(my_vals[1]), 'num balls': int(my_vals[2]), 'num': my_line[0],
                            'density': float(my_vals[3]), 'val': float(my_line[plotting_choice])})
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

for dp in my_data:
    # Check if the data has been added before
    if dp['rad std'] in lists:
        if dp['density'] in lists[dp['rad std']]:
            lists[dp['rad std']][dp['density']].append(dp['val'])
        else:
            lists[dp['rad std']][dp['density']] = [dp['val']]
    else:
        lists[dp['rad std']] = {
            dp['density']: [dp['val']]}

for _ in lists:
    print(_, [__ for __ in lists[_]])

my_densities = [round(_, 3) for _ in np.linspace(0.05, 0.5, 10)]
my_sds = [round(_, 3) for _ in np.linspace(0.05, 0.5, 10)]

# Initialize lists using list comprehensions
datavvm, datavsm, datapvm, datapsm, = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvms, datavsms, datapvms, datapsms, = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvps, datavsps, datapvps, datapsps, = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]


# Iterate over my_sds and my_densities using nested loops
for i, num in enumerate(my_sds):
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
            continue
        if plotting_choice in {0, 1, 2, 3, 4, 5, 6, 11, 12, 17, 18}:
            mean = np.mean(curr_data)
            sd = np.std(curr_data)
        else:
            mean = 100 * np.mean(curr_data)
            sd = 100 * np.std(curr_data)

        # Append means and mean +/- SD to respective lists
        datavvm[i].append(mean); datavvms[i].append(mean - sd); datavvps[i].append(mean + sd)


# Coefficient of Variation (CV) and Density values
cmap = plt.cm.rainbow  # Choose a colormap that does not have yellow and works well in grayscale
norm = Normalize(vmin=min(my_sds), vmax=max(my_sds))
sm = ScalarMappable(norm=norm, cmap=cmap)
fig, ax = plt.subplots(figsize=(8, 6))

for i, sd in enumerate(my_sds):
    # Colors for each line based on 'sd' which is used as an index into the colormap
    color = cmap(norm(sd))
    ax.plot(my_densities, datavvm[i], color=color)
    ax.fill_between(my_densities, datavvms[i], datavvps[i], color=color, alpha=0.2)

# Adding a color bar that uses the created ScalarMappable
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Coefficient of Variation (CV)', fontdict=dict(size=25))


# Set plot titles and labels
ax.set_xticks(np.arange(my_densities[0] + 0.05, my_densities[-1] + 0.05, 0.1))
# ax.set_ylim([0, 30])
ax.set_title(title_dict[plotting_choice], fontsize=20)

ax.set_xlabel('Density', fontsize=25)
ax.set_ylabel(y_label_dict[plotting_choice], fontsize=25)
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
