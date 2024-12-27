import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename
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
              7: 'Avg Abs % Diff\nNeighbor # (AW base)',
              8: 'STD Abs % Diff\nNeighbor # (AW base)',
              9: 'Avg Abs % Diff\nNeighbor # (Pow base)',
              10: 'STD Abs % Diff\nNeighbor # (Pow base)',
              11: 'Average Sphericity (AW)', 12: 'Average Sphericity (Pow)',
              13: 'Avg Abs % Diff\nSphericity (AW base)',
              14: 'STD Abs % Diff\nSphericity (AW base)',
              15: 'Avg Abs % Diff\nSphericity (Pow base)',
              16: 'STD Abs % Diff\nSphericity (Pow base)',
              17: 'Avg Max Spike Dist (AW)', 18: 'Avg Max Spike Dist (Pow)',
              19: 'Avg Abs % Diff\nMax Spike Dist (AW base)',
              20: 'STD Abs % Diff\nMax Spike Dist (AW base)',
              21: 'Avg Abs % Diff\nMax Spike Dist (Pow base)',
              22: 'STD Abs % Diff\nMax Spike Dist (Pow base)'}
y_label_dict = {0: 'File', 1: '# of Overlaps', 2: 'STD # of Overlaps', 3: '# of Overlaps',
                4: '# of Overlaps', 5: 'Neighbors', 6: 'Neighbors',
                7: 'Abs % Diff', 8: 'STD Abs % Diff',
                9: 'Abs % Diff',
                10: 'STD Abs % Diff',
                11: 'Sphericity', 12: 'Sphericity', 13: 'Abs % Diff',
                14: 'STD Abs % Diff',
                15: 'Abs % Diff',
                16: 'STD Abs % Diff',
                17: 'Spike Distance', 18: 'Spike Distance', 19: 'Abs % Diff',
                20: 'STD Abs % Diff',
                21: 'Abs % Diff',
                22: 'STD Abs % Diff'}

y_ranges = {0: [None, None], 1: [0, 10], 2: [0, 5], 3: [0, 3], 4: [0, 30], 5: [11.5, 16], 6: [11.5, 16], 7: [0, 50],
            8: [0, 50], 9: [0, 50], 10: [0, 50], 11: [0.80, 0.91], 12: [0.80, 0.91], 13: [0, 10], 14: [0, 12], 15: [0, 12],
            16: [0, 12], 17: [0, 7.5], 18: [0, 7.5], 19: [0, 50], 20: [0, 39], 21: [0, 34], 22: [0, 39]}
adjustments = {}
good_plots = []
for j in range(23):
    plotting_choice = j
    if j in good_plots:
        continue
    print('\n', j, title_dict[j])
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
                my_data.append({'rad std': float(my_vals[1]), 'num balls': int(my_vals[2]), 'mean': my_line[0],
                                'density': float(my_vals[3]), 'olap': float(my_vals[4]), 'dist': my_vals[5],
                                'PBC': my_vals[7], 'val': float(my_line[plotting_choice])})
            except ValueError:
                continue
            except IndexError:
                continue
    try:
        plot_type = 'gamma'
        cell_type = my_data[-1]['olap']
    except IndexError:
        continue

    # print(plot_type, cell_type)

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
                lists[dp['rad std']][dp['density']].append(dp['val'])
            else:
                lists[dp['rad std']][dp['density']] = [dp['val']]
        else:
            lists[dp['rad std']] = {
                dp['density']: [dp['val']]}
    # for _ in lists:
    #     print(_, [__ for __ in lists[_]])

    my_densities.sort()
    my_sds.sort()

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

            for data1 in curr_data:
                # Calculate the Z-scores
                z_scores = np.abs((data1 - np.mean(data1)) / np.std(data1))

                # Set a Z-score threshold (e.g., 3)
                z_score_threshold = 1

                # Exclude outliers based on the Z-score
                filtered_data = np.array(data1)[z_scores < z_score_threshold]
                means.append(100 * np.mean(filtered_data))
                sds.append(np.std([_ * 100 for _ in filtered_data]))

            # Append means and mean +/- SD to respective lists
            datavvm[i].append(mean); datavvms[i].append(mean - sd); datavvps[i].append(mean + sd)


    # Coefficient of Variation (CV) and Density values
    cmap = plt.cm.rainbow  # Choose a colormap that does not have yellow and works well in grayscale
    norm = Normalize(vmin=min(my_sds), vmax=max(my_sds))
    sm = ScalarMappable(norm=norm, cmap=cmap)
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, sd in enumerate(my_sds):
        # Skip the no overlap low cv average sphericity
        if cell_type == 0.0 and plotting_choice in {11, 12} and sd == 0.05:
            continue
        # Colors for each line based on 'sd' which is used as an index into the colormap
        color = cmap(norm(sd))
        if cell_type == 0.0 and plotting_choice in {1, 2, 3, 4}:
            ax.plot(my_densities, [0 for _ in my_densities], color=color)
        else:
            ax.plot(my_densities, datavvm[i], color=color)
            ax.fill_between(my_densities, datavvms[i], datavvps[i], color=color, alpha=0.2)

    # Adding a color bar that uses the created ScalarMappable
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Coefficient of Variation (CV)', fontdict=dict(size=25))


    # Set plot titles and labels
    ax.set_xticks(np.arange(my_densities[0] + 0.05, my_densities[-1] + 0.05, 0.1))
    ax.set_ylim(y_ranges[plotting_choice])
    ax.set_title(title_dict[plotting_choice], fontsize=30)

    ax.set_xlabel('Density', fontsize=25)
    ax.set_ylabel(y_label_dict[plotting_choice], fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=12)

    cbar.ax.tick_params(labelsize=20, size=10, width=2, length=12)
    plt.tight_layout()

    # Show the plot
    plt.pause(0.1)

    # Ask the user whether to save the figure
    counter = 0
    while True:
        save_response = input('Save figure{}? (yes/fix notes): '.format('' if counter == 0 else ' again')).strip().lower()
        if save_response in {'yes', 'y'}:
            root = tk.Tk()
            root.withdraw()
            root.wm_attributes('-topmost', 1)
            file_path = asksaveasfilename(defaultextension=".png",
                                          filetypes=[("PNG files", "*.png"),
                                                     ("SVG files", "*.svg")])
            if file_path:
                try:
                    fig.savefig(file_path)  # Save the figure to the chosen path
                    print(f'Figure saved: {file_path}')
                except FileNotFoundError:
                    print("Something went wrong when trying to save your file:\n  {}".format(file_path))
            good_plots.append(plotting_choice)
        elif save_response in {'n', 'no'}:
            good_plots.append(plotting_choice)
            break
        else:
            print(f'{save_response}')
            if counter == 0 or 'Notes' in save_response:
                adjustments[plotting_choice] = save_response
            break
        counter += 1

    plt.close(fig)  # Close the figure

print("My plot adjustments:  \n")
for _ in adjustments:
    print(_, title_dict[_])
    print(f'\n   {adjustments[_]}\n\n')

print('These plots are done:\n\n', *good_plots)
