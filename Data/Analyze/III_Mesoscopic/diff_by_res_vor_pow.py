import csv
import os.path
import tkinter as tk
from tkinter import filedialog
import numpy as np

from Data.Analyze.tools.compare.compare_files import compare_files
import matplotlib.pyplot as plt


def percent_diff(pdb_log_dir, file_name, sa_vol='vol', martini=False):
    # Set up the martini variable so that if it is needed we can add it
    martini_list, martini_name = [], []
    if martini:
        martini_list = ['_martini']
        martini_name = ['Martini']
    # Create the list of names
    names = ['', '_ad', '_ncap', '_scbb_ad', '_scbb_ncap'] + martini_list
    # Create the directories and duplicate them for the vor and power logs
    pdb_files = [pdb_log_dir + file_name + name + '.pdb' for name in names for _ in range(2)]
    # Create the log files lists
    log_files = [pdb_log_dir + file_name + name + _ + '_logs.csv' for name in names for _ in ['_vor', '_pow']]
    # Check to see if the data has been processed yet
    if not os.path.exists(pdb_log_dir + file_name + '_residue_data.csv'):
        # Get the data from the log files and sort it
        my_info = compare_files(pdb_files, log_files, avg_distros=True, by_residues=True)
        # Put the data into a csv file for later access
        with open(pdb_log_dir + file_name + '_residue_data.csv', 'w') as res_file:
            res_fl = csv.writer(res_file)
            res_fl.writerow(['file', 'residue type', 'name', 'residue', 'volume', 'surface area'])
            for file in my_info['residues']:
                for res_type in my_info['residues'][file]:
                    for res_name in my_info['residues'][file][res_type]:
                        for res in my_info['residues'][file][res_type][res_name]:
                            res_fl.writerow([file, res_type, res_name, res] + [my_info['residues'][file][res_type][res_name][res][_] for _ in my_info['residues'][file][res_type][res_name][res]])
    # Start re-sorting the data
    vols, sas = {}, {}
    with open(pdb_log_dir + file_name + '_residue_data.csv', 'r') as res_file:
        res_reader = csv.reader(res_file)
        for i, line in enumerate(res_reader):
            if i == 0:
                continue
            if len(line) == 0 or line[1] == 'other':
                continue
            # File
            if line[0] not in vols:
                vols[line[0]] = {}
                sas[line[0]] = {}
            # Residue Type
            if line[1] not in vols[line[0]]:
                vols[line[0]][line[1]] = {}
                sas[line[0]][line[1]] = {}
            # Residue Class
            if line[2] not in vols[line[0]][line[1]]:
                vols[line[0]][line[1]][line[2]] = []
                sas[line[0]][line[1]][line[2]] = []
            vols[line[0]][line[1]][line[2]].append(line[4])
            sas[line[0]][line[1]][line[2]].append(line[5])

    vol_res_data = {}
    sa_res_data = {}
    res_names = []
    for file in vols:
        vol_res_data[file] = {}
        sa_res_data[file] = {}
        for res_type in vols[file]:
            for res_name in vols[file][res_type]:
                res_names.append(res_name)
                vol_by_type = [float(_) for _ in vols[file][res_type][res_name]]
                vol_res_data[file][res_name] = {'avg': np.mean(vol_by_type), 'max': max(vol_by_type),
                                                'min': min(vol_by_type), 'sd': np.std(vol_by_type), 'data': vol_by_type}
                sa_by_type = [float(_) for _ in sas[file][res_type][res_name]]
                sa_res_data[file][res_name] = {'avg': np.mean(sa_by_type), 'max': max(sa_by_type),
                                               'min': min(sa_by_type), 'sd': np.std(sa_by_type), 'data': sa_by_type}

    # Sample data
    labels = ['Atoms', 'Avg Dist', 'Encapsulate', 'SC/BB AD', 'SC/BB Encap.'] + martini_name

    # Get the percent difference for each residue and average
    means, std_errs = [0], [0]
    for file in sas:
        if file == file_name:
            continue
        my_perc_diff = []
        try:
            for res_name in vols[file]['nucs']:
                for i in range(len(vols[file]['nucs'][res_name])):
                    my_perc_diff.append(abs(float(vols[file_name]['nucs'][res_name][i]) - float(
                        vols[file]['nucs'][res_name][i])) / float(vols[file_name]['nucs'][res_name][i]))
        except KeyError:
            for res_name in vols[file]['amines']:
                for i in range(len(vols[file]['amines'][res_name])):
                    my_perc_diff.append(abs(float(vols[file_name]['amines'][res_name][i]) - float(
                        vols[file]['amines'][res_name][i])) / float(vols[file_name]['amines'][res_name][i]))
        means.append(np.mean(my_perc_diff) * 100)
        std_errs.append(np.std(my_perc_diff) / np.sqrt(len(my_perc_diff)) * 100)

    # Your existing code
    std1s = std_errs[::2]
    std2s = std_errs[1::2]
    mean1s = means[::2]
    mean2s = means[1::2]
    ymax = max(means)
    bar_width = 0.35
    x = range(len(labels))

    # Create the bar graph
    bar1 = plt.bar(x, mean1s, width=bar_width, label='Additively Weighted', color='skyblue', edgecolor='black')
    bar2 = plt.bar([i + bar_width for i in x], mean2s, width=bar_width, label='Power', color='orange',
                   edgecolor='black')

    # Add error bars with custom style
    for i, bar in enumerate(bar1):
        plt.errorbar(bar.get_x() + bar.get_width() / 2, mean1s[i], yerr=std1s[i], capsize=5, capthick=2,
                     color='black', alpha=0.8)
    for i, bar in enumerate(bar2):
        plt.errorbar(bar.get_x() + bar.get_width() / 2, mean2s[i], yerr=std2s[i], capsize=5, capthick=2,
                     color='black', alpha=0.8)

    value = 'Volume'
    if sa_vol == 'sa':
        value = 'Surface Area'
    # Add labels and title
    plt.ylabel('% Difference ' + value, fontdict=dict(size=15))
    plt.title(file_name + 'Average Residue Variance ({})'.format(sa_vol.upper()), fontdict=dict(size=20))

    # Angle the labels and add values at the top of the bars
    plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right', font=dict(size=10))
    for i, v in enumerate(mean1s):
        my_height = v / 2
        if i in {0, 1, 3, 4}:
            my_height = v + 25
        plt.text(i, my_height, str(round(v, 2)) + ' %', ha='center', va='center', rotation=90, font=dict(size=12))
    for i, v in enumerate(mean2s):
        my_height = v / 2
        if i in {0, 1, 3, 4}:
            my_height = v + 25
        plt.text(i + bar_width, my_height, str(round(v, 2)) + ' %', ha='center', va='center', rotation=90,
                 fontdict=dict(size=12))
    plt.ylim(0, 1.3 * ymax)

    # Add legend with appropriate layout
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2)

    # Show the plot
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    my_pdb_folder = filedialog.askdirectory()
    percent_diff(my_pdb_folder + '/', '1BNA', 'vol', False)
