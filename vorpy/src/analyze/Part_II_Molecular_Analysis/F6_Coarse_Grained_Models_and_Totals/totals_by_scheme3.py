import csv
import os.path
import tkinter as tk
from tkinter import filedialog
import numpy as np

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.compare_files import compare_files2
from vorpy.src.analyze.tools.plot_templates.bar import bar

def percent_diff(residue_file, file_name):
    # Start re-sorting the data
    vols, sas = {}, {}
    with open(residue_file, 'r') as res_file:
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

    # Get the percent difference for each residue and average
    vol_means, vol_std_errs = [], []
    sa_means, sa_std_errs = [], []
    for file in sas:
        # Set up the percent difference list to get the mean and std error from later
        my_perc_diff = []
        # We already have the information for the atomic resolution
        if file == file_name:
            try:
                for res_name in vols[file]['aminos']:
                    for i in range(len(vols[file]['aminos'][res_name])):
                        my_perc_diff.append((float(vols[file_name]['aminos'][res_name][i])))
            except KeyError:
                pass
            # try:
            #     for res_name in vols[file]['nucs']:
            #         for i in range(len(vols[file]['nucs'][res_name])):
            #             my_perc_diff.append(float(vols[file_name]['nucs'][res_name][i]))
            # except KeyError:
            #     pass
            vol_means.append(round(sum(my_perc_diff), 3))
            vol_std_errs.append(round(np.std(my_perc_diff) / np.sqrt(len(my_perc_diff)) * 100, 3))
            # Set up the percent difference list to get the mean and std error from later
            my_perc_diff = []
            try:
                for res_name in sas[file]['aminos']:
                    for i in range(len(sas[file]['aminos'][res_name])):
                        my_perc_diff.append((float(sas[file_name]['aminos'][res_name][i])))
            except KeyError:
                pass
            # try:
            #     for res_name in vols[file]['nucs']:
            #         for i in range(len(vols[file]['nucs'][res_name])):
            #             my_perc_diff.append(float(vols[file_name]['nucs'][res_name][i]))
            # except KeyError:
            #     pass
            sa_means.append(round(sum(my_perc_diff) * 100, 3))
            sa_std_errs.append(round(np.std(my_perc_diff) / np.sqrt(len(my_perc_diff)) * 100, 3))

            continue
        # Volume data

        # try:
        #     for res_name in vols[file]['nucs']:
        #         for i in range(len(vols[file]['nucs'][res_name])):
        #             my_perc_diff.append((float(vols[file]['nucs'][res_name][i]) - float(
        #                 vols[file_name]['nucs'][res_name][i])))
        # except KeyError:
        #     pass
        try:
            for res_name in vols[file]['aminos']:
                for i in range(len(vols[file]['aminos'][res_name])):
                    my_perc_diff.append((float(vols[file]['aminos'][res_name][i]) - float(
                        vols[file_name]['aminos'][res_name][i])))
        except KeyError:
            pass
        vol_means.append(round(sum(my_perc_diff), 3))
        vol_std_errs.append(round(np.std(my_perc_diff) / np.sqrt(len(my_perc_diff)) * 100, 3))

        # Surface Area data
        my_perc_diff = []
        # try:
        #     for res_name in sas[file]['nucs']:
        #         for i in range(len(sas[file]['nucs'][res_name])):
        #             try:
        #                 my_perc_diff.append((float(sas[file]['nucs'][res_name][i]) - float(
        #                     sas[file_name]['nucs'][res_name][i])))
        #             except ZeroDivisionError:
        #                 continue
        # except KeyError:
        #     pass
        try:
            for res_name in sas[file]['aminos']:
                for i in range(len(sas[file]['aminos'][res_name])):
                    my_perc_diff.append((float(sas[file]['aminos'][res_name][i]) - float(
                        sas[file_name]['aminos'][res_name][i])))
        except KeyError:
            pass
        sa_means.append(round(sum(my_perc_diff) * 100, 3))
        sa_std_errs.append(round(np.std(my_perc_diff) / np.sqrt(len(my_perc_diff)) * 100, 3))

    # Your existing code
    vol_data = (vol_means[::2], vol_means[1::2], vol_std_errs[::2], vol_std_errs[1::2])
    sa_data = (sa_means[::2], sa_means[1::2], sa_std_errs[::2], sa_std_errs[1::2])

    return vol_data, sa_data


if __name__ == '__main__':
    # Go to the logs and pdbs folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    logs_pdb_folder = filedialog.askdirectory()
    # Get the model name
    my_model_name = ''
    log_files = []
    for file in os.listdir(logs_pdb_folder):
        filename = os.fsdecode(file)
        log_files.append(os.path.join(logs_pdb_folder, filename))
    # Get the data from the log files and sort it
    my_info = compare_files2(log_files)
    # Set up the dictionary for sorting

    order = ['Atomic', 'Encapsulation', 'Encapsulation Split', 'Average Distance', 'Average Distance Split', 'Average Distance Mass Weighted', 'Average Distance Mass Weighted Split']
    # Set up the 8 different tyupes in a dictionary:
    my_types = {'Atomic': {}, 'Encapsulation': {}, 'Encapsulation Split': {}, 'Average Distance': {}, 
                'Average Distance Split': {}, 'Average Distance Mass Weighted': {}, 'Average Distance Mass Weighted Split': {}}

    for file in my_info:
        # If the file has atoms in it, add it to the Atomic dictionary
        if 'atom' in file:
            my_types['Atomic'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        # If the file has encap and split in it, add it to the Encapsulation dictionary
        elif 'encap' in file and 'sr' in file:
            my_types['Encapsulation Split'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        # If the file has encap in it, add it to the Encapsulation Split dictionary
        elif 'encap' in file:
            my_types['Encapsulation'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        # If the file has ad, mw, and sr in it, add it to the Average Distance Mass Weighted Split dictionary
        elif 'ad' in file and 'mw' in file and 'sr' in file:
            my_types['Average Distance Mass Weighted Split'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        # If the file has ad, mw, and sr in it, add it to the Average Distance Mass Weighted dictionary
        elif 'ad' in file and 'mw' in file:
            my_types['Average Distance Mass Weighted'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        # If the file has ad and sr in it, add it to the Average Distance Split dictionary
        elif 'ad' in file and 'sr' in file:
            my_types['Average Distance Split'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        # If the file has ad in it, add it to the Average Distance dictionary
        elif 'ad' in file:
            my_types['Average Distance'][file[-1][:-4]] = {'vol': my_info[file]['group data']['Volume'], 'sa': my_info[file]['group data']['Surface Area']}
        else:
            print(f'{file} is not a valid file name')
            continue

    volx1, volx2, sax1, sax2 = [], [], [], []
    # set up the standard values for the additively weighted atomic values
    aw_vol = my_types['Atomic']['aw']['vol']
    aw_sa = my_types['Atomic']['aw']['sa']
    for my_key in order:
        volx1.append(my_types[my_key]['aw']['vol'])
        volx2.append(my_types[my_key]['pow']['vol'])
        sax1.append(my_types[my_key]['aw']['sa'])
        sax2.append(my_types[my_key]['pow']['sa'])
    print(volx1)
    print(volx2)
    print(sax1)
    print(sax2)

    # Volume bar Plot
    bar([[0] + [round(100 * (_ - volx1[0])/volx1[0], 3) for _ in volx1[1:]], [round(100 * (_ - volx1[0])/volx1[0], 3) for _ in volx2]], None, ["1", "2", "3", "4", "5", "6", "7"], ['Additively Weighted', 'Power'], my_model_name + ' Residue Volume % Diff',
        'CG Schemes', '% Difference', Show=True, print_vals_on_bars=True, unit='%', xtick_label_size=30, xlabel_size=30,
        ylabel_size=30, ytick_label_size=30, tick_width=2, tick_length=12)
    # Surface Area Plot
    bar([[0] + [round(100 * (_ - sax1[0])/sax1[0], 3) for _ in sax1[1:]], [round(100 * (_ - sax1[0])/sax1[0], 3) for _ in sax2]], None, ["1", "2", "3", "4", "5", "6", "7"], ['Additively Weighted', 'Power'], my_model_name + ' Residue SA % Diff',
        'CG Schemes', '% Difference', Show=True, print_vals_on_bars=True, unit='%', xtick_label_size=30, xlabel_size=30,
        ylabel_size=30, ytick_label_size=30, tick_width=2, tick_length=12)
