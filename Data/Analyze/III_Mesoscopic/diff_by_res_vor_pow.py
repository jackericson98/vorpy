import csv
import os.path
import tkinter as tk
from tkinter import filedialog
import numpy as np
from Data.Analyze.tools.compare.compare_files import compare_files
from Data.Analyze.tools.plot_templates.bar import bar


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
    vol_means, vol_std_errs = [0], [0]
    sa_means, sa_std_errs = [0], [0]
    for file in sas:
        # We already have the information for the atomic resolution
        if file == file_name:
            continue
        # Volume data
        # Set up the percent difference list to get the mean and std error from later
        my_perc_diff = []
        try:
            for res_name in vols[file]['nucs']:
                for i in range(len(vols[file]['nucs'][res_name])):
                    my_perc_diff.append(abs(float(vols[file_name]['nucs'][res_name][i]) - float(
                        vols[file]['nucs'][res_name][i])) / float(vols[file_name]['nucs'][res_name][i]))
        except KeyError:
            pass
        try:
            for res_name in vols[file]['aminos']:
                for i in range(len(vols[file]['aminos'][res_name])):
                    my_perc_diff.append(abs(float(vols[file_name]['aminos'][res_name][i]) - float(
                        vols[file]['aminos'][res_name][i])) / float(vols[file_name]['aminos'][res_name][i]))
        except KeyError:
            pass
        vol_means.append(round(np.mean(my_perc_diff) * 100, 3))
        vol_std_errs.append(round(np.std(my_perc_diff) / np.sqrt(len(my_perc_diff)) * 100, 3))

        # Surface Area data
        my_perc_diff = []
        try:
            for res_name in sas[file]['nucs']:
                for i in range(len(sas[file]['nucs'][res_name])):
                    try:
                        my_perc_diff.append(abs(float(sas[file_name]['nucs'][res_name][i]) - float(
                            sas[file]['nucs'][res_name][i])) / float(sas[file_name]['nucs'][res_name][i]))
                    except ZeroDivisionError:
                        continue
        except KeyError:
            pass
        try:
            for res_name in sas[file]['aminos']:
                for i in range(len(sas[file]['aminos'][res_name])):
                    my_perc_diff.append(abs(float(sas[file_name]['aminos'][res_name][i]) - float(
                        sas[file]['aminos'][res_name][i])) / float(sas[file_name]['aminos'][res_name][i]))
        except KeyError:
            pass
        sa_means.append(round(np.mean(my_perc_diff) * 100, 3))
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
    log_files, pdb_files = [], []
    for file in os.listdir(logs_pdb_folder):
        filename = os.fsdecode(file)
        if filename.endswith('a.pdb'):
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
            # Get the name of the model
            my_model_name = filename[:-6]
        elif filename.endswith('.csv') and '_a_logs' in filename:
            log_files.append(os.path.join(logs_pdb_folder, filename))

    # Go through the files in the folder sorting them
    for file in os.listdir(logs_pdb_folder):
        filename = os.fsdecode(file)
        if '_a.pdb' in filename or '_a_logs' in filename:
            continue
        if filename.endswith('.pdb'):
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
        elif filename.endswith('.csv') and 'logs' in filename:
            log_files.append(os.path.join(logs_pdb_folder, filename))
    # Check to see if the data has been processed yet
    if not os.path.exists(logs_pdb_folder + '/residue_data.csv'):
        # Get the data from the log files and sort it
        my_info = compare_files(pdb_files, log_files, avg_distros=True, by_residues=True)
        # Put the data into a csv file for later access
        with open(logs_pdb_folder + '/residue_data.csv', 'w') as res_file:
            res_fl = csv.writer(res_file)
            res_fl.writerow(['file', 'residue type', 'name', 'residue', 'volume', 'surface area'])
            for file in my_info['residues']:
                for res_type in my_info['residues'][file]:
                    for res_name in my_info['residues'][file][res_type]:
                        for res in my_info['residues'][file][res_type][res_name]:
                            res_fl.writerow(
                                [file, res_type, res_name, res] + [my_info['residues'][file][res_type][res_name][res][_]
                                                                   for _ in
                                                                   my_info['residues'][file][res_type][res_name][res]])
    vol_dat, sa_dat  = percent_diff(logs_pdb_folder + '/residue_data.csv', file_name=my_model_name + '_a')
    # Sample data
    labels_dict = {'a': 'Atoms', 'ad_mw': 'Avg Dist (mw)', 'ad': 'Avg Dist', 'ncap': 'Encapsulate',
                   'scbb_ad': 'SC/BB AD', 'scbb_ncap': 'SC/BB Encap.', 'scbb_ad_mw': 'SC/BB AD (mw)',
                   'martini': 'Martini', }
    labels = []
    for file in pdb_files[::2]:
        labels.append(labels_dict[file[len(logs_pdb_folder) + len(my_model_name) + 2:-4]])
    # Volume bar Plot
    bar(vol_dat[:2], vol_dat[2:], labels, ['Additively Weighted', 'Power'], my_model_name + ' Residue Volume % Diff',
        'CG Schemes', '% Difference', Show=True, print_vals_on_bars=True)
    # Surface Area Plot
    bar(sa_dat[:2], sa_dat[2:], labels, ['Additively Weighted', 'Power'], my_model_name + ' Residue SA % Diff',
        'CG Schemes', '% Difference', Show=True, print_vals_on_bars=True)
