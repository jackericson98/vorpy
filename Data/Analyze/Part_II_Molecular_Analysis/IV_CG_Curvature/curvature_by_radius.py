import tkinter as tk
import os
from tkinter import filedialog

from Data.Analyze.tools.compare.read_logs import read_logs
from System.system import System
from Data.Analyze.tools.plot_templates.scatter import scatter
from Data.Analyze.tools.calcs import sort_lists


# Get the logs and pdbs folder

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
            # Get the name of the model
            my_model_name = filename[:-6]
        elif filename.endswith('.csv') and '_a_logs' in filename and '_pow' not in filename:
            log_files.append(os.path.join(logs_pdb_folder, filename))

    # Go through the files in the folder sorting them
    for file in os.listdir(logs_pdb_folder):
        filename = os.fsdecode(file)
        if '_a.pdb' in filename or '_a_logs' in filename:
            continue
        if filename.endswith('.pdb') and 'martini' not in filename:
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
        elif filename.endswith('.csv') and 'logs' in filename and 'aw' in filename and 'martini' not in filename:
            log_files.append(os.path.join(logs_pdb_folder, filename))

    amino_dict = {
        "ala": "A",  # Alanine
        "arg": "R",  # Arginine
        "asn": "N",  # Asparagine
        "asp": "D",  # Aspartic acid
        "cys": "C",  # Cysteine
        "gln": "Q",  # Glutamine
        "glu": "E",  # Glutamic acid
        "gly": "G",  # Glycine
        "his": "H",  # Histidine
        "ile": "I",  # Isoleucine
        "leu": "L",  # Leucine
        "lys": "K",  # Lysine
        "met": "M",  # Methionine
        "phe": "F",  # Phenylalanine
        "pro": "P",  # Proline
        "ser": "S",  # Serine
        "thr": "T",  # Threonine
        "trp": "W",  # Tryptophan
        "tyr": "Y",  # Tyrosine
        "val": "V"   # Valine
    }
    nucs_dict = {
        'dc': 'C',   # Cytosine
        'dg': 'G',   # Guanine
        'dt': 'T',   # Tyrosine
        'da': 'A',   # Adenine
        'du': 'U',   # Uracil
    }
    # Now that we have a list of logs and pdbs for each style, we need to read the logs and the pdbs
    my_logs_list, my_systems, rad_datax, rad_datay, rad_datares = [], [], [], [], []
    max_x, max_y = 0, 0
    for i, file in enumerate(pdb_files):
        my_system = System(file=file)
        my_logs = read_logs(log_files[i], no_sol=True)
        my_systems.append(my_system)
        my_logs_list.append(my_logs)
        # Get the data from curvature and sidechain size
        rad_dic = {}
        for j, atom in my_system.atoms.iterrows():
            rad_dic[atom['num']] = atom['rad']

        rad_datax.append([])
        rad_datay.append([])
        rad_datares.append([])

        for j, atom in my_logs['atoms'].iterrows():
            sys_atom = my_system.atoms.loc[my_system.atoms['num'] == atom['num']].to_dict(orient='records')[0]

            max_x, max_y = max(max_x, sys_atom['rad']), max(max_y, atom['max curv'])
            try:
                if sys_atom['res_name'].lower().strip() in {'sol', 'na'}:
                    continue
            except IndexError:
                continue
            except TypeError:
                print(sys_atom)
                continue
            # If the system is a sidechain backbone system assign the appropriate modifier
            if 'scbb' in my_system.name:
                # Proteins
                if sys_atom['res_name'].lower().strip() in amino_dict:
                    try:
                        rad_datax[-1].append(rad_dic[atom['num']])
                        rad_datay[-1].append(atom['max curv'])
                    except KeyError:
                        print("CUM")
                        continue

                    if sys_atom['element'].lower().strip() == 'pb':
                        rad_datares[-1].append('$' + amino_dict[sys_atom['res_name'].lower().strip()] + 'S$')
                    else:
                        rad_datares[-1].append('$' + amino_dict[sys_atom['res_name'].lower().strip()] + 'B$')

                # Nucleics
                else:
                    continue
                    try:
                        rad_datax[-1].append(rad_dic[atom['num']])
                        rad_datay[-1].append(atom['max curv'])
                    except KeyError:
                        continue
                    if sys_atom['element'].strip().lower() == 'pb':
                        rad_datares[-1].append('$' + nucs_dict[sys_atom['res_name'].lower().strip()] + 'S$')
                    elif sys_atom['element'].strip().lower() == 'bi':
                        rad_datares[-1].append('$' + nucs_dict[sys_atom['res_name'].lower().strip()] + 'N$')
                    else:
                        rad_datares[-1].append('$' + nucs_dict[sys_atom['res_name'].lower().strip()] + 'P$')
            else:
                try:
                    rad_datax[-1].append(rad_dic[atom['num']])
                    rad_datay[-1].append(atom['max curv'])
                except KeyError:
                    continue
                if sys_atom['res_name'].lower().strip() in amino_dict:
                    rad_datares[-1].append('$' + amino_dict[sys_atom['res_name'].lower().strip()] + '$')
                else:
                    continue
                    rad_datares[-1].append('$' + nucs_dict[sys_atom['res_name'].lower().strip()] + '$')


    # Sample data
    # Get the labels

    labels_dict = {'a': '1', 'ad_mw': '6', 'ad': '4', 'ncap': '2',
                   'scbb_ad': '5', 'scbb_ncap': '3', 'scbb_ad_mw': '7',
                   'martini': '8'}

    labels = []
    for file in pdb_files:
        labels.append(labels_dict[file[len(logs_pdb_folder) + len(my_model_name) + 2:-4]])


    labels, rad_datax, rad_datay, rad_datares = sort_lists(labels, rad_datax, rad_datay, rad_datares)

    # # Bin the data by grid
    # xs, ys, alphas, marker_sizes, markers = [], [], [], [], []
    # x_num_bins, y_num_bins
    # for i in range(len(rad_datax)):


    print(len(labels), labels)
    print(len(rad_datax), rad_datax)
    print(len(rad_datay), rad_datay)
    scatter(rad_datax, rad_datay, Show=True, y_range=[0, 1], x_axis_title='Ball Radius', y_axis_title='Curvature',
            title='{} Curvature Map'.format(my_model_name), alpha=0.1, markers=rad_datares)


