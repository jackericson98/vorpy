import tkinter as tk
import os
from tkinter import filedialog

from Data.Analyze.tools.compare.read_logs import read_logs
from System.system import System
from Data.Analyze.tools.plot_templates.scatter import scatter


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
        elif filename.endswith('.csv') and '_a_logs' in filename:
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

    # Now that we have a list of logs and pdbs for each style, we need to read the logs and the pdbs
    my_logs_list, my_systems, rad_datax, rad_datay = [], [], [], []
    for i, file in enumerate(pdb_files):
        my_system = System(file=file)
        my_logs = read_logs(log_files[i])
        my_systems.append(my_system)
        my_logs_list.append(my_logs)

        # Get the data from curvature and sidechain size
        rad_dic = {}
        for j, atom in my_system.atoms.iterrows():
            if atom['res'].name == 'SOL':
                continue
            rad_dic[atom['num']] = atom['rad']

        rad_datax.append([])
        rad_datay.append([])

        for j, atom in my_logs['atoms'].iterrows():
            try:
                rad_datax[-1].append(rad_dic[atom['num']])
                rad_datay[-1].append(atom['max curv'])
            except KeyError:
                continue

    # Sample data
    # Get the labels

    labels_dict = {'a': '1', 'ad_mw': '6', 'ad': '4', 'ncap': '2',
                   'scbb_ad': '5', 'scbb_ncap': '3', 'scbb_ad_mw': '7',
                   'martini': '8'}
    labels = []
    for file in pdb_files:
        labels.append(labels_dict[file[len(logs_pdb_folder) + len(my_model_name) + 2:-4]])


    def sort_3_lists(lista, listb, listc):
        # Zipping lists together and sorting by the first list
        sorted_lists = sorted(zip(lista, listb, listc), key=lambda x: x[0])

        # Unpacking the sorted lists
        lista, listb, listc = zip(*sorted_lists)

        # Converting tuples back to lists if needed
        lista = list(lista)
        listb = list(listb)
        listc = list(listc)

        # Return the lists
        return lista, listb, listc


    labels, rad_datax, rad_datay = sort_3_lists(labels, rad_datax, rad_datay)

    scatter(rad_datax, rad_datay, labels=labels, Show=True, y_range=[0, 1], legend_orientation='vertical',
            x_axis_title='Ball Radius', y_axis_title='Curvature', title='{} Curvature Map'.format(my_model_name),
            alpha=0.5, legend_title='Scheme', legend_title_size=15)


