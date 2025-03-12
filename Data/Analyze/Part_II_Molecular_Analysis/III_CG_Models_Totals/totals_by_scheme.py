import os
import matplotlib.pyplot as plt
from Data.Analyze.tools.compare.compare_files import compare_files
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.plot_templates.bar import bar

"""
Plotting the totals for different schemes

Conventions:
1. All logs and pdbs must be in the same folder.
2. One pdb per pair of logs
3. Type Names - atom, ad, ncap, scbb_ad, scbb_ncap, martini
4. Pdbs = model_name + '_' + type + '.pdb'
5. Logs = model_name + '_' + type + scheme (aw or pow) + '_logs.csv'

Choose a metric below (sa or vol):
"""

metric = 'vol'

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
    my_info = compare_files(pdb_files=pdb_files,
                            log_files=log_files, totals=True)

    # Sample data
    # Get the labels

    labels_dict = {'a': '1', 'ad_mw': '6', 'ad': '4', 'ncap': '2',
                   'scbb_ad': '5', 'scbb_ncap': '3', 'scbb_ad_mw': '7',
                   'martini': '8'}
    labels = []
    for file in pdb_files[::2]:
        labels.append(labels_dict[file[len(logs_pdb_folder) + len(my_model_name) + 2:-4]])

    data = [round(my_info['totals'][_][metric], 2) for _ in my_info['totals']]  # Sample data for the first set

    data1 = data[::2]
    data2 = data[1::2]
    ymax = max(data)
    # Bar width
    bar_width = 0.35

    # Index for the x-axis
    x = range(len(labels))

    # Set the title
    # Add labels and title
    if metric == 'sa':
        unit = ' \u212B\u00B2'
        ylabel = 'Surface Area' + unit
        title = '{} Total Surface Area'.format(my_model_name.capitalize())
    elif metric == 'vol':
        unit = ' \u212B\u00B3'
        ylabel = 'Volume' + unit
        title = '{} Total Volume'.format(my_model_name.capitalize())

    # Calculate the % difference
    diff_data1 = [round(100 * (data1[i] - data1[0]) / data1[0], 3) for i in range(len(data1))]
    diff_data2 = [round(100 * (data2[i] - data2[0]) / data2[0], 3) for i in range(len(data2))]

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


    labels1, diff_data1, diff_data2 = sort_3_lists(labels, diff_data1, diff_data2)

    # Plot the difference data
    bar([diff_data1, diff_data2], x_names=labels1, legend_names=['Additively Weighted', 'Power'], x_axis_title='Scheme',
        y_axis_title='% Difference', Show=True, print_vals_on_bars=True, unit='%', xtick_label_size=30, xlabel_size=30,
        ylabel_size=30, ytick_label_size=30, tick_width=2, tick_length=12, bar_width=1)

    # Plot the bar graph
    bar([data1, data2], x_names=labels1, legend_names=['Additively Weighted', 'Power'],
        x_axis_title='Scheme', y_axis_title=ylabel, Show=True, print_vals_on_bars=True, unit=unit, y_range=[0, None],
        xtick_label_size=30, xlabel_size=30, ylabel_size=30, ytick_label_size=30, tick_width=2, tick_length=12)
