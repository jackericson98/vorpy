import os
import matplotlib.pyplot as plt
from Data.Analyze.tools.compare.compare_files import compare_files
import tkinter as tk
from tkinter import filedialog

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

metric = 'sa'

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
        if filename.endswith('atom.pdb'):
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
            pdb_files.append(os.path.join(logs_pdb_folder, filename))
            # Get the name of the model
            my_model_name = filename[:-9]
        elif filename.endswith('.csv') and 'atom' in filename:
            log_files.append(os.path.join(logs_pdb_folder, filename))

    # Go through the files in the folder sorting them
    for file in os.listdir(logs_pdb_folder):
        filename = os.fsdecode(file)
        if 'atom' in filename:
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
    labels_dict = {'atom': 'Atoms', 'ad': 'Avg Dist', 'ncap': 'Encapsulate', 'scbb_ad': 'SC/BB AD', 'scbb_ncap': 'SC/BB Encap.', 'martini': 'Martini'}
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

    # Create the bar graph
    plt.bar(x, data1, width=bar_width, label='Additively Weighted', color='skyblue', edgecolor='black')
    plt.bar([i + bar_width for i in x], data2, width=bar_width, label='Power', color='orange', edgecolor='black')

    # Add labels and title
    if metric == 'sa':
        plt.ylabel('Surface Area', fontdict=dict(size=15))
        plt.title('{} Surface Area by Scheme'.format(my_model_name), fontdict=dict(size=20))
        unit = ' \u212B\u00B2'
    elif metric == 'vol':
        plt.ylabel('Volume', fontdict=dict(size=15))
        plt.title('{} Volume by Scheme'.format(my_model_name), fontdict=dict(size=20))
        unit = ' \u212B\u00B3'

    # Angle the labels and add values at the top of the bars
    plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right')
    for i, v in enumerate(data1):
        plt.text(i, v / 2, str(v) + unit, ha='center', va='center', rotation=90)
    for i, v in enumerate(data2):
        plt.text(i + bar_width, v / 2, str(v) + unit, ha='center', va='center', rotation=90)
    plt.ylim(0, 1.25 * ymax)
    # Add legend with appropriate layout
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2)

    # Show the plot
    plt.tight_layout()
    plt.show()
