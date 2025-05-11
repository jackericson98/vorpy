import matplotlib.pyplot as plt
from Data.Analyze.tools.compare.compare_files import compare_files
import tkinter as tk
from tkinter import filedialog
import os


if __name__ == '__main__':
    # Choose between max and avg
    value = 'avg'
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
                            log_files=log_files, totals=True, curv=True)

    # Sample data
    labels_dict = {'a': 'Atoms', 'ad_mw': 'Avg Dist (mw)', 'ad': 'Avg Dist', 'ncap': 'Encapsulate',
                   'scbb_ad': 'SC/BB AD', 'scbb_ncap': 'SC/BB Encap.', 'scbb_ad_mw': 'SC/BB AD (mw)',
                   'martini': 'Martini', }
    labels = []
    for file in pdb_files[::2]:
        labels.append(labels_dict[file[len(logs_pdb_folder) + len(my_model_name) + 2:-4]])

    data = [round(my_info['totals'][_]['{} curv'.format(value)], 3) for _ in my_info['totals']]  # Sample data for the first set
    data1 = data[::2]
    data2 = data[1::2]
    max_height = max(data1)
    # Bar width
    bar_width = 0.35

    # Index for the x-axis
    x = range(len(labels))

    # Create the bar graph
    plt.bar(x, data1, width=bar_width, label='Additively Weighted')

    # Add labels and title
    if value == 'avg':
        plt.ylabel('Average Curvature (Gaussian)')
        plt.title('{} Average Curvature by Scheme'.format(my_model_name.capitalize()))
    else:
        plt.ylabel('Average Curvature (Gaussian)')
        plt.title('{} Max Curvature by Scheme'.format(my_model_name.capitalize()))

    # Angle the labels and add values at the top of the bars
    plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right')
    for i, v in enumerate(data1):
        plt.text(i, max_height / 2, str(v), ha='center', va='center', rotation=90)

    # Show the plot
    plt.tight_layout()
    plt.show()
