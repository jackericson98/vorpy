import os
import sys
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.system.system import System
from vorpy.src.group.group import Group
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.plot_templates.scatter import scatter

from vorpy.src.analyze.tools.batch.get_files import get_all_files



def plot_vols():
    """Plots the average percentage differences for the given systems"""
    folder = tk.Tk()
    folder.withdraw()
    folder = filedialog.askdirectory()
    # get the aw, pow, and prm logs
    aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
    pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
    prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)
    # Get the title
    title = folder.split('/')[-1][2:] + " Volume Comparison"
    color_dict = {'C': 'grey', 'O': 'r', 'N': 'b', 'P': 'darkorange', 'H': 'pink', 'S': 'y', 'Se': 'sandybrown'}
    # Create the lists
    aw_vols, pow_vols, prm_vols, colors, labels = [], [], [], [], []
    # Loop through the atoms and get the volume differences
    for i, atom in aw_logs['atoms'].iterrows():
        # Get the power atom
        pow_atom = pow_logs['atoms'].loc[pow_logs['atoms']['Index'] == atom['Index']].to_dict(orient='records')[0]
        # Get the primitive atom
        prm_atom = prm_logs['atoms'].loc[prm_logs['atoms']['Index'] == atom['Index']].to_dict(orient='records')[0]
        # Add the volumes to the lists
        aw_vols.append(atom['Volume'])
        pow_vols.append(pow_atom['Volume'])
        prm_vols.append(prm_atom['Volume'])
        if atom['Name'] in color_dict:
            element = atom['Name']
        elif atom['Name'][:2].lower() == 'se':
            element = 'Se'
        else:
            element = atom['Name'][0]
        labels.append(element)
        colors.append(color_dict[element])
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot([5, 20], [5, 20], color='black', linestyle='--', linewidth=3, alpha=0.7)

    # Create legend handles that show the color and marker for each unique element
    import matplotlib.lines as mlines
    unique_elements = []
    unique_colors = []
    for l, c in zip(labels, colors):
        if l not in unique_elements:
            unique_elements.append(l)
            unique_colors.append(c)
    # All markers are 'o', so use Line2D with marker='o' for legend
    legend_handles = [
        mlines.Line2D([], [], color=c, marker='o', linestyle='None', markersize=8, label=l)
        for l, c in zip(unique_elements, unique_colors)
    ]
    legend_labels = unique_elements

    # Plot the data
    scatter(
        xs=[aw_vols],
        ys=[pow_vols], 
        title=title, 
        Show=True, 
        colors=[colors], 
        x_axis_title='AW Volume',
        y_axis_title='Pow Volume', 
        x_range=[3, 22],
        y_range=[3, 22],
        ax=ax, 
        fig=fig, 
        legend_title='Element', 
        legend_labels=legend_labels,
        legend_handles=legend_handles,
        alpha=0.5, 
        marker_size=100,
        x_tick_labels=[5, 10, 15, 20], 
        y_tick_labels=[5, 10, 15, 20],
        y_tick_label_locs=[5, 10, 15, 20],
        x_tick_label_locs=[5, 10, 15, 20],
        xtick_label_size=25, 
        ytick_label_size=25,
        xlabel_size=25,
        ylabel_size=25,
        title_size=25,
        legend_entry_size=20,
        legend_title_size=20, 
        legend_bbox_to_anchor=(1.5, 0.97),
        axis_line_thickness=2, 
        tick_width=3
    )

if __name__ == "__main__":
    plot_vols()
    