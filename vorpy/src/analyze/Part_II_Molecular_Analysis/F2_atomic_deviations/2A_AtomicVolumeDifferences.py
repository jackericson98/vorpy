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


def plot_vols(by_element=False, by_curvature=False):
    """Plots the average percentage differences for the given systems"""
    folder = tk.Tk()
    folder.withdraw()
    folder = filedialog.askdirectory()
    # get the aw, pow, and prm logs
    try:
        aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)
    except FileNotFoundError:
        aw_logs = read_logs2(os.path.join(folder, 'aw/aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow/pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm/prm_logs.csv'), all_=False, balls=True)

    # Get the title
    title = folder.split('/')[-1][2:] + " Volume Comparison"
    color_dict = {'C': 'grey', 'O': 'r', 'N': 'b', 'P': 'darkorange', 'H': 'pink', 'S': 'y', 'Se': 'sandybrown'}
    # Create the lists
    aw_vols, pow_vols, prm_vols, colors, labels, curv_list = [], [], [], [], [], []
    # Loop through the atoms and get the volume differences
    for i, atom in aw_logs['atoms'].iterrows():
        # Get the power atom
        pow_atom = pow_logs['atoms'].loc[pow_logs['atoms']['Index'] == atom['Index']].to_dict(orient='records')[0]
        # Get the primitive atom
        prm_atom = prm_logs['atoms'].loc[prm_logs['atoms']['Index'] == atom['Index']].to_dict(orient='records')[0]
        if (atom['Volume'] < 3 or atom['Volume'] > 22 or
            pow_atom['Volume'] < 3 or pow_atom['Volume'] > 22 or
            prm_atom['Volume'] < 3 or prm_atom['Volume'] > 22):
            continue
        # Add the volumes to the lists
        aw_vols.append(atom['Volume'])
        pow_vols.append(pow_atom['Volume'])
        prm_vols.append(prm_atom['Volume'])
        curv_list.append(atom['Maximum Mean Curvature'])
        if atom['Name'] in color_dict:
            element = atom['Name']
        elif atom['Name'][:2].lower() == 'se':
            element = 'Se'
        else:
            element = atom['Name'][0]
        labels.append(element)
        colors.append(color_dict[element])
    fig, ax = plt.subplots(figsize=(8, 6))
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import matplotlib.lines as mlines

    # Create legend handles that show the color and marker for each unique element
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

    # Determine coloring scheme
    if by_element:
        # Color by element (as before)
        plot_colors = [colors]
        legend_handles_to_use = legend_handles
        legend_labels_to_use = legend_labels
        add_colorbar = False
    elif by_curvature:
        # Color by curvature - pass curvature values directly for color mapping
        plot_colors = [curv_list]  # Pass curvature values for color mapping
        legend_handles_to_use = None
        legend_labels_to_use = None
        add_colorbar = True
        # Store normalization and colormap for colorbar
        norm = mcolors.Normalize(vmin=min(curv_list), vmax=max(curv_list))
        cmap = cm.viridis
    else:
        # Default to element coloring
        plot_colors = [colors]
        legend_handles_to_use = legend_handles
        legend_labels_to_use = legend_labels
        add_colorbar = False

    ax.plot([5, 20], [5, 20], color='black', linestyle='--', linewidth=3, alpha=0.7)

    # Plot the data
    if by_curvature:
        # For curvature coloring, use matplotlib's scatter directly with color mapping
        scatter_plot = ax.scatter(aw_vols, pow_vols, c=curv_list, cmap=cmap, norm=norm, 
                                 alpha=0.5, s=30, edgecolors='none')
        
        # Set up the plot formatting
        ax.set_xlabel('AW Volume', fontsize=25)
        ax.set_ylabel('Pow Volume', fontsize=25)
        ax.set_title(title, fontsize=25)
        ax.set_xlim(3, 22)
        ax.set_ylim(3, 22)
        ax.set_xticks([5, 10, 15, 20])
        ax.set_yticks([5, 10, 15, 20])
        ax.tick_params(axis='both', which='major', labelsize=25, width=3, length=12)
        
        # Set axis line thickness to 2pt
        for spine in ax.spines.values():
            spine.set_linewidth(2)
        
        # Add colorbar
        # cbar = plt.colorbar(scatter_plot, ax=ax, pad=0.02)
        # cbar.set_label('Mean Curvature', fontsize=25)
        # cbar.ax.tick_params(labelsize=25, width=2, length=12)
        
        # # Increase number of ticks on colorbar
        # cbar.locator = plt.MaxNLocator(nbins=4)
        # cbar.update_ticks()
        
        # Make plot fit in frame
        plt.tight_layout()
        plt.show()
    else:
        # Use the original scatter function for element coloring
        scatter(
            xs=[aw_vols],
            ys=[pow_vols], 
            title=title, 
            Show=True, 
            colors=plot_colors, 
            x_axis_title='AW Volume',
            y_axis_title='Pow Volume', 
            x_range=[3, 22],
            y_range=[3, 22],
            ax=ax, 
            fig=fig, 
            # legend_title='Element',
            # legend_labels=legend_labels_to_use,
            # legend_handles=legend_handles_to_use,
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
            # legend_entry_size=20,
            # legend_title_size=20,
            # legend_bbox_to_anchor=(1.5, 0.97),
            axis_line_thickness=2, 
            tick_width=3
        )


if __name__ == "__main__":
    plot_vols(by_element=True)
    