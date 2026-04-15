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


def plot_vols(by_element=False, by_curvature=False, alpha=0.2):
    """Plot combined AW vs Pow volumes from all selected folders on one figure."""
    folders = []

    while True:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory(title="Pick A Folder")
        root.destroy()
        print(folder)

        if folder == '' or folder is None:
            break

        folders.append(folder)

    if not folders:
        print("No folders selected.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    color_dict = {
        'C': 'grey',
        'O': 'r',
        'N': 'b',
        'P': 'darkorange',
        'H': 'pink',
        'S': 'y',
        'Se': 'sandybrown'
    }

    # Master combined lists across all folders
    aw_vols = []
    pow_vols = []
    prm_vols = []
    colors = []
    labels = []
    curv_list = []

    for folder in folders:
        try:
            aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
            pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
            prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)
        except FileNotFoundError:
            aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True)
            pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True)
            prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True)

        # Optional: speed up matching by indexing once
        pow_atoms = pow_logs['atoms'].set_index('Index')
        prm_atoms = prm_logs['atoms'].set_index('Index')

        for _, atom in aw_logs['atoms'].iterrows():
            idx = atom['Index']

            if idx not in pow_atoms.index or idx not in prm_atoms.index:
                continue

            pow_atom = pow_atoms.loc[idx]
            prm_atom = prm_atoms.loc[idx]

            if (
                atom['Volume'] < 3 or atom['Volume'] > 22 or
                pow_atom['Volume'] < 3 or pow_atom['Volume'] > 22 or
                prm_atom['Volume'] < 3 or prm_atom['Volume'] > 22
            ):
                continue

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

    if not aw_vols:
        print("No valid data found to plot.")
        return

    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import matplotlib.lines as mlines

    unique_elements = []
    unique_colors = []

    for label, color in zip(labels, colors):
        if label not in unique_elements:
            unique_elements.append(label)
            unique_colors.append(color)

    legend_handles = [
        mlines.Line2D([], [], color=c, marker='o', linestyle='None', markersize=8, label=l)
        for l, c in zip(unique_elements, unique_colors)
    ]

    title = "Combined Volume Comparison"

    ax.plot([5, 20], [5, 20], color='black', linestyle='--', linewidth=3, alpha=0.7)

    if by_curvature:
        norm = mcolors.Normalize(vmin=min(curv_list), vmax=max(curv_list))
        cmap = cm.viridis

        scatter_plot = ax.scatter(
            aw_vols,
            pow_vols,
            c=curv_list,
            cmap=cmap,
            norm=norm,
            alpha=0.5,
            s=30,
            edgecolors='none'
        )

        ax.set_xlabel('AW Volume', fontsize=25)
        ax.set_ylabel('Pow Volume', fontsize=25)
        ax.set_title(title, fontsize=25)
        ax.set_xlim(3, 22)
        ax.set_ylim(3, 22)
        ax.set_xticks([5, 10, 15, 20])
        ax.set_yticks([5, 10, 15, 20])
        ax.tick_params(axis='both', which='major', labelsize=25, width=3, length=12)

        for spine in ax.spines.values():
            spine.set_linewidth(2)

        # Uncomment if you want the colorbar back
        # cbar = plt.colorbar(scatter_plot, ax=ax, pad=0.02)
        # cbar.set_label('Mean Curvature', fontsize=25)
        # cbar.ax.tick_params(labelsize=25, width=2, length=12)

        plt.tight_layout()

    else:
        scatter(
            xs=[aw_vols],
            ys=[pow_vols],
            title=title,
            Show=False,
            colors=[colors],
            x_axis_title='AW Volume',
            y_axis_title='Pow Volume',
            x_range=[3, 22],
            y_range=[3, 22],
            ax=ax,
            fig=fig,
            alpha=alpha,
            marker_size=50,
            x_tick_labels=[5, 10, 15, 20],
            y_tick_labels=[5, 10, 15, 20],
            y_tick_label_locs=[5, 10, 15, 20],
            x_tick_label_locs=[5, 10, 15, 20],
            xtick_label_size=25,
            ytick_label_size=25,
            xlabel_size=25,
            ylabel_size=25,
            title_size=25,
            axis_line_thickness=2,
            tick_width=3
        )

    plt.show()


if __name__ == "__main__":
    plot_vols(by_element=True, alpha=0.8)
    