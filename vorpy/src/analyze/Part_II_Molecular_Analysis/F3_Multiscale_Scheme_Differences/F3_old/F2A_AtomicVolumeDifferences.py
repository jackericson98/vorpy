import os
import sys
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.system.system import System
from vorpy.src.group.group import Group
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.plot_templates.scatter import scatter

from vorpy.src.analyze.tools.batch.get_files import get_all_files


def add_outside_element_list(fig, unique_elements, unique_colors):
    x_marker = 0.80
    x_name = 0.835
    y_start = 0.90
    y_step = 0.035

    for i, (element, color) in enumerate(zip(unique_elements, unique_colors)):
        y = y_start - i * y_step

        if y < 0.06:
            break

        fig.text(
            x_marker,
            y,
            '●',
            color=color,
            fontsize=14,
            fontweight='bold',
            ha='left',
            va='center'
        )

        fig.text(
            x_name,
            y,
            str(element),
            color='black',
            fontsize=12,
            ha='left',
            va='center'
        )


def apply_2c_axis_style(
        fig,
        ax,
        title,
        volume_range,
        unique_elements=None,
        unique_colors=None,
        save_png=None,
        save_svg=None,
        show=True
):
    vmin, vmax = volume_range

    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ticks = [5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_xticklabels([f"{t:.1f}" for t in ticks])
    ax.set_yticklabels([f"{t:.1f}" for t in ticks])

    ax.set_xlabel('AW Volume', fontsize=24)
    ax.set_ylabel('Pow Volume', fontsize=24)
    ax.set_title(title, fontsize=22)

    ax.tick_params(axis='both', which='major', labelsize=20, width=2.5, length=10)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    ax.set_aspect('equal', adjustable='box')

    fig.subplots_adjust(left=0.12, right=0.77, bottom=0.12, top=0.90)

    if unique_elements is not None and unique_colors is not None:
        add_outside_element_list(fig, unique_elements, unique_colors)

    if save_png is not None:
        plt.savefig(save_png, dpi=300, bbox_inches='tight')

    if save_svg is not None:
        plt.savefig(save_svg, bbox_inches='tight')

    if show:
        plt.show()

    plt.close(fig)


def plot_vols(
        folders=None,
        by_element=False,
        by_curvature=False,
        alpha=0.2,
        molecule_class='protein',
        volume_range=(3, 22),
        downsample_fraction=0.1,
        random_seed=42,
        show_legend=True,
        show_colorbar=True,
        point_size=10,
        save_png=None,
        save_svg=None,
        show=True
):
    """Plot combined AW vs Pow volumes from all selected folders on one figure.

    Filtering is intentionally matched to the 2C clustered plots:
    - Same default volume range: 3 <= volume <= 22 for AW, Pow, and PRM.
    - Same protein point downsampling: 10% of valid protein atoms, random_state=42.
    """
    if folders is None:
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

    fig, ax = plt.subplots(figsize=(12, 9))

    color_dict = {
        'H': '#1f77b4',  # blue
        'C': '#ff7f0e',  # orange
        'N': '#2ca02c',  # green
        'O': '#d62728',  # red
        'P': '#9467bd',  # purple (new, consistent extension)
        'S': '#8c564b',  # brown
        'Se': '#e377c2'  # pink-purple
    }

    # Master combined rows across all folders.
    # Keep these as rows first so any later downsampling removes the same atoms
    # from every plotted attribute.
    records = []

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

            if volume_range is not None:
                vmin, vmax = volume_range

                if (
                        atom['Volume'] < vmin or atom['Volume'] > vmax or
                        pow_atom['Volume'] < vmin or pow_atom['Volume'] > vmax or
                        prm_atom['Volume'] < vmin or prm_atom['Volume'] > vmax
                ):
                    continue

            if atom['Name'] in color_dict:
                element = atom['Name']
            elif atom['Name'][:2].lower() == 'se':
                element = 'Se'
            else:
                element = atom['Name'][0]

            records.append({
                'AW': float(atom['Volume']),
                'Pow': float(pow_atom['Volume']),
                'Prm': float(prm_atom['Volume']),
                'Curvature': float(atom['Maximum Mean Curvature']),
                'Element': element,
                'Color': color_dict[element],
            })

    if not records:
        print("No valid data found to plot.")
        return

    original_count = len(records)

    rng = np.random.default_rng(random_seed)
    keep_count = max(1, int(round(original_count * downsample_fraction)))
    keep_indices = set(rng.choice(original_count, size=keep_count, replace=False).tolist())
    records = [row for i, row in enumerate(records) if i in keep_indices]

    print(
        f"Protein downsampling matched to 2C: kept {len(records)} of "
        f"{original_count} valid points ({downsample_fraction:.0%})."
    )

    aw_vols = [row['AW'] for row in records]
    pow_vols = [row['Pow'] for row in records]
    prm_vols = [row['Prm'] for row in records]
    curv_list = [row['Curvature'] for row in records]
    labels = [row['Element'] for row in records]
    colors = [row['Color'] for row in records]

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
        mlines.Line2D([], [], color=c, marker='o', linestyle='None', markersize=7, label=l)
        for l, c in zip(unique_elements, unique_colors)
    ]

    pretty_class = {
        'small_molecule': 'Small Molecule',
        'rna': 'RNA',
        'dna': 'DNA',
        'protein': 'Protein'
    }.get(molecule_class, molecule_class)

    title = (
        f"{pretty_class} Volume Comparison "
        f"(Downsample: {downsample_fraction:.0%})"
    )

    if volume_range is not None:
        vmin, vmax = volume_range
    else:
        vmin = min(min(aw_vols), min(pow_vols))
        vmax = max(max(aw_vols), max(pow_vols))

    ax.plot(
        [3, 22],
        [3, 22],
        linestyle='--',
        linewidth=3.5,
        color='black',
        alpha=0.9,
        zorder=0
    )

    if by_curvature:
        norm = mcolors.Normalize(vmin=min(curv_list), vmax=max(curv_list))
        cmap = cm.viridis

        scatter_plot = ax.scatter(
            aw_vols,
            pow_vols,
            c=curv_list,
            cmap=cmap,
            norm=norm,
            alpha=alpha,
            s=point_size,
            edgecolors='none'
        )

        if show_colorbar:
            cbar = plt.colorbar(scatter_plot, ax=ax, pad=0.02)
            cbar.set_label('Maximum Mean Curvature', fontsize=20)
            cbar.ax.tick_params(labelsize=18, width=2, length=8)

        if show_legend:
            add_outside_element_list(fig, unique_elements, unique_colors)

        apply_2c_axis_style(
            fig=fig,
            ax=ax,
            title=title,
            volume_range=(vmin, vmax),
            unique_elements=None,
            unique_colors=None,
            save_png=save_png,
            save_svg=save_svg,
            show=show
        )

    else:
        ax.scatter(
            aw_vols,
            pow_vols,
            c=colors,
            alpha=alpha,
            s=point_size,
            edgecolors='none',
            zorder=1
        )

        apply_2c_axis_style(
            fig=fig,
            ax=ax,
            title=title,
            volume_range=(vmin, vmax),
            unique_elements=unique_elements if show_legend else None,
            unique_colors=unique_colors if show_legend else None,
            save_png=save_png,
            save_svg=save_svg,
            show=show
        )


if __name__ == "__main__":
    bif = 'E:/Molecular'
    bof = (
        'E:/OneDrive - Georgia State University/GSU NSC/Manuscripts'
        '/Ericson Voronoi DNA/P2/fig2_atomic_level_scheme_deviations/2A_Full_Plots'
    )

    SMOL_SETTS = {
        'class': 'small_molecule',
        'folders': [os.path.join(bif, _) for _ in ['B_EDTA', 'C_DB1976', 'I_T4LP/JZ4']],
        'output': os.path.join(bof, 'small molecule'),
        'downsample_fraction': 1.0,
        'point_alpha': 0.75,
        'point_size': 18,
    }

    RNA_SETTS = {
        'class': 'rna',
        'folders': [os.path.join(bif, _) for _ in ['G_Hammerhead']],
        'output': os.path.join(bof, 'rna'),
        'downsample_fraction': 0.75,
        'point_alpha': 0.35,
        'point_size': 14,
    }

    DNA_SETTS = {
        'class': 'dna',
        'folders': [os.path.join(bif, _) for _ in ['D_Hairpin', 'F_BDNA', 'K_NCP_DNA']],
        'output': os.path.join(bof, 'dna'),
        'downsample_fraction': 0.25,
        'point_alpha': 0.25,
        'point_size': 12,
    }

    PROT_SETTS = {
        'class': 'protein',
        'folders': [os.path.join(bif, _) for _ in ['E_Cambrin', 'H_p53tet', 'I_T4LP', 'J_Streptavidin', 'L_BSA', 'm_NCP_Protein']],
        'output': os.path.join(bof, 'protein'),
        'downsample_fraction': 0.10,
        'point_alpha': 0.20,
        'point_size': 10,
    }

    settings = {
        'smol': SMOL_SETTS,
        'rna': RNA_SETTS,
        'dna': DNA_SETTS,
        'prot': PROT_SETTS,
    }

    current = 'prot'
    cfg = settings[current]

    os.makedirs(cfg['output'], exist_ok=True)

    plot_vols(
        folders=cfg['folders'],
        by_element=False,
        by_curvature=True,
        alpha=cfg['point_alpha'],
        point_size=cfg['point_size'],
        molecule_class=cfg['class'],
        downsample_fraction=cfg['downsample_fraction'],
        save_png=os.path.join(cfg['output'], f"{current}_element_volume_plot.png"),
        save_svg=os.path.join(cfg['output'], f"{current}_element_volume_plot.svg"),
        show=True
    )
