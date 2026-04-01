import os
import sys
import tkinter as tk
from tkinter import filedialog
from typing import List, Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
import pandas as pd


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


FINAL_ATOM_ALIASES = {

    # =========================
    # TERMINAL HYDROGENS
    # =========================
    'H1': 'HT',
    'H2': 'HT',
    'H3': 'HT',

    # =========================
    # ALPHA HYDROGENS
    # =========================
    'HA1': 'HA',
    'HA2': 'HA',

    # =========================
    # BETA HYDROGENS
    # =========================
    'HB1': 'HB',
    'HB2': 'HB',
    'HB3': 'HB',

    # =========================
    # GAMMA HYDROGENS
    # =========================
    'HG1': 'HG',
    'HG2': 'HG',
    'HG11': 'HG',
    'HG12': 'HG',
    'HG13': 'HG',
    'HG21': 'HG',
    'HG22': 'HG',
    'HG23': 'HG',

    # =========================
    # DELTA HYDROGENS
    # =========================
    'HD1': 'HD',
    'HD2': 'HD',
    'HD3': 'HD',
    'HD11': 'HD',
    'HD12': 'HD',
    'HD13': 'HD',
    'HD21': 'HD',
    'HD22': 'HD',
    'HD23': 'HD',

    # =========================
    # EPSILON HYDROGENS (KEEP HE2 SEPARATE)
    # =========================
    'HE1': 'HE',
    'HE3': 'HE',

    'HE21': 'HE2',
    'HE22': 'HE2',

    # =========================
    # ZETA HYDROGENS
    # =========================
    'HZ2': 'HZ',
    'HZ3': 'HZ',

    # =========================
    # ETA HYDROGENS
    # =========================
    'HH1': 'HH',
    'HH2': 'HH',
    'HH11': 'HH',
    'HH12': 'HH',
    'HH21': 'HH',
    'HH22': 'HH',

    # =========================
    # CARBON SYMMETRY
    # =========================
    'CG1': 'CG',
    'CG2': 'CG',

    'CD1': 'CD',
    'CD2': 'CD',

    'CE2': 'CE',
    'CE3': 'CE',

    'CZ2': 'CZ',
    'CZ3': 'CZ',

    # =========================
    # NITROGEN SYMMETRY
    # =========================
    'NH1': 'NH',
    'NH2': 'NH',

    # =========================
    # OXYGEN GROUPING (CLUSTER-INFORMED)
    # =========================

    # backbone / hydroxyl merged
    'O':  'O_backbone',
    'OH': 'O_backbone',

    # carboxyl split (position-sensitive)
    'OE1': 'O_carboxyl_1',
    'OD1': 'O_carboxyl_1',

    'OE2': 'O_carboxyl_2',
    'OD2': 'O_carboxyl_2',

    # terminal oxygens grouped with carboxyl_2
    'OC1': 'O_carboxyl_2',
    'OC2': 'O_carboxyl_2',

    'OT1': 'O_carboxyl_2',
    'OT2': 'O_carboxyl_2',

    # =========================
    # OPTIONAL CARBON MERGES (from clustering)
    # =========================
    'CH2': 'C_aromatic_outer',
    'CZ':  'C_aromatic_outer',
    'C6':  'C_aromatic_outer',

    # optional CB neighborhood merge
    'C2': 'CB',
    'C3': 'CB',
}

DIRECT_LABEL_GROUPS = {
    'O_backbone',
    'O_carboxyl_1',
    'O_carboxyl_2',
    'CB',
    'C_aromatic_outer',
}


def canonicalize_atom_name(atom_name: str, alias_dict=None) -> str:
    if alias_dict is None:
        alias_dict = FINAL_ATOM_ALIASES

    name = str(atom_name).strip().upper()
    return alias_dict.get(name, name)


def get_atom_class(name: str) -> str:
    name = name.upper()

    if name.startswith('H'):
        return 'H'

    if name.startswith('C'):
        if name == 'CA':
            return 'CA'
        if name == 'CB':
            return 'CB'
        return 'C'

    if name.startswith('N'):
        if name in ['N', 'HN']:
            return 'N_backbone'
        if name in ['NZ', 'NH']:
            return 'N_terminal'
        return 'N_side'

    if name.startswith('O'):
        return 'O'

    return 'Other'


def _read_scheme_logs(folder: str):
    """
    Reads aw/pow/prm logs from either:
      folder/{aw_logs.csv,pow_logs.csv,prm_logs.csv}
    or
      folder/{aw/aw_logs.csv, pow/pow_logs.csv, prm/prm_logs.csv}
    """
    try:
        aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)

    except FileNotFoundError:
        aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True)

    return aw_logs, pow_logs, prm_logs


def select_folders_multi(title: str = "Select a folder (Cancel to finish)"):
    """
    Tkinter doesn’t provide a native multi-directory picker, so this prompts repeatedly.
    Select folders one at a time; press Cancel when you’re done.
    """
    root = tk.Tk()
    root.withdraw()

    folders = []
    while True:
        folder = filedialog.askdirectory(title=title)
        if not folder:
            break
        folders.append(folder)

    root.destroy()
    return folders


def collect_atom_volume_points(
    folders: List[str],
    atom_name_field: str = 'Name',
    volume_range: Optional[tuple] = None
) -> pd.DataFrame:
    """
    Collect pooled AW vs Pow atom volumes across selected folders.
    Returns dataframe with:
      Folder, Index, AtomName, AW, Pow
    """
    records: List[Dict[str, object]] = []

    for folder in folders:
        aw_logs, pow_logs, _ = _read_scheme_logs(folder)

        aw_atoms = aw_logs['atoms']
        pow_atoms = pow_logs['atoms']

        pow_lookup = {
            int(row['Index']): row
            for _, row in pow_atoms.iterrows()
        }

        for _, atom in aw_atoms.iterrows():
            idx = int(atom['Index'])

            if idx not in pow_lookup:
                continue

            pow_atom = pow_lookup[idx]

            aw_v = float(atom['Volume'])
            pow_v = float(pow_atom['Volume'])

            if volume_range is not None:
                vmin, vmax = volume_range
                if aw_v < vmin or aw_v > vmax or pow_v < vmin or pow_v > vmax:
                    continue

            raw_name = atom.get(atom_name_field, '')
            atom_name = str(raw_name).strip().upper()

            if atom_name == '':
                continue

            canonical_name = canonicalize_atom_name(atom_name)

            records.append({
                'Folder': folder,
                'Index': idx,
                'AtomName': atom_name,
                'CanonicalName': canonical_name,
                'AW': aw_v,
                'Pow': pow_v
            })

    return pd.DataFrame(records)


def cluster_groups_constrained(stats_df, eps=0.5, min_samples=2):
    stats_df = stats_df.copy()
    stats_df['Class'] = stats_df['GroupName'].apply(get_atom_class)

    cluster_labels = []

    for cls in stats_df['Class'].unique():
        sub = stats_df[stats_df['Class'] == cls]

        coords = sub[['Mean_AW', 'Mean_Pow']].values

        if len(coords) < 2:
            labels = [-1] * len(sub)
        else:
            from sklearn.cluster import DBSCAN
            labels = DBSCAN(eps=eps, min_samples=min_samples).fit(coords).labels_

        cluster_labels.extend(zip(sub.index, labels))

    cluster_map = dict(cluster_labels)
    stats_df['Cluster'] = stats_df.index.map(cluster_map)

    return stats_df


def cluster_groups(stats_df, eps=0.35, min_samples=2):
    """
    Cluster canonical groups based on their (Mean_AW, Mean_Pow).

    Returns:
        stats_df with new column 'Cluster'
    """

    coords = stats_df[['Mean_AW', 'Mean_Pow']].values

    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)

    stats_df = stats_df.copy()
    stats_df['Cluster'] = clustering.labels_

    return stats_df


def compute_name_volume_stats(df: pd.DataFrame, group_col: str = 'CanonicalName') -> pd.DataFrame:
    """
    Compute center and spread in AW/Pow space for each plotted group.
    group_col should usually be 'CanonicalName'.
    """
    if len(df) == 0:
        return pd.DataFrame(
            columns=['GroupName', 'Count', 'Mean_AW', 'Mean_Pow', 'SD_AW', 'SD_Pow', 'Members']
        )

    grouped_rows = []

    for group_name, subdf in df.groupby(group_col):
        member_names = sorted(subdf['AtomName'].unique().tolist())

        grouped_rows.append({
            'GroupName': group_name,
            'Count': len(subdf),
            'Mean_AW': subdf['AW'].mean(),
            'Mean_Pow': subdf['Pow'].mean(),
            'SD_AW': subdf['AW'].std(ddof=1) if len(subdf) > 1 else 0.0,
            'SD_Pow': subdf['Pow'].std(ddof=1) if len(subdf) > 1 else 0.0,
            'Members': ', '.join(member_names)
        })

    stats_df = pd.DataFrame(grouped_rows)

    stats_df['SD_AW'] = stats_df['SD_AW'].fillna(0.0)
    stats_df['SD_Pow'] = stats_df['SD_Pow'].fillna(0.0)

    stats_df = stats_df.sort_values(
        ['Count', 'GroupName'],
        ascending=[False, True]
    ).reset_index(drop=True)

    return stats_df


def build_cluster_alias_dict(stats_df):
    """
    Suggest merges: groups in same cluster → same super group
    """
    cluster_alias = {}

    for cluster_id in stats_df['Cluster'].unique():
        if cluster_id == -1:
            continue  # noise

        sub = stats_df[stats_df['Cluster'] == cluster_id]

        names = sorted(sub['GroupName'].tolist())

        # choose representative name (first or centroid-like)
        rep = names[0]

        for name in names:
            cluster_alias[name] = rep

    return cluster_alias


def print_clusters(stats_df):
    print("\nCluster assignments:\n")

    for cluster_id in sorted(stats_df['Cluster'].unique()):
        sub = stats_df[stats_df['Cluster'] == cluster_id]

        print(f"Cluster {cluster_id}:")
        for _, row in sub.iterrows():
            print(f"  {row['GroupName']:>6s} | AW={row['Mean_AW']:.2f}, Pow={row['Mean_Pow']:.2f}")
        print()


def print_name_volume_stats(stats_df: pd.DataFrame):
    print("\nGrouped atom-name centers and standard deviations in AW vs Pow volume space:\n")

    for _, row in stats_df.iterrows():
        print(
            f"{row['GroupName']:>6s} | "
            f"n={int(row['Count']):>6d} | "
            f"AW={row['Mean_AW']:>8.3f} ± {row['SD_AW']:.3f} | "
            f"Pow={row['Mean_Pow']:>8.3f} ± {row['SD_Pow']:.3f} | "
            f"Members: {row['Members']}"
        )


def build_name_color_map(atom_names: List[str], cmap_name: str = 'tab20'):
    """
    Assign a color to each atom name.
    If there are more names than the colormap length, it cycles.
    """
    cmap = plt.get_cmap(cmap_name)
    n_colors = getattr(cmap, 'N', len(atom_names))

    unique_names = list(atom_names)
    color_map = {}

    for i, name in enumerate(unique_names):
        color_map[name] = cmap(i % n_colors)

    return color_map


def get_plot_color(name):
    if name.startswith('H'):
        return '#1f77b4'  # blue

    if name.startswith('C'):
        return '#ff7f0e'  # orange

    if name.startswith('N'):
        return '#2ca02c'  # green

    if name.startswith('O'):
        return '#d62728'  # red

    return 'gray'


def plot_name_volume_groups(
    atom_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    title: str,
    volume_range: Optional[tuple] = None,
    show_points: bool = True,
    show_centers: bool = True,
    show_errorbars: bool = False,
    annotate: bool = True,
    min_count_for_label: int = 2,
    label_top_n: Optional[int] = None,
    plot_min_count: int = 2,
    point_size: float = 18,
    point_alpha: float = 0.25,
    center_size: float = 220,
    label_fontsize: int = 11,
    save: Optional[str] = None,
    show: bool = True
):
    fig, ax = plt.subplots(figsize=(12, 9))

    # Only groups above threshold get plotted / centered / labeled / legend entries
    plot_stats_df = stats_df[stats_df['Count'] >= plot_min_count].copy()
    plot_stats_df = plot_stats_df.sort_values('Mean_AW').reset_index(drop=True)
    plot_stats_df['PlotNumber'] = np.arange(1, len(plot_stats_df) + 1)

    group_names = plot_stats_df['GroupName'].tolist()
    # color_map = build_name_color_map(group_names, cmap_name='tab20')

    for _, row in plot_stats_df.iterrows():
        name = row['GroupName']
        color = get_plot_color(name)

        group_df = atom_df[atom_df['CanonicalName'] == name]

        if show_points:
            ax.scatter(
                group_df['AW'],
                group_df['Pow'],
                s=point_size,
                alpha=point_alpha,
                color=color,
                label=name
            )

        if show_errorbars:
            ax.errorbar(
                row['Mean_AW'],
                row['Mean_Pow'],
                xerr=row['SD_AW'],
                yerr=row['SD_Pow'],
                fmt='none',
                ecolor=color,
                elinewidth=2.0,
                capsize=4,
                alpha=0.95,
                zorder=3
            )

        if show_centers:
            ax.text(
                row['Mean_AW'],
                row['Mean_Pow'],
                str(int(row['PlotNumber'])),
                color=color,
                fontsize=14,
                fontweight='bold',
                ha='center',
                va='center',
                zorder=5
            )

    if annotate:
        label_df = plot_stats_df.copy()

        if label_top_n is not None:
            label_df = label_df.head(label_top_n)

        for _, row in label_df.iterrows():
            if int(row['Count']) < min_count_for_label:
                continue

            if row['GroupName'] in DIRECT_LABEL_GROUPS:
                ax.text(
                    row['Mean_AW'] + 0.18,
                    row['Mean_Pow'] + 0.18,
                    row['GroupName'],
                    fontsize=11,
                    color='black',
                    ha='left',
                    va='bottom',
                    zorder=6
                )

    if volume_range is not None:
        vmin, vmax = volume_range
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)

        ax.plot(
            [vmin, vmax],
            [vmin, vmax],
            linestyle='--',
            linewidth=3.5,
            color='black',
            alpha=0.9
        )

    else:
        if len(atom_df) > 0:
            xmin = min(atom_df['AW'].min(), atom_df['Pow'].min())
            xmax = max(atom_df['AW'].max(), atom_df['Pow'].max())
        else:
            xmin, xmax = 0.0, 1.0

        pad = 0.05 * (xmax - xmin if xmax > xmin else 1.0)

        ax.set_xlim(xmin - pad, xmax + pad)
        ax.set_ylim(xmin - pad, xmax + pad)

        ax.plot(
            [xmin - pad, xmax + pad],
            [xmin - pad, xmax + pad],
            linestyle='--',
            linewidth=2.5,
            color='black',
            alpha=0.7
        )

    ax.set_xlabel('AW Volume', fontsize=24)
    ax.set_ylabel('Pow Volume', fontsize=24)
    ax.set_title(title, fontsize=22)

    ax.tick_params(axis='both', which='major', labelsize=20, width=2.5, length=10)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    unique_handles = []
    unique_labels = []

    for h, lab in zip(handles, labels):
        if lab not in seen:
            unique_handles.append(h)
            unique_labels.append(lab)
            seen.add(lab)

    if len(unique_labels) > 0:
        ax.legend(
            unique_handles,
            unique_labels,
            fontsize=10,
            frameon=True,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5)
        )

    ax.set_aspect('equal', adjustable='box')
    fig.subplots_adjust(left=0.12, right=0.78, bottom=0.12, top=0.90)

    if save is not None:
        plt.savefig(save, dpi=300)

    if show:
        plt.show()

    plt.close(fig)


def main(
    atom_name_field: str = 'Name',
    volume_range: tuple = (3, 22),
    save_csv: bool = True,
    save_plot: bool = False,
    show_points: bool = True,
    show_centers: bool = True,
    show_errorbars: bool = True,
    annotate: bool = True,
    min_count_for_label: int = 200,
    label_top_n: Optional[int] = None,
    plot_min_count: int = 2,
    use_auto_clustering: bool = True
):
    folders = select_folders_multi()
    if len(folders) == 0:
        print("No folders selected.")
        return

    atom_df = collect_atom_volume_points(
        folders=folders,
        atom_name_field=atom_name_field,
        volume_range=volume_range
    )

    if len(atom_df) == 0:
        print("No matching AW/Pow atom volume pairs were found.")
        return

    stats_df = compute_name_volume_stats(atom_df, group_col='CanonicalName')

    print_name_volume_stats(stats_df)

    if use_auto_clustering:
        stats_df = cluster_groups(stats_df, eps=0.35, min_samples=2)
        print_clusters(stats_df)

    out_dir = folders[0]

    if save_csv:
        csv_path = os.path.join(out_dir, 'atom_group_volume_stats.csv')
        stats_df.to_csv(csv_path, index=False)
        print(f"\nSaved CSV -> {csv_path}")

    plot_path = None
    if save_plot:
        plot_path = os.path.join(out_dir, 'atom_group_volume_stats.png')

    title = (
        f"Grouped atom-name groupings in AW vs Pow space "
        f"(n_folders={len(folders)}, n_atoms={len(atom_df)})"
    )

    plot_name_volume_groups(
        atom_df=atom_df,
        stats_df=stats_df,
        title=title,
        volume_range=volume_range,
        show_points=show_points,
        show_centers=show_centers,
        show_errorbars=show_errorbars,
        annotate=annotate,
        min_count_for_label=min_count_for_label,
        label_top_n=label_top_n,
        plot_min_count=plot_min_count,
        save=plot_path,
        show=True
    )

    if plot_path is not None:
        print(f"Saved plot -> {plot_path}")


if __name__ == "__main__":
    main(
        atom_name_field='Name',
        volume_range=(3, 22),
        save_csv=True,
        save_plot=False,
        show_points=True,
        show_centers=True,
        show_errorbars=False,
        annotate=True,
        min_count_for_label=2,
        label_top_n=None,
        plot_min_count=2
    )