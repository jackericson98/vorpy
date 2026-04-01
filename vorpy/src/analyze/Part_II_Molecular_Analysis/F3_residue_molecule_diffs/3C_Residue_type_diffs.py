import os
import sys
import math
import warnings
import tkinter as tk

from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np


warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    module=r"PIL|matplotlib\.backends\._backend_tk",
)

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


AA_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
]

NA_ORDER = ["DA", "DC", "DG", "DT", "DU"]

RES_ORDER = AA_ORDER + NA_ORDER


# -----------------------------
# USER SETTINGS
# -----------------------------
EXCLUDE_KEYS = ['A', 'B', 'C']

MAX_PERCENT_DIFF = None
# Example:
# MAX_PERCENT_DIFF = 150

SORT_BY = None
# Options:
# None           -> use RES_ORDER
# 'pow_mean'     -> sort by mean power deviation
# 'prm_mean'     -> sort by mean primitive deviation
# 'combined_mean' -> sort by average of power/prm means

SHOW_POINTS = True
SHOW_MEANS = True
SHOW_COUNTS = True

FIGSIZE = (24, 10)
DPI = 300

POWER_COLOR = '#d62728'
PRIMITIVE_COLOR = '#7f3fbf'

POWER_ALPHA = 0.45
PRIMITIVE_ALPHA = 0.45

POINT_ALPHA = 0.20
POINT_SIZE = 18

VIOLIN_WIDTH = 0.34
OFFSET = 0.20

Y_RANGE = None
# Example:
# Y_RANGE = [0, 80]

TITLE = 'Figure 3C — Residue-Type Volume Deviation Distributions Across All Models'
Y_LABEL = 'Absolute % Volume Deviation from AW'
X_LABEL = 'Residue Type'

SAVE_NAME = 'figure_3c_residue_type_deviation_distribution.png'


def normalize_residue_label(res: str) -> str:
    r = str(res).strip().upper()

    if r in {"A", "C", "G", "U", "T"}:
        if r == "A":
            return "DA"
        if r == "C":
            return "DC"
        if r == "G":
            return "DG"
        if r == "T":
            return "DT"
        if r == "U":
            return "DU"

    if r in {"DA", "DC", "DG", "DT", "DU"}:
        return r

    if r in AA_ORDER:
        return r

    return r


def pick_folder(folder=None):
    if folder is not None:
        return folder

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    selected = filedialog.askdirectory()

    root.destroy()

    return selected


def safe_sem(values):
    if len(values) <= 1:
        return 0.0

    return float(np.std(values, ddof=1) / np.sqrt(len(values)))


def get_matching_atom_by_index(atom_df, idx):
    match = atom_df.loc[atom_df['Index'] == idx]

    if len(match) == 0:
        return None

    return match.iloc[0]


def collect_residue_type_diffs(folder=None, exclude_keys=None, max_percent_diff=None):
    """
    Collect raw residue-level absolute percent volume deviations from AW
    across all systems, grouped by residue type.

    Returns
    -------
    residue_diffs : dict
        {
            'ALA': {'pow': [...], 'prm': [...]},
            ...
        }
    """
    if exclude_keys is None:
        exclude_keys = []

    folder = pick_folder(folder)

    residue_diffs = {
        res: {
            'pow': [],
            'prm': [],
        }
        for res in RES_ORDER
    }

    missing_systems = []
    processed_systems = []

    for subfolder in os.listdir(folder):
        sub_path = os.path.join(folder, subfolder)

        if not os.path.isdir(sub_path):
            continue

        sys_key = subfolder.split('_')[0]

        if sys_key in exclude_keys:
            continue

        aw_path = os.path.join(sub_path, 'aw', 'aw_logs.csv')
        pow_path = os.path.join(sub_path, 'pow', 'pow_logs.csv')
        prm_path = os.path.join(sub_path, 'prm', 'prm_logs.csv')

        if not (os.path.exists(aw_path) and os.path.exists(pow_path) and os.path.exists(prm_path)):
            missing_systems.append(subfolder)
            continue

        print(f'Processing {subfolder}...')

        aw_logs = read_logs2(aw_path, all_=False, balls=True, surfs=False)
        pow_logs = read_logs2(pow_path, all_=False, balls=True, surfs=False)
        prm_logs = read_logs2(prm_path, all_=False, balls=True, surfs=False)

        aw_atoms = aw_logs['atoms']
        pow_atoms = pow_logs['atoms']
        prm_atoms = prm_logs['atoms']

        res_data = {}
        skipped_atoms = 0

        for _, atom in aw_atoms.iterrows():
            idx = atom['Index']

            pow_atom = get_matching_atom_by_index(pow_atoms, idx)
            prm_atom = get_matching_atom_by_index(prm_atoms, idx)

            if pow_atom is None or prm_atom is None:
                skipped_atoms += 1
                continue

            res_name = normalize_residue_label(atom['Residue'])
            res_seq = atom['Residue Sequence']
            res_key = f'{res_name}_{res_seq}'

            if res_key not in res_data:
                res_data[res_key] = {
                    'res_name': res_name,
                    'aw_vol': 0.0,
                    'pow_vol': 0.0,
                    'prm_vol': 0.0,
                }

            res_data[res_key]['aw_vol'] += float(atom['Volume'])
            res_data[res_key]['pow_vol'] += float(pow_atom['Volume'])
            res_data[res_key]['prm_vol'] += float(prm_atom['Volume'])

        filtered_pow = 0
        filtered_prm = 0
        kept_pow = 0
        kept_prm = 0

        for _, vals in res_data.items():
            res_name = vals['res_name']
            aw_vol = vals['aw_vol']
            pow_vol = vals['pow_vol']
            prm_vol = vals['prm_vol']

            if res_name not in residue_diffs:
                continue

            if aw_vol <= 0.0:
                continue

            pow_diff = abs((pow_vol - aw_vol) / aw_vol) * 100.0
            prm_diff = abs((prm_vol - aw_vol) / aw_vol) * 100.0

            if (max_percent_diff is None) or (pow_diff <= max_percent_diff):
                residue_diffs[res_name]['pow'].append(pow_diff)
                kept_pow += 1
            else:
                filtered_pow += 1

            if (max_percent_diff is None) or (prm_diff <= max_percent_diff):
                residue_diffs[res_name]['prm'].append(prm_diff)
                kept_prm += 1
            else:
                filtered_prm += 1

        print(
            f'  residues kept: pow={kept_pow}, prm={kept_prm} | '
            f'filtered: pow={filtered_pow}, prm={filtered_prm} | '
            f'skipped atoms={skipped_atoms}'
        )

        processed_systems.append(subfolder)

    print('\nFinished collecting residue deviations.')
    print(f'Processed systems: {len(processed_systems)}')

    if len(missing_systems) > 0:
        print(f'Skipped systems with missing logs: {len(missing_systems)}')
        for name in missing_systems:
            print(f'  - {name}')

    return residue_diffs


def get_residue_order(residue_diffs, sort_by=None):
    if sort_by is None:
        return [r for r in RES_ORDER if r in residue_diffs]

    stats = []

    for res in RES_ORDER:
        pow_vals = residue_diffs[res]['pow']
        prm_vals = residue_diffs[res]['prm']

        pow_mean = float(np.mean(pow_vals)) if len(pow_vals) > 0 else -1.0
        prm_mean = float(np.mean(prm_vals)) if len(prm_vals) > 0 else -1.0
        combined_mean = float(np.mean(pow_vals + prm_vals)) if len(pow_vals + prm_vals) > 0 else -1.0

        stats.append((res, pow_mean, prm_mean, combined_mean))

    if sort_by == 'pow_mean':
        stats.sort(key=lambda x: x[1], reverse=True)
    elif sort_by == 'prm_mean':
        stats.sort(key=lambda x: x[2], reverse=True)
    elif sort_by == 'combined_mean':
        stats.sort(key=lambda x: x[3], reverse=True)
    else:
        return [r for r in RES_ORDER if r in residue_diffs]

    return [x[0] for x in stats]


def add_jittered_points(ax, x_center, values, color, jitter_width=0.06, alpha=0.2, size=15):
    if len(values) == 0:
        return

    xs = np.random.uniform(
        low=x_center - jitter_width,
        high=x_center + jitter_width,
        size=len(values),
    )

    ax.scatter(
        xs,
        values,
        s=size,
        alpha=alpha,
        color=color,
        linewidths=0,
        zorder=3,
    )


def style_violin_body(violin_dict, color, alpha):
    for body in violin_dict['bodies']:
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(alpha)
        body.set_linewidth(1.0)

    for key in ['cbars', 'cmins', 'cmaxes']:
        if key in violin_dict:
            violin_dict[key].set_color(color)
            violin_dict[key].set_linewidth(1.2)


def plot_residue_deviation_distributions(
    residue_diffs,
    save_path=None,
    sort_by=None,
    show_points=True,
    show_means=True,
    show_counts=True,
    y_range=None,
):
    ordered_residues = get_residue_order(residue_diffs, sort_by=sort_by)

    pow_data = [residue_diffs[res]['pow'] for res in ordered_residues]
    prm_data = [residue_diffs[res]['prm'] for res in ordered_residues]

    x_positions = np.arange(len(ordered_residues), dtype=float)
    pow_positions = x_positions - OFFSET
    prm_positions = x_positions + OFFSET

    fig, ax = plt.subplots(figsize=FIGSIZE)

    safe_pow_data = [vals if len(vals) > 0 else [np.nan] for vals in pow_data]
    safe_prm_data = [vals if len(vals) > 0 else [np.nan] for vals in prm_data]

    pow_violin = ax.violinplot(
        safe_pow_data,
        positions=pow_positions,
        widths=VIOLIN_WIDTH,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )

    prm_violin = ax.violinplot(
        safe_prm_data,
        positions=prm_positions,
        widths=VIOLIN_WIDTH,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )

    style_violin_body(pow_violin, POWER_COLOR, POWER_ALPHA)
    style_violin_body(prm_violin, PRIMITIVE_COLOR, PRIMITIVE_ALPHA)

    max_y = 0.0

    for i, res in enumerate(ordered_residues):
        pow_vals = residue_diffs[res]['pow']
        prm_vals = residue_diffs[res]['prm']

        local_max = 0.0

        if len(pow_vals) > 0:
            local_max = max(local_max, max(pow_vals))

        if len(prm_vals) > 0:
            local_max = max(local_max, max(prm_vals))

        max_y = max(max_y, local_max)

        if show_points:
            add_jittered_points(
                ax=ax,
                x_center=pow_positions[i],
                values=pow_vals,
                color=POWER_COLOR,
                alpha=POINT_ALPHA,
                size=POINT_SIZE,
            )

            add_jittered_points(
                ax=ax,
                x_center=prm_positions[i],
                values=prm_vals,
                color=PRIMITIVE_COLOR,
                alpha=POINT_ALPHA,
                size=POINT_SIZE,
            )

        if show_means:
            if len(pow_vals) > 0:
                pow_mean = float(np.mean(pow_vals))
                ax.scatter(
                    pow_positions[i],
                    pow_mean,
                    marker='_',
                    s=250,
                    linewidths=2.2,
                    color=POWER_COLOR,
                    zorder=5,
                )

            if len(prm_vals) > 0:
                prm_mean = float(np.mean(prm_vals))
                ax.scatter(
                    prm_positions[i],
                    prm_mean,
                    marker='_',
                    s=250,
                    linewidths=2.2,
                    color=PRIMITIVE_COLOR,
                    zorder=5,
                )

    if y_range is not None:
        ax.set_ylim(y_range[0], y_range[1])
        top_y = y_range[1]
    else:
        top_y = max_y * 1.15 if max_y > 0 else 1.0
        ax.set_ylim(0, top_y)

    if show_counts:
        count_y = ax.get_ylim()[1] * 0.98

        for i, res in enumerate(ordered_residues):
            pow_n = len(residue_diffs[res]['pow'])
            prm_n = len(residue_diffs[res]['prm'])

            ax.text(
                x_positions[i],
                count_y,
                f'{pow_n}',
                ha='center',
                va='top',
                fontsize=10,
                color='black',
                rotation=0,
            )

    ax.set_title(TITLE, fontsize=24, pad=18)
    ax.set_xlabel(X_LABEL, fontsize=20, labelpad=12)
    ax.set_ylabel(Y_LABEL, fontsize=20, labelpad=12)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(ordered_residues, rotation=45, ha='right', fontsize=14)

    ax.tick_params(axis='y', labelsize=14, width=2, length=8)
    ax.tick_params(axis='x', width=2, length=8)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    legend_handles = [
        plt.Line2D([0], [0], color=POWER_COLOR, lw=8, alpha=POWER_ALPHA, label='Power vs AW'),
        plt.Line2D([0], [0], color=PRIMITIVE_COLOR, lw=8, alpha=PRIMITIVE_ALPHA, label='Primitive vs AW'),
    ]

    ax.legend(handles=legend_handles, fontsize=14, frameon=False, loc='upper right')

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=DPI, bbox_inches='tight')
        print(f'\nSaved figure to:\n{save_path}')

    plt.show()


def print_summary_table(residue_diffs, ordered_residues):
    print('\nResidue summary:')
    print(
        f'{"Residue":<8} '
        f'{"n_pow":>8} {"mean_pow":>12} {"sem_pow":>12} '
        f'{"n_prm":>8} {"mean_prm":>12} {"sem_prm":>12}'
    )

    for res in ordered_residues:
        pow_vals = residue_diffs[res]['pow']
        prm_vals = residue_diffs[res]['prm']

        pow_mean = float(np.mean(pow_vals)) if len(pow_vals) > 0 else float('nan')
        prm_mean = float(np.mean(prm_vals)) if len(prm_vals) > 0 else float('nan')

        pow_sem = safe_sem(pow_vals) if len(pow_vals) > 0 else float('nan')
        prm_sem = safe_sem(prm_vals) if len(prm_vals) > 0 else float('nan')

        print(
            f'{res:<8} '
            f'{len(pow_vals):>8} {pow_mean:>12.4f} {pow_sem:>12.4f} '
            f'{len(prm_vals):>8} {prm_mean:>12.4f} {prm_sem:>12.4f}'
        )


def main():
    folder = pick_folder()

    if folder is None or folder == '':
        print('No folder selected.')
        return

    residue_diffs = collect_residue_type_diffs(
        folder=folder,
        exclude_keys=EXCLUDE_KEYS,
        max_percent_diff=MAX_PERCENT_DIFF,
    )

    ordered_residues = get_residue_order(residue_diffs, sort_by=SORT_BY)

    print_summary_table(residue_diffs, ordered_residues)

    save_path = os.path.join(folder, SAVE_NAME)

    plot_residue_deviation_distributions(
        residue_diffs=residue_diffs,
        save_path=save_path,
        sort_by=SORT_BY,
        show_points=SHOW_POINTS,
        show_means=SHOW_MEANS,
        show_counts=SHOW_COUNTS,
        y_range=Y_RANGE,
    )


if __name__ == '__main__':
    main()