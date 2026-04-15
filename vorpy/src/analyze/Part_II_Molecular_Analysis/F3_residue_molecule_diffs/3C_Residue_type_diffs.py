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
SHOW_POW = False
SHOW_PRM = True


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
USE_PERCENTILE_LIMIT = True
PERCENTILE_LIMIT = 99.5

USE_MANUAL_Y_RANGE = True
MANUAL_Y_RANGE = [-10.0, 10.0]

SYMMETRIC_Y = False

TITLE = 'Figure 3C — Residue-Type Signed Volume Differences Across All Models'
Y_LABEL = 'Signed % Volume Difference from AW'
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


def check_index_uniqueness(atom_df, label='atoms'):
    counts = atom_df['Index'].value_counts()
    dupes = counts[counts > 1]

    print(f'\n[{label}] total atoms = {len(atom_df)}')
    print(f'[{label}] unique indices = {atom_df["Index"].nunique()}')
    print(f'[{label}] duplicated index count = {len(dupes)}')

    if len(dupes) > 0:
        print(f'[{label}] first duplicated indices:')
        print(dupes.head(20))


def build_atom_lookup(atom_df):
    lookup = {}

    for _, atom in atom_df.iterrows():
        key = (
            int(atom['Index']),
            str(atom['Residue']).strip().upper(),
            int(atom['Residue Sequence']),
            str(atom['Atom']).strip().upper() if 'Atom' in atom_df.columns else None,
        )

        if key in lookup:
            if not isinstance(lookup[key], list):
                lookup[key] = [lookup[key]]
            lookup[key].append(atom)
        else:
            lookup[key] = atom

    return lookup


def get_matching_atom(atom, lookup):
    key = (
        int(atom['Index']),
        str(atom['Residue']).strip().upper(),
        int(atom['Residue Sequence']),
        str(atom['Atom']).strip().upper() if 'Atom' in atom.index else None,
    )

    match = lookup.get(key, None)

    if isinstance(match, list):
        return match[0]

    return match


def get_residue_key(atom):
    res_name = normalize_residue_label(atom['Residue'])
    res_seq = int(atom['Residue Sequence'])

    chain = ''
    for candidate in ['Chain', 'Chain ID', 'ChainID', 'Subunit']:
        if candidate in atom.index:
            chain = str(atom[candidate]).strip()
            break

    return f'{chain}_{res_name}_{res_seq}'


def clip_values_by_percentile(values, percentile_limit=99.0):
    if len(values) == 0:
        return values

    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return []

    lower_pct = (100.0 - percentile_limit) / 2.0
    upper_pct = 100.0 - lower_pct

    lo = np.percentile(vals, lower_pct)
    hi = np.percentile(vals, upper_pct)

    clipped = [v for v in values if lo <= v <= hi]

    print(f"Clipping [{percentile_limit}%]: kept {len(clipped)}/{len(values)} "
          f"(range: {lo:.3f} → {hi:.3f})")

    return clipped


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

        pow_lookup = build_atom_lookup(pow_atoms)
        prm_lookup = build_atom_lookup(prm_atoms)

        # check_index_uniqueness(aw_atoms, 'AW')
        # check_index_uniqueness(pow_atoms, 'POW')
        # check_index_uniqueness(prm_atoms, 'PRM')

        res_data = {}
        skipped_atoms = 0

        for _, atom in aw_atoms.iterrows():
            idx = atom['Index']

            pow_atom = get_matching_atom(atom, pow_lookup)
            prm_atom = get_matching_atom(atom, prm_lookup)

            if pow_atom is None or prm_atom is None:
                skipped_atoms += 1
                continue

            res_name = normalize_residue_label(atom['Residue'])
            res_key = get_residue_key(atom)

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

            if pow_atom is not None:
                if str(pow_atom['Residue']).strip().upper() != str(atom['Residue']).strip().upper():
                    print('POW residue mismatch:')
                    print('  AW :', atom[['Index', 'Residue', 'Residue Sequence']].to_dict())
                    print('  POW:', pow_atom[['Index', 'Residue', 'Residue Sequence']].to_dict())
                    break

        filtered_pow = 0
        filtered_prm = 0
        kept_pow = 0
        kept_prm = 0

        for _, vals in res_data.items():
            res_name = vals['res_name']
            aw_vol = vals['aw_vol']
            pow_vol = vals['pow_vol']
            prm_vol = vals['prm_vol']
            if aw_vol > pow_vol:
                print(f"residue volumes: aw = {aw_vol}, pow = {pow_vol}")

            if res_name not in residue_diffs:
                continue

            if aw_vol <= 0.0:
                continue

            pow_diff = ((pow_vol - aw_vol) / aw_vol) * 100.0
            prm_diff = ((prm_vol - aw_vol) / aw_vol) * 100.0

            if pow_diff < 0:
                print(
                    f"NEGATIVE RESIDUE DIFF: {res_name} | pow_diff = {pow_diff:.4f} | aw = {aw_vol:.4f} | pow = {pow_vol:.4f}")

            if (max_percent_diff is None) or (abs(pow_diff) <= max_percent_diff):
                residue_diffs[res_name]['pow'].append(pow_diff)
                kept_pow += 1
            else:
                filtered_pow += 1

            if (max_percent_diff is None) or (abs(prm_diff) <= max_percent_diff):
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


def get_plot_y_limits(values, use_percentile_limit=False, percentile_limit=99.0,
                      use_manual_y_range=False, manual_y_range=None,
                      symmetric_y=False):
    if use_manual_y_range and manual_y_range is not None:
        return manual_y_range[0], manual_y_range[1]

    if len(values) == 0:
        return -1.0, 1.0

    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return -1.0, 1.0

    if use_percentile_limit:
        lower_pct = (100.0 - percentile_limit) / 2.0
        upper_pct = 100.0 - lower_pct

        y_min = float(np.percentile(vals, lower_pct))
        y_max = float(np.percentile(vals, upper_pct))
    else:
        y_min = float(np.min(vals))
        y_max = float(np.max(vals))

    if symmetric_y:
        max_abs = max(abs(y_min), abs(y_max))
        pad = max_abs * 0.08 if max_abs > 0 else 1.0
        return -(max_abs + pad), (max_abs + pad)

    span = y_max - y_min

    if span <= 0:
        pad = 1.0
    else:
        pad = span * 0.08

    return y_min - pad, y_max + pad


def debug_dataset_summary(label, data_lists):
    flat = [v for vals in data_lists for v in vals if np.isfinite(v)]

    print(f"\n--- {label} SUMMARY ---")
    print(f"groups: {len(data_lists)}")
    print(f"total values: {len(flat)}")

    if len(flat) == 0:
        print("min: None")
        print("max: None")
        print("mean: None")
        return

    print(f"min: {np.min(flat):.6f}")
    print(f"max: {np.max(flat):.6f}")
    print(f"mean: {np.mean(flat):.6f}")


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

    if not SHOW_POW and not SHOW_PRM:
        print("Nothing to plot: both SHOW_POW and SHOW_PRM are False.")
        return

    pow_data = [residue_diffs[res]['pow'] for res in ordered_residues]
    prm_data = [residue_diffs[res]['prm'] for res in ordered_residues]

    print("Initial POW counts:", [len(v) for v in pow_data[:10]])
    print("Initial PRM counts:", [len(v) for v in prm_data[:10]])
    print("Total POW values:", sum(len(v) for v in pow_data))
    print("Total PRM values:", sum(len(v) for v in prm_data))

    debug_dataset_summary("POW BEFORE CLIP", pow_data)
    debug_dataset_summary("PRM BEFORE CLIP", prm_data)

    if USE_PERCENTILE_LIMIT:
        lower_pct = (100.0 - PERCENTILE_LIMIT) / 2.0
        upper_pct = 100.0 - lower_pct

        if SHOW_POW:
            all_pow_for_clip = [v for vals in pow_data for v in vals if np.isfinite(v)]

            print(f"POW values available for clipping: {len(all_pow_for_clip)}")

            if len(all_pow_for_clip) == 0:
                print("Skipping POW percentile clipping: no finite POW values found.")
            else:
                pow_lo = np.percentile(all_pow_for_clip, lower_pct)
                pow_hi = np.percentile(all_pow_for_clip, upper_pct)

                print(f"\nGLOBAL POW CLIP [{PERCENTILE_LIMIT}%]")
                print(f"lower_pct = {lower_pct:.4f}")
                print(f"upper_pct = {upper_pct:.4f}")
                print(f"bounds = {pow_lo:.6f} -> {pow_hi:.6f}")

                before_pow_counts = [len(vals) for vals in pow_data]

                pow_data = [
                    [v for v in vals if pow_lo <= v <= pow_hi]
                    for vals in pow_data
                ]

                after_pow_counts = [len(vals) for vals in pow_data]

                print(f"POW total before = {sum(before_pow_counts)}")
                print(f"POW total after  = {sum(after_pow_counts)}")
                print(f"POW removed      = {sum(before_pow_counts) - sum(after_pow_counts)}")

        if SHOW_PRM:
            all_prm_for_clip = [v for vals in prm_data for v in vals if np.isfinite(v)]

            print(f"PRM values available for clipping: {len(all_prm_for_clip)}")

            if len(all_prm_for_clip) == 0:
                print("Skipping PRM percentile clipping: no finite PRM values found.")
            else:
                prm_lo = np.percentile(all_prm_for_clip, lower_pct)
                prm_hi = np.percentile(all_prm_for_clip, upper_pct)

                print(f"\nGLOBAL PRM CLIP [{PERCENTILE_LIMIT}%]")
                print(f"lower_pct = {lower_pct:.4f}")
                print(f"upper_pct = {upper_pct:.4f}")
                print(f"bounds = {prm_lo:.6f} -> {prm_hi:.6f}")

                before_prm_counts = [len(vals) for vals in prm_data]

                prm_data = [
                    [v for v in vals if prm_lo <= v <= prm_hi]
                    for vals in prm_data
                ]

                after_prm_counts = [len(vals) for vals in prm_data]

                print(f"PRM total before = {sum(before_prm_counts)}")
                print(f"PRM total after  = {sum(after_prm_counts)}")
                print(f"PRM removed      = {sum(before_prm_counts) - sum(after_prm_counts)}")

    debug_dataset_summary("POW AFTER CLIP", pow_data)
    debug_dataset_summary("PRM AFTER CLIP", prm_data)

    x_positions = np.arange(len(ordered_residues), dtype=float)
    if SHOW_POW and SHOW_PRM:
        pow_positions = x_positions - OFFSET
        prm_positions = x_positions + OFFSET
    elif SHOW_POW:
        pow_positions = x_positions
        prm_positions = x_positions + OFFSET
    elif SHOW_PRM:
        pow_positions = x_positions - OFFSET
        prm_positions = x_positions

    fig, ax = plt.subplots(figsize=FIGSIZE)

    safe_pow_data = [vals if len(vals) > 0 else [np.nan] for vals in pow_data]
    safe_prm_data = [vals if len(vals) > 0 else [np.nan] for vals in prm_data]


    if SHOW_POW:
        pow_violin = ax.violinplot(
            safe_pow_data,
            positions=pow_positions,
            widths=VIOLIN_WIDTH,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        style_violin_body(pow_violin, POWER_COLOR, POWER_ALPHA)

    if SHOW_PRM:
        prm_violin = ax.violinplot(
            safe_prm_data,
            positions=prm_positions,
            widths=VIOLIN_WIDTH,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        style_violin_body(prm_violin, PRIMITIVE_COLOR, PRIMITIVE_ALPHA)

    for i, res in enumerate(ordered_residues):
        pow_vals = pow_data[i]
        prm_vals = prm_data[i]

        if show_points:
            if SHOW_POW:
                add_jittered_points(
                    ax=ax,
                    x_center=pow_positions[i],
                    values=pow_vals,
                    color=POWER_COLOR,
                    alpha=POINT_ALPHA,
                    size=POINT_SIZE,
                )

            if SHOW_PRM:
                add_jittered_points(
                    ax=ax,
                    x_center=prm_positions[i],
                    values=prm_vals,
                    color=PRIMITIVE_COLOR,
                    alpha=POINT_ALPHA,
                    size=POINT_SIZE,
                )

        if show_means:
            if SHOW_POW and len(pow_vals) > 0:
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

            if SHOW_PRM and len(prm_vals) > 0:
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

    all_vals = []

    if SHOW_POW:
        for vals in pow_data:
            all_vals.extend(vals)

    if SHOW_PRM:
        for vals in prm_data:
            all_vals.extend(vals)

    print("\n" + "=" * 80)
    print("DEBUG: VALUES ACTUALLY BEING PLOTTED")
    print("=" * 80)

    all_pow_vals = []
    for vals in pow_data:
        all_pow_vals.extend(vals)

    all_prm_vals = []
    for vals in prm_data:
        all_prm_vals.extend(vals)

    print("\nPOW stats:")
    print("min:", np.min(all_pow_vals) if len(all_pow_vals) > 0 else None)
    print("max:", np.max(all_pow_vals) if len(all_pow_vals) > 0 else None)
    print("mean:", np.mean(all_pow_vals) if len(all_pow_vals) > 0 else None)

    print("\nPRM stats:")
    print("min:", np.min(all_prm_vals) if len(all_prm_vals) > 0 else None)
    print("max:", np.max(all_prm_vals) if len(all_prm_vals) > 0 else None)
    print("mean:", np.mean(all_prm_vals) if len(all_prm_vals) > 0 else None)

    neg_pow = sum(v < 0 for v in all_pow_vals)
    total_pow = len(all_pow_vals)
    print(f"\nPOW negative fraction: {neg_pow}/{total_pow} = {neg_pow / total_pow if total_pow else 0:.4f}")

    print("=" * 80 + "\n")

    y_min, y_max = get_plot_y_limits(
        values=all_vals,
        use_percentile_limit=USE_PERCENTILE_LIMIT,
        percentile_limit=PERCENTILE_LIMIT,
        use_manual_y_range=USE_MANUAL_Y_RANGE if y_range is None else False,
        manual_y_range=MANUAL_Y_RANGE if y_range is None else None,
        symmetric_y=SYMMETRIC_Y,
    )

    if y_range is not None:
        ax.set_ylim(y_range[0], y_range[1])
    else:
        ax.set_ylim(y_min, y_max)

    print(f"Plot y-limits: {ax.get_ylim()}")
    print(f"USE_PERCENTILE_LIMIT={USE_PERCENTILE_LIMIT}, "
          f"PERCENTILE_LIMIT={PERCENTILE_LIMIT}, "
          f"USE_MANUAL_Y_RANGE={USE_MANUAL_Y_RANGE}, "
          f"MANUAL_Y_RANGE={MANUAL_Y_RANGE}, "
          f"SYMMETRIC_Y={SYMMETRIC_Y}")

    if show_counts:
        y0, y1 = ax.get_ylim()
        count_y = y1 - 0.02 * (y1 - y0)

        for i, res in enumerate(ordered_residues):
            pow_n = len(pow_data[i])
            prm_n = len(prm_data[i])

            ax.text(
                x_positions[i],
                count_y,
                f'{pow_n}',
                ha='center',
                va='top',
                fontsize=14,
                color='black',
                rotation=0,
            )

    ax.set_title(TITLE, fontsize=48, pad=18)

    ax.set_xlabel(X_LABEL, fontsize=40, labelpad=12)
    ax.set_ylabel(Y_LABEL, fontsize=40, labelpad=12)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(ordered_residues, rotation=45, ha='right', fontsize=28)

    ax.tick_params(axis='y', labelsize=28, width=2, length=8)
    ax.tick_params(axis='x', labelsize=28, width=2, length=8)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    legend_handles = []

    if SHOW_POW:
        legend_handles.append(
            plt.Line2D([0], [0], color=POWER_COLOR, lw=8, alpha=POWER_ALPHA, label='Power vs AW')
        )

    if SHOW_PRM:
        legend_handles.append(
            plt.Line2D([0], [0], color=PRIMITIVE_COLOR, lw=8, alpha=PRIMITIVE_ALPHA, label='Primitive vs AW')
        )

    ax.axhline(0, color='black', linewidth=1.5, linestyle='--', zorder=1)
    ax.axvline(19.5, color='black', linewidth=1.5, linestyle='--', zorder=1)
    ax.legend(
        handles=legend_handles,
        fontsize=12,
        frameon=False,
        loc='upper right',
        bbox_to_anchor=(1.0, 0.92)  # ↓ lower from 1.0 → 0.92
    )

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
