import os
import sys
import numpy as np
import tkinter as tk
from tkinter import filedialog

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '', '../..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.plot_templates.bar import bar
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2

POWER_COLOR = '#d62728'
PRIMITIVE_COLOR = '#7f3fbf'


def _res_key(atom):
    return f"{atom['Residue']}_{atom['Residue Sequence']}"


def _stats(values):
    if not values:
        return 0.0, 0.0
    avg = float(np.mean(values))
    se = float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0
    return avg, se


def get_res_sa_neighbors(logs):
    """Return outer SA and unique touching-residue count for every residue."""
    idx_to_res = {int(atom['Index']): _res_key(atom) for _, atom in logs['atoms'].iterrows()}
    res_sa = {res: 0.0 for res in set(idx_to_res.values())}
    res_neighbors = {res: set() for res in set(idx_to_res.values())}

    for _, surf in logs['surfs'].iterrows():
        b1, b2 = (int(_) for _ in surf['Balls'])
        r1, r2 = idx_to_res.get(b1), idx_to_res.get(b2)
        if r1 == r2:
            continue
        area = float(surf['Surface Area'])
        if r1 is not None:
            res_sa[r1] += area
        if r2 is not None:
            res_sa[r2] += area
        if r1 is not None and r2 is not None:
            res_neighbors[r1].add(r2)
            res_neighbors[r2].add(r1)

    return res_sa, {res: len(neighbors) for res, neighbors in res_neighbors.items()}


def get_res_data(folder=None, exclude_keys=None, max_percent_diff=None, absolute=True):
    """Aggregate matched atoms into corresponding residues, then average residue-level % differences."""
    if exclude_keys is None:
        exclude_keys = []
    if folder is None:
        root = tk.Tk(); root.withdraw(); root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()

    sys_res_data = {}
    for subfolder in os.listdir(folder):
        my_key = subfolder.split('_')[0]
        if my_key in exclude_keys:
            continue
        sub_path = os.path.join(folder, subfolder)
        try:
            aw = read_logs2(os.path.join(sub_path, 'aw_logs.csv'), all_=False, balls=True, surfs=True)
            pow_ = read_logs2(os.path.join(sub_path, 'pow_logs.csv'), all_=False, balls=True, surfs=True)
            prm = read_logs2(os.path.join(sub_path, 'prm_logs.csv'), all_=False, balls=True, surfs=True)
        except FileNotFoundError:
            print(f"{sub_path} logs not found")
            continue

        aw_sa, aw_neigh = get_res_sa_neighbors(aw)
        pow_sa, pow_neigh = get_res_sa_neighbors(pow_)
        prm_sa, prm_neigh = get_res_sa_neighbors(prm)
        pow_atoms = pow_['atoms'].set_index('Index')
        prm_atoms = prm['atoms'].set_index('Index')
        res_vols = {}

        for _, atom in aw['atoms'].iterrows():
            idx = atom['Index']
            if idx not in pow_atoms.index or idx not in prm_atoms.index:
                continue
            key = _res_key(atom)
            if key not in res_vols:
                res_vols[key] = {'aw': 0.0, 'pow': 0.0, 'prm': 0.0}
            res_vols[key]['aw'] += float(atom['Volume'])
            res_vols[key]['pow'] += float(pow_atoms.loc[idx]['Volume'])
            res_vols[key]['prm'] += float(prm_atoms.loc[idx]['Volume'])

        diffs = {'pow': {'vol': [], 'sa': [], 'neighbors': []}, 'prm': {'vol': [], 'sa': [], 'neighbors': []}}
        for key, vols in res_vols.items():
            vals = {
                'vol': (vols['aw'], vols['pow'], vols['prm']),
                'sa': (aw_sa.get(key, 0.0), pow_sa.get(key, 0.0), prm_sa.get(key, 0.0)),
                'neighbors': (aw_neigh.get(key, 0), pow_neigh.get(key, 0), prm_neigh.get(key, 0))
            }
            for metric, (aw_val, pow_val, prm_val) in vals.items():
                if aw_val <= 0:
                    continue
                for scheme, value in [('pow', pow_val), ('prm', prm_val)]:
                    diff = (value - aw_val) / aw_val * 100.0
                    if absolute:
                        diff = abs(diff)
                    if max_percent_diff is None or abs(diff) <= max_percent_diff:
                        diffs[scheme][metric].append(diff)

        sys_res_data[my_key] = {}
        for scheme in ['pow', 'prm']:
            sys_res_data[my_key][scheme] = {}
            for metric in ['vol', 'sa', 'neighbors']:
                avg, se = _stats(diffs[scheme][metric])
                sys_res_data[my_key][scheme][metric] = {'avg': avg, 'se': se}

    return dict(sorted(sys_res_data.items()))


def plot_data(data, ylim=None, absolute=True, plot_scheme='both'):
    x_names = list(data.keys())
    prefix = 'Avg Abs %' if absolute else 'Avg %'
    titles = {
        'vol': 'Average Residue Volume Difference',
        'sa': 'Average Residue Surface Area Difference',
        'neighbors': 'Average Residue Neighbor Count Difference'
    }

    for metric, ylabel in [('vol', f'{prefix} Volume Diff'), ('sa', f'{prefix} Surface Area Diff'),
                           ('neighbors', f'{prefix} Neighbor Count Diff')]:
        title = titles[metric]
        pow_avg = [data[k]['pow'][metric]['avg'] for k in x_names]
        prm_avg = [data[k]['prm'][metric]['avg'] for k in x_names]
        pow_se = [data[k]['pow'][metric]['se'] for k in x_names]
        prm_se = [data[k]['prm'][metric]['se'] for k in x_names]

        if absolute and ylim is None:
            max_y = max([a + e for a, e in zip(pow_avg, pow_se)] + [a + e for a, e in zip(prm_avg, prm_se)], default=0.0)
            metric_ylim = [0, max_y * 1.1 if max_y > 0 else 1.0]
        else:
            metric_ylim = ylim

        if plot_scheme == 'both':
            bar([pow_avg, prm_avg], x_names=x_names, Show=True, y_axis_title=ylabel, x_axis_title='Model',
                errors=[pow_se, prm_se], y_range=metric_ylim, xtick_label_size=25, ytick_label_size=25,
                ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2, title=title,
                colors=[POWER_COLOR, PRIMITIVE_COLOR], legend_names=['Pow vs AW', 'Prm vs AW'])
        elif plot_scheme == 'pow':
            bar([pow_avg], x_names=x_names, Show=True, y_axis_title=ylabel, x_axis_title='Model',
                errors=[pow_se], y_range=metric_ylim, xtick_label_size=25, ytick_label_size=25,
                ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2, title=title + 'pow',
                colors=[POWER_COLOR], legend_names=['Pow vs AW'])
        elif plot_scheme == 'prm':
            bar([prm_avg], x_names=x_names, Show=True, y_axis_title=ylabel, x_axis_title='Model',
                errors=[prm_se], y_range=metric_ylim, xtick_label_size=25, ytick_label_size=25,
                ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2, title=title + 'prm',
                colors=[PRIMITIVE_COLOR], legend_names=['Prm vs AW'])
        elif plot_scheme == 'separate':
            bar([pow_avg], x_names=x_names, Show=True, y_axis_title=ylabel, x_axis_title='Model',
                errors=[pow_se], y_range=metric_ylim, xtick_label_size=25, ytick_label_size=25,
                ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2, title=title + 'pow',
                colors=[POWER_COLOR], legend_names=['Pow vs AW'])
            bar([prm_avg], x_names=x_names, Show=True, y_axis_title=ylabel, x_axis_title='Model',
                errors=[prm_se], y_range=metric_ylim, xtick_label_size=25, ytick_label_size=25,
                ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2, title=title + 'prm',
                colors=[PRIMITIVE_COLOR], legend_names=['Prm vs AW'])
        else:
            raise ValueError("plot_scheme must be 'both', 'pow', 'prm', or 'separate'")


if __name__ == '__main__':
    absolute = False
    plot_scheme = 'separate'
    res_data = get_res_data(exclude_keys=['A', 'B', 'C'], max_percent_diff=200.0, absolute=absolute)
    print(res_data)
    plot_data(res_data, ylim=None, absolute=absolute, plot_scheme=plot_scheme)
