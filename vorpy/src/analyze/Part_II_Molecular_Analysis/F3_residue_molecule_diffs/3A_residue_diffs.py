import os
import sys
import numpy as np
import tkinter as tk
from tkinter import filedialog



# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)



from vorpy.src.system.system import System
from vorpy.src.group.group import Group
from vorpy.src.analyze.tools.plot_templates.bar import bar
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



def get_res_sa(logs):
    """
    Aggregate surface area by residue from logs['surfs'], counting only
    outer surfaces (where exactly one ball belongs to the residue).
    """
    # Create the dictionary for the residue data
    res_data = {}

    # Create a dictionary to find the residue key for a given surface
    res_key_dict = {}

    # Loop through and get all of the residues
    for i, atom in logs['atoms'].iterrows():
        # Create the residue key
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"

        # Check to see if the residue key is in the dictionary
        if res_key not in res_data:
            res_data[res_key] = [atom['Index']]
        else:
            res_data[res_key].append(atom['Index'])

        # Add the residue key to the dictionary
        res_key_dict[atom['Index']] = res_key

    # Create the residue surface area dictionary
    res_sa_dict = {}

    # Loop through the surfaces and add the surface area if they are outer surfaces
    for i, surf in logs['surfs'].iterrows():
        # Get the first ball in the surface
        ball_1, ball_2 = surf['Balls']

        # Get the residue key for one of the balls (if present)
        if ball_1 in res_key_dict:
            res_key = res_key_dict[ball_1]
        elif ball_2 in res_key_dict:
            res_key = res_key_dict[ball_2]
        else:
            # No residue for either ball
            continue

        # Initialize residue SA if needed
        if res_key not in res_sa_dict:
            res_sa_dict[res_key] = 0.0

        # Only count surfaces where exactly one ball is in this residue
        both_in_same_res = (
            ball_1 in res_data[res_key] and
            ball_2 in res_data[res_key]
        )
        if both_in_same_res:
            continue

        res_sa_dict[res_key] += surf['Surface Area']

    return res_sa_dict


def get_res_data(folder=None, exclude_keys=None, get_sa=False, max_percent_diff=None):
    """
    For each system (subfolder), compute residue-level volumes and SAs
    for AW, Power, and Primitive (PRM), then compute the average absolute
    % differences vs AW across residues.

    We store:
        - Power vs AW: avg abs % vol diff, avg abs % SA diff, and SEs
        - Primitive vs AW: same

    Optionally exclude residue-level outliers whose absolute % difference
    exceeds max_percent_diff (applied separately for Power and Primitive,
    for volume and SA).

    Parameters
    ----------
    folder : str or None
        Path to parent folder containing system subfolders.
        If None, a folder dialog will be opened.
    exclude_keys : list of str or None
        Subfolder prefixes to skip (e.g. ['A', 'B', 'C']).
    get_sa : bool
        Whether to compute surface-area-based differences.
    max_percent_diff : float or None
        If not None, residues with abs(% diff) > max_percent_diff
        are excluded from the averages (for each scheme separately).

    Returns
    -------
    sys_res_data : dict
        {
          sys_key: {
            'pow_vs_aw': {
                'avg vol diff': float,
                'avg sa diff': float,
                'se vol diff': float,
                'se sa diff': float
            },
            'prm_vs_aw': {
                'avg vol diff': float,
                'avg sa diff': float,
                'se vol diff': float,
                'se sa diff': float
            }
          },
          ...
        }
    """
    if exclude_keys is None:
        exclude_keys = []

    if folder is None:
        # Get the dropbox (or root) folder via GUI
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()

    # Create the dictionary for the system residue averages
    sys_res_data = {}

    # Go through the folder and get the systems
    for subfolder in os.listdir(folder):
        my_key = subfolder.split('_')[0]
        if my_key in exclude_keys:
            continue

        sub_path = os.path.join(folder, subfolder)

        # Read logs for AW, POW, PRM
        try:
            aw_logs = read_logs2(os.path.join(sub_path, 'aw_logs.csv'), all_=False, balls=True, surfs=True)
            pow_logs = read_logs2(os.path.join(sub_path, 'pow_logs.csv'), all_=False, balls=True, surfs=True)
            prm_logs = read_logs2(os.path.join(sub_path, 'prm_logs.csv'), all_=False, balls=True, surfs=True)
        except FileNotFoundError:
            aw_logs = read_logs2(os.path.join(sub_path, 'aw', 'aw_logs.csv'), all_=False, balls=True, surfs=True)
            pow_logs = read_logs2(os.path.join(sub_path, 'pow', 'pow_logs.csv'), all_=False, balls=True, surfs=True)
            prm_logs = read_logs2(os.path.join(sub_path, 'prm', 'prm_logs.csv'), all_=False, balls=True, surfs=True)

        # Get residue-level surface areas if requested
        if get_sa:
            aw_res_sa = get_res_sa(aw_logs)
            pow_res_sa = get_res_sa(pow_logs)
            prm_res_sa = get_res_sa(prm_logs)
        else:
            aw_res_sa = {}
            pow_res_sa = {}
            prm_res_sa = {}

        # Create the dictionary for the residue data
        # res_data[res_key] = {
        #   'aw':  {'vol': ..., 'sa': ...},
        #   'pow': {'vol': ..., 'sa': ...},
        #   'prm': {'vol': ..., 'sa': ...}
        # }
        res_data = {}

        # Aggregate volumes by residue
        for i, atom in aw_logs['atoms'].iterrows():
            idx = atom['Index']

            # Get the corresponding POW and PRM atoms by index
            pow_atom = pow_logs['atoms'].loc[pow_logs['atoms']['Index'] == idx].to_dict(orient='records')
            prm_atom = prm_logs['atoms'].loc[prm_logs['atoms']['Index'] == idx].to_dict(orient='records')

            if not pow_atom or not prm_atom:
                # If we can't find matching atoms, skip
                continue

            pow_atom = pow_atom[0]
            prm_atom = prm_atom[0]

            # Create the residue key
            res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"

            # Initialize residue entry if needed
            if res_key not in res_data:
                res_data[res_key] = {
                    'aw':  {'vol': 0.0, 'sa': 0.0},
                    'pow': {'vol': 0.0, 'sa': 0.0},
                    'prm': {'vol': 0.0, 'sa': 0.0},
                }

            # Sum volumes for each scheme
            res_data[res_key]['aw']['vol']  += atom['Volume']
            res_data[res_key]['pow']['vol'] += pow_atom['Volume']
            res_data[res_key]['prm']['vol'] += prm_atom['Volume']

        # Add the surface areas to the residues (if available)
        for res_key, vals in res_data.items():
            if get_sa:
                vals['aw']['sa']  = aw_res_sa.get(res_key, 0.0)
                vals['pow']['sa'] = pow_res_sa.get(res_key, 0.0)
                vals['prm']['sa'] = prm_res_sa.get(res_key, 0.0)
            else:
                vals['aw']['sa']  = 0.0
                vals['pow']['sa'] = 0.0
                vals['prm']['sa'] = 0.0

        # Compute absolute % differences per residue (Power vs AW, Primitive vs AW)
        pow_vol_diffs = []
        pow_sa_diffs = []
        prm_vol_diffs = []
        prm_sa_diffs = []

        # (Optional) track how many got filtered out for debugging
        pow_vol_filtered = 0
        pow_sa_filtered = 0
        prm_vol_filtered = 0
        prm_sa_filtered = 0

        for res_key, vals in res_data.items():
            aw_vol = vals['aw']['vol']
            pow_vol = vals['pow']['vol']
            prm_vol = vals['prm']['vol']
            aw_sa = vals['aw']['sa']
            pow_sa = vals['pow']['sa']
            prm_sa = vals['prm']['sa']

            # Volume % diffs vs AW
            if aw_vol > 0.0:
                pow_vol_diff = abs((pow_vol - aw_vol) / aw_vol) * 100.0
                prm_vol_diff = abs((prm_vol - aw_vol) / aw_vol) * 100.0

                if (max_percent_diff is None) or (pow_vol_diff <= max_percent_diff):
                    pow_vol_diffs.append(pow_vol_diff)
                else:
                    pow_vol_filtered += 1

                if (max_percent_diff is None) or (prm_vol_diff <= max_percent_diff):
                    prm_vol_diffs.append(prm_vol_diff)
                else:
                    prm_vol_filtered += 1

            # Surface area % diffs vs AW
            if get_sa and aw_sa > 0.0:
                pow_sa_diff = abs((pow_sa - aw_sa) / aw_sa) * 100.0
                prm_sa_diff = abs((prm_sa - aw_sa) / aw_sa) * 100.0

                if (max_percent_diff is None) or (pow_sa_diff <= max_percent_diff):
                    pow_sa_diffs.append(pow_sa_diff)
                else:
                    pow_sa_filtered += 1

                if (max_percent_diff is None) or (prm_sa_diff <= max_percent_diff):
                    prm_sa_diffs.append(prm_sa_diff)
                else:
                    prm_sa_filtered += 1

        # Optional debug print
        if max_percent_diff is not None:
            print(
                f"{my_key}: "
                f"filtered pow_vol={pow_vol_filtered}, prm_vol={prm_vol_filtered}, "
                f"pow_sa={pow_sa_filtered}, prm_sa={prm_sa_filtered} above {max_percent_diff:.1f}%"
            )

        # Safeguards & stats for Power
        if len(pow_vol_diffs) == 0:
            pow_avg_vol_diff = 0.0
            pow_se_vol_diff = 0.0
        else:
            pow_avg_vol_diff = float(np.mean(pow_vol_diffs))
            pow_se_vol_diff = (
                float(np.std(pow_vol_diffs, ddof=1) / np.sqrt(len(pow_vol_diffs)))
                if len(pow_vol_diffs) > 1 else
                0.0
            )

        if len(pow_sa_diffs) == 0:
            pow_avg_sa_diff = 0.0
            pow_se_sa_diff = 0.0
        else:
            pow_avg_sa_diff = float(np.mean(pow_sa_diffs))
            pow_se_sa_diff = (
                float(np.std(pow_sa_diffs, ddof=1) / np.sqrt(len(pow_sa_diffs)))
                if len(pow_sa_diffs) > 1 else
                0.0
            )

        # Safeguards & stats for Primitive
        if len(prm_vol_diffs) == 0:
            prm_avg_vol_diff = 0.0
            prm_se_vol_diff = 0.0
        else:
            prm_avg_vol_diff = float(np.mean(prm_vol_diffs))
            prm_se_vol_diff = (
                float(np.std(prm_vol_diffs, ddof=1) / np.sqrt(len(prm_vol_diffs)))
                if len(prm_vol_diffs) > 1 else
                0.0
            )

        if len(prm_sa_diffs) == 0:
            prm_avg_sa_diff = 0.0
            prm_se_sa_diff = 0.0
        else:
            prm_avg_sa_diff = float(np.mean(prm_sa_diffs))
            prm_se_sa_diff = (
                float(np.std(prm_sa_diffs, ddof=1) / np.sqrt(len(prm_sa_diffs)))
                if len(prm_sa_diffs) > 1 else
                0.0
            )

        # Store per-system data
        sys_res_data[my_key] = {
            'pow_vs_aw': {
                'avg vol diff': pow_avg_vol_diff,
                'avg sa diff': pow_avg_sa_diff,
                'se vol diff': pow_se_vol_diff,
                'se sa diff': pow_se_sa_diff,
            },
            'prm_vs_aw': {
                'avg vol diff': prm_avg_vol_diff,
                'avg sa diff': prm_avg_sa_diff,
                'se vol diff': prm_se_vol_diff,
                'se sa diff': prm_se_sa_diff,
            }
        }

    return sys_res_data


def plot_data(sys_res_data):
    """
    Make TWO grouped bar plots:

    1) Volume:
        - Avg absolute % volume difference (Power vs AW)
        - Avg absolute % volume difference (Primitive vs AW)

    2) Surface Area:
        - Avg absolute % SA difference (Power vs AW)
        - Avg absolute % SA difference (Primitive vs AW)
    """
    x_names = list(sys_res_data.keys())

    # Volume stats
    pow_vol_avg = [sys_res_data[k]['pow_vs_aw']['avg vol diff'] for k in x_names]
    prm_vol_avg = [sys_res_data[k]['prm_vs_aw']['avg vol diff'] for k in x_names]

    pow_vol_se = [sys_res_data[k]['pow_vs_aw']['se vol diff'] for k in x_names]
    prm_vol_se = [sys_res_data[k]['prm_vs_aw']['se vol diff'] for k in x_names]

    # SA stats
    pow_sa_avg = [sys_res_data[k]['pow_vs_aw']['avg sa diff'] for k in x_names]
    prm_sa_avg = [sys_res_data[k]['prm_vs_aw']['avg sa diff'] for k in x_names]

    pow_sa_se = [sys_res_data[k]['pow_vs_aw']['se sa diff'] for k in x_names]
    prm_sa_se = [sys_res_data[k]['prm_vs_aw']['se sa diff'] for k in x_names]

    # Plot 1: Volume
    bar(
        [pow_vol_avg, prm_vol_avg],
        x_names=x_names,
        legend_names=['Power vs AW', 'Primitive vs AW'],
        Show=True,
        y_axis_title='Avg Absolute % Difference',
        x_axis_title='Model',
        title='Residue-Level Absolute Volume Percent Differences (vs AW)',
        errors=[pow_vol_se, prm_vol_se],
        y_range=[0, None],
        xtick_label_size=25,
        ytick_label_size=25,
        ylabel_size=30,
        xlabel_size=30,
        tick_length=12,
        tick_width=2
    )

    # Plot 2: Surface Area
    bar(
        [pow_sa_avg, prm_sa_avg],
        x_names=x_names,
        legend_names=['Power vs AW', 'Primitive vs AW'],
        Show=True,
        y_axis_title='Avg Absolute % Difference',
        x_axis_title='Model',
        title='Residue-Level Absolute Surface Area Percent Differences (vs AW)',
        errors=[pow_sa_se, prm_sa_se],
        y_range=[0, None],
        xtick_label_size=25,
        ytick_label_size=25,
        ylabel_size=30,
        xlabel_size=30,
        tick_length=12,
        tick_width=2
    )


if __name__ == '__main__':
    # Example: drop systems A/B/C and clip residue outliers above 200% diff
    sys_res_data = get_res_data(
        get_sa=True,
        exclude_keys=['A', 'B', 'C'],
        max_percent_diff=200.0   # set to None for no outlier filtering
    )

    plot_data(sys_res_data)
