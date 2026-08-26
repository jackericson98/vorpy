import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def get_res_vol(logs):
    """Calculate total volume for each residue."""
    res_vol = {}

    for _, atom in logs['atoms'].iterrows():
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
        res_vol.setdefault(res_key, 0.0)
        res_vol[res_key] += atom['Volume']

    return res_vol


def get_res_sa(logs):
    """Calculate total boundary surface area for each residue."""
    res_key_dict = {}

    for _, atom in logs['atoms'].iterrows():
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
        res_key_dict[atom['Index']] = res_key

    res_sa = {key: 0.0 for key in set(res_key_dict.values())}

    for _, surf in logs['surfs'].iterrows():
        ball_1, ball_2 = surf['Balls']
        res_1 = res_key_dict.get(ball_1)
        res_2 = res_key_dict.get(ball_2)

        # Surface is internal to a residue
        if res_1 is not None and res_1 == res_2:
            continue

        # Add interface to each residue bordering the surface
        if res_1 is not None:
            res_sa[res_1] += surf['Surface Area']
        if res_2 is not None:
            res_sa[res_2] += surf['Surface Area']

    return res_sa


def plot_residue_data(aw_logs, pow_logs, prm_logs):
    aw_vol, pow_vol, prm_vol = get_res_vol(aw_logs), get_res_vol(pow_logs), get_res_vol(prm_logs)
    aw_sa, pow_sa, prm_sa = get_res_sa(aw_logs), get_res_sa(pow_logs), get_res_sa(prm_logs)

    residues = [key for key in aw_vol if key in pow_vol and key in prm_vol]
    x = np.arange(len(residues))

    # Plot settings
    figsize = (8, 6)
    title_size = 24
    label_size = 22
    tick_size = 16
    legend_size = 18
    line_width = 2.5

    # ----------------------------
    # Residue Volume
    # ----------------------------
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(x, [aw_vol[r] for r in residues], label='AW', color='red', linewidth=line_width)
    ax.plot(x, [pow_vol[r] for r in residues], label='Power', color='blue', linewidth=line_width)
    ax.plot(x, [prm_vol[r] for r in residues], label='Primitive', color='green', linewidth=line_width)

    ax.set_xlabel('Residue', fontsize=label_size)
    ax.set_ylabel(r'Volume ($\AA^3$)', fontsize=label_size)
    ax.set_title('Residue Volumes', fontsize=title_size)
    ax.tick_params(axis='both', labelsize=tick_size)
    ax.legend(fontsize=legend_size)
    ax.grid(axis='y', alpha=0.25)

    # Don't label every residue if there are many
    step = max(1, len(residues) // 10)
    ax.set_xticks(x[::step])
    ax.set_xticklabels(residues[::step], rotation=45, ha='right')

    plt.tight_layout()
    plt.show()

    # ----------------------------
    # Residue Surface Area
    # ----------------------------
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(x, [aw_sa[r] for r in residues], label='AW', color='red', linewidth=line_width)
    ax.plot(x, [pow_sa[r] for r in residues], label='Power', color='blue', linewidth=line_width)
    ax.plot(x, [prm_sa[r] for r in residues], label='Primitive', color='green', linewidth=line_width)

    ax.set_xlabel('Residue', fontsize=label_size)
    ax.set_ylabel(r'Surface Area ($\AA^2$)', fontsize=label_size)
    ax.set_title('Residue Surface Areas', fontsize=title_size)
    ax.tick_params(axis='both', labelsize=tick_size)
    ax.legend(fontsize=legend_size)
    ax.grid(axis='y', alpha=0.25)

    ax.set_xticks(x[::step])
    ax.set_xticklabels(residues[::step], rotation=45, ha='right')

    plt.tight_layout()
    plt.show()


def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    aw_file = filedialog.askopenfilename(title='Select aw logs')
    pow_file = filedialog.askopenfilename(title='Select pow logs')
    prm_file = filedialog.askopenfilename(title='Select prm logs')


    aw_logs = read_logs2(aw_file, all_=False, balls=True, surfs=True)
    pow_logs = read_logs2(pow_file, all_=False, balls=True, surfs=True)
    prm_logs = read_logs2(prm_file, all_=False, balls=True, surfs=True)

    plot_residue_data(aw_logs, pow_logs, prm_logs)


if __name__ == '__main__':
    main()