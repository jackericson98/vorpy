import os
import sys
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


CHAINS = ['A', 'B', 'C', 'D']

FIGSIZE = (8, 6)
TITLE_SIZE = 24
LABEL_SIZE = 22
TICK_SIZE = 18


def get_chain_data(logs):
    atoms = logs['atoms']
    surfs = logs['surfs']

    atom_dict = atoms.set_index('Index').to_dict(orient='index')

    chain_volumes = {chain: 0.0 for chain in CHAINS}
    chain_sa = {chain: 0.0 for chain in CHAINS}

    # -------------------------
    # Chain volumes
    # -------------------------
    for _, atom in atoms.iterrows():
        chain = str(atom['Chain']).upper()

        if chain in chain_volumes:
            chain_volumes[chain] += float(atom['Volume'])

    # -------------------------
    # Chain surface areas
    # -------------------------
    for _, surf in surfs.iterrows():
        ball_1, ball_2 = surf['Balls']

        if ball_1 not in atom_dict or ball_2 not in atom_dict:
            continue

        chain_1 = str(atom_dict[ball_1]['Chain']).upper()
        chain_2 = str(atom_dict[ball_2]['Chain']).upper()

        area = float(surf['Surface Area'])

        # Internal surface within the same chain
        if chain_1 == chain_2:
            continue

        # Surface belongs to the boundary of each chain it touches
        if chain_1 in chain_sa:
            chain_sa[chain_1] += area

        if chain_2 in chain_sa:
            chain_sa[chain_2] += area

    return chain_volumes, chain_sa


def print_data(chain_volumes, chain_sa):
    print('\nChain Geometry')
    print('-' * 45)
    print(f"{'Chain':<10}{'Volume (Å³)':>17}{'SA (Å²)':>17}")
    print('-' * 45)

    for chain in CHAINS:
        print(f"{chain:<10}{chain_volumes[chain]:>17.3f}{chain_sa[chain]:>17.3f}")

    volumes = np.array([chain_volumes[c] for c in CHAINS])
    sas = np.array([chain_sa[c] for c in CHAINS])

    print('\nVolume')
    print(f'Mean: {np.mean(volumes):.3f} Å³')
    print(f'SD:   {np.std(volumes, ddof=1):.3f} Å³')
    print(f'CV:   {100 * np.std(volumes, ddof=1) / np.mean(volumes):.3f}%')

    print('\nSurface Area')
    print(f'Mean: {np.mean(sas):.3f} Å²')
    print(f'SD:   {np.std(sas, ddof=1):.3f} Å²')
    print(f'CV:   {100 * np.std(sas, ddof=1) / np.mean(sas):.3f}%')


def plot_data(chain_volumes, chain_sa):
    volumes = np.array([chain_volumes[c] for c in CHAINS])
    sas = np.array([3624, 3558, 3596, 3631])

    # Percent deviation from the four-chain mean
    vol_dev = 100 * (volumes - np.mean(volumes)) / np.mean(volumes)
    sa_dev = 100 * (sas - np.mean(sas)) / np.mean(sas)

    x = np.arange(len(CHAINS))
    offset = 0.10

    fig, ax = plt.subplots(figsize=(8, 6))

    # Zero = four-chain mean
    ax.axhline(0, color='black', linewidth=1.5, linestyle='--', alpha=0.7)

    # Volume and SA
    ax.scatter(x - offset, vol_dev, s=140, color='red', label='Volume', zorder=3)
    ax.scatter(x + offset, sa_dev, s=140, color='blue', marker='s', label='Surface Area', zorder=3)

    # Connect equivalent measurements for each chain
    for i in range(len(CHAINS)):
        ax.plot([x[i] - offset, x[i] + offset], [vol_dev[i], sa_dev[i]], color='gray', linewidth=1, alpha=0.5)

    ax.set_title('Geometric Equivalency Across Chains', fontsize=24)
    ax.set_xlabel('Chain', fontsize=22)
    ax.set_ylabel('Deviation from Mean (%)', fontsize=22)

    ax.set_xticks(x)
    ax.set_xticklabels(CHAINS)
    ax.tick_params(axis='both', labelsize=18)

    ax.legend(fontsize=16)
    ax.grid(axis='y', alpha=0.2)

    # Symmetric y-axis
    limit = max(abs(np.concatenate([vol_dev, sa_dev]))) * 1.3
    ax.set_ylim(-limit, limit)

    plt.tight_layout()
    plt.show()


def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    log_file = filedialog.askopenfilename(title='Select p53tet whole-system logs')
    if not log_file:
        return

    logs = read_logs2(log_file, all_=False, balls=True, surfs=True)

    chain_volumes, chain_sa = get_chain_data(logs)

    print_data(chain_volumes, chain_sa)
    plot_data(chain_volumes, chain_sa)


if __name__ == '__main__':
    main()