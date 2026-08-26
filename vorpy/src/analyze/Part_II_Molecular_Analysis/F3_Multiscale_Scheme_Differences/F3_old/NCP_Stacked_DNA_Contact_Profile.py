import os
import sys
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


# ============================================================
# SETTINGS - CHANGE THESE TO YOUR NCP CHAINS
# ============================================================

DNA_CHAINS = ['a', 'b']

# Histone number -> chain ID
HISTONE_CHAINS = {
    '2': 'c',
    '3': 'd',
    '4': 'e',
    '5': 'f',
    '6': 'g',
    '7': 'h',
    '8': 'i',
    '9': 'j'
}

HISTONE_COLORS = {
    '2': '#e41a1c',
    '3': '#377eb8',
    '4': '#4daf4a',
    '5': '#984ea3',
    '6': '#ff7f00',
    '7': '#00bfc4',
    '8': '#f781bf',
    '9': '#a65628'
}

FIGSIZE = (10, 6)
TITLE_SIZE = 24
LABEL_SIZE = 22
TICK_SIZE = 16
LEGEND_SIZE = 14


def get_histone_contact_profile(logs):
    atoms = logs['atoms']
    surfs = logs['surfs']

    atom_dict = atoms.set_index('Index').to_dict(orient='index')
    chain_to_histone = {str(chain).lower(): histone for histone, chain in HISTONE_CHAINS.items()}

    # profile[dna_chain][residue][histone] = area
    profile = {}

    # Initialize every nucleotide
    for _, atom in atoms.iterrows():
        dna_chain = str(atom['Chain']).lower()

        if dna_chain not in DNA_CHAINS:
            continue

        res_seq = int(atom['Residue Sequence'])
        res_name = atom['Residue']

        profile.setdefault(dna_chain, {})
        profile[dna_chain].setdefault(res_seq, {
            'residue': res_name,
            'histones': {histone: 0.0 for histone in HISTONE_CHAINS}
        })

    # Read DNA-histone surfaces
    for _, surf in surfs.iterrows():
        ball_1, ball_2 = surf['Balls']

        if ball_1 not in atom_dict or ball_2 not in atom_dict:
            continue

        atom_1 = atom_dict[ball_1]
        atom_2 = atom_dict[ball_2]

        chain_1 = str(atom_1['Chain']).lower()
        chain_2 = str(atom_2['Chain']).lower()

        # Determine which side is DNA and which side is histone.
        if chain_1 in DNA_CHAINS and chain_2 in chain_to_histone:
            dna_atom = atom_1
            histone_chain = chain_2
        elif chain_2 in DNA_CHAINS and chain_1 in chain_to_histone:
            dna_atom = atom_2
            histone_chain = chain_1
        else:
            continue

        dna_chain = str(dna_atom['Chain']).lower()
        res_seq = int(dna_atom['Residue Sequence'])
        histone = chain_to_histone[histone_chain]

        profile[dna_chain][res_seq]['histones'][histone] += float(surf['Surface Area'])

    return profile


def plot_stacked_profile(profile):
    for dna_chain, residues in profile.items():

        positions = sorted(residues)
        bottom = np.zeros(len(positions))

        fig, ax = plt.subplots(figsize=FIGSIZE)

        for histone in HISTONE_CHAINS:
            areas = np.array([
                residues[pos]['histones'][histone]
                for pos in positions
            ])

            ax.bar(
                positions,
                areas,
                bottom=bottom,
                width=1.0,
                color=HISTONE_COLORS[histone],
                label=f'Histone {histone}'
            )

            bottom += areas

        ax.set_title(f'DNA–Histone Contact Profile — DNA Chain {dna_chain}', fontsize=TITLE_SIZE)
        ax.set_xlabel('Nucleotide Position', fontsize=LABEL_SIZE)
        ax.set_ylabel(r'Interface Area ($\AA^2$)', fontsize=LABEL_SIZE)

        ax.tick_params(axis='both', labelsize=TICK_SIZE)
        ax.legend(fontsize=LEGEND_SIZE, ncol=2)

        ax.set_xlim(min(positions) - 1, max(positions) + 1)

        plt.tight_layout()
        plt.show()


def print_histone_totals(profile):
    totals = {histone: 0.0 for histone in HISTONE_CHAINS}

    for dna_chain, residues in profile.items():
        for res_seq, data in residues.items():
            for histone, area in data['histones'].items():
                totals[histone] += area

    print('\nDNA-Histone Interface Totals')
    print('-' * 35)

    for histone, area in totals.items():
        print(f'Histone {histone}: {area:8.3f} Å²')

    print('-' * 35)
    print(f'Total:     {sum(totals.values()):8.3f} Å²')


def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    log_file = filedialog.askopenfilename(title='Select whole-system logs')

    if not log_file:
        return

    logs = read_logs2(log_file, all_=False, balls=True, surfs=True)

    profile = get_histone_contact_profile(logs)

    print_histone_totals(profile)
    plot_stacked_profile(profile)


if __name__ == '__main__':
    main()