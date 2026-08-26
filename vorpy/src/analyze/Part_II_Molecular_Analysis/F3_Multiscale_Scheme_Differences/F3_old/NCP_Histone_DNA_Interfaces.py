import os
import sys
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


# ============================================================
# SETTINGS
# ============================================================

DNA_CHAINS = ['a', 'b']

PROTEIN_CHAINS = ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

DNA_1_COLOR = 'blue'
DNA_2_COLOR = 'red'

FIGSIZE = (10, 6)
TITLE_SIZE = 24
LABEL_SIZE = 22
TICK_SIZE = 16
LEGEND_SIZE = 18
LINE_WIDTH = 2.5


def get_dna_contact_profile(logs, dna_chains):
    """
    Sum DNA-protein Voronoi interface area for every nucleotide.

    Returns:
        {
            chain: {
                residue_sequence: {
                    'residue': residue name,
                    'area': total DNA-protein interface area
                }
            }
        }
    """

    atoms = logs['atoms']
    surfs = logs['surfs']

    # Fast atom index -> atom information lookup
    atom_dict = atoms.set_index('Index').to_dict(orient='index')

    # Initialize every DNA residue at zero so non-contacting nucleotides
    # still appear in the profile
    profile = {}

    for _, atom in atoms.iterrows():
        chain = str(atom['Chain']).lower()

        if chain not in dna_chains:
            continue

        res_seq = int(atom['Residue Sequence'])
        res_name = atom['Residue']

        profile.setdefault(chain, {})
        profile[chain].setdefault(res_seq, {'residue': res_name, 'area': 0.0})

    # Go through every Voronoi surface
    for _, surf in surfs.iterrows():
        ball_1, ball_2 = surf['Balls']

        if ball_1 not in atom_dict or ball_2 not in atom_dict:
            continue

        atom_1 = atom_dict[ball_1]
        atom_2 = atom_dict[ball_2]

        chain_1 = str(atom_1['Chain']).lower()
        chain_2 = str(atom_2['Chain']).lower()

        # Keep only DNA-histone surfaces.
        if chain_1 in dna_chains and chain_2 in PROTEIN_CHAINS:
            dna_atom = atom_1
        elif chain_2 in dna_chains and chain_1 in PROTEIN_CHAINS:
            dna_atom = atom_2
        else:
            continue

        chain = str(dna_atom['Chain']).lower()
        res_seq = int(dna_atom['Residue Sequence'])

        profile[chain][res_seq]['area'] += float(surf['Surface Area'])

    return profile


def print_profile(profile):
    """Print nucleotide-level DNA-protein interface areas."""

    for chain, residues in profile.items():
        print(f"\nDNA Chain {chain}")
        print("-" * 40)

        for res_seq in sorted(residues):
            data = residues[res_seq]
            print(f"{data['residue']:>3} {res_seq:>4}: {data['area']:8.3f} Å²")

        total = sum(data['area'] for data in residues.values())
        print(f"\nTotal DNA-protein interface: {total:.3f} Å²")


def plot_profile(profile):
    """Plot DNA-protein interface area by nucleotide."""

    fig, ax = plt.subplots(figsize=FIGSIZE)

    colors = [DNA_1_COLOR, DNA_2_COLOR]

    for i, (chain, residues) in enumerate(profile.items()):
        positions = sorted(residues)
        areas = [residues[pos]['area'] for pos in positions]

        ax.plot(positions, areas, linewidth=LINE_WIDTH, color=colors[i % len(colors)], label=f'DNA Chain {chain}')

    ax.set_title('DNA–Histone Contact Profile', fontsize=TITLE_SIZE)
    ax.set_xlabel('Nucleotide Position', fontsize=LABEL_SIZE)
    ax.set_ylabel(r'DNA–Histone Interface Area ($\AA^2$)', fontsize=LABEL_SIZE)

    ax.tick_params(axis='both', labelsize=TICK_SIZE)
    ax.legend(fontsize=LEGEND_SIZE)
    ax.grid(axis='y', alpha=0.25)

    plt.tight_layout()
    plt.show()


def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    log_file = filedialog.askopenfilename(title='Select whole-system logs')

    if not log_file:
        return

    # We only need atoms + surfaces
    logs = read_logs2(log_file, all_=False, balls=True, surfs=True)

    profile = get_dna_contact_profile(logs, DNA_CHAINS)

    print_profile(profile)
    plot_profile(profile)


if __name__ == '__main__':
    main()