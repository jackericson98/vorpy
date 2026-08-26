import os
import sys
import numpy as np
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


# ============================================================
# SETTINGS
# ============================================================

LIGAND_RESIDUE = 'JZ4'
LIGAND_SEQUENCE = 164
LIGAND_CHAIN = 'B'

FIGSIZE = (8, 6)
TITLE_SIZE = 24
LABEL_SIZE = 22
TICK_SIZE = 16


def get_ligand_residue_interfaces(logs):
    atoms = logs['atoms']
    surfs = logs['surfs']

    atom_dict = atoms.set_index('Index').to_dict(orient='index')

    # Get ligand atom indices
    ligand_atoms = set(
        atoms[
            (atoms['Residue'].astype(str).str.upper() == LIGAND_RESIDUE.upper()) &
            (atoms['Residue Sequence'] == LIGAND_SEQUENCE) &
            (atoms['Chain'].astype(str).str.upper() == LIGAND_CHAIN.upper())
        ]['Index']
    )

    print(f'Ligand atoms found: {len(ligand_atoms)}')

    if not ligand_atoms:
        raise ValueError('No ligand atoms found. Check ligand residue, sequence, and chain.')

    # residue_interfaces[key] = total ligand-residue interface area
    residue_interfaces = {}

    for _, surf in surfs.iterrows():
        ball_1, ball_2 = surf['Balls']

        ball_1_lig = ball_1 in ligand_atoms
        ball_2_lig = ball_2 in ligand_atoms

        # Need exactly one ligand atom
        if ball_1_lig == ball_2_lig:
            continue

        other_ball = ball_2 if ball_1_lig else ball_1

        if other_ball not in atom_dict:
            continue

        atom = atom_dict[other_ball]

        # Ignore solvent
        if str(atom['Residue']).upper() == 'SOL':
            continue

        res_key = f"{atom['Residue']}{int(atom['Residue Sequence'])}"
        residue_interfaces.setdefault(res_key, 0.0)
        residue_interfaces[res_key] += float(surf['Surface Area'])

    return residue_interfaces


def print_interfaces(residue_interfaces):
    print('\nLigand–Residue Interface Areas')
    print('-' * 35)

    for residue, area in sorted(residue_interfaces.items(), key=lambda x: x[1], reverse=True):
        print(f'{residue:<12}{area:8.3f} Å²')

    print('-' * 35)
    print(f'Contacting residues: {len(residue_interfaces)}')
    print(f'Total interface:     {sum(residue_interfaces.values()):.3f} Å²')


def plot_interfaces(residue_interfaces):
    # Residue classes
    residue_types = {
        'hydrophobic': {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO', 'GLY'},
        'polar': {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS'},
        'positive': {'LYS', 'ARG', 'HIS'},
        'negative': {'ASP', 'GLU'}
    }

    type_colors = {
        'hydrophobic': '#4daf4a',
        'polar': '#377eb8',
        'positive': '#984ea3',
        'negative': '#e41a1c',
        'other': '#999999'
    }

    # Sort largest -> smallest
    data = sorted(residue_interfaces.items(), key=lambda x: x[1], reverse=True)

    residues = [x[0] for x in data]
    areas = [x[1] for x in data]

    # Determine residue type from first 3 characters
    def get_residue_type(residue):
        res_name = residue[:3].upper()

        for res_type, names in residue_types.items():
            if res_name in names:
                return res_type

        return 'other'

    res_types = [get_residue_type(res) for res in residues]
    colors = [type_colors[t] for t in res_types]

    mean_area = np.mean(areas)

    fig, ax = plt.subplots(figsize=(10, 6))

    bars = ax.bar(residues, areas, color=colors, edgecolor='black', linewidth=0.7)

    # Mean line
    ax.axhline(
        mean_area,
        color='black',
        linestyle='--',
        linewidth=2,
        label=f'Mean = {mean_area:.1f} Å²'
    )

    # Values above bars
    for bar, area in zip(bars, areas):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(areas) * 0.015,
            f'{area:.1f}',
            ha='center',
            va='bottom',
            fontsize=11,
            rotation=0
        )

    # Legend entries for residue classes
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor=type_colors['hydrophobic'], edgecolor='black', label='Hydrophobic'),
        Patch(facecolor=type_colors['polar'], edgecolor='black', label='Polar'),
        Patch(facecolor=type_colors['positive'], edgecolor='black', label='Positive'),
        Patch(facecolor=type_colors['negative'], edgecolor='black', label='Negative')
    ]

    mean_handle = plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2, label=f'Mean = {mean_area:.1f} Å²')

    ax.legend(handles=legend_handles + [mean_handle], fontsize=13, ncol=2)

    ax.set_title('Ligand–Residue Interface Area', fontsize=24)
    ax.set_xlabel('Binding-Site Residue', fontsize=22)
    ax.set_ylabel(r'Interface Area ($\AA^2$)', fontsize=22)

    ax.tick_params(axis='y', labelsize=16)
    ax.tick_params(axis='x', labelsize=14, rotation=45)

    # Give labels above bars some room
    ax.set_ylim(0, max(areas) * 1.18)

    ax.grid(axis='y', alpha=0.2)

    plt.tight_layout()
    plt.show()

def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    log_file = filedialog.askopenfilename(title='Select whole-system logs')

    if not log_file:
        return

    logs = read_logs2(log_file, all_=False, balls=True, surfs=True)

    residue_interfaces = get_ligand_residue_interfaces(logs)

    print_interfaces(residue_interfaces)
    plot_interfaces(residue_interfaces)


if __name__ == '__main__':
    main()