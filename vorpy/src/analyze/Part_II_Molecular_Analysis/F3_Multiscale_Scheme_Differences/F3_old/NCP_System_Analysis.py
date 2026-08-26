import os
import sys
import tkinter as tk
from tkinter import filedialog

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


DNA_CHAINS = ['a', 'b']
PROTEIN_CHAINS = ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']


def get_group_stats(atoms, chains):
    chains = [c.lower() for c in chains]
    mask = atoms['Chain'].astype(str).str.lower().isin(chains)
    my_atoms = atoms[mask]

    return {
        'count': len(my_atoms),
        'volume': my_atoms['Volume'].sum()
    }


def get_surface_areas(logs):
    atoms = logs['atoms']
    surfs = logs['surfs']

    atom_dict = atoms.set_index('Index').to_dict(orient='index')

    whole_sa = 0.0
    protein_sa = 0.0
    dna_sa = 0.0

    for _, surf in surfs.iterrows():
        ball_1, ball_2 = surf['Balls']

        if ball_1 not in atom_dict or ball_2 not in atom_dict:
            continue

        atom_1 = atom_dict[ball_1]
        atom_2 = atom_dict[ball_2]

        chain_1 = str(atom_1['Chain']).lower()
        chain_2 = str(atom_2['Chain']).lower()

        res_1 = str(atom_1['Residue']).upper()
        res_2 = str(atom_2['Residue']).upper()

        area = float(surf['Surface Area'])

        type_1 = 'dna' if chain_1 in DNA_CHAINS else 'protein' if chain_1 in PROTEIN_CHAINS else 'sol' if res_1 == 'SOL' else 'other'
        type_2 = 'dna' if chain_2 in DNA_CHAINS else 'protein' if chain_2 in PROTEIN_CHAINS else 'sol' if res_2 == 'SOL' else 'other'

        pair = {type_1, type_2}

        # Protein-SOL:
        # contributes to protein SA and total external system SA
        if pair == {'protein', 'sol'}:
            protein_sa += area
            whole_sa += area

        # DNA-SOL:
        # contributes to DNA SA and total external system SA
        elif pair == {'dna', 'sol'}:
            dna_sa += area
            whole_sa += area

        # Protein-DNA:
        # contributes to both group SAs, but is internal to the whole system
        elif pair == {'protein', 'dna'}:
            protein_sa += area
            dna_sa += area

    return whole_sa, protein_sa, dna_sa


def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    log_file = filedialog.askopenfilename(title='Select whole-system logs')
    if not log_file:
        return

    logs = read_logs2(log_file, all_=False, balls=True, surfs=True)

    atoms = logs['atoms']

    protein = get_group_stats(atoms, PROTEIN_CHAINS)
    dna = get_group_stats(atoms, DNA_CHAINS)

    whole = {
        'count': len(atoms),
        'volume': atoms['Volume'].sum()
    }

    whole_sa, protein_sa, dna_sa = get_surface_areas(logs)

    protein['surface_area'] = protein_sa
    dna['surface_area'] = dna_sa
    whole['surface_area'] = whole_sa

    print('\nSystem Statistics')
    print('-' * 65)
    print(f"{'Group':<15}{'Atoms':>10}{'Volume (Å³)':>18}{'SA (Å²)':>18}")
    print('-' * 65)

    for name, data in [('Protein', protein), ('DNA', dna), ('Whole System', whole)]:
        print(f"{name:<15}{data['count']:>10}{data['volume']:>18.2f}{data['surface_area']:>18.2f}")


if __name__ == '__main__':
    main()