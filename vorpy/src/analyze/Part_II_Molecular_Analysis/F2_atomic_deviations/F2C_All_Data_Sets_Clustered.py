import os
import sys
import tkinter as tk
from tkinter import filedialog
from typing import Dict, List, Optional, Set

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from collections import defaultdict


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.compare.sol_facing_percentage import compute_sol_facing_percent
from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import compare_manual_vs_ml
from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import make_sol_facing_binary
from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import print_cluster_summary
from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import run_clustering
from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import summarize_clusters


SMALL_MOL_ATOM_ALIASES = {
    # -----------------
    # HYDROGENS
    # -----------------
    # Do not split H by atom number here.
    # Hydrogens will be reassigned later using AW/Pow volume regime.
    'H': 'H',

    **{f'H{i}': 'H' for i in range(1, 40)},

    # -----------------
    # OXYGENS
    # -----------------
    # Keep the current coherent split, but give it chemical language.
    # These are still empirical small-molecule oxygen classes.
    'O1': 'O_acceptor_low',
    'O2': 'O_acceptor_low',
    'O5': 'O_acceptor_low',
    'O7': 'O_acceptor_low',

    'O3': 'O_acceptor_high',
    'O4': 'O_acceptor_high',
    'O6': 'O_acceptor_high',
    'O10': 'O_acceptor_high',

    'O': 'O_other',

    # -----------------
    # NITROGENS
    # -----------------
    'N': 'N',
    'N1': 'N',

    # -----------------
    # CARBONS
    # -----------------
    'C1': 'C_small',
    'C2': 'C_small',
    'C4': 'C_small',

    'C3': 'C_medium',
    'C8': 'C_medium',

    'C5': 'C_large',
    'C6': 'C_large',
    'C7': 'C_large',
    'C9': 'C_large',

    'C': 'C_generic',
}

SMALL_MOL_MANUAL_GROUPS = {
    # -----------------
    # EDTA
    # -----------------
    'EDTA': {
        # Oxygens
        'O3': 'O_Mg_neighbor',
        'O4': 'O_Mg_neighbor',
        'O6': 'O_Mg_neighbor',
        'O1': 'O_Mg_neighbor',

        'O': 'O_sol_facing',
        'O2': 'O_sol_facing',
        'O5': 'O_sol_facing',
        'O7': 'O_sol_facing',

        # Carbons
        'C3': 'C_carboxyl_linker',
        'C5': 'C_carboxyl_linker',
        'C7': 'C_carboxyl_linker',
        'C9': 'C_carboxyl_linker',

        'C': 'C_amine_methylene',
        'C1': 'C_amine_methylene',
        'C2': 'C_amine_methylene',
        'C4': 'C_amine_methylene',
        'C6': 'C_amine_methylene',
        'C8': 'C_amine_methylene',

        # Hydrogens
        'H4': 'H_amine_acetate_methylene',
        'H5': 'H_amine_acetate_methylene',
        'H6': 'H_amine_acetate_methylene',
        'H7': 'H_amine_acetate_methylene',
        'H8': 'H_amine_acetate_methylene',
        'H9': 'H_amine_acetate_methylene',
        'H10': 'H_amine_acetate_methylene',
        'H11': 'H_amine_acetate_methylene',

        'H': 'H_amine_ethylene',
        'H1': 'H_amine_ethylene',
        'H2': 'H_amine_ethylene',
        'H3': 'H_amine_ethylene',
    },

    # Your logs may label EDTA as MOL
    'MOL': {
        'O3': 'O_Mg_neighbor',
        'O4': 'O_Mg_neighbor',
        'O6': 'O_Mg_neighbor',
        'O1': 'O_Mg_neighbor',

        'O': 'O_sol_facing',
        'O2': 'O_sol_facing',
        'O5': 'O_sol_facing',
        'O7': 'O_sol_facing',

        'C3': 'C_carboxyl_linker',
        'C5': 'C_carboxyl_linker',
        'C7': 'C_carboxyl_linker',
        'C9': 'C_carboxyl_linker',

        # EDTA ethylene bridge carbons: bonded to each other + N + H/H
        'C': 'C_amine_ethylene',
        'C1': 'C_amine_ethylene',

        # EDTA acetate methylene carbons: bonded to N + carboxyl carbon + H/H
        'C2': 'C_amine_acetate_methylene',
        'C4': 'C_amine_acetate_methylene',
        'C6': 'C_amine_acetate_methylene',
        'C8': 'C_amine_acetate_methylene',

        # Hydrogens
        'H4': 'H_amine_acetate_methylene',
        'H5': 'H_amine_acetate_methylene',
        'H6': 'H_amine_acetate_methylene',
        'H7': 'H_amine_acetate_methylene',
        'H8': 'H_amine_acetate_methylene',
        'H9': 'H_amine_acetate_methylene',
        'H10': 'H_amine_acetate_methylene',
        'H11': 'H_amine_acetate_methylene',

        'H': 'H_amine_ethylene',
        'H1': 'H_amine_ethylene',
        'H2': 'H_amine_ethylene',
        'H3': 'H_amine_ethylene',
    },

    'JZ4': {
        # Aromatic carbons: 2 carbons, 1 hydrogen
        'C5': 'C_aromatic_C2_H',
        'C6': 'C_aromatic_C2_H',
        'C7': 'C_aromatic_C2_H',
        'C8': 'C_aromatic_C2_H',

        # Saturated/alkyl-like carbons: 2 carbons, 2 hydrogens
        'C2': 'C_C2_H2',
        'C3': 'C_C2_H2',

        # Methyl-like carbon: 1 carbon, 3 hydrogens
        'C1': 'C_C1_H3',

        # Aromatic bridge/substituted carbon: 3 carbons
        'C4': 'C_aromatic_C3',

        # Aromatic carbon bonded to oxygen
        'C9': 'C_aromatic_C2_O',

        # Hydrogens
        'H1': 'H_terminal_methyl',
        'H2': 'H_terminal_methyl',
        'H3': 'H_terminal_methyl',

        'H4': 'H_chain_methylene',
        'H5': 'H_chain_methylene',
        'H6': 'H_chain_methylene',
        'H7': 'H_chain_methylene',

        'H21': 'H_aromatic_CH',
        'H18': 'H_aromatic_CH',
        'H19': 'H_aromatic_CH',
        'H20': 'H_aromatic_CH',

        'H22': 'H_aromatic_OH',
    },
}

DB1976_INDEX_GROUPS = {
    # Carbons
    0: 'C_CN2',
    9: 'C_CN2',
    14: 'C_CN2',
    23: 'C_CN2',

    5: 'C_C3',
    19: 'C_C3',

    12: 'C_C2_Se',
    27: 'C_C2_Se',

    3: 'C_C2',
    6: 'C_C2',
    7: 'C_C2',
    17: 'C_C2',
    20: 'C_C2',
    21: 'C_C2',

    13: 'C_C2_H',
    26: 'C_C2_H',

    2: 'C_C2_N',
    4: 'C_C2_N',
    16: 'C_C2_N',
    18: 'C_C2_N',

    # Nitrogens
    10: 'N_H2_C',
    11: 'N_H2_C',
    24: 'N_H2_C',
    25: 'N_H2_C',

    1: 'N_C2_H',
    15: 'N_C2_H',

    8: 'N_C2',
    22: 'N_C2',

    # Hydrogens
    30: 'H_terminal_NH2',
    31: 'H_terminal_NH2',
    32: 'H_terminal_NH2',
    33: 'H_terminal_NH2',
    36: 'H_terminal_NH2',
    37: 'H_terminal_NH2',
    38: 'H_terminal_NH2',
    39: 'H_terminal_NH2',

    34: 'H_ring_C_H',
    40: 'H_ring_C_H',

    29: 'H_ring_N_H',
    35: 'H_ring_N_H',
}


SMALL_MOL_GROUP_ORDER = [
    # Carbons first
    'C_amine_ethylene',
    'C_amine_acetate_methylene',
    'C_carboxyl_linker',
    'C_C2_H2',
    'C_C1_H3',
    'C_aromatic_C3',
    'C_aromatic_C2_H',
    'C_aromatic_C2_O',
    'C_CN2',
    'C_C3',
    'C_C2_Se',
    'C_C2',
    'C_C2_H',
    'C_C2_N',

    # Oxygens second
    'O_Mg_neighbor',
    'O_sol_facing',

    # Nitrogens third
    'N',
    'N_C2_H',
    'N_H2_C',
    'N_C2',

    # Other
    'SE',
]


def sort_plot_groups(plot_stats_df: pd.DataFrame, molecule_class: str) -> pd.DataFrame:
    df = plot_stats_df.copy()
    mol = str(molecule_class).strip().lower().replace(' ', '_')

    if mol == 'small_molecule':
        group_order = SMALL_MOL_GROUP_ORDER
    elif mol == 'rna':
        group_order = RNA_GROUP_ORDER
    elif mol == 'dna':
        group_order = DNA_GROUP_ORDER
    elif mol == 'protein':
        group_order = PROTEIN_GROUP_ORDER
    else:
        group_order = None

    if group_order is None:
        df = df.sort_values(['Mean_AW', 'Mean_Pow']).reset_index(drop=True)
        df['PlotNumber'] = np.arange(1, len(df) + 1)
        return df

    order_lookup = {
        name: i
        for i, name in enumerate(group_order)
    }

    df['SortOrder'] = df['GroupName'].map(order_lookup).fillna(9999)
    df = df.sort_values(['SortOrder', 'Mean_AW', 'Mean_Pow']).reset_index(drop=True)
    df['PlotNumber'] = np.arange(1, len(df) + 1)
    df = df.drop(columns=['SortOrder'])

    return df


PROTEIN_ATOM_ALIASES = {
    # terminal hydrogens
    'H1': 'H',
    'H2': 'H',
    'H3': 'H',

    # alpha hydrogens
    'HA1': 'HA',
    'HA2': 'HA',

    # beta hydrogens
    'HB1': 'HB',
    'HB2': 'HB',
    'HB3': 'HB',

    # gamma hydrogens
    'HG1': 'HG',
    'HG2': 'HG',
    'HG11': 'HG',
    'HG12': 'HG',
    'HG13': 'HG',
    'HG21': 'HG',
    'HG22': 'HG',
    'HG23': 'HG',

    # delta hydrogens
    'HD1': 'HD',
    'HD2': 'HD',
    'HD3': 'HD',
    'HD11': 'HD',
    'HD12': 'HD',
    'HD13': 'HD',
    'HD21': 'HD',
    'HD22': 'HD',
    'HD23': 'HD',

    # epsilon hydrogens
    'HE1': 'HE',
    'HE2': 'HE',
    'HE3': 'HE',
    'HE21': 'HE',
    'HE22': 'HE',

    # zeta hydrogens
    'HZ1': 'HZ',
    'HZ2': 'HZ',
    'HZ3': 'HZ',

    # eta hydrogens
    'HH1': 'HH',
    'HH2': 'HH',
    'HH11': 'HH',
    'HH12': 'HH',
    'HH21': 'HH',
    'HH22': 'HH',
}


CHEM_CLASSES = {
    # Backbone
    'C': 'backbone_C',
    'CA': 'backbone_CA',
    'N': 'backbone_N',
    'O': 'backbone_O',

    # Aliphatic carbons
    'CB': 'aliphatic_C',
    'CG': 'aliphatic_C',
    'CD': 'aliphatic_C',
    'CE': 'aliphatic_C',

    # Aromatic carbons
    'CD1': 'aromatic_C',
    'CD2': 'aromatic_C',
    'CE1': 'aromatic_C',
    'CE2': 'aromatic_C',
    'CZ': 'aromatic_C',

    # Oxygens
    'OD1': 'carbonyl_O',
    'OD2': 'carbonyl_O',
    'OE1': 'carbonyl_O',
    'OE2': 'carbonyl_O',
    'OH': 'hydroxyl_O',

    # Nitrogens
    'ND1': 'aromatic_N',
    'NE2': 'aromatic_N',
    'NZ': 'charged_N',
}


def get_chem_class(group_name):
    name = str(group_name).upper()

    # Strip residue suffix after underscore
    base = name.split('_')[0]

    # Handle wildcard names like CD*_FY
    base = base.replace('*', '')

    # Backbone carbons stay separate
    if base == 'C':
        return 'backbone_C'

    if base == 'CA':
        return 'alpha_C'

    # Aliphatic side-chain carbons
    if base in {'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE'}:
        # Aromatic-specific labels should not merge with aliphatic labels
        if any(tag in name for tag in ['FY', 'FHWY', 'W', 'HIS', 'TYR', 'PHE', 'TRP']):
            return f'aromatic_{base}'

        return f'aliphatic_{base}'

    # Aromatic ring extension carbons
    if base in {'CE1', 'CE2', 'CE3', 'CZ', 'CZ2', 'CZ3', 'CH2'}:
        return f'aromatic_{base}'

    # Nucleic/small-molecule carbons should not merge with protein carbons
    if base in {'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'}:
        return f'other_{base}'

    # Oxygens
    if base in {'O', 'OT1', 'OT2'}:
        return 'backbone_O'

    if base in {'OD1', 'OD2', 'OE1', 'OE2'}:
        return 'carboxyl_O'

    if base in {'OG', 'OG1', 'OH'}:
        return 'hydroxyl_O'

    # Nitrogens
    if base == 'N':
        return 'backbone_N'

    if base in {'ND1', 'NE1', 'NE2'}:
        return f'aromatic_or_amide_{base}'

    if base in {'NZ', 'NH1', 'NH2', 'NE'}:
        return f'charged_N_{base}'

    # Hydrogens: keep prefix families separate
    if base.startswith('H'):
        return base

    return base


PROTEIN_ATOM_RESIDUE_ALIASES = {
    # Residue-aware overrides go here first.
    # Exact match: (atom_name, residue_name)
    # Fallback for all other residues: (atom_name, '*')
    ('CA', 'GLY'): 'CA_G',
    ('CA', '*'): 'CA_*',
    ('CB', 'SER'): 'CB_S',
    ('CB', 'VAL'): 'CB_IV',
    ('CB', 'ILE'): 'CB_IV',
    ('CB', 'THR'): 'CB_T',
    ('CB', 'ALA'): 'CB_A',
    ('CB', '*'): 'CB_*',
    ('CD', 'LYS'): 'CD_K',
    ('CD', 'ARG'): 'CD_PR',
    ('CD', 'PRO'): 'CD_PR',
    ('CD', 'GLU'): 'CD_EQ',
    ('CD', 'GLN'): 'CD_EQ',
    ('CD', 'ILE'): 'CD_LI',
    ('CD1', 'LEU'): 'CD_LI',
    ('CD2', 'LEU'): 'CD_LI',
    ('CD1', 'PHE'): 'CD*_FY',
    ('CD1', 'TYR'): 'CD*_FY',
    ('CD2', 'PHE'): 'CD*_FY',
    ('CD2', 'TYR'): 'CD*_FY',
    ('CD1', 'TRP'): 'CD1_W',
    ('CD1', '*'): 'CD1_*',
    ('CD2', 'TRP'): 'CD2_W',
    ('CE', 'LYS'): 'CE_K',
    ('CE', 'MET'): 'CE_M',
    ('CE', '*'): 'CE_*',
    ('CE1', 'HIS'): 'CE1_H',
    ('CE1', 'PHE'): 'CE*_FY',
    ('CE1', 'TYR'): 'CE*_FY',
    ('CE2', 'PHE'): 'CE*_FY',
    ('CE2', 'TYR'): 'CE*_FY',
    ('CE2', 'TRP'): 'CE*_W',
    ('CE3', 'TRP'): 'CE*_W',
    ('CG', 'LEU'): 'CG_L',
    ('CG', 'HIS'): 'CG_FHWY',
    ('CG', 'TYR'): 'CG_FHWY',
    ('CG', 'PHE'): 'CG_FHWY',
    ('CG', 'TRP'): 'CG_FHWY',
    ('CG', 'ASN'): 'CG_DN',
    ('CG', 'ASP'): 'CG_DN',
    ('CG', 'ARG'): 'CG_EKPQR',
    ('CG', 'GLU'): 'CG_EKPQR',
    ('CG', 'PRO'): 'CG_EKPQR',
    ('CG', 'GLN'): 'CG_EKPQR',
    ('CG', 'LYS'): 'CG_EKPQR',
    ('CG', '*'): 'CG_*',
    ('CG1', 'VAL'): 'CG_V',
    ('CG1', 'ILE'): 'CG_I',
    ('CG2', 'ILE'): 'CG_ITV',
    ('CG2', 'VAL'): 'CG_ITV',
    ('CG2', 'THR'): 'CG_ITV',
    ('H1', '*'): 'H*',
    ('H2', '*'): 'H*',
    ('H3', '*'): 'H*',
    ('HA', '*'): 'HA*',
    ('HA1', '*'): 'HA*',
    ('HA2', '*'): 'HA*',
    ('HB', '*'): 'HB*',
    ('HB1', '*'): 'HB*',
    ('HB2', '*'): 'HB*',
    ('HB3', '*'): 'HB*',
    ('HD', '*'): 'HD*',
    ('HD11', '*'): 'HD*',
    ('HD12', '*'): 'HD*',
    ('HD13', '*'): 'HD*',
    ('HD2', '*'): 'HD*',
    ('HD21', '*'): 'HD*',
    ('HD22', '*'): 'HD*',
    ('HD23', '*'): 'HD*',
    ('HD3', '*'): 'HD*',
    ('HE', '*'): 'HE*',
    ('HE1', '*'): 'HE*',
    ('HE2', '*'): 'HE*',
    ('HE21', '*'): 'HE*',
    ('HE22', '*'): 'HE*',
    ('HE3', '*'): 'HE*',
    ('HG', '*'): 'HG*',
    ('HG1', '*'): 'HG*',
    ('HG11', '*'): 'HG*',
    ('HG12', '*'): 'HG*',
    ('HG13', '*'): 'HG*',
    ('HG2', '*'): 'HG*',
    ('HG21', '*'): 'HG*',
    ('HG22', '*'): 'HG*',
    ('HG23', '*'): 'HG*',
    ('HH', '*'): 'HH*',
    ('HH11', '*'): 'HH*',
    ('HH12', '*'): 'HH*',
    ('HH2', '*'): 'HH*',
    ('HH21', '*'): 'HH*',
    ('HH22', '*'): 'HH*',
    ('HN', '*'): 'HN',
    ('HZ', '*'): 'HZ*',
    ('HZ1', '*'): 'HZ*',
    ('HZ2', '*'): 'HZ*',
    ('HZ3', '*'): 'HZ*',
    ('N', '*'): 'N',
    ('N', 'PRO'): 'NP',
    ('ND1', '*'): 'ND1',
    ('ND2', '*'): 'ND2',
    ('NE', '*'): 'NE*',
    ('NE1', '*'): 'NE*',
    ('NE2', 'HIS'): 'NE*',
    ('NE2', 'GLN'): 'NEQ',
    ('NH1', '*'): 'NH*',
    ('NH2', '*'): 'NH*',
    ('O', '*'): 'O',
    ('OH', '*'): 'OH',
    # ('OC1', '*'): 'OC*',
    # ('OC2', '*'): 'OC*',
    ('OD1', 'ASP'): 'ODA',
    ('OD1', 'ASN'): 'ODN',
    ('OD2', 'ASP'): 'ODA',
    ('OE1', 'GLU'): 'OEE',
    ('OE2', 'GLU'): 'OEE',
    ('OE1', 'GLN'): 'OEQ',
    ('OG', '*'): 'OG',
    ('OG1', '*'): 'OG'
}

PROTEIN_GROUP_ORDER = [
    # -----------------
    # CARBONS
    # backbone first
    # -----------------
    'C',
    'CA_*',

    # beta carbons
    'CB_A',
    'CB_IV',
    'CB_T',
    'CB_*',

    # gamma carbons
    'CG_L',
    'CG_EKPQR',
    'CG_DN',
    'CG_FHWY',
    'CG_ITV',

    # delta carbons
    'CD_K',
    'CD_PR',
    'CD_EQ',
    'CD_LI',
    'CD*_FY',
    'CD2',

    # epsilon / aromatic extension carbons
    'CE_K',
    'CE*_FY',
    'CZ',

    # -----------------
    # OXYGENS
    # backbone, carboxylate, hydroxyl
    # -----------------
    'O',
    'ODA',
    'OEE',
    'OG',

    # -----------------
    # NITROGENS
    # backbone, charged, side-chain
    # -----------------
    'N',
    'NZ',
    'NE*',
    'NH*',

    # -----------------
    # HYDROGENS
    # backbone/alpha first, then side-chain depth
    # -----------------
    'H',
    'HN',
    'HA*',
    'HB*',
    'HG*',
    'HD',
    'HD*',
    'HE*',
    'HH*',
    'HZ*',
]


RNA_ATOM_ALIASES = {
    # phosphate
    'P': 'P',
    'OP1': 'O_phosphate',
    'OP2': 'O_phosphate',
    'O1P': 'O_phosphate',
    'O2P': 'O_phosphate',

    'O5\'': 'O_backbone',
    'O3\'': 'O_backbone',

    # sugar
    'C1\'': 'C_sugar',
    'C2\'': 'C_sugar',
    'C3\'': 'C_sugar',
    'C4\'': 'C_sugar',
    # 'C5\'': 'C_sugar',

    'O2\'': 'O_sugar',
    'O4\'': 'O_sugar',

    # base carbons
    # 'C2': 'C_base',
    # 'C4': 'C_base',
    # 'C5': 'C_base',
    # 'C6': 'C_base',
    # 'C8': 'C_base',

    # base nitrogens
    # 'N1': 'N_base',
    'N2': 'N_base1',
    'N3': 'N_base2',
    'N4': 'N_base1',
    'N6': 'N_base1',
    'N7': 'N_base2',
    # 'N9': 'N_base',

    # base oxygens
    'O2': 'O_base',
    'O4': 'O_base',
    'O6': 'O_base',

    # sugar hydrogens
    'H1\'': 'H_sugar',
    'H2\'': 'H_sugar',
    'H2\'\'': 'H_sugar',
    'H2\'1': 'H_sugar',
    'H2\'2': 'H_sugar',
    'H3\'': 'H_sugar',
    'H4\'': 'H_sugar',
    'H5\'': 'H_sugar',
    'H5\'\'': 'H_sugar',
    'H5\'1': 'H_sugar',
    'H5\'2': 'H_sugar',

    # base hydrogens
    'H1': 'H_base',
    'H2': 'H_base',
    'H3': 'H_base',
    'H5': 'H_base',
    'H6': 'H_base',
    'H8': 'H_base',

    # exocyclic amino hydrogens
    'H21': 'H_base_exocyclic',
    'H22': 'H_base_exocyclic',
    'H41': 'H_base_exocyclic',
    'H42': 'H_base_exocyclic',
    'H61': 'H_base_exocyclic',
    'H62': 'H_base_exocyclic',

    # rare terminal hydroxyl hydrogens
    'HO5\'': 'H_backbone',
    'HO3\'': 'H_backbone',
}

RNA_ATOM_RESIDUE_ALIASES = {
    # C5 split
    ('C5', 'U'): 'C5_UC',
    ('C5', 'C'): 'C5_UC',
    ('C5', 'G'): 'C5_AG',
    ('C5', 'A'): 'C5_AG',

    # N1 split
    ('N1', 'A'): 'N1_A',
    ('N1', 'G'): 'N1_G',
    ('N1', 'U'): 'N1_UC',
    ('N1', 'C'): 'N1_UC',
}

RNA_GROUP_ORDER = [
    # Carbons first: backbone/sugar, then base
    "C_sugar",
    "C1'",
    "C2'",
    "C3'",
    "C4'",
    "C5'",
    "C2",
    "C4",
    "C5_UC",
    "C5_AG",
    "C6",
    "C8",

    # Oxygens second: phosphate/backbone/sugar, then base
    "O_phosphate",
    "O_backbone",
    "O_sugar",
    "O_base",
    "O_base",
    "O_base",

    # Nitrogens third
    "N_base1",
    "N_base2",
    "N1_A",
    "N1_G",
    "N1_UC",
    "N2",
    "N3",
    "N4",
    "N6",
    "N7",
    "N9",

    # Phosphorus / other
    "P",
]


DNA_ATOM_ALIASES = {
    # phosphate
    'P': 'P',
    'OP1': 'O_phosphate',
    'OP2': 'O_phosphate',
    'O1P': 'O_phosphate',
    'O2P': 'O_phosphate',

    'O5\'': 'O_backbone',
    'O3\'': 'O_backbone',

    # common terminal naming
    'H5T': 'H_backbone',
    'H3T': 'H_backbone',

    # sugar
    'C1\'': 'C_sugar1',
    'C2\'': 'C_sugar2',
    'C3\'': 'C_sugar3',
    'C4\'': 'C_sugar3',
    'C5\'': 'C_sugar4',

    'O4\'': 'O_sugar',

    # base carbons
    # 'C2': 'C_base',
    # 'C4': 'C_base',
    # 'C5': 'C_base',
    # 'C6': 'C_base',
    # 'C8': 'C_base',

    # thymine methyl
    # 'C7': 'C_base_methyl',
    'C5M': 'C_base_methyl',

    # base nitrogens
    'N1': 'N_base1',
    'N2': 'N_base2',
    # 'N3': 'N_base',
    'N4': 'N_base2',
    'N6': 'N_base2',
    # 'N7': 'N_base',
    # 'N9': 'N_base',

    # base oxygens
    'O2': 'O_base',
    'O4': 'O_base',
    'O6': 'O_base',

    # sugar hydrogens
    'H1\'': 'H_sugar',
    'H2\'': 'H_sugar',
    'H2\'\'': 'H_sugar',
    'H2\'1': 'H_sugar',
    'H2\'2': 'H_sugar',
    'H3\'': 'H_sugar',
    'H4\'': 'H_sugar',
    'H5\'': 'H_sugar',
    'H5\'\'': 'H_sugar',
    'H5\'1': 'H_sugar',
    'H5\'2': 'H_sugar',

    # base hydrogens
    # 'H1': 'H_base',
    # 'H2': 'H_base',
    # 'H3': 'H_base',
    # 'H5': 'H_base',
    # 'H6': 'H_base',
    # 'H8': 'H_base',

    # exocyclic amino hydrogens
    'H21': 'H_base_exocyclic',
    'H22': 'H_base_exocyclic',
    'H41': 'H_base_exocyclic',
    'H42': 'H_base_exocyclic',
    'H61': 'H_base_exocyclic',
    'H62': 'H_base_exocyclic',

    # thymine methyl hydrogens
    'H71': 'H_base_methyl',
    'H72': 'H_base_methyl',
    'H73': 'H_base_methyl',
    'HM1': 'H_base_methyl',
    'HM2': 'H_base_methyl',
    'HM3': 'H_base_methyl',

    # rare terminal hydroxyl hydrogens
    'HO5\'': 'H_backbone',
    'HO3\'': 'H_backbone',
}

DNA_ATOM_RESIDUE_ALIASES = {
    # Split N_base1 by residue class
    ('N1', 'DA'): 'N1_A',
    ('N1', 'A'): 'N1_A',

    ('N1', 'DG'): 'N1_G',
    ('N1', 'G'): 'N1_G',

    ('N1', 'DC'): 'N1_CT',
    ('N1', 'DT'): 'N1_CT',
    ('N1', 'C'): 'N1_CT',
    ('N1', 'T'): 'N1_CT',

    # ('C5', 'DA'): 'C5_A',
    # ('C5', 'DC'): 'C'
}

DNA_GROUP_ORDER = [
    # --- Carbons (backbone → sugar → base) ---
    "C_backbone",
    "C_sugar1",
    "C_sugar2",
    "C_sugar3",
    "C_sugar4",
    "C1'",
    "C2'",
    "C3'",
    "C4'",
    "C5'",
    "C2",
    "C4",
    "C5",
    "C6",
    "C7",   # thymine methyl
    "C8",

    # --- Oxygens ---
    "O_phosphate",
    "O_backbone",
    "O_sugar",
    "O_base",
    "O_base",
    "O_base",

    # --- Nitrogens ---
    "N_base1",
    "N1_A",
    "N1_G",
    "N1_CT",
    "N_base2",
    "N1",
    "N2",
    "N3",
    "N4",
    "N6",
    "N7",
    "N9",

    # --- Phosphorus ---
    "P",
]


DIRECT_LABEL_GROUPS: Set[str] = {
    'O_backbone',
    'O_carboxyl_1',
    'O_carboxyl_2',
}


FORCE_INCLUDE_GROUPS: Set[str] = set()
#     'O_backbone',
#     'O_carboxyl_1',
#     'O_carboxyl_2',
#     'CA',
#     'CB',
#     'CG',
#     'CD',
#     'CE',
#     'C',
#     'CZ',
#     'CE2',
#     'CD1',
#     'CD2',
#     'CE1',
# }


SMALL_MOL_FORCE_INCLUDE_GROUPS: Set[str] = set()


COV_RADII = {
    'H': 0.31,
    'C': 0.76,
    'N': 0.71,
    'O': 0.66,
}


def is_bonded(a, b):
    r1 = COV_RADII.get(a['element'], 0.7)
    r2 = COV_RADII.get(b['element'], 0.7)
    cutoff = 1.25 * (r1 + r2)

    dist = np.linalg.norm(a['coord'] - b['coord'])
    return dist < cutoff


def build_bonds(atoms):
    bonds = {i: [] for i in range(len(atoms))}

    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            if is_bonded(atoms[i], atoms[j]):
                bonds[i].append(j)
                bonds[j].append(i)

    return bonds


def classify_carbon(i, atoms, bonds):
    neighbors = bonds[i]
    neighbor_elements = [atoms[j]['element'] for j in neighbors]

    n = len(neighbors)

    if 'O' in neighbor_elements:
        return 'C_carbonyl'

    if n == 4:
        return 'C_sp3'

    if n == 3:
        return 'C_sp2'

    return 'C_other'


def classify_nitrogen(i, atoms, bonds):
    neighbors = bonds[i]
    neighbor_elements = [atoms[j]['element'] for j in neighbors]

    if neighbor_elements.count('H') >= 2:
        return 'N_amine'

    if 'C' in neighbor_elements and len(neighbors) == 3:
        return 'N_amide'

    return 'N_other'


def classify_oxygen(i, atoms, bonds):
    neighbors = bonds[i]
    neighbor_elements = [atoms[j]['element'] for j in neighbors]

    if 'H' in neighbor_elements:
        return 'O_hydroxyl'

    if neighbor_elements.count('C') == 1:
        return 'O_carbonyl'

    return 'O_other'


def classify_hydrogen(i, atoms, bonds):
    neighbors = bonds[i]

    if not neighbors:
        return 'H'

    parent = atoms[neighbors[0]]['element']

    if parent == 'O':
        return 'H_on_O'
    elif parent == 'N':
        return 'H_on_N'
    else:
        return 'H_on_C'


def get_force_include_groups(molecule_class: str) -> Set[str]:
    mol = str(molecule_class).strip().lower().replace(' ', '_')

    if mol == 'small_molecule':
        return SMALL_MOL_FORCE_INCLUDE_GROUPS

    return FORCE_INCLUDE_GROUPS


DNA_RESIDUES = {'DA', 'DC', 'DG', 'DT', 'DI'}
RNA_RESIDUES = {'RA', 'RC', 'RG', 'RU', 'A', 'C', 'G', 'U'}
SMALL_MOL_RESIDUES = {'EDTA', 'LIG', 'UNK', 'MOL'}  # extend as needed

DATA_ROOT = r'E:\Molecular'

OUTPUT_ROOT = (
    r'E:\OneDrive - Georgia State University\GSU NSC\Manuscripts'
    r'\Ericson Voronoi DNA\P2\fig2_atomic_level_scheme_deviations\2C_Full_Plots_Clustered'
)

MOLECULE_GROUPS = {
    'small molecule': ['B_EDTA', 'C_DB1976', r'I_T4LP\JZ4'],
    'dna': ['F_BDNA', 'D_Hairpin', 'K_NCP_DNA'],
    'rna': ['G_Hammerhead'],
    'protein': ['E_Cambrin', 'H_p53tet', 'I_T4LP', 'J_Streptavidin', 'L_BSA', 'm_NCP_Protein'],
    'complex': ['N_NCP'],
}

SETTINGS_COMBINATIONS = [
    # Manual baseline
    {
        'mode': 'manual',
        'use_ml_clustering': False,
        'ml_method': 'manual',
        'n_clusters': None,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },

    # KMeans
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kmeans',
        'n_clusters': 8,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kmeans',
        'n_clusters': 12,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kmeans',
        'n_clusters': 16,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kmeans',
        'n_clusters': 12,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': ['AtomName', 'ResidueName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },

    # Agglomerative
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'agglomerative',
        'n_clusters': 8,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'agglomerative',
        'n_clusters': 12,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'agglomerative',
        'n_clusters': 16,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'agglomerative',
        'n_clusters': 20,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': ['AtomName', 'ResidueName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },

    # KPrototypes
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kprototypes',
        'n_clusters': 8,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kprototypes',
        'n_clusters': 12,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kprototypes',
        'n_clusters': 16,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName', 'ResidueName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },

    # DBSCAN
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'dbscan',
        'n_clusters': None,
        'eps': 0.6,
        'min_samples': 8,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'dbscan',
        'n_clusters': None,
        'eps': 0.8,
        'min_samples': 10,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'dbscan',
        'n_clusters': None,
        'eps': 1.0,
        'min_samples': 12,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },

    # HDBSCAN
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'hdbscan',
        'n_clusters': None,
        'eps': None,
        'min_samples': 8,
        'min_cluster_size': 20,
        'numerical_cols': ['AW', 'Pow'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'hdbscan',
        'n_clusters': None,
        'eps': None,
        'min_samples': 10,
        'min_cluster_size': 25,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': [],
        'boolean_cols': [],
        'use_sol_binary': False,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'hdbscan',
        'n_clusters': None,
        'eps': None,
        'min_samples': 12,
        'min_cluster_size': 35,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': [],
        'use_sol_binary': False,
    },

    # A couple with solvent binary, if available later
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'agglomerative',
        'n_clusters': 12,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName'],
        'boolean_cols': ['SolFacingBinary'],
        'use_sol_binary': True,
    },
    {
        'mode': 'ml',
        'use_ml_clustering': True,
        'ml_method': 'kprototypes',
        'n_clusters': 12,
        'eps': None,
        'min_samples': None,
        'min_cluster_size': None,
        'numerical_cols': ['AW', 'Pow', 'DeltaV'],
        'categorical_cols': ['AtomName', 'ResidueName'],
        'boolean_cols': ['SolFacingBinary'],
        'use_sol_binary': True,
    },
]


def groups_are_similar(g1, g2):
    d_aw = abs(g1['Mean_AW'] - g2['Mean_AW'])
    d_pow = abs(g1['Mean_Pow'] - g2['Mean_Pow'])

    distance = np.sqrt(d_aw**2 + d_pow**2)

    spread_ok = max(g1['Spread'], g2['Spread']) < 2.0

    return distance < 0.6 and spread_ok


def can_merge(g1, g2):
    chem1 = get_chem_class(g1['GroupName'])
    chem2 = get_chem_class(g2['GroupName'])

    chem_ok = (chem1 == chem2)

    return chem_ok and groups_are_similar(g1, g2)


def resolve_group_folders(data_root: str, relative_folders: List[str]) -> List[str]:
    folders = []

    for rel_path in relative_folders:
        full_path = os.path.join(data_root, rel_path)
        if os.path.isdir(full_path):
            folders.append(full_path)
        else:
            print(f"WARNING: folder not found -> {full_path}")

    return folders


def find_pdb_file(folder: str) -> Optional[str]:
    for root, _, files in os.walk(folder):
        for file in files:
            if file.lower().endswith('.pdb'):
                return os.path.join(root, file)

    return None


def parse_pdb_atoms(pdb_path: str) -> List[Dict[str, object]]:
    atoms = []

    with open(pdb_path, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            name = line[12:16].strip().upper()
            residue = line[17:20].strip().upper()
            chain = line[21].strip()
            res_seq = line[22:26].strip()

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            element = line[76:78].strip().upper()
            if not element:
                element = ''.join([c for c in name if c.isalpha()])[0].upper()

            atoms.append({
                'PDBOrderIndex': len(atoms),
                'Name': name,
                'ResidueName': residue,
                'Chain': chain,
                'Residue Sequence': res_seq,
                'element': element,
                'coord': np.array([x, y, z], dtype=float),
            })

    return atoms


def classify_small_molecule_atom(pdb_i: int, pdb_atoms: List[Dict[str, object]], bonds: Dict[int, List[int]]) -> str:
    atom = pdb_atoms[pdb_i]
    element = str(atom['element']).upper()
    name = str(atom['Name']).upper()

    if element == 'H':
        return 'H'

    if element == 'C':
        base_class = classify_carbon(pdb_i, pdb_atoms, bonds)
        return f"{base_class}_C"

    if element == 'N':
        base_class = classify_nitrogen(pdb_i, pdb_atoms, bonds)
        return f"{base_class}_N"

    if element == 'O':
        base_class = classify_oxygen(pdb_i, pdb_atoms, bonds)
        return f"{base_class}_O"

    if element == 'SE':
        return 'SE'

    return element


def get_sol_environment_label(sol_pct: float) -> str:
    if pd.isna(sol_pct):
        return 'env_unknown'

    if sol_pct < 20.0:
        return 'buried'

    if sol_pct < 60.0:
        return 'partial'

    return 'exposed'


def canonicalize_atom_name(
    atom_name: str,
    molecule_class: str,
    residue_name: str = ''
) -> str:
    atom = str(atom_name).strip().upper()
    residue = str(residue_name).strip().upper()
    mol = str(molecule_class).strip().lower()

    if mol == 'protein':
        if (atom, residue) in PROTEIN_ATOM_RESIDUE_ALIASES:
            return PROTEIN_ATOM_RESIDUE_ALIASES[(atom, residue)]

        if (atom, '*') in PROTEIN_ATOM_RESIDUE_ALIASES:
            return PROTEIN_ATOM_RESIDUE_ALIASES[(atom, '*')]

        return PROTEIN_ATOM_ALIASES.get(atom, atom)

    if mol == 'dna':
        if (atom, residue) in DNA_ATOM_RESIDUE_ALIASES:
            return DNA_ATOM_RESIDUE_ALIASES[(atom, residue)]

        return DNA_ATOM_ALIASES.get(atom, atom)

    if mol == 'rna':
        if (atom, residue) in RNA_ATOM_RESIDUE_ALIASES:
            return RNA_ATOM_RESIDUE_ALIASES[(atom, residue)]

        return RNA_ATOM_ALIASES.get(atom, atom)

    if mol == 'small_molecule':
        if residue in SMALL_MOL_MANUAL_GROUPS:
            if atom in SMALL_MOL_MANUAL_GROUPS[residue]:
                return SMALL_MOL_MANUAL_GROUPS[residue][atom]

        return SMALL_MOL_ATOM_ALIASES.get(atom, atom)

    return atom


def apply_small_molecule_index_groups(
    canonical_name: str,
    folder: str,
    atom_index: int
) -> str:
    folder_name = os.path.basename(os.path.normpath(folder)).upper()

    # DB1976-specific index-based grouping
    if 'DB1976' in folder.upper() or folder_name == 'C_DB1976':
        return DB1976_INDEX_GROUPS.get(int(atom_index), canonical_name)

    return canonical_name


def refine_small_molecule_group_by_volume(
    canonical_name: str,
    atom_name: str,
    aw_v: float,
    pow_v: float
) -> str:
    atom = str(atom_name).strip().upper()
    group = str(canonical_name).strip()

    mean_v = 0.5 * (float(aw_v) + float(pow_v))

    # Hydrogens were fully overlapping when split by atom number.
    # Split them by their actual geometric/volume regime instead.
    if group == 'H':
        if mean_v < 7.0:
            return 'H_compact'

        if mean_v < 10.5:
            return 'H_standard'

        return 'H_expanded'

    return group


def tag_small_molecule_environmental_outliers(atom_df: pd.DataFrame) -> pd.DataFrame:
    df = atom_df.copy()

    if 'CanonicalName' not in df.columns:
        return df

    h_mask = df['CanonicalName'].astype(str).str.startswith('H_')

    if h_mask.sum() < 5:
        return df

    h_df = df[h_mask].copy()

    h_center = h_df[['AW', 'Pow']].median()
    h_mad = (h_df[['AW', 'Pow']] - h_center).abs().median()
    h_mad = h_mad.replace(0, np.nan)

    z_aw = 0.6745 * (h_df['AW'] - h_center['AW']) / h_mad['AW']
    z_pow = 0.6745 * (h_df['Pow'] - h_center['Pow']) / h_mad['Pow']

    robust_dist = np.sqrt(z_aw**2 + z_pow**2)

    outlier_indices = h_df.index[robust_dist > 4.0]

    df.loc[outlier_indices, 'CanonicalName'] = 'H_environmental_outlier'

    return df


def infer_molecule_class(atoms_df: pd.DataFrame) -> str:
    residue_names = set(
        str(x).strip().upper()
        for x in atoms_df['Residue'].dropna().unique()
    )

    dna_residues = {'DA', 'DC', 'DG', 'DT', 'DI'}
    rna_residues = {'A', 'C', 'G', 'U', 'RA', 'RC', 'RG', 'RU'}
    protein_markers = {'ALA', 'GLY', 'VAL', 'LEU', 'SER', 'THR', 'ASP', 'GLU', 'LYS', 'ARG'}

    if residue_names & protein_markers:
        return 'protein'

    if residue_names & dna_residues:
        return 'dna'

    if residue_names & rna_residues:
        return 'rna'

    return 'small_molecule'


def _read_scheme_logs(folder: str):
    """
    Reads aw/pow/prm logs from either:
      folder/{aw_logs.csv,pow_logs.csv,prm_logs.csv}
    or
      folder/{aw/aw_logs.csv, pow/pow_logs.csv, prm/prm_logs.csv}
    """
    try:
        aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True, surfs=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True, surfs=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True, surfs=True)

    except FileNotFoundError:
        aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True, surfs=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True, surfs=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True, surfs=True)

    return aw_logs, pow_logs, prm_logs


def build_residue_order_lookup(pdb_atoms):
    lookup = defaultdict(list)

    for i, atom in enumerate(pdb_atoms):
        key = (atom['ResidueName'], atom['Chain'], atom['Residue Sequence'])
        lookup[key].append(i)

    return lookup


def build_pdb_name_lookup(pdb_atoms: List[Dict[str, object]]) -> Dict[str, int]:
    lookup = {}

    for i, atom in enumerate(pdb_atoms):
        name = str(atom['Name']).strip().upper()
        lookup[name] = i

    return lookup


def select_folders_multi(title: str = "Select a folder (Cancel to finish)") -> List[str]:
    root = tk.Tk()
    root.withdraw()

    folders = []

    while True:
        folder = filedialog.askdirectory(title=title)
        if not folder:
            break

        print(f"Loaded: {os.path.basename(folder)}")
        folders.append(folder)

    root.destroy()
    return folders


def select_output_folder(title: str = "Select output folder for saved plot data") -> Optional[str]:
    root = tk.Tk()
    root.withdraw()
    try:
        root.attributes('-topmost', True)
    except Exception:
        pass
    root.update_idletasks()
    root.lift()
    folder = filedialog.askdirectory(parent=root, title=title)
    root.destroy()

    if not folder:
        return None

    return folder


def attach_sol_facing_pct(atoms_df: pd.DataFrame, surfs_df: pd.DataFrame, out_col: str) -> pd.DataFrame:
    atoms_df = atoms_df.copy()

    if surfs_df is None or len(surfs_df) == 0:
        atoms_df[out_col] = np.nan
        return atoms_df

    surfs_df = surfs_df.copy()

    # Expand list-valued Balls column into Ball 1 / Ball 2
    if 'Balls' in surfs_df.columns:
        surfs_df['Ball 1'] = surfs_df['Balls'].apply(
            lambda x: x[0] if isinstance(x, (list, tuple)) and len(x) > 0 else np.nan
        )
        surfs_df['Ball 2'] = surfs_df['Balls'].apply(
            lambda x: x[1] if isinstance(x, (list, tuple)) and len(x) > 1 else np.nan
        )

    # Normalize area column name if needed
    if 'Surface Area' in surfs_df.columns:
        pass
    elif 'Area' in surfs_df.columns:
        surfs_df = surfs_df.rename(columns={'Area': 'Surface Area'})

    sol_df = compute_sol_facing_percent(atoms_df.copy(), surfs_df.copy())

    if not isinstance(sol_df, pd.DataFrame):
        atoms_df[out_col] = np.nan
        return atoms_df

    if 'Index' not in sol_df.columns or 'SolFacingPct' not in sol_df.columns:
        atoms_df[out_col] = np.nan
        return atoms_df

    sol_df = sol_df[['Index', 'SolFacingPct']].copy()
    sol_df = sol_df.rename(columns={'SolFacingPct': out_col})

    atoms_df = atoms_df.merge(sol_df, on='Index', how='left')

    return atoms_df


def collect_atom_volume_points(
    folders: List[str],
    atom_name_field: str = 'Name',
    molecule_class: str = 'protein',
    volume_range: Optional[tuple] = None
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []

    mol = str(molecule_class).strip().lower().replace(' ', '_')

    def get_first_present(row, candidates, default=np.nan):
        for col in candidates:
            if col in row.index:
                return row[col]
        return default

    for folder in folders:
        aw_logs, pow_logs, _ = _read_scheme_logs(folder)

        aw_atoms = aw_logs['atoms'].copy()
        pow_atoms = pow_logs['atoms'].copy()

        aw_surfs = aw_logs['surfs'] if 'surfs' in aw_logs else None
        pow_surfs = pow_logs['surfs'] if 'surfs' in pow_logs else None

        aw_atoms = attach_sol_facing_pct(aw_atoms, aw_surfs, 'AWSolFacingPct')
        pow_atoms = attach_sol_facing_pct(pow_atoms, pow_surfs, 'PowSolFacingPct')

        pdb_atoms = None
        pdb_bonds = None
        pdb_by_order = None

        if mol == 'small_molecule':
            pdb_path = find_pdb_file(folder)

            if pdb_path is not None:
                pdb_atoms = parse_pdb_atoms(pdb_path)
                pdb_bonds = build_bonds(pdb_atoms)
                pdb_by_order = list(range(len(pdb_atoms)))

                print(f"\nLoaded small-molecule PDB for classification: {pdb_path}")
                print(f"PDB atoms found: {len(pdb_atoms)}")
                print(f"AW atoms found: {len(aw_atoms)}")
            else:
                print(f"\nWARNING: no PDB found in {folder}; using fallback canonicalization.")

        pow_lookup = {
            int(row['Index']): row
            for _, row in pow_atoms.iterrows()
        }

        small_mol_counter = 0

        for _, atom in aw_atoms.iterrows():
            idx = int(atom['Index'])

            if idx not in pow_lookup:
                continue

            pow_atom = pow_lookup[idx]

            aw_v = float(atom['Volume'])
            pow_v = float(pow_atom['Volume'])

            if volume_range is not None:
                vmin, vmax = volume_range
                if aw_v < vmin or aw_v > vmax or pow_v < vmin or pow_v > vmax:
                    continue

            atom_name = str(atom.get(atom_name_field, '')).strip().upper()

            # if molecule_class == 'small_molecule' and atom_name.upper().startswith('H'):
            #     continue

            if not atom_name:
                continue

            residue_name = str(
                atom.get('Residue', atom.get('Residue Name', atom.get('ResName', '')))
            ).strip().upper()

            if should_exclude_from_grouping(residue_name):
                continue

            chain = get_first_present(atom, ['Chain', 'chain'], default='')
            res_seq = get_first_present(atom, ['Residue Sequence', 'ResSeq', 'res_seq'], default='')

            canonical_name = canonicalize_atom_name(
                atom_name=atom_name,
                molecule_class=molecule_class,
                residue_name=residue_name
            )
            if mol == 'small_molecule':
                canonical_name = apply_small_molecule_index_groups(
                    canonical_name=canonical_name,
                    folder=folder,
                    atom_index=idx
                )

            coord_x = get_first_present(atom, ['x', 'X', 'coord_x', 'CoordX'])
            coord_y = get_first_present(atom, ['y', 'Y', 'coord_y', 'CoordY'])
            coord_z = get_first_present(atom, ['z', 'Z', 'coord_z', 'CoordZ'])

            aw_sol_pct = (
                float(atom['AWSolFacingPct'])
                if 'AWSolFacingPct' in atom.index and pd.notna(atom['AWSolFacingPct'])
                else np.nan
            )

            pow_sol_pct = (
                float(pow_atom['PowSolFacingPct'])
                if 'PowSolFacingPct' in pow_atom.index and pd.notna(pow_atom['PowSolFacingPct'])
                else np.nan
            )

            delta_sol_pct = (
                pow_sol_pct - aw_sol_pct
                if pd.notna(aw_sol_pct) and pd.notna(pow_sol_pct)
                else np.nan
            )

            # if mol == 'small_molecule':
            #     sol_env = get_sol_environment_label(aw_sol_pct)
            #     canonical_name = f"{canonical_name}_{sol_env}"

            records.append({
                'Folder': folder,
                'Index': idx,
                'Chain': chain,
                'ResidueName': residue_name,
                'Residue Sequence': res_seq,
                'AtomName': atom_name,
                'CanonicalName': canonical_name,
                'AW': aw_v,
                'Pow': pow_v,
                'DeltaV': pow_v - aw_v,
                'CoordX': float(coord_x) if pd.notna(coord_x) else np.nan,
                'CoordY': float(coord_y) if pd.notna(coord_y) else np.nan,
                'CoordZ': float(coord_z) if pd.notna(coord_z) else np.nan,
                'AWSolFacingPct': aw_sol_pct,
                'PowSolFacingPct': pow_sol_pct,
                'DeltaSolFacingPct': delta_sol_pct,
            })

    return pd.DataFrame(records)


def analyze_group_outliers(df, min_group_size=3, z_thresh=2.5, top_k=20):
    """
    Identify:
    - Small / ungrouped clusters
    - Within-group outliers
    - Global outliers

    Returns:
        dict of DataFrames
    """

    results = {}

    # -----------------------
    # 1. Group stats
    # -----------------------
    group_stats = df.groupby('CanonicalName').agg(
        n=('AW', 'count'),
        aw_mean=('AW', 'mean'),
        pow_mean=('Pow', 'mean'),
        aw_std=('AW', 'std'),
        pow_std=('Pow', 'std')
    ).reset_index()

    results['group_stats'] = group_stats

    # -----------------------
    # 2. Small / ungrouped groups
    # -----------------------
    small_groups = group_stats[group_stats['n'] <= min_group_size].copy()
    results['small_groups'] = small_groups.sort_values('n')

    # -----------------------
    # 3. Within-group z-score outliers
    # -----------------------
    df = df.copy()

    df = df.merge(group_stats, on='CanonicalName', how='left')

    df['z_aw'] = (df['AW'] - df['aw_mean']) / (df['aw_std'].replace(0, np.nan))
    df['z_pow'] = (df['Pow'] - df['pow_mean']) / (df['pow_std'].replace(0, np.nan))

    df['z_total'] = np.sqrt(df['z_aw']**2 + df['z_pow']**2)

    within_outliers = df[df['z_total'] > z_thresh].copy()
    results['within_group_outliers'] = within_outliers.sort_values('z_total', ascending=False)

    global_mean = df[['AW', 'Pow']].mean().values

    df['global_dist'] = np.sqrt(
        (df['AW'] - global_mean[0])**2 +
        (df['Pow'] - global_mean[1])**2
    )

    global_outliers = df.sort_values('global_dist', ascending=False).head(top_k)
    results['global_outliers'] = global_outliers

    results['annotated_df'] = df

    return results


def find_ungrouped_atoms(
    atom_df: pd.DataFrame,
    aw_range: Optional[tuple] = None,
    pow_range: Optional[tuple] = None,
    top_k: int = 50
):
    """
    Find atoms whose canonical name is unchanged, meaning they were not grouped
    by the alias dictionary.

    Optional AW / Pow window lets you inspect a suspicious patch.
    """
    df = atom_df.copy()

    # ungrouped = canonicalization did nothing
    df = df[df['CanonicalName'] == df['AtomName']].copy()

    if aw_range is not None:
        df = df[(df['AW'] >= aw_range[0]) & (df['AW'] <= aw_range[1])]

    if pow_range is not None:
        df = df[(df['Pow'] >= pow_range[0]) & (df['Pow'] <= pow_range[1])]

    if len(df) == 0:
        print("\nNo ungrouped atoms found in this selection.")
        return df, pd.DataFrame(), pd.DataFrame()

    summary_by_atom = (
        df.groupby('AtomName')
        .agg(
            Count=('AtomName', 'size'),
            Mean_AW=('AW', 'mean'),
            Mean_Pow=('Pow', 'mean'),
            Residues=('ResidueName', lambda x: ', '.join(sorted(set(map(str, x))))),
        )
        .reset_index()
        .sort_values(['Count', 'Mean_AW', 'Mean_Pow'], ascending=[False, True, True])
    )

    summary_by_atom_residue = (
        df.groupby(['AtomName', 'ResidueName'])
        .agg(
            Count=('AtomName', 'size'),
            Mean_AW=('AW', 'mean'),
            Mean_Pow=('Pow', 'mean'),
        )
        .reset_index()
        .sort_values(['Count', 'Mean_AW', 'Mean_Pow'], ascending=[False, True, True])
    )

    print("\n=== UNGROUPED ATOMS: BY ATOM NAME ===")
    print(summary_by_atom.head(top_k).to_string(index=False))

    print("\n=== UNGROUPED ATOMS: BY ATOM NAME + RESIDUE ===")
    print(summary_by_atom_residue.head(top_k).to_string(index=False))

    return df, summary_by_atom, summary_by_atom_residue


def extract_patch_atoms(
    atom_df: pd.DataFrame,
    aw_range: tuple,
    pow_range: tuple,
    min_count: int = 1
):
    """
    Extract atoms in a specific AW/Pow region and print
    PyMOL-ready index selections.
    """

    df = atom_df.copy()

    # filter by region
    df = df[
        (df['AW'] >= aw_range[0]) & (df['AW'] <= aw_range[1]) &
        (df['Pow'] >= pow_range[0]) & (df['Pow'] <= pow_range[1])
    ].copy()

    if len(df) == 0:
        print("No atoms found in patch.")
        return df

    # summary
    summary = (
        df.groupby(['AtomName', 'ResidueName'])
        .size()
        .reset_index(name='Count')
        .sort_values('Count', ascending=False)
    )

    summary = summary[summary['Count'] >= min_count]

    print("\n=== PATCH SUMMARY ===")
    print(summary.to_string(index=False))

    # print raw atoms
    print("\n=== PATCH ATOMS ===")
    cols = ['Index', 'AtomName', 'ResidueName', 'Chain', 'Residue Sequence', 'AW', 'Pow']
    cols = [c for c in cols if c in df.columns]

    print(df[cols].sort_values(['ResidueName', 'Residue Sequence']).to_string(index=False))

    # build PyMOL selections
    print("\n=== PYMOL SELECTIONS ===")

    # 1. Precise atom-level selection using chain / resi / name
    required_cols = {'AtomName', 'ResidueName', 'Index'}
    has_chain = 'Chain' in df.columns
    has_resi = 'Residue Sequence' in df.columns

    if has_chain and has_resi and 'AtomName' in df.columns:
        selection_terms = []

        for _, row in df.sort_values(['Chain', 'Residue Sequence', 'AtomName']).iterrows():
            chain = str(row['Chain']).strip()
            resi = str(row['Residue Sequence']).strip()
            atom_name = str(row['AtomName']).strip()

            # PyMOL atom selection: chain X and resi 12 and name C4
            term = f"(chain {chain} and resi {resi} and name {atom_name})"
            selection_terms.append(term)

        selection_str = " or ".join(selection_terms)

        print("\nPyMOL precise selection:")
        print(f"select patch_atoms, ({selection_str})")

    # 2. Optional debug print using your internal Index shifted to PyMOL-style +1
    elif 'Index' in df.columns:
        pymol_indices = sorted((df['Index'] + 1).unique())
        index_str = " or ".join([f"index {i}" for i in pymol_indices])

        print("\nPyMOL index-only selection (+1 corrected):")
        print(f"select patch_atoms_idx, ({index_str})")

    return df


def compute_name_volume_stats(df: pd.DataFrame, group_col: str = 'CanonicalName') -> pd.DataFrame:
    if len(df) == 0:
        return pd.DataFrame(
            columns=[
                'GroupName', 'Count', 'Mean_AW', 'Mean_Pow',
                'SD_AW', 'SD_Pow', 'Var_AW', 'Var_Pow', 'Spread', 'Members'
            ]
        )

    grouped_rows = []

    for group_name, subdf in df.groupby(group_col):
        member_names = sorted(subdf['AtomName'].unique().tolist())

        var_aw = subdf['AW'].var(ddof=1) if len(subdf) > 1 else 0.0
        var_pow = subdf['Pow'].var(ddof=1) if len(subdf) > 1 else 0.0

        grouped_rows.append({
            'GroupName': group_name,
            'Count': len(subdf),
            'Mean_AW': subdf['AW'].mean(),
            'Mean_Pow': subdf['Pow'].mean(),
            'SD_AW': subdf['AW'].std(ddof=1) if len(subdf) > 1 else 0.0,
            'SD_Pow': subdf['Pow'].std(ddof=1) if len(subdf) > 1 else 0.0,
            'Var_AW': var_aw,
            'Var_Pow': var_pow,
            'Spread': np.sqrt(var_aw + var_pow),
            'Members': ', '.join(member_names),
        })

    stats_df = pd.DataFrame(grouped_rows)

    stats_df['SD_AW'] = stats_df['SD_AW'].fillna(0.0)
    stats_df['SD_Pow'] = stats_df['SD_Pow'].fillna(0.0)
    stats_df['Spread'] = stats_df['Spread'].fillna(0.0)

    stats_df = stats_df.sort_values(['Mean_AW', 'Mean_Pow']).reset_index(drop=True)

    return stats_df


def filter_plot_groups(
    stats_df: pd.DataFrame,
    plot_min_count: int = 50,
    max_spread: Optional[float] = None,
    force_include: Optional[Set[str]] = None
) -> pd.DataFrame:
    if force_include is None:
        force_include = set()

    mask = stats_df['Count'] >= plot_min_count

    if max_spread is not None:
        mask &= stats_df['Spread'] <= max_spread

    mask |= stats_df['GroupName'].isin(force_include)

    plot_stats_df = stats_df[mask].copy()
    plot_stats_df = plot_stats_df.sort_values(['Mean_AW', 'Mean_Pow']).reset_index(drop=True)
    plot_stats_df['PlotNumber'] = np.arange(1, len(plot_stats_df) + 1)

    return plot_stats_df


def is_hydrogen_group(group_name: str) -> bool:
    name = str(group_name).upper()
    return name == 'H' or name.startswith('H_')


EXCLUDE_RESIDUES_FROM_GROUPING = {'GDP', 'OMC'}


def should_exclude_from_grouping(residue_name: str) -> bool:
    return str(residue_name).strip().upper() in EXCLUDE_RESIDUES_FROM_GROUPING


def build_merged_group_map(
    stats_df: pd.DataFrame,
    distance_cutoff: float = 0.60,
    spread_cutoff: float = 2.50,
    min_count_to_merge: int = 20
) -> Dict[str, str]:
    """
    Build a conservative merge map for protein groups.

    Groups are merged only when:
    1. They share the same chemistry class.
    2. Their AW/Pow centers are close.
    3. Their spreads are not too large.
    4. At least one group is small enough that merging is useful.
    """
    merge_map = {}

    stats_df = stats_df.copy()
    stats_df['ChemClass'] = stats_df['GroupName'].apply(get_chem_class)

    for chem_class, chem_df in stats_df.groupby('ChemClass'):
        chem_df = chem_df.sort_values(['Mean_AW', 'Mean_Pow']).reset_index(drop=True)

        used = set()

        for i, row_i in chem_df.iterrows():
            group_i = row_i['GroupName']

            if group_i in used:
                continue

            merge_members = [group_i]
            used.add(group_i)

            for j, row_j in chem_df.iterrows():
                group_j = row_j['GroupName']

                if group_j in used:
                    continue

                d_aw = row_i['Mean_AW'] - row_j['Mean_AW']
                d_pow = row_i['Mean_Pow'] - row_j['Mean_Pow']
                distance = float(np.sqrt(d_aw ** 2 + d_pow ** 2))

                count_ok = (
                    row_i['Count'] < min_count_to_merge and
                    row_j['Count'] < min_count_to_merge
                )

                spread_ok = (
                    row_i['Spread'] <= spread_cutoff and
                    row_j['Spread'] <= spread_cutoff
                )

                if distance <= distance_cutoff and spread_ok and count_ok:
                    merge_members.append(group_j)
                    used.add(group_j)

            if len(merge_members) == 1:
                merge_name = group_i
            else:
                merge_name = f"{chem_class}_merged"

            for member in merge_members:
                merge_map[member] = merge_name

    return merge_map


def print_name_volume_stats(stats_df: pd.DataFrame):
    print("\nGrouped atom-name centers and standard deviations in AW vs Pow volume space:\n")

    for _, row in stats_df.iterrows():
        print(
            f"{str(row['GroupName']):>14s} | "
            f"AW={row['Mean_AW']:>6.3f} ± {row['SD_AW']:.3f} | "
            f"Pow={row['Mean_Pow']:>6.3f} ± {row['SD_Pow']:.3f} | "
            f"Spread={row['Spread']:.3f} | "
            f"n={int(row['Count']):>5d} | "
            f"Members: {row['Members']}"
        )


def print_group_plot_table(plot_stats_df: pd.DataFrame):
    print("\nPlotted groups:\n")

    for _, row in plot_stats_df.iterrows():
        print(
            f"{int(row['PlotNumber']):>2d} | "
            f"{str(row['GroupName']):<18s} | "
            f"AW={row['Mean_AW']:>6.3f}, "
            f"Pow={row['Mean_Pow']:>6.3f} | "
            f"Spread={row['Spread']:.3f} | "
            f"n={int(row['Count']):>5d} | "
            f"Members: {row['Members']}"
        )


def save_group_plot_table(plot_stats_df: pd.DataFrame, out_path: str):
    cols = ['PlotNumber', 'GroupName', 'Mean_AW', 'Mean_Pow', 'Spread', 'Count', 'Members']
    plot_stats_df[cols].to_csv(out_path, index=False)
    print(f"\nSaved plotted group table -> {out_path}")


def get_plot_color(name: str) -> str:
    name = str(name).upper()
    base = name.split('_')[0]

    if base.startswith('H'):
        return '#1f77b4'

    if base.startswith('C'):
        return '#ff7f0e'

    if base.startswith('N'):
        return '#2ca02c'

    if base.startswith('O'):
        return '#d62728'

    if base.startswith('P'):
        return '#9467bd'

    if base.startswith('S') or base.startswith('SE'):
        return '#8c564b'

    return 'gray'


def add_covariance_ellipse(
    ax,
    x,
    y,
    color,
    n_std: float = 1.5,
    face_alpha: float = 0.10,
    edge_alpha: float = 0.75,
    linewidth: float = 2.0,
    zorder: float = 2
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) < 3:
        return

    cov = np.cov(x, y)

    if np.any(~np.isfinite(cov)):
        return

    vals, vecs = np.linalg.eigh(cov)

    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    if np.any(vals < 0):
        return

    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width = 2.0 * n_std * np.sqrt(vals[0])
    height = 2.0 * n_std * np.sqrt(vals[1])

    fill = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        edgecolor=color,
        linewidth=linewidth,
        alpha=face_alpha,
        zorder=zorder,
    )
    ax.add_patch(fill)

    edge = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=angle,
        facecolor='none',
        edgecolor=color,
        linewidth=linewidth,
        alpha=edge_alpha,
        zorder=zorder + 0.1
    )

    edge.set_path_effects([
        pe.Stroke(linewidth=linewidth + 2.5, foreground='white'),
        pe.Normal()
    ])

    ax.add_patch(edge)


def add_outside_group_list(fig, plot_stats_df: pd.DataFrame):
    x_num = 0.80
    x_name = 0.835
    y_start = 0.90
    y_step = 0.024

    for i, (_, row) in enumerate(plot_stats_df.iterrows()):
        y = y_start - i * y_step
        name = row['GroupName']
        if str(name).isdigit():
            color = plt.cm.tab20(int(name) % 20)
        else:
            color = get_plot_color(str(name))

        if y < 0.06:
            break

        fig.text(
            x_num,
            y,
            f"{int(row['PlotNumber']):>2d}",
            color=color,
            fontsize=10,
            fontweight='bold',
            ha='left',
            va='center',
            path_effects=[
                pe.Stroke(linewidth=2.5, foreground='white'),
                pe.Normal()
            ]
        )

        fig.text(
            x_name,
            y,
            row['GroupName'],
            color='black',
            fontsize=10,
            ha='left',
            va='center'
        )
        # --- Add hydrogen entry ---
        # --- Add hydrogen entry ---
        fig.text(
            0.80,  # adjust if needed
            0.10,  # place at bottom of legend block
            u"\u25CF  Hydrogens",
            fontsize=16,
            color='#1f77b4'
        )


def plot_ml_clusters_aw_pow(
    df: pd.DataFrame,
    title: str,
    save: Optional[str] = None,
    show: bool = True
):
    fig, ax = plt.subplots(figsize=(10, 8))



    cluster_ids = sorted(df['MLCluster'].dropna().unique().tolist())

    for cluster_id in cluster_ids:
        subdf = df[df['MLCluster'] == cluster_id]

        label = f"Cluster {cluster_id}"
        if cluster_id == -1:
            label = "Noise"

        ax.scatter(
            subdf['AW'],
            subdf['Pow'],
            s=16,
            alpha=0.35,
            label=label
        )

    xmin = min(df['AW'].min(), df['Pow'].min())
    xmax = max(df['AW'].max(), df['Pow'].max())
    pad = 0.05 * (xmax - xmin if xmax > xmin else 1.0)

    ax.plot(
        [xmin - pad, xmax + pad],
        [xmin - pad, xmax + pad],
        linestyle='--',
        linewidth=2.5,
        color='black'
    )

    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(xmin - pad, xmax + pad)

    ax.set_xlabel('AW Volume', fontsize=22)
    ax.set_ylabel('Pow Volume', fontsize=22)
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2.0, length=8)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    ax.set_aspect('equal', adjustable='box')
    # ax.legend(fontsize=10, loc='best')

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')

    if show:
        plt.show()

    plt.close(fig)


def make_settings_name(
    mode: str,
    ml_method: str,
    n_clusters: Optional[int],
    eps: Optional[float],
    min_samples: Optional[int],
    min_cluster_size: Optional[int],
    numeric_cols: List[str],
    categorical_cols: List[str],
    use_sol_binary: bool,
    point_alpha: float,
    ellipse_n_std: float
) -> str:
    num_tag = '-'.join(numeric_cols) if numeric_cols else 'none'
    cat_tag = '-'.join(categorical_cols) if categorical_cols else 'none'

    parts = [
        mode,
        ml_method,
        f"k{n_clusters}" if n_clusters is not None else None,
        f"eps{eps}" if eps is not None else None,
        f"minsamp{min_samples}" if min_samples is not None else None,
        f"minclust{min_cluster_size}" if min_cluster_size is not None else None,
        f"num_{num_tag}",
        f"cat_{cat_tag}",
        f"solbin_{int(use_sol_binary)}",
        f"alpha_{point_alpha}",
        f"ell_{ellipse_n_std}",
    ]

    safe = "_".join([str(p) for p in parts if p is not None])
    safe = safe.replace("'", "").replace(" ", "")
    return safe


def plot_name_volume_groups(
    atom_df: pd.DataFrame,
    plot_stats_df: pd.DataFrame,
    title: str,
    group_col: str = 'CanonicalName',
    volume_range: Optional[tuple] = None,
    show_points: bool = True,
    show_numbers: bool = True,
    annotate_direct_groups: bool = True,
    point_size: float = 18,
    point_alpha: float = 0.10,
    number_fontsize: int = 14,
    direct_label_fontsize: int = 11,
    show_ellipses: bool = True,
    ellipse_alpha: float = 0.5,
    ellipse_n_std: float = 1.5,
    ellipse_min_count: int = 50,
    ellipse_max_spread: Optional[float] = None,
    save_png: Optional[str] = None,
    save_svg: Optional[str] = None,
    downsample_fraction: float = 1.0,
    show: bool = True,
    group_hydrogens: bool = False
):
    fig, ax = plt.subplots(figsize=(12, 9))

    # Plot hydrogens as background only when they are not treated as groups.
    if not group_hydrogens and group_col in atom_df.columns:
        h_df = atom_df[
            atom_df[group_col].astype(str).apply(is_hydrogen_group)
        ].copy()

        if len(h_df) > 0:
            if downsample_fraction < 1.0:
                h_df = h_df.sample(frac=downsample_fraction, random_state=42)

            ax.scatter(
                h_df['AW'],
                h_df['Pow'],
                s=point_size,
                alpha=point_alpha,
                color='#1f77b4',
                edgecolors='none',
                zorder=2
            )

    print_group_plot_table(plot_stats_df)

    for _, row in plot_stats_df.iterrows():
        name = row['GroupName']
        if is_hydrogen_group(name) and not group_hydrogens:
            continue
        if group_col == 'MLCluster':
            color = plt.cm.tab20(int(name) % 20)
        else:
            color = get_plot_color(str(name))

        # --- DOWNSAMPLE ---
        full_group_df = atom_df[atom_df[group_col].astype(str) == str(name)]

        # --- DOWNSAMPLE (scatter only) ---
        plot_group_df = full_group_df

        if downsample_fraction < 1.0 and len(full_group_df) > 0:
            plot_group_df = full_group_df.sample(
                frac=downsample_fraction,
                random_state=42
            )

        print(f"\nGroupName from plot_stats_df: {name}")
        print(f"Matching rows in atom_df['CanonicalName'] == {name}: {len(full_group_df)}")
        if len(full_group_df) == 0:
            print("WARNING: empty group_df")
            print("Sample CanonicalName values:")
            print(atom_df['CanonicalName'].dropna().astype(str).unique()[:20])

        # --- Draw connecting line for 2-point groups ---
        if len(full_group_df) == 2:
            xs = full_group_df['AW'].values
            ys = full_group_df['Pow'].values

            ax.plot(
                xs,
                ys,
                linewidth=2.5,
                color=get_plot_color(name),
                alpha=0.9,
                zorder=0  # behind points
            )

        if show_points:
            ax.scatter(
                plot_group_df['AW'],
                plot_group_df['Pow'],
                s=point_size,
                alpha=point_alpha,
                color=color,
                zorder=1
            )

        can_draw_ellipse = (
                show_ellipses and
                (group_hydrogens or not is_hydrogen_group(name)) and
                len(full_group_df) >= ellipse_min_count and
                (ellipse_max_spread is None or row['Spread'] <= ellipse_max_spread)
        )

        if can_draw_ellipse:
            add_covariance_ellipse(
                ax=ax,
                x=full_group_df['AW'].to_numpy(),
                y=full_group_df['Pow'].to_numpy(),
                color=color,
                n_std=ellipse_n_std,
                face_alpha=0.10,
                edge_alpha=ellipse_alpha,
                linewidth=1.0,
                zorder=3
            )

        if show_numbers:
            ax.text(
                row['Mean_AW'],
                row['Mean_Pow'],
                str(int(row['PlotNumber'])),
                color=color,
                fontsize=number_fontsize,
                fontweight='bold',
                ha='center',
                va='center',
                zorder=6,
                path_effects=[
                    pe.Stroke(linewidth=3.0, foreground='white'),
                    pe.Normal()
                ]
            )

        if annotate_direct_groups and row['GroupName'] in DIRECT_LABEL_GROUPS:
            ax.text(
                row['Mean_AW'] + 0.18,
                row['Mean_Pow'] + 0.18,
                row['GroupName'],
                fontsize=direct_label_fontsize,
                color='black',
                ha='left',
                va='bottom',
                zorder=6
            )

    if volume_range is not None:
        vmin, vmax = volume_range

        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)

        ax.plot(
            [vmin, vmax],
            [vmin, vmax],
            linestyle='--',
            linewidth=3.5,
            color='black',
            alpha=0.9,
            zorder=0
        )

    else:
        xmin = min(atom_df['AW'].min(), atom_df['Pow'].min())
        xmax = max(atom_df['AW'].max(), atom_df['Pow'].max())
        pad = 0.05 * (xmax - xmin if xmax > xmin else 1.0)

        ax.set_xlim(xmin - pad, xmax + pad)
        ax.set_ylim(xmin - pad, xmax + pad)

        ax.plot(
            [xmin - pad, xmax + pad],
            [xmin - pad, xmax + pad],
            linestyle='--',
            linewidth=3.5,
            color='black',
            alpha=0.9,
            zorder=0
        )

    ax.set_xlabel('AW Volume', fontsize=24)
    ax.set_ylabel('Pow Volume', fontsize=24)
    ax.set_title(title, fontsize=22)

    ax.tick_params(axis='both', which='major', labelsize=20, width=2.5, length=10)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    ax.set_aspect('equal', adjustable='box')

    fig.subplots_adjust(left=0.12, right=0.77, bottom=0.12, top=0.90)

    add_outside_group_list(fig, plot_stats_df)

    if save_png is not None:
        plt.savefig(save_png, dpi=300, bbox_inches='tight')

    if save_svg is not None:
        plt.savefig(save_svg, bbox_inches='tight')

    if show:
        plt.show()

    plt.close(fig)


def run_batch_cluster_plots():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    for molecule_label, rel_paths in MOLECULE_GROUPS.items():
        folders = resolve_group_folders(DATA_ROOT, rel_paths)

        if len(folders) == 0:
            print(f"Skipping {molecule_label}: no folders found.")
            continue

        print(f"\n=== RUNNING GROUP: {molecule_label} ===")
        print(f"Folders: {folders}")

        manual_dir = os.path.join(OUTPUT_ROOT, molecule_label, 'manual')
        os.makedirs(manual_dir, exist_ok=True)

        for cfg in SETTINGS_COMBINATIONS:
            mode = cfg['mode']
            method = cfg['ml_method']

            out_subdir = os.path.join(OUTPUT_ROOT, molecule_label, method)
            os.makedirs(out_subdir, exist_ok=True)

            settings_name = make_settings_name(
                mode=mode,
                ml_method=method,
                n_clusters=cfg['n_clusters'],
                eps=cfg['eps'],
                min_samples=cfg['min_samples'],
                min_cluster_size=cfg['min_cluster_size'],
                numeric_cols=cfg['numerical_cols'],
                categorical_cols=cfg['categorical_cols'],
                use_sol_binary=cfg['use_sol_binary'],
                point_alpha=0.2,
                ellipse_n_std=1.2
            )

            png_path = os.path.join(out_subdir, f"{settings_name}.png")
            svg_path = os.path.join(out_subdir, f"{settings_name}.svg")

            print(f"\nRunning: {molecule_label} | {settings_name}")

            try:
                main(
                    folders=folders,
                    atom_name_field='Name',
                    molecule_class='small_molecule' if molecule_label == 'small molecule' else (
                        'dna' if molecule_label == 'dna' else (
                            'rna' if molecule_label == 'rna' else (
                                'protein' if molecule_label == 'protein' else 'protein'
                            )
                        )
                    ),
                    volume_range=(3, 22),
                    save_csv=False,
                    save_plot=True,
                    show_points=True,
                    show_numbers=True,
                    annotate_direct_groups=False,
                    plot_min_count=50,
                    max_spread=3.5,
                    ellipse_min_count=50,
                    ellipse_max_spread=None,
                    ellipse_n_std=1.2,
                    point_alpha=0.2,
                    use_ml_clustering=cfg['use_ml_clustering'],
                    ml_method=method if method != 'manual' else 'kmeans',
                    use_sol_binary=cfg['use_sol_binary'],
                    n_clusters=cfg['n_clusters'] if cfg['n_clusters'] is not None else 12,
                    min_samples=cfg['min_samples'] if cfg['min_samples'] is not None else 10,
                    eps=cfg['eps'] if cfg['eps'] is not None else 1.0,
                    min_cluster_size=cfg['min_cluster_size'] if cfg['min_cluster_size'] is not None else 25,
                    numerical_cols=cfg['numerical_cols'],
                    categorical_cols=cfg['categorical_cols'],
                    boolean_cols=cfg['boolean_cols'],
                    output_base=out_subdir,
                    output_name=settings_name,
                )
            except Exception as e:
                print(f"\nFAILED: {molecule_label} | {settings_name}")
                print(f"Reason: {e}")

                fail_log = os.path.join(out_subdir, "_failed_runs.txt")
                with open(fail_log, "a", encoding="utf-8") as f:
                    f.write(f"{molecule_label} | {settings_name}\n")
                    f.write(f"{type(e).__name__}: {e}\n\n")

                continue


def main(
    atom_name_field: str = 'Name',
    folders: Optional[List[str]] = None,
    output_base: Optional[str] = None,
    output_name: Optional[str] = None,
    volume_range: tuple = (3, 22),
    molecule_class: str = 'rna',
    save_csv: bool = True,
    save_plot: bool = False,
    show_points: bool = True,
    show_numbers: bool = True,
    annotate_direct_groups: bool = True,
    plot_min_count: int = 10,
    max_spread: Optional[float] = 2.35,
    ellipse_min_count: int = 10,
    ellipse_max_spread: Optional[float] = 2.15,
    ellipse_n_std: float = 2,
    point_alpha: float = 0.2,
    use_ml_clustering: bool = True,
    ml_method: str = 'kprototypes',  # kmeans, agglomerative, dbscan, hdbscan, kprototypes
    use_sol_binary: bool = True,
    sol_threshold: float = 20.0,
    n_clusters: int = 12,
    min_samples: int = 10,
    eps: float = 1.0,
    min_cluster_size: int = 25,
    numerical_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    show_plot: bool = True,
    prompt_save_after_plot: bool = False,
    downsample_fraction: Optional[float] = None,
    group_hydrogens: bool = False
):
    if numerical_cols is None:
        numerical_cols = ['AW', 'Pow', 'x', 'y']

    if folders is None:
        folders = select_folders_multi()

    if len(folders) == 0:
        print("No folders selected.")
        return

    if len(folders) == 0:
        print("No folders selected.")
        return

    atom_df = collect_atom_volume_points(
        folders=folders,
        atom_name_field=atom_name_field,
        molecule_class=molecule_class,
        volume_range=volume_range
    )

    print("\n=== HYDROGEN CHECK AFTER COLLECTION ===")
    print(atom_df[atom_df['AtomName'].astype(str).str.upper().str.startswith('H')][
              ['Folder', 'Index', 'AtomName', 'CanonicalName', 'AW', 'Pow']
          ].to_string(index=False))

    if str(molecule_class).strip().lower().replace(' ', '_') == 'small_molecule':
        atom_df = tag_small_molecule_environmental_outliers(atom_df)

    if downsample_fraction is None:
        downsample_fraction = 1.0

    print("\n=== UNIQUE CANONICAL NAMES ===")
    print(sorted(atom_df['CanonicalName'].unique()))

    # Temporary test mode: skip solvent-facing calculation until surfs_df is wired in
    if 'SolFacingPct' not in atom_df.columns:
        atom_df['SolFacingPct'] = np.nan

    if len(atom_df) == 0:
        print("No matching AW/Pow atom volume pairs were found.")
        return

    if use_ml_clustering:
        ml_df = atom_df.copy()

        print("\n=== ENTERED ML CLUSTERING BLOCK ===")
        print(f"use_ml_clustering = {use_ml_clustering}")
        print(f"ml_method = {ml_method}")
        print(f"starting rows in ml_df = {len(ml_df)}")
        print(f"columns in ml_df = {list(ml_df.columns)}")

        ml_df['DeltaV'] = ml_df['Pow'] - ml_df['AW']

        if use_sol_binary and 'SolFacingPct' in ml_df.columns:
            ml_df = make_sol_facing_binary(
                ml_df,
                source_col='SolFacingPct',
                out_col='SolFacingBinary',
                threshold=sol_threshold
            )

        default_numeric_cols = ['AW', 'Pow', 'x', 'y']
        default_categorical_cols = ['AtomName', 'ResidueName']
        default_boolean_cols = ['SolFacingBinary'] if use_sol_binary and 'SolFacingBinary' in ml_df.columns else []

        numeric_cols = numerical_cols if numerical_cols is not None else default_numeric_cols
        categorical_cols = categorical_cols if categorical_cols is not None else default_categorical_cols
        boolean_cols = boolean_cols if boolean_cols is not None else default_boolean_cols

        numeric_cols = [col for col in numeric_cols if col in ml_df.columns and ml_df[col].notna().any()]
        categorical_cols = [col for col in categorical_cols if col in ml_df.columns]
        boolean_cols = [col for col in boolean_cols if col in ml_df.columns]
        # numeric_cols = ['AW', 'Pow', 'DeltaV']
        # categorical_cols = []
        print("\n=== FEATURE SELECTION BEFORE DROPNA ===")
        print(f"boolean_cols = {boolean_cols}")
        print(f"numeric_cols = {numeric_cols}")
        print(f"categorical_cols candidate = {[col for col in ['AtomName', 'ResidueName'] if col in ml_df.columns]}")
        print("\nNon-null counts:")
        print(ml_df[numeric_cols + [col for col in ['AtomName', 'ResidueName'] if col in ml_df.columns]].notna().sum())

        ml_df = ml_df.dropna(subset=numeric_cols + categorical_cols).copy()

        if len(ml_df) > 0:
            if ml_method in ['kmeans', 'agglomerative']:
                ml_df = run_clustering(
                    df=ml_df,
                    method=ml_method,
                    numeric_cols=numeric_cols,
                    categorical_cols=categorical_cols,
                    boolean_cols=boolean_cols,
                    n_clusters=n_clusters
                )

            elif ml_method == 'dbscan':
                ml_df = run_clustering(
                    df=ml_df,
                    method=ml_method,
                    numeric_cols=numeric_cols,
                    categorical_cols=categorical_cols,
                    boolean_cols=boolean_cols,
                    eps=eps,
                    min_samples=min_samples
                )

            elif ml_method == 'hdbscan':
                ml_df = run_clustering(
                    df=ml_df,
                    method=ml_method,
                    numeric_cols=numeric_cols,
                    categorical_cols=categorical_cols,
                    boolean_cols=boolean_cols,
                    min_cluster_size=min_cluster_size,
                    min_samples=min_samples
                )

            elif ml_method == 'kprototypes':
                ml_df = run_clustering(
                    df=ml_df,
                    method=ml_method,
                    numeric_cols=numeric_cols,
                    categorical_cols=categorical_cols,
                    boolean_cols=boolean_cols,
                    n_clusters=n_clusters
                )
            print("\n=== CLUSTERING RETURNED ===")
            print(f"returned columns = {list(ml_df.columns)}")

            if 'MLCluster' in ml_df.columns:
                print("MLCluster value counts:")
                print(ml_df['MLCluster'].value_counts(dropna=False).sort_index())

                print("\nSample rows with MLCluster:")
                cols_to_show = [col for col in ['AtomName', 'ResidueName', 'CanonicalName', 'AW', 'Pow', 'x', 'y', 'SolFacingPct', 'SolFacingBinary', 'MLCluster'] if col in ml_df.columns]
                print(ml_df[cols_to_show].head(25).to_string(index=False))
            else:
                print("WARNING: MLCluster column was not added by run_clustering.")

            cluster_summary = summarize_clusters(ml_df)
            print_cluster_summary(cluster_summary)
            print("\n=== CHECKING MANUAL VS ML ===")
            if 'CanonicalName' in ml_df.columns and 'MLCluster' in ml_df.columns:
                print("Unique CanonicalName count:", ml_df['CanonicalName'].nunique())
                print("Unique MLCluster count:", ml_df['MLCluster'].nunique(dropna=True))
            else:
                print("Missing CanonicalName or MLCluster for comparison.")
            comparison_table = compare_manual_vs_ml(
                ml_df,
                manual_col='CanonicalName',
                ml_col='MLCluster'
            )

            print("\n=== MANUAL VS ML CLUSTER CROSS-TAB ===")
            print(comparison_table.to_string())

            out_dir = output_base if output_base is not None else folders[0]
            os.makedirs(out_dir, exist_ok=True)
            ml_df.to_csv(os.path.join(out_dir, f'ml_clusters_{ml_method}.csv'), index=False)
            cluster_summary.to_csv(os.path.join(out_dir, f'ml_cluster_summary_{ml_method}.csv'), index=False)
            comparison_table.to_csv(os.path.join(out_dir, f'manual_vs_ml_{ml_method}.csv'))

            print(f"\nSaved ML clustering outputs for method: {ml_method}")

            print("\n=== RETURNING TO MANUAL GROUP ANALYSIS ===")
            print("The remaining stats/plots below still use CanonicalName, not MLCluster.")

    if use_ml_clustering and 'MLCluster' in ml_df.columns:
        ml_df['MLCluster'] = ml_df['MLCluster'].astype(str)
        stats_df = compute_name_volume_stats(ml_df, group_col='MLCluster')
        plot_df = ml_df
        plot_group_col = 'MLCluster'

    else:
        stats_df = compute_name_volume_stats(atom_df, group_col='CanonicalName')

        mol = str(molecule_class).strip().lower().replace(' ', '_')

        if mol == 'protein':
            merge_map = build_merged_group_map(
                stats_df=stats_df,
                distance_cutoff=0.40,
                spread_cutoff=1.50,
                min_count_to_merge=00
            )

            atom_df['MergedCanonicalName'] = atom_df['CanonicalName'].map(merge_map).fillna(atom_df['CanonicalName'])

            print("\n=== PROTEIN MERGE MAP ===")
            for old_name, new_name in sorted(merge_map.items()):
                if old_name != new_name:
                    print(f"{old_name} -> {new_name}")

            stats_df = compute_name_volume_stats(atom_df, group_col='MergedCanonicalName')
            plot_df = atom_df
            plot_group_col = 'MergedCanonicalName'

        else:
            plot_df = atom_df
            plot_group_col = 'CanonicalName'

    print_name_volume_stats(stats_df)

    # if str(molecule_class).strip().lower().replace(' ', '_') == 'small_molecule':
    #     bad_groups = stats_df[
    #         stats_df['Members'].astype(str).str.contains(',') &
    #         stats_df['GroupName'].astype(str).str.contains('C_')
    #         ]
    #
    #     if len(bad_groups) > 0:
    #         print("\nWARNING: Some small-molecule carbon groups still contain multiple atom names:")
    #         print(bad_groups[['GroupName', 'Count', 'Spread', 'Members']].to_string(index=False))

    plot_stats_df = filter_plot_groups(
        stats_df=stats_df,
        plot_min_count=plot_min_count,
        max_spread=max_spread,
        force_include=set() if str(molecule_class).lower().replace(' ', '_') == 'small_molecule'
                            else get_force_include_groups(molecule_class)
        )

    if not group_hydrogens:
        plot_stats_df = plot_stats_df[
            ~plot_stats_df['GroupName'].apply(is_hydrogen_group)
        ].copy()

    plot_stats_df = sort_plot_groups(plot_stats_df, molecule_class)

    if str(molecule_class).strip().lower().replace(' ', '_') == 'small_molecule':
        plot_stats_df = plot_stats_df[
            plot_stats_df['GroupName'] != 'H_environmental_outlier'
            ].copy()

        plot_stats_df['PlotNumber'] = np.arange(1, len(plot_stats_df) + 1)

    results = analyze_group_outliers(atom_df)
    annotated_df = results['annotated_df']

    print("\n=== SMALL GROUPS ===")
    print(results['small_groups'])

    print("\n=== WITHIN-GROUP OUTLIERS ===")
    print(results['within_group_outliers'][[
        'CanonicalName', 'AtomName', 'ResidueName', 'AW', 'Pow', 'z_total'
    ]].head(20))

    print("\n=== GLOBAL OUTLIERS ===")
    print(results['global_outliers'][[
        'CanonicalName', 'AtomName', 'ResidueName', 'AW', 'Pow', 'global_dist'
    ]])

    for name, group in annotated_df.groupby('CanonicalName'):
        group = group.sort_values('z_total', ascending=False)
        worst = group.iloc[0]

        if pd.notna(worst['z_total']) and worst['z_total'] > 2.5:
            residue_name = worst['ResidueName'] if 'ResidueName' in worst.index else 'UNKNOWN'
            print(f"{name}: worst z = {worst['z_total']:.2f} | {worst['AtomName']} ({residue_name})")

    if output_name is None:
        output_name = 'atom_group_volume_plot'

    title = f"Grouped atom volumes in AW vs Pow space (n_folders={len(folders)}, n_atoms={len(atom_df)})"

    plot_name_volume_groups(
        atom_df=plot_df,
        plot_stats_df=plot_stats_df,
        title=title,
        group_col=plot_group_col,
        volume_range=volume_range,
        show_points=show_points,
        show_numbers=show_numbers,
        annotate_direct_groups=annotate_direct_groups,
        point_size=10,
        point_alpha=point_alpha,
        number_fontsize=10,
        direct_label_fontsize=11,
        show_ellipses=True,
        ellipse_n_std=ellipse_n_std,
        ellipse_min_count=ellipse_min_count,
        ellipse_max_spread=ellipse_max_spread,
        save_png=None,
        save_svg=None,
        downsample_fraction=downsample_fraction,
        show=show_plot,
        group_hydrogens=group_hydrogens
    )

    if prompt_save_after_plot:
        out_dir = select_output_folder("Select folder to save plot data and files")
        if out_dir is None:
            print("\nNo output folder selected. Skipping save step.")
            return
    else:
        if output_base is None:
            out_dir = folders[0]
        else:
            out_dir = output_base

    os.makedirs(out_dir, exist_ok=True)

    if save_csv:
        stats_csv_path = os.path.join(out_dir, f'{output_name}_stats.csv')
        key_csv_path = os.path.join(out_dir, f'{output_name}_group_key.csv')
        plot_points_csv_path = os.path.join(out_dir, f'{output_name}_plot_points.csv')

        stats_df.to_csv(stats_csv_path, index=False)
        save_group_plot_table(plot_stats_df, key_csv_path)
        plot_df.to_csv(plot_points_csv_path, index=False)

        print(f"\nSaved full stats CSV -> {stats_csv_path}")
        print(f"Saved plotted points CSV -> {plot_points_csv_path}")

    if save_plot:
        png_path = os.path.join(out_dir, f'{output_name}.png')
        svg_path = os.path.join(out_dir, f'{output_name}.svg')

        plot_name_volume_groups(
            atom_df=plot_df,
            plot_stats_df=plot_stats_df,
            title=title,
            group_col=plot_group_col,
            volume_range=volume_range,
            show_points=show_points,
            show_numbers=show_numbers,
            annotate_direct_groups=annotate_direct_groups,
            point_size=10,
            point_alpha=point_alpha,
            number_fontsize=10,
            direct_label_fontsize=11,
            show_ellipses=True,
            ellipse_n_std=ellipse_n_std,
            ellipse_min_count=ellipse_min_count,
            ellipse_max_spread=ellipse_max_spread,
            save_png=png_path,
            save_svg=svg_path,
            downsample_fraction=downsample_fraction,
            show=False,
            group_hydrogens=group_hydrogens
        )

        print(f"Saved PNG -> {png_path}")
        print(f"Saved SVG -> {svg_path}")


if __name__ == "__main__":
    bif = 'E:/Molecular'
    bof = (
        'E:/OneDrive - Georgia State University/GSU NSC/Manuscripts'
        '/Ericson Voronoi DNA/P2/fig2_atomic_level_scheme_deviations/2C_Full_Plots_Clustered'
    )

    SMOL_SETTS = {
        'class': 'small_molecule',
        'folders': [os.path.join(bif, _) for _ in ['B_EDTA', 'C_DB1976', 'I_T4LP/JZ4']],
        'output': os.path.join(bof, 'small molecule'),
        'output_name': 'small_molecule_manual',
        'downsample_fraction': 1.0,
        'point_alpha': 0.75,
        'plot_min_count': 2,
        'ellipse_min_count': 3,
        'max_spread': 4.5,
        'ellipse_max_spread': 4.5,
        'ellipse_n_std': 1.5,
        'min_samples': 2,
        'min_cluster_size': 2,
        'eps': 1.5,
        'n_clusters': 6,
        'group_hydrogens': False
    }

    RNA_SETTS = {
        'class': 'rna',
        'folders': [os.path.join(bif, _) for _ in ['G_Hammerhead']],
        'output': os.path.join(bof, 'rna'),
        'output_name': 'rna_manual',
        'downsample_fraction': 1.0,
        'point_alpha': 0.35,
        'plot_min_count': 10,
        'ellipse_min_count': 10,
        'max_spread': 4,
        'ellipse_max_spread': 4,
        'ellipse_n_std': 1.5,
        'min_samples': 5,
        'min_cluster_size': 10,
        'eps': 2.0,
        'n_clusters': 10,
        'group_hydrogens': True
    }

    DNA_SETTS = {
        'class': 'dna',
        'folders': [os.path.join(bif, _) for _ in ['D_Hairpin', 'F_BDNA', 'K_NCP_DNA']],
        'output': os.path.join(bof, 'dna'),
        'output_name': 'dna_manual',
        'downsample_fraction': 0.25,
        'point_alpha': 0.25,
        'plot_min_count': 75,
        'ellipse_min_count': 75,
        'max_spread': 4,
        'ellipse_max_spread': 4,
        'ellipse_n_std': 1.6,
        'min_samples': 10,
        'min_cluster_size': 25,
        'eps': 2.0,
        'n_clusters': 12,
        'group_hydrogens': False
    }

    PROT_SETTS = {
        'class': 'protein',
        'folders': [os.path.join(bif, _) for _ in ['E_Cambrin', 'H_p53tet', 'I_T4LP', 'J_Streptavidin', 'L_BSA', 'm_NCP_Protein']],
        'output': os.path.join(bof, 'protein'),
        'output_name': 'protein_manual',
        'downsample_fraction': 0.05,
        'point_alpha': 0.20,
        'plot_min_count': 300,
        'ellipse_min_count': 300,
        'max_spread': 6,
        'ellipse_max_spread': 6,
        'ellipse_n_std': 1.3,
        'min_samples': 15,
        'min_cluster_size': 300,
        'eps': 5,
        'n_clusters': 16,
        'group_hydrogens': False
    }

    settings = {
        'smol': SMOL_SETTS,
        'rna': RNA_SETTS,
        'dna': DNA_SETTS,
        'prot': PROT_SETTS,
    }

    current = 'rna'
    cfg = settings[current]

    main(
        atom_name_field='Name',
        folders=cfg['folders'],
        output_base=cfg['output'],
        output_name=cfg['output_name'],
        volume_range=(3, 22),
        molecule_class=cfg['class'],
        save_csv=True,
        save_plot=True,
        show_points=True,
        show_numbers=True,
        annotate_direct_groups=False,
        plot_min_count=cfg['plot_min_count'],
        max_spread=cfg['max_spread'],
        ellipse_min_count=cfg['ellipse_min_count'],
        ellipse_max_spread=cfg['ellipse_max_spread'],
        ellipse_n_std=cfg['ellipse_n_std'],
        point_alpha=cfg['point_alpha'],
        downsample_fraction=cfg['downsample_fraction'],
        use_ml_clustering=False,
        ml_method='kprototypes',
        use_sol_binary=True,
        sol_threshold=20.0,
        n_clusters=cfg['n_clusters'],
        min_samples=cfg['min_samples'],
        eps=cfg['eps'],
        min_cluster_size=cfg['min_cluster_size'],
        numerical_cols=None,
        categorical_cols=None,
        boolean_cols=None,
        prompt_save_after_plot=False,
        group_hydrogens=cfg['group_hydrogens']
    )
