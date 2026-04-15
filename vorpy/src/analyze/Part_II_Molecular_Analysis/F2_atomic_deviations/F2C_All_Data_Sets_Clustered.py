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
    # Hydrogens
    'H': 'H', 'H1': 'H', 'H2': 'H', 'H3': 'H',

    # Carbons
    'C': 'C',
    'C1': 'C', 'C2': 'C', 'C3': 'C', 'C4': 'C',
    'C5': 'C', 'C6': 'C', 'C7': 'C', 'C8': 'C', 'C9': 'C',

    # Nitrogen
    'N': 'N', 'N1': 'N', 'N2': 'N', 'N3': 'N',

    # Oxygen (first pass: collapse all)
    'O': 'O',
    'O1': 'O', 'O2': 'O', 'O3': 'O', 'O4': 'O',
    'O5': 'O', 'O6': 'O', 'O7': 'O',
}

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

    # carbon symmetry
    # 'CG1': 'CG',
    # 'CG2': 'CG',
    # 'CE3': 'CE2',
    # 'CH2': 'CZ',
    # 'CZ2': 'CZ',
    # 'CZ3': 'CZ',

    # nitrogen symmetry
    # 'NH1': 'NH',
    # 'NH2': 'NH',
    # 'NE1': 'NE',
    # 'ND1': 'NE',

    # oxygen grouping
    # 'O': 'O_backbone',
    # 'OH': 'O_backbone',

    # 'OE1': 'O_carboxyl_1',
    # 'OD1': 'O_carboxyl_1',
    #
    # 'OE2': 'O_carboxyl_2',
    # 'OD2': 'O_carboxyl_2',
    #
    # 'OC1': 'O_carboxyl_2',
    # 'OC2': 'O_carboxyl_2',
    # 'OT1': 'O_carboxyl_2',
    # 'OT2': 'O_carboxyl_2',

    # optional CB neighborhood merge
    # 'C2': 'CB',
    # 'C3': 'CB',
}

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
    'O2': 'O_base_2',
    'O4': 'O_base_4',
    'O6': 'O_base_6',

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


DIRECT_LABEL_GROUPS: Set[str] = {
    'O_backbone',
    'O_carboxyl_1',
    'O_carboxyl_2',
}


FORCE_INCLUDE_GROUPS: Set[str] = {
    'O_backbone',
    'O_carboxyl_1',
    'O_carboxyl_2',
    'CA',
    'CB',
    'CG',
    'CD',
    'CE',
    'C',
    'CZ',
    'CE2',
    'CD1',
    'CD2',
    'CE1',
}


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


def resolve_group_folders(data_root: str, relative_folders: List[str]) -> List[str]:
    folders = []

    for rel_path in relative_folders:
        full_path = os.path.join(data_root, rel_path)
        if os.path.isdir(full_path):
            folders.append(full_path)
        else:
            print(f"WARNING: folder not found -> {full_path}")

    return folders


def canonicalize_atom_name(
    atom_name: str,
    molecule_class: str,
    residue_name: str = ''
) -> str:
    atom = str(atom_name).strip().upper()
    mol = str(molecule_class).strip().lower()

    if mol == 'protein':
        return PROTEIN_ATOM_ALIASES.get(atom, atom)

    if mol == 'dna':
        return DNA_ATOM_ALIASES.get(atom, atom)

    if mol == 'rna':
        return RNA_ATOM_ALIASES.get(atom, atom)

    if mol == 'small_molecule':
        return SMALL_MOL_ATOM_ALIASES.get(atom, atom)

    return atom


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
        aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)

    except FileNotFoundError:
        aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True)

    return aw_logs, pow_logs, prm_logs


def select_folders_multi(title: str = "Select a folder (Cancel to finish)") -> List[str]:
    root = tk.Tk()
    root.withdraw()

    folders = []

    while True:
        folder = filedialog.askdirectory(title=title)
        if not folder:
            break
        folders.append(folder)

    root.destroy()
    return folders


def collect_atom_volume_points(
    folders: List[str],
    atom_name_field: str = 'Name',
    molecule_class: str = 'protein',
    volume_range: Optional[tuple] = None
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []

    for folder in folders:
        aw_logs, pow_logs, _ = _read_scheme_logs(folder)

        aw_atoms = aw_logs['atoms']
        pow_atoms = pow_logs['atoms']

        pow_lookup = {
            int(row['Index']): row
            for _, row in pow_atoms.iterrows()
        }

        for _, atom in aw_atoms.iterrows():
            idx = int(atom['Index'])
            chain = atom['Chain']
            res_seq = atom['Residue Sequence']

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

            if not atom_name:
                continue

            residue_name = str(atom.get('Residue', atom.get('Residue Name', atom.get('ResName', '')))).strip().upper()
            canonical_name = canonicalize_atom_name(
                atom_name=atom_name,
                molecule_class=molecule_class,
                residue_name=residue_name
            )

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
                'x': float(atom['x']) if 'x' in atom.index else np.nan,
                'y': float(atom['y']) if 'y' in atom.index else np.nan,
                'SolFacingPct': float(atom['SolFacingPct']) if 'SolFacingPct' in atom.index else np.nan,
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
    if name.startswith('H'):
        return '#1f77b4'

    if name.startswith('C'):
        return '#ff7f0e'

    if name.startswith('N'):
        return '#2ca02c'

    if name.startswith('O'):
        return '#d62728'

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
        zorder=zorder
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
    ellipse_n_std: float = 1.5,
    ellipse_min_count: int = 50,
    ellipse_max_spread: Optional[float] = None,
    save_png: Optional[str] = None,
    save_svg: Optional[str] = None,
    show: bool = True
):
    fig, ax = plt.subplots(figsize=(12, 9))

    print_group_plot_table(plot_stats_df)

    for _, row in plot_stats_df.iterrows():
        name = row['GroupName']
        if group_col == 'MLCluster':
            color = plt.cm.tab20(int(name) % 20)
        else:
            color = get_plot_color(str(name))

        group_df = atom_df[atom_df[group_col].astype(str) == str(name)]
        print(f"\nGroupName from plot_stats_df: {name}")
        print(f"Matching rows in atom_df['CanonicalName'] == {name}: {len(group_df)}")
        if len(group_df) == 0:
            print("WARNING: empty group_df")
            print("Sample CanonicalName values:")
            print(atom_df['CanonicalName'].dropna().astype(str).unique()[:20])

        if show_points:
            ax.scatter(
                group_df['AW'],
                group_df['Pow'],
                s=point_size,
                alpha=point_alpha,
                color=color,
                zorder=1
            )

        can_draw_ellipse = (
            show_ellipses and
            len(group_df) >= ellipse_min_count and
            (ellipse_max_spread is None or row['Spread'] <= ellipse_max_spread)
        )

        if can_draw_ellipse:
            add_covariance_ellipse(
                ax=ax,
                x=group_df['AW'].to_numpy(),
                y=group_df['Pow'].to_numpy(),
                color=color,
                n_std=ellipse_n_std,
                face_alpha=0.20,
                edge_alpha=1.0,
                linewidth=1.5,
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
    molecule_class: str = 'protein',
    save_csv: bool = True,
    save_plot: bool = False,
    show_points: bool = True,
    show_numbers: bool = True,
    annotate_direct_groups: bool = True,
    plot_min_count: int = 100,
    max_spread: Optional[float] = 2.35,
    ellipse_min_count: int = 150,
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
    boolean_cols=None
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

    # ungrouped_df, ungrouped_summary, ungrouped_atom_res_summary = find_ungrouped_atoms(
    #     atom_df,
    #     aw_range=(11.0, 13.5),
    #     pow_range=(13.0, 15.0),
    #     top_k=100
    # )
    #
    # if len(ungrouped_df) > 0:
    #     print("\n=== RAW UNGROUPED ATOMS IN PATCH ===")
    #     cols = [c for c in ['AtomName', 'ResidueName', 'AW', 'Pow', 'CanonicalName'] if c in ungrouped_df.columns]
    #     print(
    #         ungrouped_df[cols]
    #         .sort_values(['AtomName', 'ResidueName', 'AW', 'Pow'])
    #         .head(100)
    #         .to_string(index=False)
    #     )
    #
    # patch_df = extract_patch_atoms(
    #     atom_df,
    #     aw_range=(12.0, 13.5),
    #     pow_range=(13.0, 14.8),
    #     min_count=2
    # )

    # out_path = os.path.join(os.getcwd(), "patch_atoms.csv")
    # patch_df[['Index']].to_csv(out_path, index=False)
    #
    # print(f"\nSaved patch atoms to:\n{out_path}")

    if use_ml_clustering and 'MLCluster' in ml_df.columns:
        ml_df['MLCluster'] = ml_df['MLCluster'].astype(str)
        stats_df = compute_name_volume_stats(ml_df, group_col='MLCluster')
        plot_df = ml_df
        plot_group_col = 'MLCluster'
    else:
        stats_df = compute_name_volume_stats(atom_df, group_col='CanonicalName')
        plot_df = atom_df
        plot_group_col = 'CanonicalName'

    print_name_volume_stats(stats_df)

    plot_stats_df = filter_plot_groups(
        stats_df=stats_df,
        plot_min_count=plot_min_count,
        max_spread=max_spread,
        force_include=FORCE_INCLUDE_GROUPS
    )

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

    if output_base is None:
        out_dir = folders[0]
    else:
        out_dir = output_base

    os.makedirs(out_dir, exist_ok=True)

    if output_name is None:
        output_name = 'atom_group_volume_plot'

    if save_csv:
        stats_csv_path = os.path.join(out_dir, f'{output_name}_stats.csv')
        key_csv_path = os.path.join(out_dir, f'{output_name}_group_key.csv')

        stats_df.to_csv(stats_csv_path, index=False)
        save_group_plot_table(plot_stats_df, key_csv_path)

        print(f"\nSaved full stats CSV -> {stats_csv_path}")

    png_path = None
    svg_path = None

    if save_plot:
        png_path = os.path.join(out_dir, f'{output_name}.png')
        svg_path = os.path.join(out_dir, f'{output_name}.svg')

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
        number_fontsize=12,
        direct_label_fontsize=11,
        show_ellipses=True,
        ellipse_n_std=ellipse_n_std,
        ellipse_min_count=ellipse_min_count,
        ellipse_max_spread=ellipse_max_spread,
        save_png=png_path,
        save_svg=svg_path,
        show=False
    )

    if png_path is not None:
        print(f"Saved PNG -> {png_path}")

    if svg_path is not None:
        print(f"Saved SVG -> {svg_path}")


if __name__ == "__main__":
    run_batch_cluster_plots()
