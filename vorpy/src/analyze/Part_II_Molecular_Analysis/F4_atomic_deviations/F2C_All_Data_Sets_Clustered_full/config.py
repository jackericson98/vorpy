
from typing import Set

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
