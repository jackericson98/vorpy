
"""
This module defines atomic radii for various elements and amino acid residues used in Voronoi tessellation calculations.

The module contains two main dictionaries:

1. element_radii: A dictionary mapping element symbols to their atomic radii in Angstroms.
   - Contains standard atomic radii for common elements (H, C, N, O, etc.)
   - Includes radii for metals and other special elements
   - Used as default values when specific residue/atom combinations are not found

2. special_radii: A dictionary containing detailed atomic radii for amino acid residues.
   - Keys are three-letter amino acid codes (e.g., 'ALA', 'ARG')
   - Each amino acid entry contains a sub-dictionary mapping atom names to their specific radii
   - Includes special handling for backbone atoms (N, CA, C, O) and side chain atoms
   - Provides more accurate radii for protein structure analysis

These radii are used in Voronoi tessellation calculations to determine atom-atom contacts
and molecular surface areas. The values are optimized for protein structure analysis
and molecular modeling applications.
"""


# Default radii and other radii, excluding protein main chain and side chain radii
element_radii = {
    'H': 1.30, 'HE': 1.40, 'LI': 0.76, 'BE': 0.45, 'B': 1.92, 'C': 1.80, 'N': 1.60, 'O': 1.50, 'P': 1.90,
    'S': 1.90, 'F': 1.33, 'CL': 1.81, 'BR': 1.96, 'I': 2.20, 'AL': 0.60, 'AS': 0.58, 'SE': 1.90, 'AU': 1.37, 'BA': 1.35,
    'BI': 1.03, 'CA': 1.00, 'CD': 0.95, 'CO': 0.65, 'CR': 0.73, 'CS': 1.67, 'CU': 0.73, 'FE': 0.61, 'GA': 0.62,
    'GE': 0.73, 'HG': 1.02, 'K': 1.38, 'MG': 0.72, 'MN': 0.83, 'MO': 0.69, 'NA': 1.02, 'NI': 0.69, 'PB': 1.19,
    'PD': 0.86, 'PT': 0.80, 'RB': 1.52, 'SB': 0.76, 'SC': 0.75, 'SN': 0.69, 'SR': 1.18, 'TC': 0.65, 'TI': 0.86,
    'V': 0.79, 'ZN': 0.74, 'ZR': 0.72
}

special_radii = {
    'ALA': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HB3': 1.3, 'N': 1.7, 'O': 1.49, 'OC1': 1.5, 'OC2': 1.5
    },
    'ARG': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.88, 'CG': 1.92, 'CZ': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HB3': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HD3': 1.3, 'HE': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'HH11': 1.3, 'HH12': 1.3, 'HH21': 1.3, 'HH22': 1.3, 'N': 1.7, 'NE': 1.62, 'NH1': 1.62, 'NH2': 1.67, 'O': 1.49
    },
    'ASN': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CG': 1.81, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD21': 1.3, 'HD22': 1.3, 'N': 1.7, 'ND2': 1.62, 'O': 1.49, 'OD1': 1.52
    },
    'ASP': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CG': 1.76, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'N': 1.7, 'O': 1.49, 'OD1': 1.49, 'OD2': 1.49
    },
    'CYS': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG': 1.3, 'N': 1.7, 'O': 1.49, 'S': 1.88, 'SG': 1.88
    },
    'GLN': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.81, 'CG': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HE21': 1.3, 'HE22': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'NE2': 1.62, 'O': 1.49, 'OE1': 1.52
    },
    'GLU': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.76, 'CG': 1.88, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'O': 1.49, 'OE1': 1.49, 'OE2': 1.49
    },
    'GLY': {
        'C': 1.75, 'CA': 1.9, 'H': 1.3, 'HA1': 1.3, 'HA2': 1.3, 'N': 1.7, 'O': 1.49, 'OC1': 1.5, 'OC2': 1.5
    },
    'HIS': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.74, 'CE': 1.74, 'CD2': 1.74, 'CE1': 1.74, 'CG': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'N': 1.7, 'ND1': 1.6, 'ND2': 1.6, 'NE2': 1.6, 'O': 1.49
    },
    'ILE': {
        'C': 1.75, 'CA': 1.9, 'CB': 2.01, 'CD': 1.92, 'CD1': 1.92, 'CG1': 1.92, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HD3': 1.3, 'HD11': 1.3, 'HD12': 1.3, 'HD13': 1.3, 'HG12': 1.3, 'HG13': 1.3, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.7, 'O': 1.49
    },
    'LEU': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD1': 1.92, 'CD2': 1.92, 'CG': 2.01, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD11': 1.3, 'HD12': 1.3, 'HD13': 1.3, 'HD21': 1.3, 'HD22': 1.3, 'HD23': 1.3, 'HG': 1.3, 'N': 1.7, 'O': 1.49
    },
    'LYS': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.92, 'CE': 1.88, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'HZ1': 1.3, 'HZ2': 1.3, 'HZ3': 1.3, 'N': 1.7, 'NZ': 1.67, 'O': 1.49
    },
    'MET': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CE': 1.8, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HE3': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'O': 1.49, 'SD': 1.94, 'S': 1.94
    },
    'PHE': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD1': 1.82, 'CD2': 1.82, 'CD': 1.82, 'CE1': 1.82, 'CE2': 1.82, 'CG': 1.74, 'CZ': 1.82, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HZ': 1.3, 'N': 1.7, 'O': 1.49
    },
    'PRO': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.92, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'O': 1.49
    },
    'SER': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'H1': 1.3, 'H2': 1.3, 'H3': 1.3, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG': 1.3, 'N': 1.7, 'O': 1.49, 'OG': 1.54
    },
    'THR': {
        'C': 1.75, 'CA': 1.9, 'CB': 2.01, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HG1': 1.3, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.7, 'O': 1.49, 'OG1': 1.54, 'OG': 1.54
    },
    'TRP': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD1': 1.82, 'CD2': 1.82, 'CD': 1.82, 'CE': 1.82, 'CE2': 1.74, 'CE3': 1.82, 'CG': 1.74, 'CH': 1.82, 'CH2': 1.82, 'CZ': 1.82, 'CZ1': 1.82, 'CZ2': 1.82, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HE1': 1.3, 'HE3': 1.3, 'HZ1': 1.3, 'HZ2': 1.3, 'HH2': 1.3, 'N': 1.7, 'NE1': 1.66, 'O': 1.49
    },
    'TYR': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.82, 'CD1': 1.82, 'CD2': 1.82, 'CE': 1.82, 'CE1': 1.82, 'CE2': 1.82, 'CG': 1.74, 'CZ': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HH': 1.3, 'N': 1.7, 'O': 1.49, 'OH': 1.54
    },
    'VAL': {
        'C': 1.75, 'CA': 1.9, 'CB': 2.01, 'CG1': 1.92, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HG11': 1.3, 'HG12': 1.3, 'HG13': 1.3, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.7, 'O': 1.49
    }
}


