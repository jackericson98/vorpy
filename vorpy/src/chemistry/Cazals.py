"""
Cazals.py

Cazals/Chothia-style atomic group radii for reproducing the weighted
(power/regular) alpha-complex interface calculations described for
protein-protein interfaces such as 2KAI.

This module intentionally follows the simple VorPy radii-module API:

    element_radii
    special_radii

All values are in Angstroms (A).

Base group radii:
    aliphatic carbon : 1.87
    trigonal carbon  : 1.76
    neutral nitrogen : 1.65
    charged nitrogen : 1.50
    oxygen            : 1.40
    sulfur            : 1.85
    other             : 2.00

For the Cazals solvent-accessible construction, add PROBE_RADIUS = 1.40 A
to each base radius before constructing the power diagram, i.e.

    R_i = r_i + 1.40
    weight_i = R_i ** 2

Do NOT bake the probe into special_radii if the same module is to be usable
for both base-radius inspection and the Cazals alpha=0 construction.

Notes
-----
1. Carbon and nitrogen cannot be classified correctly from element alone.
   The residue/atom-specific table below performs the classification for
   the 20 standard amino acids.

2. Hydrogen is not one of the explicitly tabulated C/N/O/S group classes.
   It therefore receives OTHER_RADIUS here. For an exact historical
   reproduction, verify whether the source structure/model includes
   explicit hydrogens; many crystallographic PDB models do not.

3. Histidine protonation and terminal/protonation variants can be chemically
   ambiguous from a bare PDB atom name. Standard HIS ring nitrogens are
   treated as neutral here. Explicit charged/protonated variants should be
   reviewed before claiming an exact reproduction.

4. The dictionary includes common PDB terminal oxygen names where practical.
"""

PROBE_RADIUS = 1.40

ALIPHATIC_C = 1.87
TRIGONAL_C = 1.76
NEUTRAL_N = 1.65
CHARGED_N = 1.50
OXYGEN = 1.40
SULFUR = 1.85
OTHER_RADIUS = 2.00

# Element-only fallbacks. C and N are necessarily ambiguous at this level,
# so these values should not be used for standard protein atoms when the
# residue and atom name are known.
element_radii = {
    'C': OTHER_RADIUS,
    'N': OTHER_RADIUS,
    'O': OXYGEN,
    'S': SULFUR,
    'H': OTHER_RADIUS,
    'P': OTHER_RADIUS,
}

# Classification helpers used to generate special_radii.
# Protein backbone:
#   CA = aliphatic carbon
#   C  = trigonal peptide carbonyl carbon
#   N  = neutral peptide nitrogen
#   O  = oxygen

_ALIPHATIC = ALIPHATIC_C
_TRIGONAL = TRIGONAL_C
_N_NEUTRAL = NEUTRAL_N
_N_CHARGED = CHARGED_N
_O = OXYGEN
_S = SULFUR

special_radii = {
    'ALA': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'OC1': _O, 'OC2': _O,
    },
    'ARG': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC, 'CD': _ALIPHATIC,
        'NE': _N_CHARGED, 'CZ': _TRIGONAL,
        'NH1': _N_CHARGED, 'NH2': _N_CHARGED,
    },
    'ASN': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _TRIGONAL, 'OD1': _O, 'ND2': _N_NEUTRAL,
    },
    'ASP': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _TRIGONAL, 'OD1': _O, 'OD2': _O,
    },
    'CYS': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'SG': _S, 'S': _S,
    },
    'GLN': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC, 'CD': _TRIGONAL,
        'OE1': _O, 'NE2': _N_NEUTRAL,
    },
    'GLU': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC, 'CD': _TRIGONAL,
        'OE1': _O, 'OE2': _O,
    },
    'GLY': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'OC1': _O, 'OC2': _O,
    },
    'HIS': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC,
        'CG': _TRIGONAL, 'ND1': _N_NEUTRAL, 'CD2': _TRIGONAL,
        'CE1': _TRIGONAL, 'NE2': _N_NEUTRAL,
    },
    'ILE': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG1': _ALIPHATIC, 'CG2': _ALIPHATIC,
        'CD': _ALIPHATIC, 'CD1': _ALIPHATIC,
    },
    'LEU': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC,
        'CD1': _ALIPHATIC, 'CD2': _ALIPHATIC,
    },
    'LYS': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC, 'CD': _ALIPHATIC,
        'CE': _ALIPHATIC, 'NZ': _N_CHARGED,
    },
    'MET': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC, 'SD': _S, 'S': _S,
        'CE': _ALIPHATIC,
    },
    'PHE': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC,
        'CG': _TRIGONAL, 'CD1': _TRIGONAL, 'CD2': _TRIGONAL,
        'CE1': _TRIGONAL, 'CE2': _TRIGONAL, 'CZ': _TRIGONAL,
    },
    'PRO': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG': _ALIPHATIC, 'CD': _ALIPHATIC,
    },
    'SER': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'OG': _O,
    },
    'THR': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'OG1': _O, 'OG': _O, 'CG2': _ALIPHATIC,
    },
    'TRP': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC,
        'CG': _TRIGONAL, 'CD1': _TRIGONAL, 'CD2': _TRIGONAL,
        'NE1': _N_NEUTRAL, 'CE2': _TRIGONAL, 'CE3': _TRIGONAL,
        'CZ2': _TRIGONAL, 'CZ3': _TRIGONAL, 'CH2': _TRIGONAL,
        # aliases appearing in some existing VorPy tables
        'CD': _TRIGONAL, 'CE': _TRIGONAL, 'CZ': _TRIGONAL,
        'CZ1': _TRIGONAL, 'CH': _TRIGONAL,
    },
    'TYR': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC,
        'CG': _TRIGONAL, 'CD1': _TRIGONAL, 'CD2': _TRIGONAL,
        'CE1': _TRIGONAL, 'CE2': _TRIGONAL, 'CZ': _TRIGONAL,
        'OH': _O,
        'CD': _TRIGONAL, 'CE': _TRIGONAL,
    },
    'VAL': {
        'N': _N_NEUTRAL, 'CA': _ALIPHATIC, 'C': _TRIGONAL, 'O': _O, 'OXT': _O,
        'CB': _ALIPHATIC, 'CG1': _ALIPHATIC, 'CG2': _ALIPHATIC,
    },
}


def base_radius(residue_name, atom_name, element=None):
    """
    Return the base Cazals/Chothia group radius.

    Standard protein residue/atom assignments are preferred. If no specific
    assignment exists, O and S can be assigned safely by element. Ambiguous
    C/N fallbacks return OTHER_RADIUS rather than silently guessing hybridization
    or charge state.
    """
    residue = str(residue_name).strip().upper()
    atom = str(atom_name).strip().upper()
    elem = "" if element is None else str(element).strip().upper()

    if residue in special_radii and atom in special_radii[residue]:
        return special_radii[residue][atom]

    if elem in element_radii:
        return element_radii[elem]

    return OTHER_RADIUS


def expanded_radius(residue_name, atom_name, element=None, probe_radius=PROBE_RADIUS):
    """Return the probe-expanded radius R_i = r_i + probe_radius."""
    return base_radius(residue_name, atom_name, element) + float(probe_radius)


def power_weight(residue_name, atom_name, element=None, probe_radius=PROBE_RADIUS):
    """Return the Cazals power weight w_i = (r_i + probe_radius)^2."""
    radius = expanded_radius(residue_name, atom_name, element, probe_radius)
    return radius * radius
