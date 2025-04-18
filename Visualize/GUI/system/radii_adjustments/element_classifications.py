

# Define basic properties of elements (symbol, atomic number, atomic mass, atomic radius)
elements = {
        'H': {'name': 'Hydrogen', 'number': 1, 'mass': 1.008, 'radius': 1.30, 'row': 1, 'column': 1,
              'group': 'Nonmetal'},
        'He': {'name': 'Helium', 'number': 2, 'mass': 4.003, 'radius': 1.40, 'row': 1, 'column': 18,
               'group': 'Noble Gas'},
        'Li': {'name': 'Lithium', 'number': 3, 'mass': 6.941, 'radius': 0.76, 'row': 2, 'column': 1,
               'group': 'Alkali Metal'},
        'Be': {'name': 'Beryllium', 'number': 4, 'mass': 9.012, 'radius': 0.45, 'row': 2, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'B': {'name': 'Boron', 'number': 5, 'mass': 10.811, 'radius': 1.92, 'row': 2, 'column': 13,
              'group': 'Metalloid'},
        'C': {'name': 'Carbon', 'number': 6, 'mass': 12.011, 'radius': 1.80, 'row': 2, 'column': 14,
              'group': 'Nonmetal'},
        'N': {'name': 'Nitrogen', 'number': 7, 'mass': 14.007, 'radius': 1.60, 'row': 2, 'column': 15,
              'group': 'Nonmetal'},
        'O': {'name': 'Oxygen', 'number': 8, 'mass': 15.999, 'radius': 1.50, 'row': 2, 'column': 16,
              'group': 'Nonmetal'},
        'F': {'name': 'Fluorine', 'number': 9, 'mass': 18.998, 'radius': 1.33, 'row': 2, 'column': 17,
              'group': 'Halogens'},
        'Ne': {'name': 'Neon', 'number': 10, 'mass': 20.180, 'radius': 1.54, 'row': 2, 'column': 18,
               'group': 'Noble Gas'},
        'Na': {'name': 'Sodium', 'number': 11, 'mass': 22.990, 'radius': 1.02, 'row': 3, 'column': 1,
               'group': 'Alkali Metal'},
        'Mg': {'name': 'Magnesium', 'number': 12, 'mass': 24.305, 'radius': 0.72, 'row': 3, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Al': {'name': 'Aluminum', 'number': 13, 'mass': 26.982, 'radius': 0.60, 'row': 3, 'column': 13,
               'group': 'Post-transition Metal'},
        'Si': {'name': 'Silicon', 'number': 14, 'mass': 28.086, 'radius': 2.10, 'row': 3, 'column': 14,
               'group': 'Metalloid'},
        'P': {'name': 'Phosphorus', 'number': 15, 'mass': 30.974, 'radius': 1.90, 'row': 3, 'column': 15,
              'group': 'Nonmetal'},
        'S': {'name': 'Sulfur', 'number': 16, 'mass': 32.066, 'radius': 1.90, 'row': 3, 'column': 16,
              'group': 'Nonmetal'},
        'Cl': {'name': 'Chlorine', 'number': 17, 'mass': 35.453, 'radius': 1.81, 'row': 3, 'column': 17,
               'group': 'Halogens'},
        'Ar': {'name': 'Argon', 'number': 18, 'mass': 39.948, 'radius': 1.88, 'row': 3, 'column': 18,
               'group': 'Noble Gas'},
        'K': {'name': 'Potassium', 'number': 19, 'mass': 39.098, 'radius': 1.38, 'row': 4, 'column': 1,
              'group': 'Alkali Metal'},
        'Ca': {'name': 'Calcium', 'number': 20, 'mass': 40.078, 'radius': 1.00, 'row': 4, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Ga': {'name': 'Gallium', 'number': 31, 'mass': 69.723, 'radius': 0.62, 'row': 4, 'column': 13,
               'group': 'Post-transition Metal'},
        'Ge': {'name': 'Germanium', 'number': 32, 'mass': 72.631, 'radius': 0.73, 'row': 4, 'column': 14,
               'group': 'Metalloid'},
        'As': {'name': 'Arsenic', 'number': 33, 'mass': 74.922, 'radius': 0.58, 'row': 4, 'column': 15,
               'group': 'Metalloid'},
        'Se': {'name': 'Selenium', 'number': 34, 'mass': 78.971, 'radius': 1.90, 'row': 4, 'column': 16,
               'group': 'Nonmetal'},
        'Br': {'name': 'Bromine', 'number': 35, 'mass': 79.904, 'radius': 1.83, 'row': 4, 'column': 17,
               'group': 'Halogens'},
        'Kr': {'name': 'Krypton', 'number': 36, 'mass': 83.798, 'radius': 2.02, 'row': 4, 'column': 18,
               'group': 'Noble Gas'},
        'Rb': {'name': 'Rubidium', 'number': 37, 'mass': 85.468, 'radius': 1.52, 'row': 5, 'column': 1,
               'group': 'Alkali Metal'},
        'Sr': {'name': 'Strontium', 'number': 38, 'mass': 87.62, 'radius': 1.18, 'row': 5, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'In': {'name': 'Indium', 'number': 49, 'mass': 114.818, 'radius': 1.93, 'row': 5, 'column': 13,
               'group': 'Post-transition Metal'},
        'Sn': {'name': 'Tin', 'number': 50, 'mass': 118.711, 'radius': 2.17, 'row': 5, 'column': 14,
               'group': 'Post-transition Metal'},
        'Sb': {'name': 'Antimony', 'number': 51, 'mass': 121.760, 'radius': 2.06, 'row': 5, 'column': 15,
               'group': 'Metalloid'},
        'Te': {'name': 'Tellurium', 'number': 52, 'mass': 127.6, 'radius': 2.06, 'row': 5, 'column': 16,
               'group': 'Metalloid'},
        'I': {'name': 'Iodine', 'number': 53, 'mass': 126.904, 'radius': 2.20, 'row': 5, 'column': 17,
              'group': 'Halogens'},
        'Xe': {'name': 'Xenon', 'number': 54, 'mass': 131.293, 'radius': 2.16, 'row': 5, 'column': 18,
               'group': 'Noble Gas'},
        'Cs': {'name': 'Cesium', 'number': 55, 'mass': 132.905, 'radius': 1.67, 'row': 6, 'column': 1,
               'group': 'Alkali Metal'},
        'Ba': {'name': 'Barium', 'number': 56, 'mass': 137.328, 'radius': 1.35, 'row': 6, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Tl': {'name': 'Thallium', 'number': 81, 'mass': 204.383, 'radius': 1.96, 'row': 6, 'column': 13,
               'group': 'Post-transition Metal'},
        'Pb': {'name': 'Lead', 'number': 82, 'mass': 207.2, 'radius': 2.02, 'row': 6, 'column': 14,
               'group': 'Post-transition Metal'},
        'Bi': {'name': 'Bismuth', 'number': 83, 'mass': 208.980, 'radius': 2.07, 'row': 6, 'column': 15,
               'group': 'Post-transition Metal'},
        'Po': {'name': 'Polonium', 'number': 84, 'mass': 208.982, 'radius': 1.97, 'row': 6, 'column': 16,
               'group': 'Metalloid'},
        'At': {'name': 'Astatine', 'number': 85, 'mass': 209.987, 'radius': 2.02, 'row': 6, 'column': 17,
               'group': 'Halogens'},
        'Rn': {'name': 'Radon', 'number': 86, 'mass': 222.018, 'radius': 2.20, 'row': 6, 'column': 18,
               'group': 'Noble Gas'},
        'Fr': {'name': 'Francium', 'number': 87, 'mass': 223.020, 'radius': 3.48, 'row': 7, 'column': 1,
               'group': 'Alkali Metal'},
        'Ra': {'name': 'Radium', 'number': 88, 'mass': 226.025, 'radius': 2.83, 'row': 7, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Zn': {'name': 'Zinc', 'number': 30, 'mass': 65.38, 'radius': 1.39, 'row': 4, 'column': 12,
               'group': 'Transition Metal'}
    }


special_radii = {
    'ALA': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HB3': 1.3, 'N': 1.7, 'O': 1.49,
        'OC1': 1.5, 'OC2': 1.5
    },
    'ARG': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.88, 'CG': 1.92, 'CZ': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3,
        'HB2': 1.3, 'HB3': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HD3': 1.3, 'HE': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'HH11': 1.3,
        'HH12': 1.3, 'HH21': 1.3, 'HH22': 1.3, 'N': 1.7, 'NE': 1.62, 'NH1': 1.62, 'NH2': 1.67, 'O': 1.49
    },
    'ASN': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CG': 1.81, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD21': 1.3,
        'HD22': 1.3, 'N': 1.7, 'ND2': 1.62, 'O': 1.49, 'OD1': 1.52
    },
    'ASP': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CG': 1.76, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'N': 1.7, 'O': 1.49,
        'OD1': 1.49, 'OD2': 1.49
    },
    'CYS': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG': 1.3, 'N': 1.7, 'O': 1.49,
        'S': 1.88, 'SG': 1.88
    },
    'GLN': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.81, 'CG': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3,
        'HE21': 1.3, 'HE22': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'NE2': 1.62, 'O': 1.49, 'OE1': 1.52
    },
    'GLU': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.76, 'CG': 1.88, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3,
        'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'O': 1.49, 'OE1': 1.49, 'OE2': 1.49
    },
    'GLY': {
        'C': 1.75, 'CA': 1.9, 'H': 1.3, 'HA1': 1.3, 'HA2': 1.3, 'N': 1.7, 'O': 1.49, 'OC1': 1.5, 'OC2': 1.5
    },
    'HIS': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.74, 'CE': 1.74, 'CD2': 1.74, 'CE1': 1.74, 'CG': 1.8, 'H': 1.3,
        'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'N': 1.7, 'ND1': 1.6, 'ND2': 1.6,
        'NE2': 1.6, 'O': 1.49
    },
    'ILE': {
        'C': 1.75, 'CA': 1.9, 'CB': 2.01, 'CD': 1.92, 'CD1': 1.92, 'CG1': 1.92, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3,
        'HB': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HD3': 1.3, 'HD11': 1.3, 'HD12': 1.3, 'HD13': 1.3, 'HG12': 1.3, 'HG13': 1.3,
        'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.7, 'O': 1.49
    },
    'LEU': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD1': 1.92, 'CD2': 1.92, 'CG': 2.01, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3,
        'HB2': 1.3, 'HD11': 1.3, 'HD12': 1.3, 'HD13': 1.3, 'HD21': 1.3, 'HD22': 1.3, 'HD23': 1.3, 'HG': 1.3, 'N': 1.7,
        'O': 1.49
    },
    'LYS': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.92, 'CE': 1.88, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3,
        'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'HZ1': 1.3, 'HZ2': 1.3,
        'HZ3': 1.3, 'N': 1.7, 'NZ': 1.67, 'O': 1.49
    },
    'MET': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CE': 1.8, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3,
        'HE1': 1.3, 'HE2': 1.3, 'HE3': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'O': 1.49, 'SD': 1.94, 'S': 1.94
    },
    'PHE': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD1': 1.82, 'CD2': 1.82, 'CD': 1.82, 'CE1': 1.82, 'CE2': 1.82, 'CG': 1.74,
        'CZ': 1.82, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3,
        'HZ': 1.3, 'N': 1.7, 'O': 1.49
    },
    'PRO': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.92, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3,
        'HD1': 1.3, 'HD2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.7, 'O': 1.49
    },
    'SER': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'H1': 1.3, 'H2': 1.3, 'H3': 1.3, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3,
        'HG': 1.3, 'N': 1.7, 'O': 1.49, 'OG': 1.54
    },
    'THR': {
        'C': 1.75, 'CA': 1.9, 'CB': 2.01, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HG1': 1.3, 'HG21': 1.3,
        'HG22': 1.3, 'HG23': 1.3, 'N': 1.7, 'O': 1.49, 'OG1': 1.54, 'OG': 1.54
    },
    'TRP': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD1': 1.82, 'CD2': 1.82, 'CD': 1.82, 'CE': 1.82, 'CE2': 1.74, 'CE3': 1.82,
        'CG': 1.74, 'CH': 1.82, 'CH2': 1.82, 'CZ': 1.82, 'CZ1': 1.82, 'CZ2': 1.82, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3,
        'HB2': 1.3, 'HD1': 1.3, 'HE1': 1.3, 'HE3': 1.3, 'HZ1': 1.3, 'HZ2': 1.3, 'HH2': 1.3, 'N': 1.7, 'NE1': 1.66,
        'O': 1.49
    },
    'TYR': {
        'C': 1.75, 'CA': 1.9, 'CB': 1.91, 'CD': 1.82, 'CD1': 1.82, 'CD2': 1.82, 'CE': 1.82, 'CE1': 1.82, 'CE2': 1.82,
        'CG': 1.74, 'CZ': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3,
        'HE2': 1.3, 'HH': 1.3, 'N': 1.7, 'O': 1.49, 'OH': 1.54
    },
    'VAL': {
        'C': 1.75, 'CA': 1.9, 'CB': 2.01, 'CG1': 1.92, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HG11': 1.3,
        'HG12': 1.3, 'HG13': 1.3, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.7, 'O': 1.49
    }, 
        'DA': {  # Adenine
        'C1\'', 'C2', 'C2\'', 'C3\'', 'C4', 'C4\'', 'C5', 'C5\'', 'C6', 'C8', 'H1\'', 'H2', 'H2\'1', 'H2\'2', 'H3\'',
        'H4\'', 'H5\'1', 'H5\'2', 'H5T', 'H61', 'H62', 'H8', 'N1', 'N3', 'N6', 'N7', 'N9', 'O1P', 'O2P', 'O3\'',
        'O4\'', 'O5\'', 'P'
    },
    'DC': {  # Cytosine
        'C1\'', 'C2', 'C2\'', 'C3\'', 'C4', 'C4\'', 'C5', 'C5\'', 'C6', 'H1\'', 'H2\'1', 'H2\'2', 'H3\'', 'H4\'', 'H41',
        'H42', 'H5', 'H5\'1', 'H5\'2', 'H6', 'N1', 'N3', 'N4', 'O1P', 'O2', 'O2P', 'O3\'', 'O4\'', 'O5\'', 'P'
    },
    'DT': {  # Thymine
        'C1\'', 'C2', 'C2\'', 'C3\'', 'C4', 'C4\'', 'C5', 'C5\'', 'C6', 'C7', 'H1\'', 'H2\'1', 'H2\'2', 'H3', 'H3\'',
        'H3T', 'H4\'', 'H5\'1', 'H5\'2', 'H6', 'H71', 'H72', 'H73', 'N1', 'N3', 'O1P', 'O2', 'O2P', 'O3\'', 'O4',
        'O4\'', 'O5\'', 'P'
    },
    'DG': {  # Guanine
        'C1\'', 'C2', 'C2\'', 'C3\'', 'C4', 'C4\'', 'C5', 'C5\'', 'C6', 'C8', 'H1\'', 'H2\'1', 'H2\'2', 'H3\'', 'H4\'',
        'H5\'1', 'H5\'2', 'H8', 'N1', 'N2', 'N3', 'N7', 'N9', 'O1P', 'O2P', 'O3\'', 'O4\'', 'O5\'', 'O6', 'P'
    },
    'U': {   # Uracil
        'C1\'', 'C2', 'C2\'', 'C3\'', 'C4', 'C4\'', 'C5', 'C5\'', 'C6', 'H1\'', 'H2\'', 'H2\'\'', 'H3', 'H3\'', 'H4\'',
        'H5', 'H5\'', 'H5\'\'', 'H6', 'N1', 'N3', 'O1P', 'O2', 'O2\'', 'O2P', 'O3\'', 'O4', 'O4\'', 'O5\'', 'P'
    }
}
