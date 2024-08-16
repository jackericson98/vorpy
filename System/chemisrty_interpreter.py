
"""
This file helps interpret different inputs for chemical classifiers
1. Residue names: lowercase general classifiers -> three letter amino acid codes


"""

residue_names = {
    **{_: 'ARG' for _ in {'r', 'arginine', 'arg', 'argi', 'argin', 'arganine'}},        # Arginine
    **{_: 'ALA' for _ in {'alanine', 'ala', 'alan'}},                                   # Alanine
    **{_: 'ASN' for _ in {'n', 'asparagine', 'asn', 'aspar', 'asparagin'}},             # Asparagine
    **{_: 'ASP' for _ in {'d', 'aspartic acid', 'asp', 'aspart', 'aspartate'}},         # Aspartic acid
    **{_: 'CYS' for _ in {'cysteine', 'cys', 'cyst'}},                                  # Cysteine
    **{_: 'GLU' for _ in {'e', 'glutamic acid', 'glu', 'glut', 'glutamate'}},           # Glutamic acid
    **{_: 'GLN' for _ in {'q', 'glutamine', 'gln', 'glutamin'}},                        # Glutamine
    **{_: 'GLY' for _ in {'glycine', 'gly', 'glycin'}},                                 # Glycine
    **{_: 'HIS' for _ in {'h', 'histidine', 'his', 'hist'}},                            # Histidine
    **{_: 'ILE' for _ in {'i', 'isoleucine', 'ile', 'isol'}},                           # Isoleucine
    **{_: 'LEU' for _ in {'l', 'leucine', 'leu', 'leuc'}},                              # Leucine
    **{_: 'LYS' for _ in {'k', 'lysine', 'lys', 'lysin'}},                              # Lysine
    **{_: 'MET' for _ in {'m', 'methionine', 'met', 'meth'}},                           # Methionine
    **{_: 'PHE' for _ in {'f', 'phenylalanine', 'phe', 'phenyl'}},                      # Phenylalanine
    **{_: 'PRO' for _ in {'p', 'proline', 'pro', 'prolin'}},                            # Proline
    **{_: 'SER' for _ in {'s', 'serine', 'ser', 'serin'}},                              # Serine
    **{_: 'THR' for _ in {'threonine', 'thr', 'threon'}},                               # Threonine
    **{_: 'TRP' for _ in {'w', 'tryptophan', 'trp', 'trypto'}},                         # Tryptophan
    **{_: 'TYR' for _ in {'y', 'tyrosine', 'tyr', 'tyros'}},                            # Tyrosine
    **{_: 'VAL' for _ in {'v', 'valine', 'val', 'valin'}},                              # Valine
    # Nucleo bases
    **{_: 'A' for _ in {'da', 'a', 'adenine', 'adenin', 'ade'}},                        # Adenine
    **{_: 'C' for _ in {'dc', 'c', 'cytosine', 'cytosin', 'cyto'}},                     # Cytosine
    **{_: 'G' for _ in {'dg', 'g', 'guanine', 'guanin', 'guan'}},                       # Guanine
    **{_: 'T' for _ in {'dt', 't', 'thymine', 'thymi', 'thym'}},                        # Thymine
    **{_: 'U' for _ in {'du', 'u', 'uracil', 'uraci', 'ura'}}                           # Uracil
}

residue_atoms = {
    'SOL': {'HW1', 'HW2', 'OW'},
    'NA': {'NA'},
    'ALA': {'C', 'CA', 'CB', 'H', 'H1', 'H2', 'H3', 'HA', 'HB1', 'HB2', 'HB3', 'N', 'O', 'OC1', 'OC2'},
    'ARG': {'C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'H', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE', 'HG1', 'HG2', 'HH11', 'HH12', 'HH21', 'HH22', 'N', 'NE', 'NH1', 'NH2', 'O'},
    'THR': {'C', 'CA', 'CB', 'CG2', 'H', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23', 'N', 'O', 'OG1'},
    'LYS': {'C', 'CA', 'CB', 'CD', 'CE', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HG1', 'HG2', 'HZ1', 'HZ2', 'HZ3', 'N', 'NZ', 'O', 'OC1', 'OC2'},
    'GLN': {'C', 'CA', 'CB', 'CD', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HE21', 'HE22', 'HG1', 'HG2', 'N', 'NE2', 'O', 'OE1'},
    'SER': {'C', 'CA', 'CB', 'H', 'H1', 'H2', 'H3', 'HA', 'HB1', 'HB2', 'HG', 'N', 'O', 'OG'},
    'GLY': {'C', 'CA', 'H', 'HA1', 'HA2', 'N', 'O', 'OC1', 'OC2'},
    'PRO': {'C', 'CA', 'CB', 'CD', 'CG', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HG1', 'HG2', 'N', 'O'},
    'LEU': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23', 'HG', 'N', 'O'},
    'VAL': {'C', 'CA', 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'N', 'O'},
    'HIS': {'C', 'CA', 'CB', 'CD2', 'CE1', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HD2', 'HE1', 'HE2', 'N', 'ND1', 'NE2', 'O'},
    'TYR': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'H', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HH', 'N', 'O', 'OH'},
    'GLU': {'C', 'CA', 'CB', 'CD', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HG1', 'HG2', 'N', 'O', 'OE1', 'OE2'},
    'ILE': {'C', 'CA', 'CB', 'CD', 'CG1', 'CG2', 'H', 'HA', 'HB', 'HD1', 'HD2', 'HD3', 'HG11', 'HG12', 'HG21', 'HG22', 'HG23', 'N', 'O'},
    'PHE': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'H', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ', 'N', 'O'},
    'ASP': {'C', 'CA', 'CB', 'CG', 'H', 'HA', 'HB1', 'HB2', 'N', 'O', 'OD1', 'OD2'},
    'MET': {'C', 'CA', 'CB', 'CE', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HE1', 'HE2', 'HE3', 'HG1', 'HG2', 'N', 'O', 'SD'},
    'ASN': {'C', 'CA', 'CB', 'CG', 'H', 'HA', 'HB1', 'HB2', 'HD21', 'HD22', 'N', 'ND2', 'O', 'OD1'},
    'CYS': {'C', 'CA', 'CB', 'H', 'HA', 'HB1', 'HB2', 'HG', 'N', 'O', 'SG'},
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


atoms = {
    **{_: 'CA' for _ in {'ca'}},
    **{_: 'CB' for _ in {'cb'}}
}


