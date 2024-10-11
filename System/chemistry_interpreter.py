
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
    **{_: 'U' for _ in {'du', 'u', 'uracil', 'uraci', 'ura'}},                          # Uracil
    # Ions
    **{_: 'NA' for _ in {'na', 'sodium', 'natrium'}},                                   # Sodium
    **{_: 'K' for _ in {'k', 'potassium', 'kalium'}},                                   # Potassium
    **{_: 'CA' for _ in {'ca', 'calcium'}},                                             # Calcium
    **{_: 'MG' for _ in {'mg', 'magnesium'}},                                           # Magnesium
    **{_: 'MN' for _ in {'mn', 'manganese'}},                                           # Manganese
    **{_: 'FE' for _ in {'fe', 'iron', 'ferrum'}},                                      # Iron
    **{_: 'CO' for _ in {'co', 'cobalt'}},                                              # Cobalt
    **{_: 'NI' for _ in {'ni', 'nickel'}},                                              # Nickel
    **{_: 'CU' for _ in {'cu', 'copper', 'cuprum'}},                                    # Copper
    **{_: 'ZN' for _ in {'zn', 'zinc'}},                                                # Zinc
    **{_: 'CD' for _ in {'cd', 'cadmium'}},                                             # Cadmium
    **{_: 'HG' for _ in {'hg', 'mercury', 'hydrargyrum'}},                              # Mercury
    **{_: 'PB' for _ in {'pb', 'lead', 'plumbum'}},                                     # Lead
    **{_: 'CL' for _ in {'cl', 'chloride', 'chlorine'}},                                # Chloride
    **{_: 'BR' for _ in {'br', 'bromide', 'bromine'}},                                  # Bromide
    **{_: 'I' for _ in {'i', 'iodide', 'iodine'}},                                      # Iodide
    **{_: 'SO4' for _ in {'so4', 'sulfate', 'sulphate'}},                               # Sulfate
    **{_: 'PO4' for _ in {'po4', 'phosphate'}},                                         # Phosphate
    **{_: 'C2H3O2' for _ in {'acetate', 'c2h3o2', 'ethanoate'}},                        # Acetate
    **{_: 'NO3' for _ in {'no3', 'nitrate'}},                                           # Nitrate
    **{_: 'NH4' for _ in {'nh4', 'ammonium'}},                                          # Ammonium
    **{_: 'H3O' for _ in {'h3o', 'hydronium'}},                                         # Hydronium
    **{_: 'MOO4' for _ in {'moo4', 'molybdate'}},                                       # Molybdate
    **{_: 'SEO4' for _ in {'seo4', 'selenate'}},                                        # Selenate
    **{_: 'VO4' for _ in {'vo4', 'vanadate'}},                                          # Vanadate
    **{_: 'WO4' for _ in {'wo4', 'tungstate'}}                                          # Tungstate
}

residue_atoms = {
    'SOL': {'HW1', 'HW2', 'OW'},
    'NA':      {'NA'},                    # Sodium
    'K':       {'K'},                     # Potassium
    'CA':      {'CA'},                    # Calcium
    'MG':      {'MG'},                    # Magnesium
    'MN':      {'MN'},                    # Manganese
    'FE':      {'FE'},                    # Iron
    'CO':      {'CO'},                    # Cobalt
    'NI':      {'NI'},                    # Nickel
    'CU':      {'CU'},                    # Copper
    'ZN':      {'ZN'},                    # Zinc
    'CD':      {'CD'},                    # Cadmium
    'HG':      {'HG'},                    # Mercury
    'PB':      {'PB'},                    # Lead
    'CL':      {'CL'},                    # Chloride
    'BR':      {'BR'},                    # Bromide
    'I':       {'I'},                     # Iodide
    'SO4':     {'S', 'O'},                # Sulfate (SO4)
    'PO4':     {'P', 'O'},                # Phosphate (PO4)
    'C2H3O2':  {'C', 'H', 'O'},           # Acetate (C2H3O2)
    'NO3':     {'N', 'O'},                # Nitrate (NO3)
    'NH4':     {'N', 'H'},                # Ammonium
    'H3O':     {'H', 'O'},                # Hydronium
    'MoO4':    {'MO', 'O'},               # Molybdate (MoO4)
    'SEO4':    {'SE', 'O'},               # Selenate (SeO4)
    'VO4':     {'V', 'O'},                # Vanadate (VO4)
    'WO4':     {'W', 'O'},                # Tungstate (WO4)
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


element_names = {
    **{_: 'H' for _ in {'h', 'hydrogen'}},
    **{_: 'HE' for _ in {'he', 'helium'}},
    **{_: 'LI' for _ in {'li', 'lithium'}},
    **{_: 'BE' for _ in {'be', 'beryllium'}},
    **{_: 'B' for _ in {'b', 'boron'}},
    **{_: 'C' for _ in {'c', 'carbon'}},
    **{_: 'N' for _ in {'n', 'nitrogen'}},
    **{_: 'O' for _ in {'o', 'oxygen'}},
    **{_: 'F' for _ in {'f', 'fluorine'}},
    **{_: 'NE' for _ in {'ne', 'neon'}},
    **{_: 'NA' for _ in {'na', 'sodium', 'natrium', 'sod'}},
    **{_: 'MG' for _ in {'mg', 'magnesium'}},
    **{_: 'AL' for _ in {'al', 'aluminium', 'aluminum'}},
    **{_: 'SI' for _ in {'si', 'silicon'}},
    **{_: 'P' for _ in {'p', 'phosphorus'}},
    **{_: 'S' for _ in {'s', 'sulfur', 'sulphur'}},
    **{_: 'CL' for _ in {'cl', 'chlorine'}},
    **{_: 'AR' for _ in {'ar', 'argon'}},
    **{_: 'K' for _ in {'k', 'potassium', 'kalium'}},
    **{_: 'CA' for _ in {'ca', 'calcium'}},
    **{_: 'SC' for _ in {'sc', 'scandium'}},
    **{_: 'TI' for _ in {'ti', 'titanium'}},
    **{_: 'V' for _ in {'v', 'vanadium'}},
    **{_: 'CR' for _ in {'cr', 'chromium'}},
    **{_: 'MN' for _ in {'mn', 'manganese'}},
    **{_: 'FE' for _ in {'fe', 'iron', 'ferrum'}},
    **{_: 'CO' for _ in {'co', 'cobalt'}},
    **{_: 'NI' for _ in {'ni', 'nickel'}},
    **{_: 'CU' for _ in {'cu', 'copper', 'cuprum'}},
    **{_: 'ZN' for _ in {'zn', 'zinc'}},
    **{_: 'GA' for _ in {'ga', 'gallium'}},
    **{_: 'GE' for _ in {'ge', 'germanium'}},
    **{_: 'AS' for _ in {'as', 'arsenic'}},
    **{_: 'SE' for _ in {'se', 'selenium'}},
    **{_: 'BR' for _ in {'br', 'bromine'}},
    **{_: 'KR' for _ in {'kr', 'krypton'}},
    **{_: 'RB' for _ in {'rb', 'rubidium'}},
    **{_: 'SR' for _ in {'sr', 'strontium'}},
    **{_: 'Y' for _ in {'y', 'yttrium'}},
    **{_: 'ZR' for _ in {'zr', 'zirconium'}},
    **{_: 'NB' for _ in {'nb', 'niobium', 'columbium'}},
    **{_: 'MO' for _ in {'mo', 'molybdenum'}},
    **{_: 'TC' for _ in {'tc', 'technetium'}},
    **{_: 'RU' for _ in {'ru', 'ruthenium'}},
    **{_: 'RH' for _ in {'rh', 'rhodium'}},
    **{_: 'PD' for _ in {'pd', 'palladium'}},
    **{_: 'AG' for _ in {'ag', 'silver', 'argentum'}},
    **{_: 'CD' for _ in {'cd', 'cadmium'}},
    **{_: 'IN' for _ in {'in', 'indium'}},
    **{_: 'SN' for _ in {'sn', 'tin', 'stannum'}},
    **{_: 'SB' for _ in {'sb', 'antimony', 'stibium'}},
    **{_: 'TE' for _ in {'te', 'tellurium'}},
    **{_: 'I' for _ in {'i', 'iodine'}},
    **{_: 'XE' for _ in {'xe', 'xenon'}},
    **{_: 'CS' for _ in {'cs', 'cesium', 'caesium'}},
    **{_: 'BA' for _ in {'ba', 'barium'}},
    **{_: 'LA' for _ in {'la', 'lanthanum'}},
    **{_: 'CE' for _ in {'ce', 'cerium'}},
    **{_: 'PR' for _ in {'pr', 'praseodymium'}},
    **{_: 'ND' for _ in {'nd', 'neodymium'}},
    **{_: 'PM' for _ in {'pm', 'promethium'}},
    **{_: 'SM' for _ in {'sm', 'samarium'}},
    **{_: 'EU' for _ in {'eu', 'europium'}},
    **{_: 'GD' for _ in {'gd', 'gadolinium'}},
    **{_: 'TB' for _ in {'tb', 'terbium'}},
    **{_: 'DY' for _ in {'dy', 'dysprosium'}},
    **{_: 'HO' for _ in {'ho', 'holmium'}},
    **{_: 'ER' for _ in {'er', 'erbium'}},
    **{_: 'TM' for _ in {'tm', 'thulium'}},
    **{_: 'YB' for _ in {'yb', 'ytterbium'}},
    **{_: 'LU' for _ in {'lu', 'lutetium'}},
    **{_: 'HF' for _ in {'hf', 'hafnium'}},
    **{_: 'TA' for _ in {'ta', 'tantalum'}},
    **{_: 'W' for _ in {'w', 'tungsten', 'wolfram'}},
    **{_: 'RE' for _ in {'re', 'rhenium'}},
    **{_: 'OS' for _ in {'os', 'osmium'}},
    **{_: 'IR' for _ in {'ir', 'iridium'}},
    **{_: 'PT' for _ in {'pt', 'platinum'}},
    **{_: 'AU' for _ in {'au', 'gold', 'aurum'}},
    **{_: 'HG' for _ in {'hg', 'mercury', 'hydrargyrum'}},
    **{_: 'TL' for _ in {'tl', 'thallium'}},
    **{_: 'PB' for _ in {'pb', 'lead', 'plumbum'}},
    **{_: 'BI' for _ in {'bi', 'bismuth'}},
    **{_: 'PO' for _ in {'po', 'polonium'}},
    **{_: 'AT' for _ in {'at', 'astatine'}},
    **{_: 'RN' for _ in {'rn', 'radon'}},
    **{_: 'FR' for _ in {'fr', 'francium'}},
    **{_: 'RA' for _ in {'ra', 'radium'}},
    **{_: 'AC' for _ in {'ac', 'actinium'}},
    **{_: 'TH' for _ in {'th', 'thorium'}},
    **{_: 'PA' for _ in {'pa', 'protactinium'}},
    **{_: 'U' for _ in {'u', 'uranium'}},
    **{_: 'NP' for _ in {'np', 'neptunium'}},
    **{_: 'PU' for _ in {'pu', 'plutonium'}},
    **{_: 'AM' for _ in {'am', 'americium'}},
    **{_: 'CM' for _ in {'cm', 'curium'}},
    **{_: 'BK' for _ in {'bk', 'berkelium'}},
    **{_: 'CF' for _ in {'cf', 'californium'}},
    **{_: 'ES' for _ in {'es', 'einsteinium'}},
    **{_: 'FM' for _ in {'fm', 'fermium'}},
    **{_: 'MD' for _ in {'md', 'mendelevium'}},
    **{_: 'NO' for _ in {'no', 'nobelium'}},
    **{_: 'LR' for _ in {'lr', 'lawrencium'}},
    **{_: 'RF' for _ in {'rf', 'rutherfordium'}},
    **{_: 'DB' for _ in {'db', 'dubnium'}},
    **{_: 'SG' for _ in {'sg', 'seaborgium'}},
    **{_: 'BH' for _ in {'bh', 'bohrium'}},
    **{_: 'HS' for _ in {'hs', 'hassium'}},
    **{_: 'MT' for _ in {'mt', 'meitnerium'}},
    **{_: 'DS' for _ in {'ds', 'darmstadtium'}},
    **{_: 'RG' for _ in {'rg', 'roentgenium'}},
    **{_: 'CN' for _ in {'cn', 'copernicium'}},
    **{_: 'NH' for _ in {'nh', 'nihonium'}},
    **{_: 'FL' for _ in {'fl', 'flerovium'}},
    **{_: 'MC' for _ in {'mc', 'moscovium'}},
    **{_: 'LV' for _ in {'lv', 'livermorium'}},
    **{_: 'TS' for _ in {'ts', 'tennessine'}},
    **{_: 'OG' for _ in {'og', 'oganesson'}}
}

my_masses = {'h': 1.008, 'he': 4.003, 'li': 6.941, 'be': 9.012, 'b': 10.811, 'c': 12.011, 'n': 14.007, 'o': 15.999,
             'f': 18.998, 'ne': 20.180, 'na': 22.990, 'mg': 24.305, 'al': 26.982, 'si': 28.086, 'p': 30.974,
             's': 32.066, 'cl': 35.453, 'ar': 39.948, 'k': 39.098, 'ca': 40.078, 'ga': 69.723, 'ge': 72.631,
             'as': 74.922, 'se': 78.971, 'br': 79.904, 'kr': 83.798, 'rb': 85.468, 'sr': 87.62, 'in': 114.818,
             'sn': 118.711, 'sb': 121.760, 'te': 27.6, 'i': 126.904, 'xe': 131.293, 'cs': 132.905, 'ba': 137.328,
             'tl': 204.383, 'pb': 207.2, 'bi': 208.980, 'po': 208.982, 'at': 209.987, 'rn': 222.018, 'fr': 223.020,
             'ra': 226.025, '': 1.80, 'W': 4.1}

