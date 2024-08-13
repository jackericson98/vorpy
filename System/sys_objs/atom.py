import numpy as np


def get_radius(atom):
    """
    Finds the radius of the ball from the symbol or vice versa
    :return: The radius of the ball from the symbol or vice versa
    """
    # Get the radius and the element from the name of the ball
    if atom['res'] is not None and atom['res'].name in special_radii:
        # Check if no ball name exists or its empty
        if atom['name'] is not None and atom['name'] != '':
            for i in range(len(atom['name'])):
                name = atom['name'][:-i]
                # Check the residue name
                if name in special_radii[atom['res'].name]:
                    atom['rad'] = special_radii[atom['res'].name][name]
    # If we have the type and just want the radius, keep scanning until we find the radius
    if atom['rad'] is None and atom['element'].lower() in element_radii:
        atom['rad'] = element_radii[atom['element'].lower()]
    # If indicated we return the symbol of ball that the radius indicates
    if atom['rad'] is None or atom['rad'] == 0:
        # Check to see if the radius is in the system
        if atom['rad'] in {element_radii[_] for _ in element_radii[1]}:
            atom['element'] = element_radii[atom['rad']]
        else:
            # Get the closest ball to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in element_radii:
                if element_radii[radius] - atom['rad'] < min_diff:
                    atom['element'] = element_radii[radius]
    return atom['rad']


class Atom:
    def __init__(self, system=None, location=None, radius=None, index='', name='', residue='', chain='', res_seq="",
                 seg_id="", element="", chn=None, res=None):

        # System groups
        self.sys = system           # System       :   Main system object
        self.res = res              # Residue      :   Residue object of which the atom is a part
        self.chn = chn              # Chain        :   Chain object of which the atom is a part

        self.loc = location         # Location     :   Set the location of the center of the sphere
        self.rad = radius           # Radius       :   Set the radius for the sphere object. Default is 1

        # Calculated Traits
        self.vol = 0                # Cell Volume  :   Volume of the voronoi cell for the atom
        self.sa = 0                 # Surface Area :   Surface area of the atom's cell
        self.curv = 0
        self.box = []               # Box          :   The grid location of the atom

        # Network objects
        self.verts = []             # Vertices     :   List of Vertex type objects
        self.surfs = []             # Surfaces     :   List of Surface type objects
        self.edges = []             # Edges        :   List of Edge type objects

        # Input traits
        self.num = index            # Number       :   The index from the initial atom file
        self.name = name            # Name         :   Name retrieved from pdb file
        self.chain = chain          # Chain        :   Molecule chain the atom is a part of
        self.residue = residue      # Residue      :   Class of molecule that the atom is a part of
        self.res_seq = res_seq      # Sequence     :   Sequence of the residue that the atom is a part of
        self.seg_id = seg_id        # Segment ID   :   Segment identifier for the atom
        self.element = element      # Symbol       :   Element of the atom

        self.rad = get_radius(self)


def make_atom(system=None, location=None, radius=None, index='', name='', residue='', chain='', chn_name='',
              res_name='', res_seq="", seg_id="", element="", chn=None, res=None):
    atom = {
        # System groups
        'sys': system,           # System       :   Main system object
        'res': res,              # Residue      :   Residue object of which the atom is a part
        'chn': chn,              # Chain        :   Chain object of which the atom is a part

        'loc': location,         # Location     :   Set the location of the center of the sphere
        'rad': radius,           # Radius       :   Set the radius for the sphere object. Default is 1

        # Calculated Traits
        'vol': 0,                # Cell Volume  :   Volume of the voronoi cell for the atom
        'sa': 0,                 # Surface Area :   Surface area of the atom's cell
        'curv': 0,
        'box': [],               # Box          :   The grid location of the atom

        # Network objects
        'averts': [],             # Vertices     :   List of Vertex type objects
        'asurfs': [],             # Surfaces     :   List of Surface type objects
        'aedges': [],             # Edges        :   List of Edge type objects

        # Input traits
        'num': index,            # Number       :   The index from the initial atom file
        'name': name,            # Name         :   Name retrieved from pdb file
        'chain': chain,          # Chain        :   Molecule chain the atom is a part of
        'chain_name': chn_name,
        'residue': residue,      # Residue      :   Class of molecule that the atom is a part of
        'res_name': res_name,
        'res_seq': res_seq,      # Sequence     :   Sequence of the residue that the atom is a part of
        'seg_id': seg_id,        # Segment ID   :   Segment identifier for the atom
        'element': element,      # Symbol       :   Element of the atom
    }
    if atom['rad'] is None:
        atom['rad'] = get_radius(atom)
    return atom


element_radii = {'h': 1.30, 'he': 1.40, 'li': 0.76, 'be': 0.45, 'b': 1.92, 'c': 1.80, 'n': 1.60, 'o': 1.50, 'f': 1.33,
                 'ne': 1.54, 'na': 1.02, 'mg': 0.72, 'al': 0.60, 'si': 2.10, 'p': 1.90, 's': 1.90, 'cl': 1.81, 'ar': 1.88,
                 'k': 1.38, 'ca': 1.00, 'ga': 0.62, 'ge': 0.73, 'as': 0.58, 'se': 1.90, 'br': 1.83, 'kr': 2.02, 'rb': 1.52,
                 'sr': 1.18, 'in': 1.93, 'sn': 2.17, 'sb': 2.06, 'te': 2.06, 'i': 2.20, 'xe': 2.16, 'cs': 1.67, 'ba': 1.35,
                 'tl': 1.96, 'pb': 2.02, 'bi': 2.07, 'po': 1.97, 'at': 2.02, 'rn': 2.20, 'fr': 3.48, 'ra': 2.83, '': 1.80,
                 'zn': 1.39}
special_radii = {''   : {'C': 1.75, 'CA': 1.90, 'N': 1.70, 'O': 1.49, 'F': 1.33, 'CL': 1.81, 'BR': 1.96, 'I': 2.20},
                 'ALA': {'CB': 1.92},
                 'ARB': {'CB': 1.91, 'CD': 1.88, 'CG': 1.92, 'CZ': 1.80, 'NE': 1.62, 'NH1': 1.62, 'NH2': 1.67},
                 'ASN': {'CB': 1.91, 'CG': 1.81, 'ND2': 1.62, 'OD1': 1.52},
                 'ASP': {'CB': 1.91, 'CG': 1.76, 'OD1': 1.49, 'OD2': 1.49},
                 'CYS': {'CB': 1.91, 'S': 1.88},
                 'GLN': {'CB': 1.91, 'CD': 1.81, 'CG': 1.80, 'NE2': 1.62, 'OE1': 1.52},
                 'GLU': {'CB': 1.91, 'CD': 1.76, 'CG': 1.88, 'OE1': 1.49, 'OE2': 1.49},
                 'HIS': {'CB': 1.91, 'CD': 1.74, 'CE': 1.74, 'CG': 1.80, 'ND1': 1.60, 'ND2': 1.60},
                 'ILE': {'CB': 2.01, 'CD1': 1.92, 'CG1': 1.92, 'CG2': 1.92},
                 'LEU': {'CB': 1.91, 'CD1': 1.92, 'CD2': 1.92, 'CG': 2.01},
                 'LYS': {'CB': 1.91, 'CD': 1.92, 'CE': 1.88, 'CG': 1.92, 'NZ': 1.67},
                 'MET': {'CB': 1.91, 'CE': 1.80, 'CG': 1.92, 'S': 1.94},
                 'PHE': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CG': 1.74, 'CZ': 1.82},
                 'PRO': {'CB': 1.91, 'CD': 1.92, 'CG': 1.92},
                 'SER': {'CB': 1.91, 'OG': 1.54},
                 'THR': {'CB': 2.01, 'CG2': 1.92, 'OG': 1.54},
                 'TRP': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CE2': 1.74, 'CG': 1.74, 'CH': 1.82, 'CZ': 1.82, 'NE1': 1.66},
                 'TYR': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CG': 1.74, 'CZ': 1.80, 'OH': 1.54},
                 'VAL': {'CB': 2.01, 'CG1': 1.92, 'CG2': 1.92}}