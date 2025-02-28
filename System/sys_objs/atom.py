import numpy as np
from radii.Standard import element_radii, special_radii


def get_element(atom):
    """
    Discerns the element from the atom name and the atom information
    """
    # Simple SOL
    if atom['res_name'].lower() in {'sol', 'hoh'}:
        if 'O' in atom['name']:
            return 'O'
        elif 'H' in atom['name']:
            return 'H'

    # Find the simple atom types:
    for name in ['C', 'O', 'H', 'N', 'P', 'S']:
        if atom['name'][0] == name:
            return name

    # Print the unknown atoms
    print(f'Unknown atom name {atom['name']}')

    # Otherwise just return h
    return 'H'


def get_radius(atom, my_radii=None):
    """
    Finds the radius of the ball from the symbol or vice versa
    :return: The radius of the ball from the symbol or vice versa
    """
    if my_radii is None:
        elements_radii, specials_radii = element_radii, special_radii
    else:
        elements_radii, specials_radii = my_radii['elements'], my_radii['specials']
    # Get the radius and the element from the name of the ball
    if atom['res'] is not None and atom['res'].name in specials_radii:
        # Check if no ball name exists or its empty
        if atom['name'] is not None and atom['name'] != '':
            for i in range(len(atom['name'])):
                name = atom['name'][:-i]
                # Check the residue name
                if name in specials_radii[atom['res'].name]:
                    atom['rad'] = specials_radii[atom['res'].name][name]
    # If we have the type and just want the radius, keep scanning until we find the radius
    if atom['rad'] is None and atom['element'].upper() in elements_radii:
        atom['rad'] = elements_radii[atom['element'].upper()]
    # If indicated we return the symbol of ball that the radius indicates
    if atom['rad'] is None or atom['rad'] == 0:
        # Check to see if the radius is in the system
        if atom['rad'] in {elements_radii[_] for _ in elements_radii[1]}:
            atom['element'] = elements_radii[atom['rad']]
        else:
            # Get the closest ball to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in elements_radii:
                if elements_radii[radius] - atom['rad'] < min_diff:
                    atom['element'] = elements_radii[radius]
    return atom['rad']


class Atom:
    def __init__(self, system=None, location=None, radius=None, index='', name='', residue='', chain='', res_seq="",
                 seg_id="", element="", chn=None, res=None, radii=None):

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

        self.rad = get_radius(self, my_radii=radii)


def make_atom(system=None, location=None, radius=None, index='', name='ball', residue='', chain='', chn_name='',
              res_name='', res_seq="0", seg_id="0", element=None, chn=None, res=None, mass=1.0, occ_choice='.',
              chn_id='', pdb_ins_code='', occupancy="", b_factor="", charge="", auth_seq_id="", auth_comp_id="",
              auth_asym_id="", auth_atom_id="", pdbx_PDB_model_num=""):
    atom = {
        # System groups
        'sys': system,           # System       :   Main system object

        'num': index,            # Number       :   The index from the initial atom file
        'loc': location,         # Location     :   Set the location of the center of the sphere
        'rad': radius,           # Radius       :   Set the radius for the sphere object. Default is 1

        # Calculated Traits
        'vol': 0,                # Cell Volume  :   Volume of the voronoi cell for the atom
        'sa': 0,                 # Surface Area :   Surface area of the atom's cell
        'curv': 0,
        'box': [],               # Box          :   The grid location of the atom

        # Network objects
        'verts': [],             # Vertices     :   List of Vertex type objects
        'surfs': [],             # Surfaces     :   List of Surface type objects
        'edges': [],             # Edges        :   List of Edge type objects

        # Molecule traits
        'name': name,            # Name         :   Name retrieved from pdb file
        'res': res,              # Residue      :   Residue object of which the atom is a part
        'chn': chn,              # Chain        :   Chain object of which the atom is a part
        'chain': chain,          # Chain        :   Molecule chain the atom is a part of
        'chain_name': chn_name,  # Chain Name   :   Name of the chain that the ball is a part of
        'chain_id': chn_id,
        'residue': residue,      # Residue      :   Class of molecule that the atom is a part of
        'res_name': res_name,    # Residue Name :   Name of the residue the ball is a part of
        'res_seq': res_seq,      # Sequence     :   Sequence of the residue that the atom is a part of
        'seg_id': seg_id,        # Segment ID   :   Segment identifier for the atom
        'element': element,      # Symbol       :   Element of the atom
        'mass': mass,             # Mass         :   Mass of the atom

        # Other Values
        'occupancy': occupancy,
        'occ_choice': occ_choice,
        'pdb_ins_code': pdb_ins_code,
        'b_factor': b_factor,
        'charge': charge,
        'auth_seq_id': auth_seq_id,
        'auth_comp_id': auth_comp_id,
        'auth_asym_id': auth_asym_id,
        'auth_atom_id': auth_atom_id,
        'pdbx_PDB_model_num': pdbx_PDB_model_num

    }
    if atom['element'] is None:
        atom['element'] = get_element(atom)
    if atom['rad'] is None:
        atom['rad'] = get_radius(atom)
    return atom
