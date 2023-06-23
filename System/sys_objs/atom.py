from System.sys_funcs.calcs.calcs import *


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

        get_radius(self)


def make_atom(system=None, location=None, radius=None, index='', name='', residue='', chain='', res_seq="", seg_id="",
              element="", chn=None, res=None):
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
        'residue': residue,      # Residue      :   Class of molecule that the atom is a part of
        'res_seq': res_seq,      # Sequence     :   Sequence of the residue that the atom is a part of
        'seg_id': seg_id,        # Segment ID   :   Segment identifier for the atom
        'element': element,      # Symbol       :   Element of the atom
    }
    return atom
