from System.calcs import *


class Atom:
    """
    Atom object class used to represent loaded atoms

    location: list
        set the location of the center of the sphere
    radius : float
        set the radius for the sphere object. Default is 1
    system : System object
        set the atom's system attribute
    element : str
        element of the atom
    chain : str
        molecule chain the atom is a part of
    residue : str
        residue of the molecule that the atom is a part of
    res_seq : int
        sequence of the residue that the atom is a part of
    name : str
        name retrieved from pdb file
    ocp : str
        Occupancy of the atom
    t_fact : str
        Temperature factor for the atom
    seg_id : str
        Segment identifier for the atom
    charge : float
        Charge of the atom
    verts : list
        Vertex objects connected to the atom
    surfs : list
        Surface objects connected to the atom
    edges : list
        Edge objects connected to the atom
    load_ndxs : list
        Holds the object indices for when the system is loaded back in
    cell_vol : float
        Volume of the voronoi cell for the atom
    box : list
        The grid location of the atom

    """
    def __init__(self, location=None, radius=None, system=None, element="", chain="", residue="", res_seq="", name="",
                 ocp="", t_fact="", seg_id="", charge=""):

        # Inherent traits
        self.loc = location     # Location     :   Set the location of the center of the sphere
        self.rad = radius       # Radius       :   Set the radius for the sphere object. Default is 1
        self.sys = system       # System       :   Set the atom's system attribute
        self.element = element  # Symbol       :   Element of the atom
        self.chain = chain      # Chain        :   Molecule chain the atom is a part of
        self.res = residue      # Residue      :   Residue of the molecule that the atom is a part of
        self.res_seq = res_seq  # Sequence     :   Sequence of the residue that the atom is a part of
        self.name = name        # Name         :   Name retrieved from pdb file
        self.occupancy = ocp    # Occupancy    :   Occupancy of the atom
        self.t_fact = t_fact    # Temp Factor  :   Temperature factor for the atom
        self.seg_id = seg_id    # Segment ID   :   Segment identifier for the atom
        self.charge = charge    # Charge       :   Charge of the atom

        # Network connections
        self.verts = []         # Vertices     :   List of Vertex type objects
        self.surfs = []         # Surfaces     :   List of Surface type objects
        self.edges = []         # Edges        :   List of Edge type objects
        self.load_ndxs = []     # Load indices :   Holds the object indices for when the system is loaded back in

        # Calculated traits
        self.cell_vol = 0       # Cell Volume  :   Volume of the voronoi cell for the atom
        self.box = []           # Box          :   The grid location of the atom


def get_radius(radius, return_symbol=False):
    """
        Finds the radius of the atom from the symbol or vice versa

    :param radius: Either the elemental symbol for the atom or it's radius
    :param return_symbol: Boolean for whether to return the symbol or not
    :return: The radius of the atom from the symbol or vice versa
    """

    # If indicated we return the symbol of atom that the radius indicates
    if return_symbol:
        # Set the atom type to nothing
        atom_type = ""
        # Check to see if the radius is in the system
        if radius in radii[1]:
            return radii[0][radii[1].index(radius)]
        else:
            # Get the closest atom to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for i in range(len(radii[1])):
                if radii[1][i] is not None and radii[1][i] - radius < min_diff:
                    atom_type = radii[0][i]
    # If we have the type and just want the radius, keep scanning until we find the radius
    else:
        radius = radius.strip()
        return radii[1][radii[0].index(radius.lower())]
    # If nothing is found to be exact return the closest atom type
    return atom_type


radii = [['h' , 'he', 'li', 'be', 'b' , 'c' , 'n' , 'o' , 'f' , 'ne', 'na', 'mg', 'al', 'si', 'p' , 's' , 'cl', 'ar',
          'k' , 'ca', 'sc', 'ti', 'v' , 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
          'rb', 'sr', 'y' , 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i' , 'xe',
          'cs', 'ba', 'la', 'hf', 'ta', 'w' , 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
          'fr', 'ra', 'ac', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og',
          'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu',
          'th', 'pa', 'u' , 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr'],
         [1.30, 1.40, 0.76, 0.45, 1.92, 1.80, 1.60, 1.50, 1.33, 1.54, 1.02, 0.72, 0.60, 2.10, 1.90, 1.90, 1.81, 1.88,
          1.38, 1.00, None, None, None, None, None, None, None, None, None, None, 0.62, 0.73, 0.58, 1.90, 1.83, 2.02,
          1.52, 1.18, None, None, None, None, None, None, None, None, None, None, 1.93, 2.17, 2.06, 2.06, 2.20, 2.16,
          1.67, 1.35, None, None, None, None, None, None, None, None, None, None, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20,
          3.48, 2.83, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
          None, None, None, None, None, None, None, None, None, None, None, None, None, None,
          None, None, None, None, None, None, None, None, None, None, None, None, None, None]]
