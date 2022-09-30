from System.calcs import *


class Atom:
    """Atom object. Created with import of file. Used to reference for building network and analyzing"""
    def __init__(self, location, radius, system=None, symbol="", chain="", res="", res_seq=""):
        # Inherent traits
        self.loc = location     # Location    :   Set the location of the center of the sphere
        self.rad = radius       # Radius      :   Set the radius for the sphere object. Default is 1
        self.sys = system       # System      :   Set the atom's system attribute
        self.element = symbol   # Symbol      :   Element of the atom
        self.chain = chain      # Chain       :   Molecule chain the atom is a part of
        self.res = res          # Residue     :   Residue of the molecule that the atom is a part of
        self.res_seq = res_seq  # Sequence    :   Sequence of the residue that the atom is a part of
        self.box = []           # Box         :   The grid location of the atom
        # Network connections
        self.verts = []         # Vertices    :   List of Vertex type objects
        self.surfs = []         # Surfaces    :   List of Surface type objects
        self.edges = []         # Edges       :   List of Edge type objects
        # Calculated traits
        self.cell_vol = 0       # Cell Volume :   Volume of the voronoi cell for the atom


# Get radius function. Goes through the bondi_radius file from voronota and gives a radius to the given atom name
def get_radius(radius, return_symbol=False):
    # If indicated we return the symbol of atom that the radius indicates
    if return_symbol:
        # Set the atom type to nothing
        atom_type = ""
        # Check to see if the radius is in the system
        if radius in radii[1]:
            return radius[0][radius[1].index(radius)]
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


radii = [['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p', 's', 'cl',
          'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as',
          'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in',
          'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb',
          'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl',
          'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk',
          'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'nh',
          'fl', 'mc', 'lv', 'ts', 'og'],
         [1.2, 1.4, 1.81, 1.53, 1.92, 1.7, 1.55, 1.52, 1.47, 1.54, 1.55, 1.73, 1.84, 2.1, 1.8, 1.8, 1.75,
          1.88, 2.75, 2.31, None, None, None, None, None, None, None, None, None, None, 1.87, 2.11, 1.85,
          1.90, 1.83, 2.02, 3.03, 2.49, None, None, None, None, None, None, None, None, None, None, 1.93,
          2.17, 2.06, 2.06, 1.98, 2.16, 3.43, 2.68, None, None, None, None, None, None, None, None, None,
          None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 1.96,
          2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, None, None, None, None, None, None, None, None, None,
          None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
          None, None, None, None, None]]