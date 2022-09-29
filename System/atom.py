from System.calcs import *


class Atom:
    """Atom object. Created with import of file. Used to reference for building network and analyzing"""
    def __init__(self, location, radius, symbol="", chain="", res="", res_seq=""):
        # Inherent traits
        self.loc = location     # Location    :   Set the location of the center of the sphere
        self.rad = radius       # Radius      :   Set the radius for the sphere object. Default is 1
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
    # Get the classifier document
    radii = open(os.getcwd() + "/Data/bondi_classifier.txt").readlines()
    atom_type = ""
    min_diff = np.inf
    # If indicated we return the symbol of atom that the radius indicates
    if return_symbol:
        # Go through each line in the classifier document to find the radius or symbol for the atom
        for line in radii:
            # Split the line
            line = line.split()
            # If the line is empty, continue
            if len(line) == 0:
                continue
            # If we get the exact radius, return it
            if line[2] == float(radius):
                return line[1]
            # Find the difference between the bondi classifier line's radius and the atom's
            new_min = abs(float(radius) - float(line[2]))
            # If the check radius is closer to the actual radius update the symbol and the minimum difference
            if new_min < min_diff:
                min_diff = new_min
                atom_type = line[1]
    # If we have the type and just want the radius, keep scanning until we find the radius
    else:
        # Go through each line in the classifier document to find the radius or symbol for the atom
        for line in radii:
            # Strip the white space around the letter(s)
            radius = radius.strip()
            # If we have found the radius return it
            if radius.lower() == line[1].lower():
                return float(line[2])
    # If nothing is found to be exact return the closest atom type
    return atom_type
