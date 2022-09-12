from System.calcs import *


class Atom:
    """Atom object. Created with import of file. Used to reference for building network and analyzing"""
    def __init__(self, location, radius, symbol="", chain="", res="", res_seq=""):
        self.loc = location  # Set the location of the center of the sphere
        self.rad = radius  # Set the radius for the sphere object. Default is 1
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.cell = True
        self.vol = 0
        self.type = symbol
        self.chain = chain
        self.res = res
        self.res_seq = res_seq
        self.box = []


# Get radius Method. Goes through the bondi_radius file from voronota and gives a radius to the given atom name
def get_radius(radius, return_symbol=False):
    # Get the classifier document
    radii = open(os.getcwd() + "/Data/bondi_classifier.txt").readlines()
    atom_type = ""
    min_diff = np.inf
    # Go through each line in the classifier document to find the radius or symbol for the atom
    for line in radii:
        # Split the line
        line = line.split()
        # If the line is empty, continue
        if len(line) == 0:
            continue
        # If indicated we return the symbol of atom that the radius indicates
        if return_symbol:
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
            if radius.lower() == line[1].lower():
                return float(line[2])
    # If nothing is found to be exact return the closest atom type
    return atom_type
