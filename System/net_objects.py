"""Network objects file"""


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, atoms):
        self.atoms = atoms  # List of Atom type objects
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.rad = 50  # Ballpark range for radius needed for the entire network.


class Atom:
    """Atom object. Created with import of file. Used to reference for building network and analyzing"""
    def __init__(self, location, radius):
        self.rad = radius  # Set the radius for the sphere object. Default is 1
        self.loc = location  # Set the location of the center of the sphere
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects


class Vertex:
    """Vertex object. Used to build the network and calculate the surfaces"""
    def __init__(self, location, radius, atoms=None):
        self.loc = location  # Location of the vertex
        self.rad = radius  # Radius of the vertex's tangential sphere
        self.atoms = atoms  # List of Atom type objects
        self.edges = []  # List of Edge type objects
        self.surfs = []  # List of Surface type objects


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms, verts, surfs=[]):
        self.atoms = atoms  # List of Atom type objects
        self.verts = verts  # List of Vertex type objects
        self.surfs = surfs
        self.loc = None
        self.rad = None
        self.dir = None
        self.points = []  # List of points on the edge. These points do not include the vertex points


class Surface:
    """Surface object. Holds the mesh data. Used to analyze."""
    def __init__(self, atoms, func, edges=None, verts=None):
        self.func = func
        self.atoms = atoms  # List of Atom type objects
        self.edges = edges  # List of Edge type objects
        self.edge_points = []
        self.verts = verts
        self.vert_points = []
        self.points = []  # List of points on the surface
