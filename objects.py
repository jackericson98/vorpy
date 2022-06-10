"""This file holds all object types needed for calculations: Molecule, Mesh, Sphere, Ray, Plane"""


class System:
    """System object. Holds everything from the import file. Used to build network for and analyze import file"""
    def __init__(self):
        self.info = {}  # Information about the system
        self.atoms = []  # List of Atom type objects
        self.net = Network(self.atoms)  # Network type object for calculations
        self.Analysis = None  # Analysis type object for data collection


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, atoms):
        self.atoms = atoms  # List of Atom type objects
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects


class Atom:
    """Atom object. Created with import of file. Used to reference for building network and analyzing"""
    def __init__(self, location, radius):
        self.rad = radius  # Set the radius for the sphere object. Default is 1
        self.loc = location  # Set the location of the center of the sphere
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects


class Vertex:
    """Voronoi S-Vertex object. Created when building the network."""
    def __init__(self, location, radius, e0, atoms):
        self.loc = location  # Location of the vertex
        self.rad = radius  # Radius of the vertex's tangential sphere
        self.atoms = atoms  # List of Atom type objects
        self.edges = [e0]  # List of Edge type objects


class Edge:
    """Voronoi S-Channel"""
    def __init__(self, atoms, v0):
        self.atoms = atoms  # List of Atom type objects
        self.verts = [v0]  # List of Vertex type objects
        self.points = []  # List of points on the edge


class Surface:
    """Defines surface object"""
    def __init__(self, func, atoms):
        self.func = func
        self.atoms = atoms  # List of Atom type objects
        self.edges = []  # List of Edge type objects
        self.points = []  # List of points on the surface

