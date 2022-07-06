"""This file holds all object types needed for calculations: Molecule, Mesh, Sphere, Ray, Plane"""
import numpy as np
import os


class System:
    """Class used to import files of all types and return a system"""
    def __init__(self, file=None):
        self.atoms = []  # List of Atom type objects
        # If no file is given, generate a random system
        if file is None:
            self.random_system()
        # If the file type is a list
        if type(file) == list:
            self.build_sys(file)
        # Grab the file
        self.file_address = file
        self.file_name = None
        self.file = None

        # Set up our
        self.name = None
        self.box = None
        self.bonds = None
        self.Analysis = None  # Analysis type object for data collection
        # Non-pertinent information
        self.header = None
        self.title = None
        self.compound = None
        self.source = None
        self.key_words = None
        self.exp_data = None
        self.author = None
        self.revisions = None
        self.journal = None
        self.remarks = None
        self.debrief = None
        self.seq_adv = None
        self.formula = None
        self.residues = None
        self.helix = None
        self.sheet = None
        self.crystal = None
        self.origin = None
        self.scale = None
        self.terminals = None
        self.het_atom = None
        self.master = None
        # Check the filetype and use the appropriate function to get it
        if file[-3:] == "pdb":
            self.get_pdb()
        elif file[-3:] == "gro":
            self.get_gro()
        elif file[-3:] == "mol":
            self.get_mol()
        self.net = Network(self.atoms)

    # Get name method. Extracts the name from the file name
    @staticmethod
    def get_name(file):
        if not file:
            return
        filename = ""
        i = -1
        # Go through each char in the path from the back and stop at the first slash
        while file[i] != "/":
            filename = filename + file[i]
            i -= 1
        # Trim the extension and the dot
        return filename[::-1][:-4]

    # Get pdb data method. Finds the lines of the file with prefixes and returns them as a list
    def get_pdb_data(self, word):

        self.file = open(self.file_address).readlines()
        self.name = self.get_name(self.file_address)
        # Split each line in the file
        for i in range(len(self.file)):
            self.file[i] = self.file[i].split()
        # Special case for Atom lines
        if word.lower() == 'atom':
            atoms = []
        # Go through each line in the file and check if the first word is the word we are looking for
            for i in range(len(self.file)):
                line = self.file[i]
                if line and line[0].lower() == 'atom':  # Check if the line starts with atom
                    atom = Atom([float(line[-7]), float(line[-5]), float(line[-4])], self.get_radius(line[-1]))
                    atoms.append(atom)
            return atoms

        # Standard case
        else:
            data = []
            for i in range(len(self.file)):
                if word == self.file[i][:len(word)]:  # check the first len(word) letters
                    # Add the split test_data to our list and remove the word at the beginning of the list
                    data.append(self.file[i].split()[1:])
        return data

    # Get pdb method. Finds the atoms and
    def get_pdb(self):
        # Open and read the file
        self.header = self.get_pdb_data('HEADER')
        self.title = self.get_pdb_data('TITLE')
        self.compound = self.get_pdb_data('COMPOUND')
        self.source = self.get_pdb_data('SOURCE')
        self.key_words = self.get_pdb_data('KEYWDS')
        self.exp_data = self.get_pdb_data('EXPDTA')
        self.author = self.get_pdb_data('AUTHOR')
        self.revisions = self.get_pdb_data('REVDAT')
        self.journal = self.get_pdb_data('JRNL')
        self.remarks = self.get_pdb_data('REMARK')
        self.debrief = self.get_pdb_data('DBREF')
        self.seq_adv = self.get_pdb_data('SEQADV')
        self.formula = self.get_pdb_data('FORMUL')
        self.residues = self.get_pdb_data('SEQRES')
        self.helix = self.get_pdb_data('HELIX')
        self.sheet = self.get_pdb_data('SHEET')
        self.crystal = self.get_pdb_data('CRYST')
        self.origin = self.get_pdb_data('ORIG')
        self.scale = self.get_pdb_data('SCALE')
        self.terminals = self.get_pdb_data('TER')
        self.het_atom = self.get_pdb_data('HETATM')
        self.master = self.get_pdb_data('MASTER')
        # Grab all the lines that start with ATOM. Creates Atom objects
        self.atoms = self.get_pdb_data('ATOM')

    # Get gro method. Finds data in a gro file
    def get_gro(self):
        self.header = self.file[0]
        self.box = self.file[-2]
        # Go through each line in the file and
        for line in self.file[2:-2]:
            atom = Atom([line[3], line[4], line[5]], self.get_radius(line[1][0]))
            if atom.rad is None:
                print(line[1])
            self.atoms.append(atom)

    # Get mol method. Finds data in a mol file
    def get_mol(self):
        for line in self.file:
            if len(line) > 6:
                self.atoms.append(Atom([line[0], line[1], line[2]], self.get_radius(line[3])))

    # Get radius Method. Goes through the bondi_radius file from voronota and gives a radius to the given atom name
    @staticmethod
    def get_radius(atom_name):
        # Get the classifier document
        radii = open(os.getcwd() + "/Data/bondi_classifier.txt").readlines()
        # Go through each line in the classifier document
        for line in radii:
            line = line.split()
            # Compare the given atom name and atom name in the line
            if atom_name.lower() == line[1].lower():
                # Get the classifier for the line (0, 1, 2, 3, 4, 5, 6, 7)
                return float(line[2])

    # Build system function. Takes in a list of coordinates and string atom names
    def build_sys(self, lr_input_list):
        # Go through each line in the input list
        for line in lr_input_list:
            # If the radius is a
            if type(line[1]) == str:
                self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], self.get_radius(line[1])))
            else:
                self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], line[1]))

    # Random system function. Creates a system with atoms placed in random locations with random radii
    def random_system(self, anums=30, dmax=15, rmax=1):
        # Create the atoms
        for i in range(anums):
            # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
            self.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))


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


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms, verts):
        self.atoms = atoms  # List of Atom type objects
        self.verts = verts  # List of Vertex type objects
        self.center = None
        self.points = []  # List of points on the edge. These points do not include the vertex points


class Surface:
    """Surface object. Holds the mesh data. Used to analyze."""
    def __init__(self, atoms, func):
        self.func = func
        self.atoms = atoms  # List of Atom type objects
        self.edges = []  # List of Edge type objects
        self.edge_points = []
        self.verts = []
        self.vert_points = []
        self.points = []  # List of points on the surface
