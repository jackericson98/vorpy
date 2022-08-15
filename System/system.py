"""This file holds all object types needed for calculations: Molecule, Mesh, Sphere, Ray, Plane"""
from System.Network.network import *
from System.analysis import *


class Atom:
    """Atom object. Created with import of file. Used to reference for building network and analyzing"""
    def __init__(self, location, radius):
        self.loc = location  # Set the location of the center of the sphere
        self.rad = radius  # Set the radius for the sphere object. Default is 1
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.cell = True
        self.cell_vol = 0


class System:
    """Class used to import files of all types and return a System"""
    def __init__(self, file=None, box_size=1.5):
        self.atoms = []  # List of Atom type objects
        # If no file is given, generate a random System
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
        self.box = self.calc_box(box_size)
        self.bonds = None
        self.Analysis = None  # Analysis type object for data collection
        # Non-pertinent information
        self.info = {}
        # Check the filetype and use the appropriate function to get it
        if file:
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

    # Calculate box function. Takes in a System and returns the dimensions of a box x times the size of the atoms
    def calc_box(self, x):
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = np.array([np.inf, np.inf, np.inf])
        max_vert = np.array([-np.inf, -np.inf, -np.inf])
        # Check each atom in the System
        for atom in self.atoms:
            # Go through x, y, z
            for i in range(3):
                # If we find that the x, y, z value is less replace the value in the mins list
                if atom.loc[i] < min_vert[i]:
                    min_vert[i] = atom.loc[i]
                # If we find that the x, y, z value is less replace the value in the mins list
                elif atom.loc[i] > max_vert[i]:
                    max_vert[i] = atom.loc[i]
        # Get the vector between the minimum and maximum vertices for the defining box
        r_box = max_vert - min_vert
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert + r_box * x, min_vert - r_box * x
        # Return the list of array turned list vertices
        return [min_vert.tolist(), max_vert.tolist()]

    # Get pdb data method. Finds the lines of the file with prefixes and returns them as a list
    def get_pdb_data(self, word):
        # Get the file information
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
                    atom = Atom([float(line[5]), float(line[6]), float(line[7])], self.get_radius(line[-1]))
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
        # .PDB file type standards.
        pdb_stds = ['HEADER', 'TITLE', 'COMPOUND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK',
                    'DBREF', 'SEQADV', 'FORMUL', 'SEQRES', 'HELIX', 'SHEET', 'CRYST', 'ORIG', 'SCALE', 'TER', 'HETATM',
                    'MASTER']
        # Define the keys
        keys = ['header', 'title', 'compound', 'source', 'key_words', 'exp_data', 'author', 'revisions', 'journal',
                'remarks', 'debrief', 'seq_adv', 'formula', 'residues', 'helix', 'sheet', 'crystal', 'origin', 'scale',
                'terminals', 'het_atom', 'master']
        # Set the keys
        for i in range(len(pdb_stds)):
            self.info[keys[i]] = self.get_pdb_data(pdb_stds[i])
        # Grab the lines that start with ATOM and create Atom objects
        self.atoms = self.get_pdb_data('ATOM')

    # Get gro method. Finds data in a gro file
    def get_gro(self):
        self.info['header'] = self.file[0]
        self.box = self.file[-2]
        # Go through each line in the file and create an atom object
        for line in self.file[2:-2]:
            self.atoms.append(Atom([line[3], line[4], line[5]], self.get_radius(line[1][0])))

    # Get mol method. Finds data in a mol file
    def get_mol(self):
        for line in self.file:
            if len(line) > 6:
                self.atoms.append(Atom([line[0], line[1], line[2]], self.get_radius(line[3])))

    # Get radius Method. Goes through the bondi_radius file from voronota and gives a radius to the given atom name
    @staticmethod
    def get_radius(atom_name):
        # Get the classifier document
        radii = open("./Data/bondi_classifier.txt").readlines()
        # Go through each line in the classifier document
        for line in radii:
            line = line.split()
            # Compare the given atom name and atom name in the line
            if atom_name.lower() == line[1].lower():
                # Get the classifier for the line (0, 1, 2, 3, 4, 5, 6, 7)
                return float(line[2])

    # Build System function. Takes in a list of coordinates and string atom names
    def build_sys(self, lr_input_list):
        # Check if atom objects are given
        if type(lr_input_list[0]) is Atom:
            self.atoms = lr_input_list
        else:
            # Go through each line in the input list
            for line in lr_input_list:
                # If the radius is a string, convert the radius using the get_radius method
                if type(line[1]) == str:
                    self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], self.get_radius(line[1])))
                else:
                    self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], line[1]))

    # Random System function. Creates a System with atoms placed in random locations with random radii
    def random_system(self, anums=30, dmax=15, rmax=1):
        # Create the atoms
        for i in range(anums):
            # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
            self.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))

    # Add Voronota data method. Takes in voronota data and adds it to the System
    def add_vta_data(self, ball_file, vert_file):
        # Set the voronota system indicator to True
        self.net.vta = True
        # Create the System and load the files
        vert_file = open(vert_file).readlines()
        ball_file = open(ball_file).readlines()
        # Interpret the balls
        balls = []
        for i in range(len(ball_file)):
            # Split the data
            data = ball_file[i].split(" ")
            # Grab the data reference for the atoms
            balls.append(self.atoms[int(data[5])])
        # Interpret the vertices
        for i in range(len(vert_file)):
            # Split the data
            data = vert_file[i].split(" ")
            # Add the vertex data
            loc, rad = [float(data[4]), float(data[5]), float(data[6])], float(data[7])
            atoms = [balls[int(data[0])], balls[int(data[1])], balls[int(data[2])], balls[int(data[3])]]
            myVert = Vertex(atoms=atoms, location=loc, radius=rad)
            self.net.verts.append(myVert)

    # Build network function. Allows user to build the network from the system object.
    def build_network(self, min_dist=0.1, surfs=True):
        # Build the network
        self.net.build(min_dist, surfs)

    # Analyze method. Finds the surface area of every surface in the system and volume of all of the cells
    def analyze(self):
        # Run analysis on the network
        analyze(self)

    # Export function. Takes the system data and creates a set of obj files in the working directory for the surfaces
    def export(self, directory=None):
        # Change to the directory indicated
        if directory:
            os.chdir(directory)
        # Set the name of the file to be created if no name exists
        if self.name is None:
            self.name = "mySystem"
        # Go through each surface and create a file for each adding the vertex points
        for i in range(len(self.net.surfs)):
            #
            file = open(str(self.name + "_surf_" + str(i) + ".obj"), 'w')
            for point in self.net.surfs[i].points:
                file.write("v " + str(round(point[0], 3)) + " " + str(round(point[1], 3)) + " " + str(round(point[2], 3)) + '\n')
        # Go through each surface opening the previously created file and add the faces
        for i in range(len(self.net.surfs)):
            file = open(str(self.name + "_surf_" + str(i) + ".obj"), 'a')
            for tri in self.net.surfs[i].tris:
                file.write("vf " + str(tri[0]) + "// " + str(tri[1]) + "// " + str(tri[2]) + "//\n")

