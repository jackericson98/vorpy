"""This file holds all object types needed for calculations: Molecule, Mesh, Sphere, Ray, Plane"""
import os
from System.Network.network import *


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


class Molecule:
    """Molecule object. Created from pdb import files or by user"""
    def __init__(self):
        self.atoms = []
        self.id = ""
        self.chain = ""


class System:
    """Class used to import files of all types and return a System"""
    def __init__(self, file=None, box_size=1.5, min_dist=1):
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
        self.name = ""
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
        self.min_dist = min_dist
        self.net = Network(self, self.atoms, box_size=box_size)
        self.mols = []

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
        # Reverse to normal and trim the extension and the dot
        return filename[::-1][:-4]

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
                    atom = Atom([float(line[5]), float(line[6]), float(line[7])], self.get_radius(line[-1]),
                                symbol=line[-1], res=line[2], chain=line[3], res_seq=line[4])
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
        keys = ['HEADER', 'TITLE', 'COMPOUND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK',
                'DBREF', 'SEQADV', 'FORMUL', 'SEQRES', 'HELIX', 'SHEET', 'CRYST', 'ORIG', 'SCALE', 'TER', 'HETATM',
                'MASTER']
        # Define the keys
        pdb_stds = ['header', 'title', 'compound', 'source', 'key_words', 'exp_data', 'author', 'revisions', 'journal',
                    'remarks', 'debrief', 'seq_adv', 'formula', 'residues', 'helix', 'sheet', 'crystal', 'origin',
                    'scale', 'terminals', 'het_atom', 'master']
        # Set the keys
        for i in range(len(pdb_stds)):
            self.info[keys[i]] = self.get_pdb_data(pdb_stds[i])
        # Grab the lines that start with ATOM and create Atom objects
        self.atoms = self.get_pdb_data('ATOM')

    # Create pdb method. Creates a pdb file type in the current working directory
    def create_pdb(self, directory=None):
        # Make sure each atom has a type before we change directories
        for atom in self.atoms:
            # Give each atom a type if not indicated
            if atom.type == "":
                atom.type = self.get_radius(atom.rad, return_symbol=True)
        # Move to the indicated directory
        if directory:
            os.chdir(directory)
        # Create the output file
        file = open(self.name + "_structure.pdb", 'w')
        # Go through each atom in the system
        for i in range(len(self.atoms)):
            a = self.atoms[i]
            loc = [str(round(a.loc[0], 3)), str(round(a.loc[1], 3)), str(round(a.loc[2], 3))]
            # Write the lines for the atom
            file.write("ATOM" + " " * (7 - len(str(i+1))) +
                       str(i + 1) + "  " +
                       a.type + " " * (4 - len(a.type)) +
                       a.res + " " * (4 - len(a.res)) +
                       a.chain + " " * (5 - len(a.chain) - len(a.res_seq)) +
                       " " * 4 + " " * (8 - len(loc[0])) +
                       loc[0] + " " * (8 - len(loc[1])) +
                       loc[1] + " " * (8 - len(loc[2])) +
                       loc[2] + " " * 2 +
                       "1.00  0.00" + " " * (12 - len(a.type)) + a.type + "\n")

    # Get gro method. Finds data in a gro file
    def get_gro(self):
        self.info['header'] = self.file[0]
        # Go through each line in the file and create an atom object
        for line in self.file[2:-2]:
            self.atoms.append(Atom([line[3], line[4], line[5]], self.get_radius(line[1][0]), symbol=line[1][0]))

    # Get mol method. Finds data in a mol file
    def get_mol(self):
        for line in self.file:
            if len(line) > 6:
                self.atoms.append(Atom([line[0], line[1], line[2]], self.get_radius(line[3]), symbol=line[3]))

    # Get radius Method. Goes through the bondi_radius file from voronota and gives a radius to the given atom name
    @staticmethod
    def get_radius(radius, return_symbol=False):
        # Get the classifier document
        radii = open("./Data/bondi_classifier.txt").readlines()
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
            myVert = Vertex(atoms=atoms, net=self.net, location=loc, radius=rad)
            self.net.verts.append(myVert)

    # Build network function. Allows user to build the network from the system object.
    def build_network(self, surfs=True):
        # Build the network
        self.net.build(surfs)

    # Analyze method. Finds the surface area of every surface in the system and volume of all the cells
    def analyze(self):
        # Run analysis on the network
        self.net.analyze()

    # Export function. Takes the system data and creates a set of obj files in the working directory for the surfaces
    def export(self, directory=None):
        # Set the name of the file to be created if no name exists
        if self.name is None:
            self.name = "mySystem"
        # Change to the directory indicated or create a directory called User_Data
        if directory:
            myDir = directory
        else:
            # If the system has a name set the data folder to it
            if len(self.name) > 0:
                myDir = os.getcwd() + "/" + self.name
                os.mkdir(myDir)
            else:
                myDir = os.getcwd() + "/User_Data"
                os.mkdir(myDir)
        # If the file is none create a pdb for the file
        if self.file is None:
            self.create_pdb(myDir)
        # Move to the new directory
        os.chdir(myDir)
        # Create a total number of operations for the percentage calculator
        tot_num = 2 * len(self.net.surfs) + len(self.atoms)
        # Set the counters to 0
        tot_verts, tot_tris = 0, 0
        # Get the total number of vertices and tris
        for surf in self.net.surfs:
            tot_verts += len(surf.points)
            tot_tris += len(surf.tris)

        # System file
        sys_file = open(self.name + "_System.off", 'w')
        sys_file.write("OFF\n" + str(tot_verts) + " " + str(tot_tris) + " 0\n\n\n")

        # Surfaces Folder
        os.mkdir(os.getcwd() + "/Surfaces")
        os.chdir("./Surfaces")
        # Go through each surface and create a file for each adding the vertex points
        surf_ndxs = []
        for i in range(len(self.net.surfs)):
            # Find the relative surface index and add it to the list
            surf_ndxs.append(str(self.atoms.index(self.net.surfs[i].atoms[0]) + 1) + "_" +
                             str(self.atoms.index(self.net.surfs[i].atoms[1]) + 1))
            # Name the file with the surface's index
            file = open(str("surf_" + surf_ndxs[i] + ".off"), 'w')
            # Start the file with "OFF", the number of points for the surface and the number of triangles
            file.write("OFF\n" +
                       str(len(self.net.surfs[i].points)) + " " + str(len(self.net.surfs[i].tris)) + " 0\n\n\n")
            # Go through the points on the surface
            for point in self.net.surfs[i].points:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(point[_], 4)) for _ in range(3)]
                sys_file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
            # Percentage printer
            per = int((i + 1) / tot_num * 100)
            print("\rExporting Files:   ", '#' * (per // 10) + ' ' * (10 - (per // 10)), per, "%", end='')
        num_verts = 0
        # Go through each surface opening the previously created file and add the faces
        for i in range(len(self.net.surfs)):
            file = open(str("surf_" + surf_ndxs[i] + ".off"), 'a')
            for tri in self.net.surfs[i].tris:
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                sys_file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " 1 0 0\n")
                file.write("3 " + str(tri[0]) + " " + str(tri[1]) + " " + str(tri[2]) + " 1 0 0\n")
            # Keep counting triangles for the system file
            num_verts += len(self.net.surfs[i].points)
            per = int((i + 1 + len(self.net.surfs)) / tot_num * 100)
            print("\rExporting Files:   ", '#' * (per // 10) + ' ' * (10 - (per // 10)), per, "%", end='')
        os.chdir("..")

        # Atoms Folder
        os.mkdir(os.getcwd() + "/Atoms")
        os.chdir("./Atoms")
        atom_info = open("atom_info.txt", 'w')
        # Add the vertices and triangles for each surface of each atom
        for i in range(len(self.atoms)):
            # Create a file for each atom
            atom_file = open(str("atom_" + str(i + 1) + "_cell.off"), 'w')
            # Calculate the number of points and triangles
            tot_verts, tot_tris = 0, 0
            for surf in self.atoms[i].surfs:
                tot_verts += len(surf.points)
                tot_tris += len(surf.tris)
            # Create the off header
            atom_file.write("OFF\n" + str(tot_verts) + " " + str(tot_tris) + " 0\n\n\n")
            # Go through each surface of the atom and add the vertices
            for j in range(len(self.atoms[i].surfs)):
                for point in self.atoms[i].surfs[j].points:
                    str_point = [str(round(point[_], 4)) for _ in range(3)]
                    atom_file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
            num_verts = 0
            # Go through each surface opening the previously created file and add the faces
            for j in range(len(self.atoms[i].surfs)):
                atom_file = open(str("atom_" + str(i + 1) + "_cell.off"), 'a')
                for tri in self.atoms[i].surfs[j].tris:
                    atom_file.write("3 " + str(tri[0] + num_verts) + " " + str(tri[1] + num_verts) + " " +
                                    str(tri[2] + num_verts) + " 1 0 0\n")
                num_verts += len(self.atoms[i].surfs[j].points)
            atom_info.write("Atom " + str(i + 1) + "\n")
            atom_info.write(" Cell Volume = {}\n".format(self.atoms[i].vol))
            atom_info.write(" Surfaces:\n")
            for j in range(len(self.atoms[i].surfs)):
                if self.atoms[i] == self.atoms[i].surfs[j].atoms[0]:
                    a1 = self.atoms[i].surfs[j].atoms[1]
                else:
                    a1 = self.atoms[i].surfs[j].atoms[0]
                atom_info.write("  Surface " + str(j + 1) + ", Made with Atom " + str(self.atoms.index(a1) + 1) +
                                ", Surface Area = " + str(self.atoms[i].surfs[j].sa) + "\n")
            atom_info.write("\n")
            per = int((i + 1 + 2 * len(self.net.surfs)) / tot_num * 100)
            print("\rExporting Files:   ", '#' * (per // 10) + ' ' * (10 - (per // 10)), per, "%", end='')

        os.chdir("..")
        os.mkdir(os.getcwd() + "/Molecules")
        os.chdir("./Molecules")

        print("\rExporting Files:    ########## 100 %")
        print("\rFiles Exported")
