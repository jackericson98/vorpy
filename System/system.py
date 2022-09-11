from System.read_write import *


class System:
    """Class used to import files of all types and return a System"""
    def __init__(self, file=None, box_size=None):
        self.atoms = []  # List of Atom type objects
        # If no file is given, generate a random System
        if file is None:
            self.random_system()
        # If the file type is a list
        if type(file) == list:
            self.build_sys(file)
            self.name = None
        else:
            self.name = get_name(file)
        # Grab the file
        self.file_address = file
        self.file_name = None
        self.file = None
        self.bonds = None
        self.Analysis = None  # Analysis type object for data collection
        # Non-pertinent information
        self.info = {}
        # Check the filetype and use the appropriate function to get it
        if file and len(file) > 4:
            if file[-3:] == "pdb":
                get_pdb(self)
            elif file[-3:] == "cif":
                get_cif(self)
            elif file[-3:] == "gro":
                get_gro(self)
            elif file[-3:] == "mol":
                get_mol(self)
        self.net = Network(self, self.atoms, box_size=box_size)
        self.mols = []
        self.output_directory = None
        self.box_size = box_size

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
                    self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], get_radius(self, line[1]),
                                      symbol=line[1], chain="None"))
                else:
                    self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], line[1],
                                           symbol=get_radius(line[1], return_symbol=True), chain="None"))

    # Random System function. Creates a System with atoms placed in random locations with random radii
    def random_system(self, anums=30, dmax=15, rmax=1):
        # Create the atoms
        for i in range(anums):
            # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
            self.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))

    # Build network function. Allows user to build the network from the system object.
    def build_network(self, get_verts=True, get_surfs=True, export_verts=False, directory=None, min_dist=None, box_size=None):
        if min_dist is not None:
            self.net.min_dist = min_dist
        if box_size is not None:
            self.net.box_size = box_size
        # Build the network
        self.net.build(get_verts=get_verts, get_surfs=get_surfs)
        # Export the vertices
        if export_verts:
            os.chdir(directory)
            file = open(os.getcwd() + "/Vertices.txt", 'w')
            for i in range(len(self.net.verts)):
                vert = self.net.verts[i]
                ndxs = [self.atoms.index(vert.atoms[i]) for i in range(4)]
                file.write(str(vert.loc[0]) + " " + str(vert.loc[1]) + " " + str(vert.loc[2]) + " " + str(vert.rad)
                           + " " + str(ndxs[0]) + " " + str(ndxs[1]) + " " + str(ndxs[2]) + " " + str(ndxs[3]) + '\n')

    # Analyze method. Finds the surface area of every surface in the system and volume of all the cells
    def analyze(self):
        # Run analysis on the network
        self.net.analyze()

    def add_verts(self, file_address):
        add_verts(self, file_address)

    def export_verts(self):
        os.chdir(self.output_directory)
        file = open(os.getcwd() + "/Vertices.txt", 'w')
        for i in range(len(self.net.verts)):
            vert = self.net.verts[i]
            ndxs = [self.atoms.index(vert.atoms[i]) for i in range(4)]
            file.write(str(vert.loc[0]) + " " + str(vert.loc[1]) + " " + str(vert.loc[2]) + " " + str(vert.rad)
                       + " " + str(ndxs[0]) + " " + str(ndxs[1]) + " " + str(ndxs[2]) + " " + str(ndxs[3]) + '\n')

    # Export method. Takes in an export type: 'Atoms', 'surfs'
    def export(self, directory=None, export_all=True, export_sys=False, export_atoms=False, export_mols=False,
               export_surfs=False, export_analysis=False, export_sys_pdb=False):
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
        self.output_directory = myDir
        os.chdir(myDir)
        # Go through all the possible user inputs and choose the correct export function
        n = 0
        lengths = [1, len(self.net.surfs), len(self.atoms), len(self.mols), 1]
        max_num_arr = [lengths[i] for i in range(len(lengths))
                       if [export_sys, export_surfs, export_atoms, export_mols, export_analysis][i] or export_all]
        max_num = sum(max_num_arr)
        # Export each
        if export_sys or export_all:
            export_mySys(self, n, max_num)
            n += 1
        if export_sys_pdb or export_all:
            create_pdb(self, os.getcwd())
        if export_surfs or export_all:
            export_mySurfs(self, n, max_num)
            n += len(self.net.surfs)
        if export_atoms or export_all:
            export_myAtoms(self, n, max_num)
            n += len(self.atoms)
        if export_mols or export_all:
            export_myMols(self, n, max_num)
            n += len(self.mols)
        if export_analysis or export_all:
            export_myAnalysis(self, n, max_num)
            n += 1
        print("\rExporting System:  ########## 100 %")
        print("\rSystem Exported")
