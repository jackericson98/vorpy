from System.read_write import *


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
        if file and len(file) > 4:
            if file[-3:] == "pdb":
                get_pdb(self)
            elif file[-3:] == "gro":
                get_gro(self)
            elif file[-3:] == "mol":
                get_mol(self)
        self.min_dist = min_dist
        self.net = Network(self, self.atoms, box_size=box_size)
        self.mols = []
        self.dir = ''

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
                    self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], get_radius(self, line[1])))
                else:
                    self.atoms.append(Atom([line[0][0], line[0][1], line[0][2]], line[1]))

    # Random System function. Creates a System with atoms placed in random locations with random radii
    def random_system(self, anums=30, dmax=15, rmax=1):
        # Create the atoms
        for i in range(anums):
            # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
            self.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))

    # Build network function. Allows user to build the network from the system object.
    def build_network(self, surfs=True):
        # Build the network
        self.net.build(surfs)

    # Analyze method. Finds the surface area of every surface in the system and volume of all the cells
    def analyze(self):
        # Run analysis on the network
        self.net.analyze()

    # Export method. Takes in an export type: 'Atoms', 'surfs'
    def export(self, directory=None, export_type=None):
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
        self.dir = myDir
        # Go through all the possible user inputs and choose the correct export function
        if export_type is None or export_type.lower() == 'all':
            export_sys(self)
            export_surfs(self)
            export_atoms(self)
            export_mols(self)
            export_analysis(self)
        # Go through the export options
        elif export_type.lower()[:5] == 'atoms':
            export_atoms(self)
        elif export_type.lower()[:5] == 'surfs':
            export_atoms(self)
        elif export_type.lower()[:3] == 'sys':
            export_sys(self)
        elif export_type.lower()[:4] == 'mols':
            export_mols(self)
        elif export_type.lower() == 'analysis':
            export_analysis(self)
