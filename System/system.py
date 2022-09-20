from System.input_system import *
from System.output_system import *


class System:
    """Class used to import files of all types and return a System"""
    def __init__(self, file=None, user_atoms=None, box_size=1.5, min_dist=0.1):

        # Set up the major system attributes
        self.info = {}
        self.atoms = []
        self.mols = []
        self.min_dist = min_dist

        # Set up the file attributes
        self.name = None
        self.file = None
        self.file_address = None

        # If a file is given read the file and set the system attributes
        if file:
            self.file = open(file).readlines()
            self.file_address = file
            self.name = get_name(file)
            # Check the file type
            if file[-3:] == "pdb":
                get_pdb(self)
            elif file[-3:] == "cif":
                get_cif(self)
            elif file[-3:] == "gro":
                get_gro(self)
            elif file[-3:] == "mol":
                get_mol(self)
            else:
                return
        # If no file is given and the user has entered atoms build system
        elif user_atoms:
            self.name = "User_Atoms"
            self.build_sys(user_atoms)
        # If no file is given, generate a random System
        else:
            self.random_system()

        # Set up the network
        self.net = Network(self, self.atoms, box_size=box_size, min_dist=min_dist)
        self.output_directory = None

    # Build System method. Takes in a list of atomic values
    def build_sys(self, user_atoms):
        # Check if the user entered Atoms into their list
        if type(user_atoms[0]) is Atom:
            self.atoms = user_atoms
            return
        # Go through each line in the input list
        for line in user_atoms:
            # If the radius is a string, convert the radius using the get_radius method
            if type(line[1]) == str:
                self.atoms.append(Atom([float(line[0][0]), float(line[0][1]), float(line[0][2])],
                                       get_radius(self, line[1]), symbol=line[1], chain="None"))
            else:
                self.atoms.append(Atom([float(line[0][0]), float(line[0][1]), float(line[0][2])], float(line[1]),
                                       symbol=get_radius(line[1], return_symbol=True), chain="None"))

    # Random System function. Creates a System with atoms placed in random locations with random radii
    def random_system(self, anums=30, dmax=15, rmax=1):
        # Create the atoms
        for i in range(anums):
            # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
            self.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))

    # Build network function. Allows user to build the network from the system object.
    def build_network(self, get_verts=True, get_surfs=True, export_verts=False, min_dist=None, box_size=None):
        if min_dist is not None:
            self.min_dist = min_dist
        if box_size is not None:
            self.net.box_size = box_size
        # Build the network
        self.net.build(get_verts=get_verts, get_surfs=get_surfs)
        # Export the vertices
        if export_verts:
            export_myVerts(self)

    # Analyze method. Finds the surface area of every surface in the system and volume of all the cells
    def analyze(self):
        # Run analysis on the network
        self.net.analyze()

    # System level add vertex method. Just a pass through to the input system file
    def add_verts(self, file_address):
        add_grant_verts(self, file_address)

    # System level Export vertices method.
    def export_verts(self):
        if self.output_directory is None:
            set_output_dir(self)
        export_myVerts(self)


    # Export method. Takes in an export type: 'Atoms', 'surfs'
    def export(self, export_all=True, export_sys=False, export_atoms=False, export_mols=False,
               export_surfs=False, export_sys_pdb=False, export_reses=False):
        if self.output_directory is None:
            set_output_dir(self)
        os.chdir(self.output_directory)
        # Go through all the possible user inputs and choose the correct export function
        n = 0
        lengths = [1, len(self.net.surfs), len(self.atoms), len(self.mols)]
        exports = [export_sys, export_surfs, export_atoms, export_mols]
        # Find the total number of things being exported for the loading bar
        max_num_arr = [lengths[i] for i in range(len(lengths)) if exports[i] or export_all]
        max_num = sum(max_num_arr)
        # Export the system
        if export_sys or export_all:
            export_mySys(self, n, max_num)
            n += 1
        # Export the pdb of the system
        if export_sys_pdb or export_all:
            create_pdb(self, os.getcwd())
        # Export the individual surfaces
        if export_surfs or export_all:
            export_mySurfs(self, n, max_num)
            n += len(self.net.surfs)
        # Export the atom cells
        if export_atoms or export_all:
            export_myAtoms(self, n, max_num)
            n += len(self.atoms)
        # Export the molecules
        if export_mols or export_all:
            export_myMols(self, n, max_num, export_residues=export_reses or export_all)
            n += len(self.mols)
        print("\rExporting System:  ########## 100 %")
        print("\rSystem Exported")
