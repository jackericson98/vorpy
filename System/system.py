import time

from System.input_system import *
from System.output_system import *


class System:
    """Class used to import files of all types and return a System"""
    def __init__(self, atoms=None, mols=None, sol=None, residues=None, data=None, name=None, base_file=None,
                 frame_files=None, nets=None, net_files=None, ndx_files=None, gui=None):

        self.atoms = atoms                    # Atoms            :    List holding the atom objects
        self.mols = mols                      # Molecules        :    List of molecules
        self.sol = sol                        # Solution         :    List of solution molecules (lists of atoms)
        self.residues = residues              # Residues         :    List of residues (lists of atoms)
        # Set up the file attributes
        self.name = name                      # Name             :    Name describing the system
        self.data = data                      # Data             :    Additional data provided by the base file
        self.net = Network(self, self.atoms)  # Network          :    Network object holding the primary network
        self.base_file = base_file            # Base file        :    Primary file address
        self.nets = nets                      # Networks         :    Different networks for other frames
        self.net_files = net_files            # Network files    :    Network files for multiple frames
        self.ndx_files = ndx_files            # Index files      :    File addresses for index file in GROMACS format
        self.frame_files = frame_files        # Frame files      :    File addresses for different frames (.gro,.pdb)
        self.output_directory = None          # Output Directory :    Output directory for the export files
        self.vorpy_directory = os.getcwd()    # Vorpy Directory  :    Directory that vorpy is running out of

        self.gui = gui                       # GUI              :    GUI Vorpy object that can be updated through sys

    # Load network method. Used to load a network that was previously calculated
    def load_net(self, net_file):
        # If no file has been loaded before, create the main network
        if self.nets is None:
            self.net_files = [net_file]
            import_net(self.net, net_file)
        else:
            self.net_files.append(net_file)
            self.nets.append(Network(self, self.atoms))
            import_net(self.nets[-1], net_file)

    # Build System method. Takes in a list of atomic values
    def build_user_atoms_sys(self, user_atoms):
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

    # Sort atoms method. Used to put atoms in their correct molecules and residues
    def sort_atoms(self):
        # Set up the chain names list
        self.mols, chain_names = [], []
        # Go through each of the atoms in the system adding the atoms to their respective chains
        for atom in self.atoms:
            # If no chain is specified, set the chain to 'None'
            if atom.chain == ' ':
                if atom.res.lower == 'sol':
                    atom.chain = 'SOL'
                else:
                    atom.chain = 'Mol'
            # If the atom's chain does not exist add it to the list of chains
            if atom.chain not in chain_names:
                self.mols.append([atom])
                chain_names.append(atom.chain)
            else:
                self.mols[chain_names.index(atom.chain)].append(atom)
        # Set up the residues names list
        self.residues, res_names = [], []
        # Set up the residues
        for atom in self.atoms:
            # Get the residue name for the atom
            res_name = [atom.res, atom.res_seq]
            # If the residue name does not exist, add it
            if res_name not in res_names:
                self.residues.append([atom])
                res_names.append(res_name)
            else:
                self.residues[res_names.index(res_name)].append(atom)

    # Load system method. Chooses the correct file type from the file provided
    def load_sys(self, file):
        # If a file is given read the file and set the system attributes
        if file:
            self.base_file = file
            self.name = get_name(file)
            # Check the file type
            if file[-3:] == "pdb":
                self.atoms, self.data = read_pdb(self)
            elif file[-3:] == "cif":
                read_cif(self)
            elif file[-3:] == "gro":
                read_gro(self)
            elif file[-3:] == "mol":
                read_mol(self)
            else:
                print("Wrong file Loser!")
                return
        self.sort_atoms()

    # Build network function. Allows user to build the network from the system object.
    def build_network(self, net_ndx=0):
        # Start the timer
        start = time.perf_counter()
        # Set the network's atoms
        self.net.atoms = self.atoms
        # Set the settings info
        self.net.min_dist, self.net.beta_val = self.gui.sys_res_flt.get(), self.gui.sys_alpha_value.get()
        self.net.box_size = self.gui.sys_box_x_flt.get()
        # Sort the atoms in the network
        self.net.sort_atoms()
        # Check to see if there are vertices loaded
        if not self.gui.use_loaded_verts.get():
            # For small systems (<= 200) run the normal algorithm
            if len(self.atoms) <= 200:
                self.net.find_verts()
            # For large systems, split the atoms into separate smaller networks top search for vertices in
            else:
                self.net.split_sys()
            self.gui.update_progress_canvas()
            print("Connecting Network")
            build(self.net)
        print("Building Surfaces")
        # Set the output directory
        set_output_dir(self, self.output_directory)
        # Export the vertices
        self.gui.update_progress_canvas()
        # Build the network
        self.net.build_surfs()
        self.gui.update_progress_canvas()
        # Analyze the network
        print("Analyzing surfaces")
        self.net.analyze()
        self.gui.update_progress_canvas()
        print("Exporting system")
        # Export the rest of the network
        self.export_net()
        # Stop the timer and measure the time
        stop = time.perf_counter()
        self.net.my_time = stop - start

    # Export network method. Exports the values calculated by the network
    def export_net(self):
        export_net(self.net)
