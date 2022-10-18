import os
import time

from System.input import *
from System.output import *
from System.Network.network import *


class System:
    """Class used to import files of all types and return a System"""
    def __init__(self, atoms=None, mols=None, sol=None, residues=None, data=None, name=None, base_file=None,
                 frame_files=None, nets=None, net_files=None, ndx_files=None, gui=None):

        # Names
        self.name = name                      # Name             :    Name describing the system
        self.atom_names = None                # Atom Names       :    List holding the names of the atoms in the system
        self.mol_names = None                 # Residue Names    :    List of molecule names
        self.res_names = None                 # Residue Names    :    List of residue names
        self.ndx_names = None                 # Index Names      :    List of names of indices corresponding to ndxs

        # Data
        self.net = Network(self, atoms)       # Network          :    Network object holding the primary network
        self.nets = nets                      # Networks         :    Different networks for other frames
        self.atoms = atoms                    # Atoms            :    List holding the atom objects
        self.mols = mols                      # Molecules        :    List of molecules
        self.residues = residues              # Residues         :    List of residues (lists of atoms)
        self.sol = sol                        # Solution         :    List of solution molecules (lists of atoms)
        self.ndxs = None                      # Indices          :    List of lists indices of atoms

        # Set up the file attributes
        self.data = data                      # Data             :    Additional data provided by the base file
        self.base_file = base_file            # Base file        :    Primary file address
        self.net_files = net_files            # Network files    :    Network files for multiple frames
        self.ndx_files = ndx_files            # Index files      :    File addresses for index file in GROMACS format
        self.frame_files = frame_files        # Frame files      :    File addresses for different frames (.gro,.pdb)
        self.output_directory = None          # Output Directory :    Output directory for the export files
        self.vorpy_directory = os.getcwd()    # Vorpy Directory  :    Directory that vorpy is running out of

        # Gui
        self.gui = gui                       # GUI               :    GUI Vorpy object that can be updated through sys

    # Load network method. Used to load a network that was previously calculated
    def load_net(self, net_file, verts_only=False):
        # If no file has been loaded before, create the main network
        if self.nets is None:
            self.net_files = [net_file]
            import_net(self.net, net_file, verts_only=verts_only)
        else:
            self.net_files.append(net_file)
            self.nets.append(Network(self, self.atoms))
            import_net(self.nets[-1], net_file, verts_only=verts_only)
        self.analysis_prep()

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
        self.mols, self.mol_names, self.atom_names = [], [], []
        # Go through each of the atoms in the system adding the atoms to their respective chains
        for atom in self.atoms:
            # Set the atom's name
            self.atom_names.append("Atom " + str(self.atoms.index(atom)) + " - " + atom.element)
            # If no chain is specified, set the chain to 'None'
            if atom.chain == ' ':
                if atom.res.lower == 'sol':
                    atom.chain = 'SOL'
                else:
                    atom.chain = 'Mol'
            # If the atom's chain does not exist add it to the list of chains
            if atom.chain not in self.mol_names:
                self.mols.append([atom])
                self.mol_names.append(atom.chain)
            else:
                self.mols[self.mol_names.index(atom.chain)].append(atom)
        # Set up the residues names list
        self.residues, self.res_names = [], []
        # Set up the residues
        for atom in self.atoms:
            # Get the residue name for the atom
            res_name = atom.res + atom.res_seq
            # If the residue name does not exist, add it
            if res_name not in self.res_names:
                self.residues.append([atom])
                self.res_names.append(res_name)
            else:
                self.residues[self.res_names.index(res_name)].append(atom)

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
        # Set the output directory
        if self.output_directory is None:
            set_output_dir(self)
        # Sort the atoms
        self.sort_atoms()

    # Build network function. Allows user to build the network from the system object.
    def build_network(self):
        # Instantiate the timer variables
        self.net.my_time, self.net.cpu_time = 0, 0
        # Start the timer
        start = time.perf_counter()
        # Set the network's atoms
        self.net.atoms = self.atoms
        # Set the settings info
        if self.gui is not None:
            self.net.min_dist, self.net.beta_val = self.gui.sys_res_flt.get(), self.gui.sys_alpha_value.get()
            self.net.box_size, self.net.sol_verts = self.gui.sys_box_x_flt.get(), self.gui.sol_verts.get()
        # Sort the atoms in the network
        self.net.sort_atoms()
        # Check to see if there are vertices loaded
        if self.gui is None or not self.gui.use_loaded_verts.get():
            # Set the main network's name to main
            self.net.name = "Main"
            # Find the vertices
            self.net.find_verts()
            # Connect the network
            self.net.connect()
        # Build the edges in the network
        self.net.build_edges()
        # Build the network
        self.net.build_surfs()
        # Analyze the network
        self.net.analyze()
        # Stop the timer and measure the time
        stop = time.perf_counter()
        self.net.my_time = stop - start
        self.analysis_prep()

    # Export network method. Exports the values calculated by the network
    def export_net(self):
        # Export the network
        export_net(self.net)

    # Export selection method. Exports either a single group's body or two group's bodies and the interface between them
    def export_selection(self, group1, group2=None, info=True):
        # Change to the designated output directory
        os.chdir(self.output_directory)
        # Export the first group's body
        export_body(group1, info_file=info)
        # Check for a second group
        if group2 is not None:
            export_body(group1, info_file=info)
            export_interface([group1, group2], info_file=info)

    # Set output directory method. Links set output directory to the system
    def set_output_directory(self):
        set_output_dir(self)

    # Analysis preparation method. Prepares the output directory and system for output. Keeps things consistent
    def analysis_prep(self):
        # Set the output directory
        self.set_output_directory()
        # Export the network
        self.export_net()
        # Export a pdb file for the system
        write_pdb(self.atoms, self.name)
        # Export a full system
        export_mySys(self)
        # Make the surfaces file
        os.mkdir(self.output_directory + "/Surfaces")
        os.chdir(self.output_directory + "/Surfaces")
        # Export the surfaces one by one
        for i in range(len(self.net.surfs)):
            # Export the
            surf = self.net.surfs[i]
            # Get a random color
            my_color = np.random.rand(3)
            # Write each of the surfaces
            write_surfs([surf], "surf_" + str(surf.ndx[0]) + "_" + str(surf.ndx[1]), my_color)

