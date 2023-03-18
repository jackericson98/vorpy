from System.sys_funcs.input import *
from System.sys_funcs.read_net import read_net
from System.sys_funcs.sort_atoms import sort_atoms
from System.sys_funcs.output import *
from Visualize.mpl_visualize import *
from numpy import seterr, random


class System:
    def __init__(self, file=None, atoms=None, verts_file=None, network_file=None, index_file=None, frame_files=None,
                 output_directory=None, gui=None, root_dir=None, print_actions=False):
        """
        Class used to import files of all types and return a System
        :param file: Base system file address
        :param atoms: List holding the atom objects
        :param verts_file: Vertex data file address in vorpy format
        :param network_file: Network data file address in vorpy format
        :param index_file: Index file address in GROMACS index format
        :param frame_files: Files for atom movements
        :param output_directory: Directory for export files to be output to
        :param gui: The GUI object (tkinter) associated with loading the system and loading/creating the network
        """

        # Names
        self.name = None                    # Name                :   Name describing the system
        self.atom_names = []                # Atom Names          :   List holding the names of the atoms in the system
        self.mol_names = []                 # Residue Names       :   List of molecule names
        self.res_names = []                 # Residue Names       :   List of residue names
        self.ndx_names = []                 # Index Names         :   List of names of indices corresponding to ndxs
        self.group_names = []               # Group Names         :   List of names of user groups for to self.groups

        # Data
        self.net = None                     # Network             :   Network object holding the primary network
        self.atoms = atoms                  # Atoms               :   List holding the atom objects
        self.user_atoms = atoms             # User Atoms          :   User provided locations and radii
        self.mols = None                    # Molecules           :   Molecule objects from the system
        self.residues = None                # Residues            :   List of residues (lists of atoms)
        self.sol = None                     # Solution            :   List of solution molecules (lists of atoms)
        self.sol_name = "None"              # Solute Name         :   Name for the solute from the atoms
        self.groups = []                    # Groups              :   List of groups in the system
        self.ndxs = []                      # Indices             :   List of lists indices of atoms
        self.radii = my_radii               # Radii               :   List of atomic radii
        self.decimals = None                # Decimals            :   Decimals setting for the whole system

        # Set up the file attributes
        self.data = None                    # Data                :   Additional data provided by the base file
        self.base_file = file               # Base file           :   Primary file address
        self.ball_file = None               # Ball file           :   Balls used to create vertices from Voronota
        self.vert_file = verts_file         # Vertex file         :   Address to the vertices of the primary system
        self.net_file = network_file        # Network files       :   Network files for multiple frames
        self.ndx_file = index_file          # Index file          :   File addresses for index file in GROMACS format
        self.frame_files = frame_files      # Frame files         :   Files storing atom movements
        self.dir = output_directory         # Output Directory    :   Output directory for the export files
        self.vpy_dir = root_dir             # Vorpy Directory     :   Directory that vorpy is running out of
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference

        # Gui
        self.gui = gui                      # GUI                 :   GUI Vorpy object that can be updated through sys
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # # Initiate the system
        self.__load_files__()

        seterr(divide='ignore', invalid='ignore')

    def __load_files__(self):
        """
        Create the system and make sure the files added in __init__ are added to the system
        :return:
        """

        # Load the system
        if self.base_file is not None:
            self.load_sys()
        elif self.user_atoms is not None:
            self.load_sys_atoms()
        else:
            return

        # Load the network
        if self.net_file is not None:
            self.load_net()

        # Load the index file
        if self.ndx_file is not None:
            self.load_ndx()

        # Get the name
        self.name = os.path.basename(self.base_file)[:-4]

    def load_sys(self, file=None):
        """
        Sets the base file for the system using one of the import file functions
        :param file: .pdb, .gro, .mol, .cif
        :return:
        """
        # If a file is given read the file and set the system attributes
        if file is not None:
            # Set the file
            self.base_file = file

        # Set the name of the system
        self.name = path.basename(self.base_file)[:-4]

        # Read PDB file
        if self.base_file[-3:] == "pdb":
            read_pdb(self)

        # Read CIF file
        elif self.base_file[-3:] == "cif":
            read_cif(self)

        # Read GRO file
        elif self.base_file[-3:] == "gro":
            read_gro(self)

        # Read MOL file
        elif self.base_file[-3:] == "mol":
            read_mol(self)

        # Name the system
        if self.name is None:
            self.name = os.path.basename(self.base_file)[:-4]

        # Sort the atoms
        sort_atoms(self)

        # If the system wants its actions printed
        if self.print_actions:
            print("{} loaded - {} atoms, {} chains, {} residues"
                  .format(self.name, len(self.atoms), len(self.mols), len(self.residues)))

    def load_verts(self, file=None, vta_ball_file=None):
        """
        Loads vorpy specific vertices file from the system level
        :param vta_ball_file: Voronota Ball file, triggers voronota reading of the verts file
        :param file: Main verts file that could be vorpy generated or voronota generated
        :return: Sets the vertex values for the network
        """
        # Check for a loaded vertex file
        if file is not None:
            self.vert_file = file

        # Check to see if the network has been created yet or not
        if self.net is None:
            self.net = Network(atoms=self.atoms, sys=self)

        # If a ball file is loaded as well, this is a voronota deal
        if vta_ball_file is None:
            read_verts(self.net, self.vert_file)
        else:
            read_vta_data(self, vert_file=file, ball_file=vta_ball_file)

        # If the system wants its actions printed
        if self.print_actions:
            print("{} vertices loaded - {} vertices".format(self.name, len(self.net.verts)))

    def load_net(self, file=None):
        """
        Used to load a network that was previously calculated
        :param file: Network file for loading
        """
        # If no file has been loaded before, create the main network
        if file is not None:
            self.net_file = file
        # Read the network file
        read_net(self, integrate=self.net is not None)

        # Print if the system requires
        if self.print_actions:
            print("\r{} network loaded - {} verts, {} surfs\n"
                  .format(self.name, len(self.net.verts), len(self.net.surfs)), end="")

    def load_ndx(self, file=None):
        """
        Reads GROMACS index files from the system level
        :return: Creates group objects for the system
        """
        # Read the ndx file
        read_ndx(self, file=file)

        # If the system wants its actions printed
        if self.print_actions:
            print("{} indices loaded - {} indices total".format(self.name, len(self.ndxs)))

    def load_sys_atoms(self):
        """
        Takes in a list of atomic values and creates atom objects for the system to interpret
        """
        # Disconnect atoms and user atoms
        self.atoms = []
        # Set the system Name
        if self.name is None:
            self.name = "User_Atoms"
        # Go through each line in the input list
        for i in range(len(self.user_atoms)):
            # Get the atom
            atom = self.user_atoms[i]
            # If the radius is a string, convert the radius using the get_radius method
            if isinstance(atom, Atom):
                self.atoms.append(atom)
            else:
                if type(atom[1]) == str:
                    self.atoms.append(Atom([float(atom[0][0]), float(atom[0][1]), float(atom[0][2])],
                                           get_radius(self, atom[1]), element=atom[1], chain="None", index=i))
                else:
                    self.atoms.append(Atom([float(atom[0][0]), float(atom[0][1]), float(atom[0][2])], float(atom[1]),
                                           element=get_radius(atom[1], system=self, return_symbol=True), chain="None",
                                           index=i))

    def random_system(self, anums=30, dmax=15, rmax=1):
        """
        Creates a System with atoms placed in random locations with random radii
        :param anums: Integer for the number of atoms in the system
        :param dmax: Maximum distance from the center for the atoms
        :param rmax: Maximum radius of an atom in the system
        :return:
        """
        # Create the atoms
        for i in range(anums):
            # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
            self.atoms.append(Atom(location=random.rand(3)*2*dmax - dmax, radius=random.rand()*rmax, index=i))

    def build_network(self, surf_res=None, max_vert=None, box_size=None, sol_verts=True, output=True,
                      calc_verts=True):
        """
        Allows user to build the network from the system object.
        :return:
        """
        # Check to see if a network exists
        if self.net is None:
            self.net = Network(self, atoms=self.atoms)
        # Build the network
        self.net.build(surf_res=surf_res, max_vert=max_vert, box_size=box_size, build_surfs=sol_verts, output=output,
                       calc_verts=calc_verts)

    def export_verts(self):
        """
        Exports the vertices after they are calculated
        :return:
        """
        export_verts(self.net)

    def export_net(self):
        """
        Exports the values calculated by the network
        :return:
        """
        # Export the network
        export_net(self.net)

    def set_output_directory(self, directory=None):
        """
        Links set output directory to the system
        :return:
        """
        set_output_dir(self, dir_name=directory)

    def exports(self, network=False, pdb=False, surfaces=False, full_network_object=False, no_sol_shell=False,
                set_atoms=False, info=False):
        """
        Prepares the output directory and system for output. Keeps things consistent
        :return:
        """
        # Export the system (/System/sys_funcs/output)
        export_sys(self, network=network, pdb=pdb, surfaces=surfaces, full_network_object=full_network_object,
                   no_sol_network_object=no_sol_shell, alter_atoms_script=set_atoms, info=info)

    def show_net(self, info=True, full_net=False, verts=False, edges=False, surfs=False, system=False):
        # Empty Network
        if self.net is None:
            print("No network constructed")
            return
        # Info section
        if info:
            # Print the header
            print("{} Network Information:\n".format(self.name))
            # Network objects information
            print("\nObject counts:\n  Vertices: {}\n  Edges   : {}\n  Surfaces: {}\n  Doublets: {}"
                  .format(len(self.net.verts), len(self.net.edges), len(self.net.surfs), len(self.net.doublets)))

            # Print the build parameters
            print("\nBuild parameters:\n  Surface Resolution: {}\n  Maximum Vertex Radius: {}\n  Retaining Box Size: {}"
                  .format(self.net.surf_res, self.net.max_vert, self.net.box_size))

        # Show section
        if full_net or verts or edges or surfs or system:
            # Set up the figure
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
            # If the full network is expected to be shown
            if full_net:
                plot_atoms(self.atoms, fig=fig, ax=ax)
                plot_verts(self.net.verts, fig=fig, ax=ax)
                plot_surfs(self.net.surfs, fig=fig, ax=ax)
                plot_edges(self.net.edges, fig=fig, ax=ax, Show=True)
            # Plot individual items
            else:
                # Plot the verts, surfs or edges
                if verts:
                    plot_verts(self.net.verts, fig=fig, ax=ax)
                if surfs:
                    plot_surfs(self.net.surfs, fig=fig, ax=ax)
                if verts:
                    plot_edges(self.net.edges, fig=fig, ax=ax)


##################################################### Atomic Radii #####################################################


my_radii = [['h' , 'he', 'li', 'be', 'b' , 'c' , 'n' , 'o' , 'f' , 'ne', 'na', 'mg', 'al', 'si', 'p' , 's' , 'cl', 'ar',
             'k' , 'ca', 'sc', 'ti', 'v' , 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
             'rb', 'sr', 'y' , 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i' , 'xe',
             'cs', 'ba', 'la', 'hf', 'ta', 'w' , 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
             'fr', 'ra', 'ac', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og',
             'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu',
             'th', 'pa', 'u' , 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr'],
            [1.30, 1.40, 0.76, 0.45, 1.92, 1.80, 1.60, 1.50, 1.33, 1.54, 1.02, 0.72, 0.60, 2.10, 1.90, 1.90, 1.81, 1.88,
             1.38, 1.00, None, None, None, None, None, None, None, None, None, None, 0.62, 0.73, 0.58, 1.90, 1.83, 2.02,
             1.52, 1.18, None, None, None, None, None, None, None, None, None, None, 1.93, 2.17, 2.06, 2.06, 2.20, 2.16,
             1.67, 1.35, None, None, None, None, None, None, None, None, None, None, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20,
             3.48, 2.83, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
             None, None, None, None, None, None, None, None, None, None, None, None, None, None,
             None, None, None, None, None, None, None, None, None, None, None, None, None, None]]
