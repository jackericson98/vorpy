import os.path

from System.sys_funcs.input import *
from System.sys_funcs.output import *
from Visualize.mpl_visualize import *
from System.sys_objs.group import Group
from System.sys_objs.molcule import Molecule
from System.sys_objs.residue import Residue
import numpy as np


class System:
    def __init__(self, file=None, atoms=None, verts_file=None, network_file=None, index_file=None, frame_files=None,
                 output_directory=None, gui=None, root_dir=None):
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
        self.user_atoms = atoms             # User Atoms          :   List of user provided locations and radii
        self.mols = None                    # Molecules           :   List of molecules
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
        self.ball_file = None
        self.vert_file = verts_file         # Vertex file         :   Address to the vertices of the primary system
        self.net_file = network_file        # Network files       :   Network files for multiple frames
        self.ndx_file = index_file          # Index file          :   File addresses for index file in GROMACS format
        self.frame_files = frame_files      # Frame files         :   Files storing atom movements
        self.dir = output_directory         # Output Directory    :   Output directory for the export files
        self.vpy_dir = root_dir             # Vorpy Directory     :   Directory that vorpy is running out of
        # Gui
        self.gui = gui                     # GUI                 :   GUI Vorpy object that can be updated through sys

        # Initiate the system
        if self.base_file is not None:
            self.load_files()
        np.seterr(divide='ignore', invalid='ignore')

    def load_files(self):
        """
        Create the system and make sure the files added in __init__ are added to the system
        :return:
        """
        # Load the system
        if self.base_file is not None:
            self.load_sys()
        elif self.user_atoms is not None:
            self.load_sys_atoms()

        # Load the network
        if self.net_file is not None:
            self.load_net()

        # Load the index file
        if self.ndx_file is not None:
            self.load_ndx()
        # Instantiate the major variables
        if self.atoms is None:
            self.atoms = []
        if self.residues is None:
            self.residues = []
        if self.mols is None:
            self.mols = []
        if self.ndxs is None:
            self.ndxs = []
        if self.data is None:
            self.data = []
        self.name = get_name(self.base_file)

    def load_sys(self, file=None):
        """
        Sets the base file for the system using one of the import file functions
        :param file: .pdb, .gro, .mol, .cif
        :return:
        """
        # If a file is given read the file and set the system attributes
        if file is None:
            # Set the file
            file = self.base_file
        else:
            self.base_file = file
        # Set the name of the system
        self.name = get_name(file)
        # Read PDB file
        if self.base_file[-3:] == "pdb":
            file_data = read_pdb(self)
            if file_data is not None:
                self.atoms, self.data = file_data
            else:
                return
        # Read CIF file
        elif self.base_file[-3:] == "cif":
            read_cif(self)
        # Read GRO file
        elif self.base_file[-3:] == "gro":
            read_gro(self)
        # Read MOL file
        elif self.base_file[-3:] == "mol":
            read_mol(self)
        else:
            return
        # Sort the atoms
        self.sort_atoms()

    def load_verts(self, file=None, vta_ball_file=None):
        """
        Loads vorpy specific vertices file from the system level
        :param vta_ball_file:
        :param file:
        :return:
        """
        if file is not None:
            self.vert_file = file
        # Check to see if the network has been created yet or not
        if self.net is None:
            self.net = Network(atoms=self.atoms, sys=self)
        if vta_ball_file is None:
            read_verts(self.net, self.vert_file)
        else:
            read_vta_data(self, vert_file=file, ball_file=vta_ball_file)

    def load_net(self, file=None):
        """
        Used to load a network that was previously calculated
        :param file:
        """
        # Set the output directory if None has been set yet
        if self.dir is None:
            self.set_output_directory()
        # If no file has been loaded before, create the main network
        if file is not None:
            self.net_file = file
        # Create the network
        self.net = Network(self, atoms=self.atoms)
        # read_net(self.net, self.net_file, verts_only=verts_only)
        read_net(self)
        print("\rnetwork loaded - {} verts, {} surfs\n".format(len(self.net.verts), len(self.net.surfs)), end="")

    def load_ndx(self, file=None):
        """
        Reads GROMACS index files from the system level
        :return:
        """
        read_ndx(self, file=file)

    # Build System method.
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
            self.atoms.append(Atom(location=np.random.rand(3)*2*dmax - dmax, radius=np.random.rand()*rmax, index=i))

    def sort_atoms(self):
        """
        Used to put atoms in their correct molecules and residues
        :return:
        """
        # Set up the chain names list
        self.mols, self.mol_names, self.atom_names, = [], [], []
        # Go through each of the atoms in the system adding the atoms to their respective chains
        for atom in self.atoms:

            # Set the atom's name
            self.atom_names.append(atom.element + str(self.atoms.index(atom)))
            # Add the solution
            if atom.mol_class.lower() == 'sol':
                if self.sol is None:
                    self.sol = Molecule(atoms=[atom])
                    atom.mol = self.sol
                else:
                    self.sol.atoms.append(atom)
                    atom.mol = self.sol
            else:
                # If no chain is specified, set the chain to 'None'
                mol_name = atom.mol_class + atom.chain
                # If the atom's chain does not exist add it to the list of chains
                if mol_name not in self.mol_names:
                    my_mol = Molecule(atoms=[atom], name=mol_name)
                    self.mols.append(my_mol)
                    self.mol_names.append(mol_name)
                    atom.mol = my_mol
                else:
                    my_mol = self.mols[self.mol_names.index(mol_name)]
                    my_mol.atoms.append(atom)
                    atom.mol = my_mol
        # Add the solution to the molecules list
        if self.sol is not None:
            self.mols.append(Molecule(atoms=self.sol))
            self.mol_names.append("SOL")
        # Set up the residues names list
        self.residues, self.res_names = [], []
        # Set up the residues
        for atom in self.atoms:
            res_name = atom.mol_class + atom.res_seq
            # If the residue name does not exist, add it
            if res_name not in self.res_names:
                my_res = Residue(atoms=[atom], sequence=atom.res_seq, seg_id=atom.seg_id, mol=atom.mol, name=res_name)
                self.residues.append(my_res)
                self.res_names.append(res_name)
                atom.res = my_res
                atom.mol.resids.append(my_res)
            else:
                my_res = self.residues[self.res_names.index(res_name)]
                my_res.atoms.append(atom)
                atom.res = my_res

    def build_network(self, surf_res=None, max_vert=None, box_size=None, sol_verts=True, output=True, flat_surfs=False,
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
                       flat_surfs=flat_surfs, calc_verts=calc_verts)

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

    def export_selection(self, group1, group2=None, info=True):
        """
        Exports either a single group's body or two group's bodies and the interface between them
        :param group1: Group object holding a set of atoms for cell analysis
        :param group2: Group object holding a set of atoms for interface with group 1 analysis
        :param info: Boolean for whether to export the information file
        :return:
        """
        # Change to the designated output directory
        os.chdir(self.dir)
        # Check for an interface request
        interface = False
        if group2 is not None:
            interface = True
            # Set the bff
            if group1.bff is not group2:
                group1.bff = group2
                group1.iface_atoms = None
            # Calculate the interface
            if group1.iface_atoms is None or len(group1.iface_atoms) == 0:
                group1.get_iface(bff=group2)
        # Export the first group's body
        group1.export(info=info, iface=interface)

    def set_output_directory(self):
        """
        Links set output directory to the system
        :return:
        """
        set_output_dir(self)
        os.chdir(self.dir)

    def exports(self, network=False, pdb=False, surfaces=False, full_network_object=False, no_sol_network_object=False,
                alter_atoms_script=False, info=False):
        """
        Prepares the output directory and system for output. Keeps things consistent
        :return:
        """
        # Check to see if the pdb directory is suitable
        if self.dir is None:
            if os.path.dirname(self.base_file)[-9:] != 'test_data':
                self.dir = os.path.dirname(self.base_file)
            else:
                self.set_output_directory()
        if network:
            os.chdir(self.dir)
            # Export the network
            self.export_net()
        if pdb:
            if not os.path.exists(self.dir + '/sys'):
                os.mkdir(self.dir + "/sys")
            os.chdir(self.dir + "/sys")
            # Export a pdb file for the system
            write_pdb(self.atoms, self.name, self)
            os.chdir(self.dir)
        if surfaces:
            if not os.path.exists(self.dir + '/surfs'):
                os.mkdir(self.dir + "/surfs")
            # Export a pdb file for the system
            for surf in self.net.surfs:
                write_surfs(surfs=[surf], file_name="_".join(surf.ndx), directory=self.dir + "/surfs")
            os.chdir(self.dir)
        if full_network_object and self.net.build_surfs:
            if not os.path.exists(self.dir + '/sys'):
                os.mkdir(self.dir + "/sys")
            os.chdir(self.dir + "/sys")
            # Export a full system
            export_mySys(self)
        # Write the alter atoms script
        if alter_atoms_script:
            if not os.path.exists(self.dir + '/sys'):
                os.mkdir(self.dir + "/sys")
            os.chdir(self.dir + "/sys")
            set_pymol_atoms(self)
        # If the user wants the surfaces of the system without the SOL
        if no_sol_network_object:
            # Create the group
            no_sol = Group(sys=self, mols=self.mols[:-1], name=self.name + "_shell")
            no_sol.exports(shell=True, info=True)
        # If the information is requested, export it
        if info:
            os.chdir(self.dir + "/sys")
            export_net_info(self.net)
        os.chdir(self.dir)


    def show_sys(self, info=True, show=False, fig=None, ax=None):  # Needs to be way more extensive
        if show:
            pass
        # Info section
        if info:
            # Print the header
            print("{} System information:".format(self.name))
            # Print the information
            print("System Object Counts:\n  Atoms    : {}\nMolecules: {}\n".format(len(self.atoms), len(self.mols)))

        # Show section
        plot_atoms(self.atoms, fig=fig, ax=ax)

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
