# import os
import os
import time

from System.input import *
from System.output import *
from System.Network.network import *
from System.group import Group


class System:
    def __init__(self, base_file=None, atoms=None, net=None, mols=None, sol=None, residues=None, data=None, name=None,
                 vert_file=None, net_file=None, index_file=None, sig_figs=3):
        """
        Class used to import files of all types and return a System
        :param name: (str) Name describing the system
        :param atoms: List holding the atom objects
        :param mols: List of molecule atom object lists
        :param sol: List of solution atom objects
        :param residues: List of residue atom object lists
        :param data: Additional data provided by the base file
        :param base_file: Base system file address
        :param vert_file: Vertex data file address in vorpy format
        :param net_file: Network data file address in vorpy format
        :param index_file: Index file address in GROMACS index format
        """

        # Names
        self.name = name                    # Name                :   Name describing the system
        self.atom_names = None              # Atom Names          :   List holding the names of the atoms in the system
        self.mol_names = None               # Residue Names       :   List of molecule names
        self.res_names = None               # Residue Names       :   List of residue names
        self.ndx_names = None               # Index Names         :   List of names of indices corresponding to ndxs

        # Data
        self.net = net                      # Network             :   Network object holding the primary network
        self.atoms = atoms                  # Atoms               :   List holding the atom objects
        self.mols = mols                    # Molecules           :   List of molecules
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.sol = sol                      # Solution            :   List of solution molecules (lists of atoms)
        self.ndxs = None                    # Indices             :   List of lists indices of atoms
        self.radii = my_radii               # Radii               :   List of atomic radii
        self.sig_figs = sig_figs            # Significant figures :   Significant figures setting for the whole system

        # Set up the file attributes
        self.data = data                    # Data                :   Additional data provided by the base file
        self.base_file = base_file          # Base file           :   Primary file address
        self.vert_file = vert_file          # Vertex file         :   Address to the vertices of the primary system
        self.net_file = net_file            # Network files       :   Network files for multiple frames
        self.ndx_file = index_file          # Index file          :   File addresses for index file in GROMACS format
        out_dir = os.getcwd() + "/Data/User_data/"
        self.output_directory = out_dir     # Output Directory    :   Output directory for the export files
        self.vorpy_directory = os.getcwd()  # Vorpy Directory     :   Directory that vorpy is running out of

        # Gui
        self.gui = None                     # GUI                 :   GUI Vorpy object that can be updated through sys
        if self.base_file is not None:
            self.load_sys(self.base_file)   # Forced loading of the base file

    def load_sys(self, file):
        """
        Sets the base file for the system using one of the import file functions
        :param file: .pdb, .gro, .mol, .cif
        :return:
        """
        # If a file is given read the file and set the system attributes
        if file:
            # Set the file
            self.base_file = file
            self.name = get_name(file)
            # Read PDB file
            if file[-3:] == "pdb":
                self.atoms, self.data = read_pdb(self)
            # Read CIF file
            elif file[-3:] == "cif":
                read_cif(self)
            # Read GRO file
            elif file[-3:] == "gro":
                read_gro(self)
            # Read MOL file
            elif file[-3:] == "mol":
                read_mol(self)
            else:
                return
        # Sort the atoms
        self.sort_atoms()

    def load_verts(self, vert_file, vta_ball_file=None):
        """
        Loads vorpy specific vertices file from the system level
        :param vta_ball_file:
        :param vert_file:
        :return:
        """
        if vta_ball_file is None:
            read_verts(self.net, vert_file)
        else:
            add_vta_data(self, vert_file=vert_file, ball_file=vta_ball_file)

    def load_net(self, net_file, verts_only=False):
        """
        Used to load a network that was previously calculated
        :param net_file:
        :param verts_only:
        """
        # If no file has been loaded before, create the main network
        self.net_file = net_file
        self.net = Network(self, atoms=self.atoms)
        read_net(self.net, net_file, verts_only=verts_only)

    def load_ndx(self):
        """
        Reads GROMACS index files from the system level
        :return:
        """
        read_ndx(self)

    # Build System method.
    def build_user_atoms_sys(self, user_atoms):
        """
        Takes in a list of atomic values and creates atom objects for the system to interpret
        :param user_atoms: List of locations and radii for the atoms in the system
        """
        # Go through each line in the input list
        for line in user_atoms:
            # If the radius is a string, convert the radius using the get_radius method
            if type(line[1]) == str:
                self.atoms.append(Atom([float(line[0][0]), float(line[0][1]), float(line[0][2])],
                                       get_radius(self, line[1]), element=line[1], chain="None"))
            else:
                self.atoms.append(Atom([float(line[0][0]), float(line[0][1]), float(line[0][2])], float(line[1]),
                                       element=get_radius(line[1], system=self, return_symbol=True), chain="None"))

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
            self.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))

    def sort_atoms(self):
        """
        Used to put atoms in their correct molecules and residues
        :return:
        """
        # Set up the chain names list
        self.mols, self.mol_names, self.atom_names, self.sol = [], [], [], []
        # Go through each of the atoms in the system adding the atoms to their respective chains
        for atom in self.atoms:
            # Set the atom's name
            self.atom_names.append("Atom " + str(self.atoms.index(atom)) + " - " + atom.element)
            # Add the solution
            if atom.res.lower() == 'sol':
                self.sol.append(atom)
            # If no chain is specified, set the chain to 'None'
            if atom.chain == ' ':
                atom.chain = 'ZZ'
            # If the atom's chain does not exist add it to the list of chains
            if atom.chain not in self.mol_names:
                self.mols.append([atom])
                self.mol_names.append(atom.chain)
            else:
                self.mols[self.mol_names.index(atom.chain)].append(atom)
        # Add the solution to the molecules list
        if self.sol is not None:
            self.mols.append(self.sol)
            self.mol_names.append("SOL")
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

    def build_network(self, surf_res=None, max_vert=None, box_size=None, sol_verts=True, output=True, flat_faces=False,
                      find_verts=True):
        """
        Allows user to build the network from the system object.
        :return:
        """
        set_output_dir(self)
        os.chdir(self.output_directory)
        # Check to see if a network exists
        if self.net is None:
            self.net = Network(self, atoms=self.atoms)
        # Build the network
        self.net.build(surf_res=surf_res, max_vert=max_vert, box_size=box_size, sol_verts=sol_verts, output=output,
                       flat_faces=flat_faces, find_verts=find_verts)

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
        os.chdir(self.output_directory)
        # Export the first group's body
        export_body(group1, info_file=info)
        # Check for a second group
        if group2 is not None:
            export_body(group1, info_file=info)
            export_iface([group1, group2], info_file=info)

    def set_output_directory(self):
        """
        Links set output directory to the system
        :return:
        """
        set_output_dir(self)

    def initial_export(self, network=True, pdb=True, surfaces=True, full_network_object=False,
                       no_sol_network_object=True, alter_atoms_script=True):
        """
        Prepares the output directory and system for output. Keeps things consistent
        :return:
        """
        os.chdir(self.output_directory)
        if network:
            # Export the network
            self.export_net()
        if pdb:
            # Export a pdb file for the system
            write_pdb(self.atoms, self.name, self)
        if full_network_object:
            # Export a full system
            export_mySys(self)
        if no_sol_network_object:
            # Export the system without the solution
            no_sol = Group(self.net, self.sol, self.name + "_no_sol")
            no_sol.get_info()
            # Export the group
            export_body(no_sol)
        if alter_atoms_script:
            # Write the alter atoms script
            set_pymol_atoms(self)
        if surfaces:
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

    def show_sys(self, info=True, show=False, fig=None, ax=None): ######## Needs to be way more extensive

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
                  .format(self.net.min_dist, self.net.max_vert, self.net.box_size))

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
