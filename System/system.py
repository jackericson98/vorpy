from os import path
from System.sys_funcs.input.pdb import read_pdb
from System.sys_funcs.input.cif import read_cif
from System.sys_funcs.input.gro import read_gro
from System.sys_funcs.input.mol import read_mol
from System.sys_funcs.input.net import read_net, read_ndx
from System.sys_funcs.input.vta import read_vta_data
from System.sys_funcs.input.verts import read_verts
from System.Network.split_net import split_net_slow, split_net
from System.sys_funcs.output.output import set_sys_dir, export_sys
from System.sys_funcs.output.net import write_verts, add_metrics
from System.Group.group import Group
from Visualize.mpl_visualize import *
from numpy import seterr
from Visualize.GUIs.periodic_table_GUI import elements


class System:
    def __init__(self, file=None, files=None, spheres=None, verts_file=None, balls_file=None, network_file=None,
                 index_file=None, frame_files=None, output_directory=None, gui=None, root_dir=None, print_actions=False,
                 atoms=None, residues=None, chains=None, segments=None, groups=None):
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
        self.chn_names = []                 # Chain Names         :   List of chain names
        self.res_names = []                 # Residue Names       :   List of residue names
        self.ndx_names = []                 # Index Names         :   List of names of indices corresponding to ndxs
        self.group_names = []               # Group Names         :   List of names of user groups for to self.groups

        # Data
        self.user_atoms = spheres             # User Atoms          :   User provided locations and radii
        self.type = 'mol'                   # Type of file        :   Holds the type of file loaded (mol, coarse, foam)
        self.foam_box = None                # Foam Retaining Box  :   Indicated in file the box that contains all balls
        self.foam_data = None               # Foam Data Info      :   Holds general information from the foam generation

        # Loadable objects
        self.spheres = spheres              # Spheres             :   List holding the atom objects
        self.atoms = atoms                  # Atoms
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.chains = chains                # Chains              :   List of the chains that make up the molecule
        self.segments = segments            # Segments            :   List of segments in the molecule
        self.sol = None                     # Solution            :   List of solution molecules (lists of atoms)

        # Settings
        self.groups = groups                # Groups              :   List of groups in the system
        self.ndxs = None                    # Indices             :   List of indices used to create groups
        self.elements = elements            # Elements            :   List of elements with mass, number, radius, group
        self.radii = my_radii               # Radii               :   List of atomic radii
        self.special_radii = special_radii  # Special Radii       :   List of special radius situations. Helpful for gro
        self.decimals = None                # Decimals            :   Decimals setting for the whole system
        self.export_type = 'large'          # Export type         :   Holds the type of objects that come out
        self.cmnds = None                   # Commands            :   All of the input commands for the sytem to be run

        # Set up the file attributes
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference
        self.files = files                  # Files               :   Files dictionary referenced for

        # Gui
        self.gui = gui                      # GUI                 :   GUI Vorpy object that can be updated through sys
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # Set the files
        self.set_files(base_file=file, ball_file=balls_file, verts_file=verts_file, ndx_file=index_file,
                       net_file=network_file, file_dir=output_directory, frame_files=frame_files, root_dir=root_dir)

        # # Initiate the system
        self.load_files()

        seterr(divide='ignore', invalid='ignore')

    def set_files(self, base_file=None, ball_file=None, verts_file=None, net_file=None, ndx_file=None, file_dir=None,
                  frame_files=None, root_dir=None):
        # Set the defaults
        defaults = {'base_file': base_file, 'ball_file': ball_file, 'verts_file': verts_file, 'net_file': net_file,
                    'ndx_file': ndx_file, 'dir': file_dir, 'frame_files': frame_files, 'root_dir': root_dir}
        # Set the files if they arent set yet
        if self.files is None:
            self.files = defaults
        # Go through the files and see if they need to be set
        for file in self.files:
            if self.files[file] is None:
                self.files[file] = defaults[file]

    def load_files(self):
        """
        Create the system and make sure the files added in __init__ are added to the system
        """

        # Load the system
        if self.files['base_file'] is not None:
            self.load_sys()
        # elif self.user_atoms is not None:
        #     self.load_sys_atoms()
        else:
            return

        # Load the network
        if self.files['net_files'] is not None:
            self.load_net()

        # Load the index file
        if self.files['ndx_file'] is not None:
            self.load_ndx()

        # Get the name
        if self.type == 'foam':
            fd = self.foam_data
            self.name = self.foam_data

        # Set the name for the system
        self.name = path.basename(self.base_file)[:-4]

    def load_sys(self, file=None):
        """
        Sets the base file for the system using one of the import file functions
        :param file: .pdb, .gro, .mol, .cif
        """
        # If a file is given read the file and set the system attributes
        if file is not None:
            # Set the file
            self.files['base_file'] = file

        # Set the name of the system
        self.name = path.basename(self.files['base_file'])[:-4]

        # Read PDB file
        if self.files['base_file'][-3:] == "pdb":
            read_pdb(self)

        # Read CIF file
        elif self.files['base_file'][-3:] == "cif":
            read_cif(self)

        # Read GRO file
        elif self.files['base_file'][-3:] == "gro":
            read_gro(self)

        # Read MOL file
        elif self.files['base_files'][-3:] == "mol":
            read_mol(self)

        # Name the system
        if self.name is None:
            self.name = path.basename(self.files['base_file'])[:-4]

        # Set the system directory
        if self.files['dir'] is None:
            self.set_output_directory()

        # If the system wants its actions printed
        if self.print_actions:
            print("{} loaded - {} atoms, {} residues, {} chain{}, ".format(self.name, len(self.atoms),
                  len(self.residues), len(self.chains), 's' if len(self.chains) > 1 else ''))

    def load_verts(self, file=None, vta_ball_file=None):
        """
        Loads vorpy specific vertices file from the system level
        :param vta_ball_file: Voronota Ball file, triggers Voronota reading of the verts file
        :param file: Main verts file that could be vorpy generated or Voronota generated
        """
        # Check for a loaded vertex file
        if file is not None:
            self.vert_file = file

        # If just verts we are loading vorpy verts
        if vta_ball_file is None:
            self.net.verts = read_verts(self.net, file)
            self.ball_file = "deez nuts"
        else:
            # If a ball file is loaded as well, this is a Voronota deal
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
        read_net(self.net, self.net_file)

        # Print if the system requires
        if self.print_actions:
            print("\r{} network loaded - {} verts, {} surfs\n"
                  .format(self.name, len(self.net.verts), len(self.net.surfs)), end="")

    def load_ndx(self, file=None):
        """
        Reads GROMACS index files from the system level
        """
        # Read the ndx file
        read_ndx(self, file=file)

        # If the system wants its actions printed
        if self.print_actions:
            print("{} indices loaded - {} indices total".format(self.name, len(self.ndxs)))

    # def load_sys_atoms(self):
    #     """
    #     Takes in a list of atomic values and creates atom objects for the system to interpret
    #     """
    #     # Disconnect atoms and user atoms
    #     self.atoms = []
    #     # Set the system Name
    #     if self.name is None:
    #         self.name = "User_Atoms"
    #     # Go through each line in the input list
    #     for i in range(len(self.user_atoms)):
    #         # Get the atom
    #         atom = self.user_atoms[i]
    #         # If the radius is a string, convert the radius using the get_radius method
    #         if isinstance(atom, Atom):
    #             self.atoms.append(atom)
    #         else:
    #             if type(atom[1]) == str:
    #                 self.atoms.append(Atom([float(atom[0][0]), float(atom[0][1]), float(atom[0][2])],
    #                                        element=atom[1], chain="None", index=i))
    #             else:
    #                 self.atoms.append(Atom([float(atom[0][0]), float(atom[0][1]), float(atom[0][2])], float(atom[1]),
    #                                        chain="None", index=i))

    # def random_system(self, anums=30, dmax=15, rmax=1):
    #     """
    #     Creates a System with atoms placed in random locations with random radii
    #     :param anums: Integer for the number of atoms in the system
    #     :param dmax: Maximum distance from the center for the atoms
    #     :param rmax: Maximum radius of an atom in the system
    #     """
    #     # Create the atoms
    #     for i in range(anums):
    #         # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
    #         self.atoms.append(Atom(location=random.rand(3)*2*dmax - dmax, radius=random.rand()*rmax, index=i))

    def print_info(self):
        atoms_var = str(len(self.spheres)) + " Atoms"
        resids_var = str(len(self.residues)) + " Residues"
        chains_var = str(len(self.chains)) + " Chains: " + ", ".join(["{} - {} atoms, {} residues"
                            .format(_.name, len(_.atoms), len(_.residues)) for _ in self.chains])
        sol_var = ""
        if self.sol is not None:
            sol_var = self.sol.name + " - " + str(len(self.sol.residues)) + " residues"
        # print(atoms_var, resids_var, chains_var, sol_var)

    def create_group(self, atoms=None, residues=None, chains=None):
        """
        Creates a group for the system
        """
        # Create the group
        print("create group")
        self.groups.append(Group(sys=self, atoms=atoms, residues=residues, chains=chains))

    def export_verts(self):
        """
        Exports the vertices after they are calculated
        """
        write_verts(self.net)

    def export_net(self):
        """
        Exports the values calculated by the network
        """
        # Export the network
        # write_net(self.net)
        pass

    def set_output_directory(self, directory=None):
        """
        Links set output directory to the system
        """
        set_sys_dir(self, dir_name=directory)

    def exports(self, all_=False, network=False, pdb=False, surfaces=False, full_network_object=False, set_atoms=False,
                info=False, logs=False, all_verts=False, all_edges=False):
        """
        Prepares the output directory and system for output. Keeps things consistent
        """
        # Export the system (/System/sys_funcs/output)
        export_sys(self, all_=all_, network=network, pdb=pdb, surfaces=surfaces, verts=all_verts, edges=all_edges,
                   full_network_object=full_network_object,
                   alter_atoms_script=set_atoms, info=info, logs=logs)

    def show_net(self, info=True, full_net=False, verts=False, edges=False, surfs=False, system=False):
        """
        Shows the network for export
        """
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
                  .format(self.net.settings['surf_res'], self.net.settings['max_vert'], self.net.settings['box_size']))

        # Show section
        if full_net or verts or edges or surfs or system:
            # Set up the figure
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
            # If the full network is expected to be shown
            if full_net:
                plot_balls(self.atoms, fig=fig, ax=ax)
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


my_radii = {'h': 1.30, 'he': 1.40, 'li': 0.76, 'be': 0.45, 'b': 1.92, 'c': 1.80, 'n': 1.60, 'o': 1.50, 'f': 1.33,
            'ne': 1.54, 'na': 1.02, 'mg': 0.72, 'al': 0.60, 'si': 2.10, 'p': 1.90, 's': 1.90, 'cl': 1.81, 'ar': 1.88,
            'k': 1.38, 'ca': 1.00, 'ga': 0.62, 'ge': 0.73, 'as': 0.58, 'se': 1.90, 'br': 1.83, 'kr': 2.02, 'rb': 1.52,
            'sr': 1.18, 'in': 1.93, 'sn': 2.17, 'sb': 2.06, 'te': 2.06, 'i': 2.20, 'xe': 2.16, 'cs': 1.67, 'ba': 1.35,
            'tl': 1.96, 'pb': 2.02, 'bi': 2.07, 'po': 1.97, 'at': 2.02, 'rn': 2.20, 'fr': 3.48, 'ra': 2.83, '': 1.80,
            'zn': 1.39}
special_radii = {''   : {'C': 1.75, 'CA': 1.90, 'N': 1.70, 'O': 1.49, 'F': 1.33, 'CL': 1.81, 'BR': 1.96, 'I': 2.20},
                 'ALA': {'CB': 1.92},
                 'ARB': {'CB': 1.91, 'CD': 1.88, 'CG': 1.92, 'CZ': 1.80, 'NE': 1.62, 'NH1': 1.62, 'NH2': 1.67},
                 'ASN': {'CB': 1.91, 'CG': 1.81, 'ND2': 1.62, 'OD1': 1.52},
                 'ASP': {'CB': 1.91, 'CG': 1.76, 'OD1': 1.49, 'OD2': 1.49},
                 'CYS': {'CB': 1.91, 'S': 1.88},
                 'GLN': {'CB': 1.91, 'CD': 1.81, 'CG': 1.80, 'NE2': 1.62, 'OE1': 1.52},
                 'GLU': {'CB': 1.91, 'CD': 1.76, 'CG': 1.88, 'OE1': 1.49, 'OE2': 1.49},
                 'HIS': {'CB': 1.91, 'CD': 1.74, 'CE': 1.74, 'CG': 1.80, 'ND1': 1.60, 'ND2': 1.60},
                 'ILE': {'CB': 2.01, 'CD1': 1.92, 'CG1': 1.92, 'CG2': 1.92},
                 'LEU': {'CB': 1.91, 'CD1': 1.92, 'CD2': 1.92, 'CG': 2.01},
                 'LYS': {'CB': 1.91, 'CD': 1.92, 'CE': 1.88, 'CG': 1.92, 'NZ': 1.67},
                 'MET': {'CB': 1.91, 'CE': 1.80, 'CG': 1.92, 'S': 1.94},
                 'PHE': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CG': 1.74, 'CZ': 1.82},
                 'PRO': {'CB': 1.91, 'CD': 1.92, 'CG': 1.92},
                 'SER': {'CB': 1.91, 'OG': 1.54},
                 'THR': {'CB': 2.01, 'CG2': 1.92, 'OG': 1.54},
                 'TRP': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CE2': 1.74, 'CG': 1.74, 'CH': 1.82, 'CZ': 1.82, 'NE1': 1.66},
                 'TYR': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CG': 1.74, 'CZ': 1.80, 'OH': 1.54},
                 'VAL': {'CB': 2.01, 'CG1': 1.92, 'CG2': 1.92}}

nucleic_acids = {'DT', 'DA', 'DG', 'DC', 'DU', 'U', 'G', 'A', 'T', 'C', 'GDP', 'OMC'}
amino_acids = {'ALA', 'ARB', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
               'THR', 'TRP', 'TYR', 'VAL', 'GLY', 'ARG'}
