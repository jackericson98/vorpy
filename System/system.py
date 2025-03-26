import os
import time
from os import path
from System.sys_funcs.input.pdb import read_pdb
from System.sys_funcs.input.cif import read_cif
from System.sys_funcs.input.gro import read_gro
from System.sys_funcs.input.mol import read_mol
from System.sys_funcs.input.balls import load_balls
from System.sys_funcs.input.txt import read_txt
from System.sys_funcs.input.net import read_net, read_ndx
from System.sys_funcs.input.vta import read_vta_data
from System.sys_funcs.input.verts import read_verts
from System.sys_funcs.output.output import set_sys_dir, export_sys
from System.Group.group import Group
from numpy import seterr
from Visualize.GUI.settings.build.periodic_table_GUI import elements
from radii.Standard import special_radii, element_radii
from System.sys_funcs.calcs.compare import compare_networks, make_interfaces


class System:
    def __init__(self, file=None, files=None, spheres=None, verts_file=None, balls_file=None, network_file=None,
                 index_file=None, frame_files=None, output_directory=None, gui=None, root_dir=None, print_actions=False,
                 atoms=None, residues=None, chains=None, segments=None, groups=None, ifaces=None, simple=False, name=None):
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

        # An initial default shell system
        self.simple = simple                # Simple System       :   Indicates the system is simple and is only a shell

        # Names
        self.name = name                    # Name                :   Name describing the system
        self.atom_names = []                # Atom Names          :   List holding the names of the atoms in the system
        self.chn_names = []                 # Chain Names         :   List of chain names
        self.res_names = []                 # Residue Names       :   List of residue names
        self.ndx_names = []                 # Index Names         :   List of names of indices corresponding to ndxs
        self.group_names = []               # Group Names         :   List of names of user groups for to self.groups

        # Data
        self.user_atoms = spheres           # User Atoms          :   User provided locations and radii
        self.type = 'mol'                   # Type of file        :   Holds the type of file loaded (mol, coarse, foam)
        self.foam_box = None                # Foam Retaining Box  :   Indicated in file the box that contains all balls
        self.foam_data = None               # Foam Data Info      :   Holds general information from the foam generation

        # Loadable objects
        self.balls = spheres                # Spheres             :   List holding the atom objects
        self.atoms = atoms                  # Atoms
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.chains = chains                # Chains              :   List of the chains that make up the molecule
        self.segments = segments            # Segments            :   List of segments in the molecule
        self.sol = None                     # Solute              :   List of solute molecules (lists of atoms)

        # Settings
        self.groups = groups                # Groups              :   List of groups in the system
        self.ifaces = ifaces                # Interfaces          :   List of interface objects between groups
        self.ndxs = None                    # Indices             :   List of indices used to create groups
        self.elements = elements            # Elements            :   List of elements with mass, number, radius, group
        self.element_radii = element_radii  # Element Radii       :   Dictionary of elements and their radii
        self.special_radii = special_radii  # Special Radii       :   Dictionary of residues and their atomic radii
        self.decimals = None                # Decimals            :   Decimals setting for the whole system
        self.export_type = 'large'          # Export type         :   Holds the type of objects that come out
        self.cmnds = None                   # Commands            :   Input commands for the system to be run

        # Set up the file attributes
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference
        self.files = files                  # Files               :   Files dictionary referenced for

        # Gui
        self.gui = gui                      # GUI                 :   GUI Vorpy object that can be updated through sys
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # Set the files
        self.set_files(base_file=file, ball_file=balls_file, verts_file=verts_file, ndx_file=index_file,
                       net_file=network_file, file_dir=output_directory, frame_files=frame_files, root_dir=root_dir)

        # Check if the System is simple
        if self.simple:
            self.make_simple()
            return

        # # Initiate the system
        self.start = time.perf_counter()
        self.load_files()

        seterr(divide='ignore', invalid='ignore')

    def make_simple(self):
        """
        Makes the system work without doing all the hubub
        """
        # Set the type first
        self.type = 'simple'
        # Set everything to None
        self.load_sys(simple=True)
        self.groups, self.atoms, self.chains, self.residues = [], [], [], []
        # Set the system name
        if self.name is None:
            self.name = 'Test'
        # Set the root directory as the working directory
        self.files['root_dir'] = os.getcwd()
        # Set the output directory
        # self.set_output_directory()

    def set_files(self, base_file=None, ball_file=None, verts_file=None, net_file=None, ndx_file=None, file_dir=None,
                  frame_files=None, root_dir=None):
        # Set the defaults
        defaults = {'base_file': base_file, 'ball_file': ball_file, 'verts_file': verts_file, 'net_file': net_file,
                    'ndx_file': ndx_file, 'dir': file_dir, 'frame_files': frame_files, 'vpy_dir': os.getcwd()}        # Set the files if they arent set yet
        
        # Get the directory two levels up from this file
        if defaults['vpy_dir'] is None or defaults['vpy_dir'][-5:] != 'vorpy':
            current_file_path = os.path.abspath(__file__)
            two_dirs_up = os.path.dirname(os.path.dirname(current_file_path))
            
            # Update the vpy_dir and root_dir in defaults
            defaults['vpy_dir'] = two_dirs_up
            defaults['root_dir'] = two_dirs_up

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
        elif self.atoms is not None:
            load_balls(self)
            self.set_output_directory()
            return

        # Load the network
        if self.files['net_file'] is not None:
            self.load_net()

        # Load the index file
        if self.files['ndx_file'] is not None:
            self.load_ndx()

        # Get the name
        if self.type == 'foam':
            fd = self.foam_data
            self.name = self.foam_data

        # Set the name for the system
        try:
            self.name = path.basename(self.files['base_file'])[:-4]
        except TypeError:
            self.name = "my_system"

    def load_sys(self, file=None, simple=False):
        """
        Sets the base file for the system using one of the import file functions
        :param file: .pdb, .gro, .mol, .cif
        """
        # If a file is given read the file and set the system attributes
        if file is not None:
            # Set the file
            self.files['base_file'] = file

        # Set the name of the system
        if self.files['base_file'] is not None:
            self.name = path.basename(self.files['base_file'])[:-4]
        else:
            self.files['base_file'] = 'No File Loaded'
            self.name = self.files['base_file']

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
        elif self.files['base_file'][-3:] == "mol":
            read_mol(self)

        # Read a txt file
        elif self.files['base_file'][-3:] == 'txt':
            read_txt(self)

        # Name the system
        if self.name is None:
            self.name = path.basename(self.files['base_file'])[:-4]

        # Set the system directory
        if self.files['dir'] is None and not simple:
            self.set_output_directory()

        # If the system wants its actions printed
        if self.print_actions and not simple:
            print("{} loaded - {} atoms, {} residues, {} chain{}, "
                  .format(self.name, len(self.atoms) if self.atoms is not None else len(self.balls),
                          len(self.residues), len(self.chains), 's' if len(self.chains) > 1 else ''))

    def set_radii(self, my_element_radii=None, my_special_radii=None):
        """
        Sets the atom radii in the spheres dataframe based on the element radii and special radii
        """
        # First check to see of the spheres actually exist
        if self.balls is None or len(self.balls) == 0 or self.type != 'mol':
            return
        # Check if the user has identified some element radii they want to assign
        if my_element_radii is not None:
            # Go through the basic elemental radii to cover all atoms
            for element in my_element_radii:
                self.balls.loc[self.balls['element'] == element, 'rad'] = my_element_radii[element]
            # Check if we need to return
            if my_special_radii is None:
                return
        # Check if the user set the special radii
        if my_special_radii is not None:
            # Go through the special radii and assign radii based on the residue and name of the atom.
            for residue in my_special_radii:
                for name in my_special_radii[residue]:
                    self.balls.loc[(self.balls['res_name'] == residue) & (self.balls['name'] == name), 'rad'] \
                        = my_special_radii[residue][name]
        # If no special or element radii were specified, call the method with the system's special and element radii
        if my_special_radii is None and my_element_radii is None:
            self.set_radii(my_element_radii=self.element_radii, my_special_radii=self.special_radii)

    def load_verts(self, file=None, vta_ball_file=None):
        """
        Loads vorpy specific vertices file from the system level
        :param vta_ball_file: Voronota Ball file, triggers Voronota reading of the verts file
        :param file: Main verts file that could be vorpy generated or Voronota generated
        """
        # Check for a loaded vertex file
        if file is not None:
            self.files['verts_file'] = file

        # If just verts we are loading vorpy verts
        if vta_ball_file is None:
            if self.groups is None:
                self.create_group()
            self.groups[0].verts = read_verts(self.groups[0], file)
        else:
            # If a ball file is loaded as well, this is a Voronota deal
            read_vta_data(self, vert_file=file, ball_file=vta_ball_file)

    def load_net(self, file=None):
        """
        Used to load a network that was previously calculated
        :param file: Network file for loading
        """
        # If no file has been loaded before, create the main network
        if file is not None:
            self.files['net_file'] = file
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

    def print_info(self):
        atoms_var = str(len(self.balls)) + " Atoms"
        resids_var = str(len(self.residues)) + " Residues"
        chains_var = str(len(self.chains)) + " Chains: " + ", ".join(["{} - {} atoms, {} residues"
                            .format(_.name, len(_.atoms), len(_.residues)) for _ in self.chains])
        sol_var = ""
        if self.sol is not None:
            sol_var = self.sol.name + " - " + str(len(self.sol.residues)) + " residues"
        print(atoms_var, resids_var, chains_var, sol_var)

    def create_group(self, atoms=None, residues=None, chains=None, make_net=False):
        """
        Creates a group for the system
        """
        # Check to see of any groups have been made
        if self.groups is None:
            self.groups = []
        # Create the group
        self.groups.append(Group(sys=self, atoms=atoms, residues=residues, chains=chains, make_net=make_net))

    def compare_networks(self, group1, group2, data_file=None):
        """

        """
        compare_networks(self, group1, group2, data_file)

    def make_interfaces(self):
        """

        """
        make_interfaces(self)

    def set_output_directory(self, directory=None):
        """
        Links set output directory to the system
        """
        print('vorpy directory', self.files['vpy_dir'])
        set_sys_dir(self, dir_name=directory)

    def exports(self, all_=False, pdb=False, set_atoms=False, info=False):
        """
        Prepares the output directory and system for output. Keeps things consistent
        """
        # Export the system (/System/sys_funcs/output)
        export_sys(self, all_=all_, pdb=pdb, alter_atoms_script=set_atoms, info=info)
