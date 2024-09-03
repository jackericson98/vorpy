import time
import os
import csv
from os import path
from System.sys_funcs.input.pdb import read_pdb
from System.sys_funcs.input.cif import read_cif
from System.sys_funcs.input.gro import read_gro
from System.sys_funcs.input.mol import read_mol
from System.sys_funcs.input.net import read_net, read_ndx
from System.sys_funcs.input.vta import read_vta_data
from System.sys_funcs.input.verts import read_verts
from System.sys_funcs.output.output import set_sys_dir, export_sys
from System.sys_funcs.output.net import write_verts
from System.Group.group import Group
from Visualize.mpl_visualize import *
from numpy import seterr
from Visualize.GUIs.periodic_table_GUI import elements
from System.radii import special_radii, element_radii


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
        self.user_atoms = spheres           # User Atoms          :   User provided locations and radii
        self.type = 'mol'                   # Type of file        :   Holds the type of file loaded (mol, coarse, foam)
        self.foam_box = None                # Foam Retaining Box  :   Indicated in file the box that contains all balls
        self.foam_data = None               # Foam Data Info      :   Holds general information from the foam generation

        # Loadable objects
        self.spheres = spheres              # Spheres             :   List holding the atom objects
        self.atoms = atoms                  # Atoms
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.chains = chains                # Chains              :   List of the chains that make up the molecule
        self.segments = segments            # Segments            :   List of segments in the molecule
        self.sol = None                     # Solute              :   List of solute molecules (lists of atoms)

        # Settings
        self.groups = groups                # Groups              :   List of groups in the system
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
        self.name = path.basename(self.files['base_file'])[:-4]

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

    def set_radii(self, my_element_radii=None, my_special_radii=None):
        """
        Sets the atom radii in the spheres dataframe based on the element radii and special radii
        """
        # First check to see of the spheres actually exist
        if self.spheres is None or len(self.spheres) == 0 or self.type != 'mol':
            return
        # Check if the user has identified some element radii they want to assign
        if my_element_radii is not None:
            # Go through the basic elemental radii to cover all atoms
            for element in my_element_radii:
                self.spheres.loc[self.spheres['element'] == element, 'rad'] = my_element_radii[element]
            # Check if we need to return
            if my_special_radii is None:
                return
        # Check if the user set the special radii
        if my_special_radii is not None:
            # Go through the special radii and assign radii based on the residue and name of the atom.
            for residue in my_special_radii:
                for name in my_special_radii[residue]:
                    self.spheres.loc[(self.spheres['res_name'] == residue) & (self.spheres['name'] == name), 'rad'] \
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
            self.files['vert_file'] = file

        # If just verts we are loading vorpy verts
        if vta_ball_file is None:
            read_verts(self.net, file)
            self.files['vert_file'] = "deez nuts"
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
        atoms_var = str(len(self.spheres)) + " Atoms"
        resids_var = str(len(self.residues)) + " Residues"
        chains_var = str(len(self.chains)) + " Chains: " + ", ".join(["{} - {} atoms, {} residues"
                            .format(_.name, len(_.atoms), len(_.residues)) for _ in self.chains])
        sol_var = ""
        if self.sol is not None:
            sol_var = self.sol.name + " - " + str(len(self.sol.residues)) + " residues"
        print(atoms_var, resids_var, chains_var, sol_var)

    def create_group(self, atoms=None, residues=None, chains=None):
        """
        Creates a group for the system
        """
        # Create the group
        self.groups.append(Group(sys=self, atoms=atoms, residues=residues, chains=chains))

    def compare_networks(self, group1, group2, data_file=None):
        """
        The goal is to take the comparison instructions and make two separate groups with their networks and compare
        their results
        """
        start = time.perf_counter()
        # Create the data storage
        data = {'vdn1': [], 'sdn1': [], 'vdn2': [], 'sdn2': [], 'rads': []}
        # Compare the networks
        for i, ball1 in group1.net.balls.iterrows():
            # Get the equivalent ball from the second group
            ball2 = group2.net.balls.iloc[i]
            # Make sure both cells are complete
            if ball1['complete'] and ball2['complete']:

                # Calculate the differences in volume and surface area for each network as the standard
                vdn1, sdn1, vdn2, sdn2, rads = ((ball2['vol'] - ball1['vol']) / ball1['vol'],
                                                (ball2['sa'] - ball1['sa']) / ball1['sa'],
                                                (ball1['vol'] - ball2['vol']) / ball2['vol'],
                                                (ball1['sa'] - ball2['sa']) / ball2['sa'], ball1['rad'])
                # Check for outliers
                if any([_ > 10 for _ in [vdn1, sdn1, vdn2, sdn2]]):
                    print('Outlier in comparison detected: {}'.format(ball1['name']))
                    continue
                # Add the data
                data['vdn1'].append(vdn1)
                data['sdn1'].append(sdn2)
                data['vdn2'].append(vdn2)
                data['sdn2'].append(sdn2)
                data['rads'].append(ball1['rad'])

                # Filter for radicals
                if any([data[_][-1] > 10 for _ in data]):
                    print(ball1['name'])
                    continue
        # Create the data line to be added to the data file
        nbs, my_line = len(data['vdn1']), []
        if nbs > 0:
            my_line = ("\r{}".format(self.files['dir']), *self.foam_data,
                       round(sum([abs(_) for _ in data['vdn1']]) / nbs, 5),  # Mean absolute difference
                       round(sum([abs(_) for _ in data['sdn1']]) / nbs, 5),  # Mean absolute difference
                       round(sum([abs(_) for _ in data['vdn2']]) / nbs, 5),  # Mean absolute difference
                       round(sum([abs(_) for _ in data['sdn2']]) / nbs, 5),  # Mean absolute difference
                       round(sum(data['vdn1']) / nbs, 5),  # Percent Difference
                       round(sum(data['sdn1']) / nbs, 5),  # Percent Difference
                       round(sum(data['vdn2']) / nbs, 5),  # Percent Difference
                       round(sum(data['sdn2']) / nbs, 5),  # Percent Difference
                       # round(np.polyfit(data['rads'], data['vdn1'], 1)[0], 5),  # Slope of the val by radius
                       # round(np.polyfit(data['rads'], data['sdn1'], 1)[0], 5),  # Slope of the val by radius
                       # round(np.polyfit(data['rads'], data['vdn2'], 1)[0], 5),  # Slope of the val by radius
                       # round(np.polyfit(data['rads'], data['sdn2'], 1)[0], 5),  # Slope of the val by radius
                       nbs, round((time.perf_counter() - start), 3))
        print(*my_line, end="")

        # Make the data file location
        if data_file is None or not path.exists(data_file):
            data_file = self.files['root_dir'] + '/Data/user_data/foam_data.csv'

        try:
            with open(data_file, 'a') as foam_file:
                foam_writer = csv.writer(foam_file)
                foam_writer.writerow(my_line)
        except PermissionError:
            with open(data_file[:-4] + '1.csv', 'a') as foam_file:
                foam_writer = csv.writer(foam_file)
                foam_writer.writerow(my_line)

    def export_verts(self):
        """
        Exports the vertices after they are calculated
        """
        write_verts(self.net)

    def set_output_directory(self, directory=None):
        """
        Links set output directory to the system
        """
        set_sys_dir(self, dir_name=directory)

    def exports(self, all_=False, pdb=False, set_atoms=False, info=False):
        """
        Prepares the output directory and system for output. Keeps things consistent
        """
        # Export the system (/System/sys_funcs/output)
        export_sys(self, all_=all_, pdb=pdb, alter_atoms_script=set_atoms, info=info)
