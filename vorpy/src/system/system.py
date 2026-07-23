import os
import time
from os import path
from numpy import seterr
from scipy.spatial import cKDTree
from collections import Counter
from vorpy.src.inputs import read_pdb
from vorpy.src.inputs import read_cif
from vorpy.src.inputs import read_gro
from vorpy.src.inputs import read_mol
from vorpy.src.inputs import read_txt
from vorpy.src.inputs import read_net
from vorpy.src.inputs import read_ndx
from vorpy.src.inputs import read_vta
from vorpy.src.inputs import read_verts
from vorpy.src.output import set_sys_dir
from vorpy.src.output import export_sys
from vorpy.src.group import Group
from vorpy.src.interface import Interface
from vorpy.src.GUI.system.radii_adjustments.periodic_table_GUI import elements
from vorpy.src.chemistry import special_radii
from vorpy.src.chemistry import element_radii
from vorpy.src.calculations import compare_networks
from vorpy.src.calculations import get_time


def infer_index_offset(index_map):
    offsets = [
        int(sys_i) - int(log_i)
        for log_i, sys_i in index_map.items()
    ]

    if not offsets:
        return 0

    offset, count = Counter(offsets).most_common(1)[0]

    if count / len(offsets) < 0.95:
        return None

    return offset


def _atom_key_from_row(row):
    return (
        str(row.get("Name", row.get("name", ""))).strip(),
        str(row.get("Residue", row.get("res_name", ""))).strip(),
        str(row.get("Residue Sequence", row.get("res_seq", ""))).strip(),
        str(row.get("Chain", row.get("chain", ""))).strip(),
    )


def _atom_key_relaxed_from_row(row):
    return (
        str(row.get("Name", row.get("name", ""))).strip(),
        str(row.get("Residue", row.get("res_name", ""))).strip(),
        str(row.get("Chain", row.get("chain", ""))).strip(),
    )


def _loc_from_row(row):
    loc = row.get("loc", None)

    if loc is not None:
        return [float(_) for _ in loc]

    if {"X", "Y", "Z"}.issubset(row.index):
        return [float(row["X"]), float(row["Y"]), float(row["Z"])]

    return None


def _dist3(a, b):
    return sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)) ** 0.5


def build_log_to_system_index_map(sys_balls, log_balls, tol=0.01, coord_fallback_tol=0.05):
    start_time = time.perf_counter()
    total = len(log_balls)

    exact_lookup = {}
    relaxed_lookup = {}

    sys_locs = []
    sys_indices = []

    for sys_i, sys_ball in sys_balls.iterrows():
        exact_key = _atom_key_from_row(sys_ball)
        relaxed_key = _atom_key_relaxed_from_row(sys_ball)

        exact_lookup.setdefault(exact_key, []).append(sys_i)
        relaxed_lookup.setdefault(relaxed_key, []).append(sys_i)

        sys_loc = _loc_from_row(sys_ball)
        if sys_loc is not None:
            sys_locs.append(sys_loc)
            sys_indices.append(sys_i)

    sys_tree = cKDTree(sys_locs) if sys_locs else None

    index_map = {}
    used_sys_indices = set()
    mismatches = []

    for count, (_, log_ball) in enumerate(log_balls.iterrows(), start=1):
        if count == 1 or count % 100 == 0 or count == total:
            my_time = time.perf_counter() - start_time
            h, m, s = get_time(my_time)
            percentage = 100.0 * count / total

            print(
                "\rRun Time = {}:{:02d}:{:2.2f} - Process: verifying log/system atom map: {} atoms - {:.2f} %"
                .format(int(h), int(m), round(s, 2), count, percentage),
                end=""
            )

        log_i = int(log_ball["Index"])
        log_key = _atom_key_from_row(log_ball)
        log_key_relaxed = _atom_key_relaxed_from_row(log_ball)
        log_loc = _loc_from_row(log_ball)

        candidate_source = "exact"
        candidates = exact_lookup.get(log_key, [])

        if not candidates:
            candidate_source = "relaxed"
            candidates = relaxed_lookup.get(log_key_relaxed, [])

        candidates = [i for i in candidates if i not in used_sys_indices]

        if not candidates and log_loc is not None and sys_tree is not None:
            candidate_source = "coordinate"

            dist, tree_i = sys_tree.query(log_loc, k=10)

            if not hasattr(dist, "__len__"):
                dist = [dist]
                tree_i = [tree_i]

            candidates = [
                sys_indices[int(i)]
                for d, i in zip(dist, tree_i)
                if int(i) < len(sys_indices)
                and d <= coord_fallback_tol
                and sys_indices[int(i)] not in used_sys_indices
            ]

        if not candidates:
            mismatches.append(
                f"Logs Ball {log_i} has no available match in system.\n"
                f"Logs exact key = {log_key}\n"
                f"Logs relaxed key = {log_key_relaxed}\n"
                f"Logs loc = {log_loc}"
            )
            continue

        best_sys_i = None
        best_dist = float("inf")

        if log_loc is not None:
            for sys_i in candidates:
                sys_loc = _loc_from_row(sys_balls.iloc[sys_i])
                if sys_loc is None:
                    continue

                d = _dist3(log_loc, sys_loc)

                if d < best_dist:
                    best_dist = d
                    best_sys_i = sys_i
        else:
            best_sys_i = candidates[0]
            best_dist = 0.0

        allowed_tol = coord_fallback_tol if candidate_source == "coordinate" else tol

        if best_sys_i is None or best_dist > allowed_tol:
            mismatches.append(
                f"Logs Ball {log_i} matched candidate by {candidate_source}, but coordinates do not align.\n"
                f"Logs exact key = {log_key}\n"
                f"Logs relaxed key = {log_key_relaxed}\n"
                f"Logs loc = {[round(_, 3) for _ in log_loc] if log_loc is not None else None}\n"
                f"Closest System Ball {best_sys_i} loc = "
                f"{[round(_, 3) for _ in _loc_from_row(sys_balls.iloc[best_sys_i])] if best_sys_i is not None else None}\n"
                f"Distance = {round(best_dist, 4)}"
            )
            continue

        index_map[log_i] = int(best_sys_i)
        used_sys_indices.add(int(best_sys_i))

    print()

    if mismatches:
        raise ValueError(
            "Logs file indices do not align with the current system and could not be remapped.\n\n"
            + "\n\n".join(mismatches[:5])
        )

    return index_map


class System:
    def __init__(self, file=None, files=None, spheres=None, verts_file=None, balls_file=None, network_file=None,
                 index_file=None, frame_files=None, output_directory=None, gui=None, root_dir=None, print_actions=False,
                 atoms=None, residues=None, chains=None, segments=None, groups=None, ifaces=None, simple=False, name=None):
        """
        Initialize a new System object for managing molecular systems and networks.

        Attributes
        ----------
        file : str, optional
            Base system file address
        files : dict, optional
            Dictionary containing file paths for various system components
        spheres : list, optional
            List of sphere objects representing atoms or particles
        verts_file : str, optional
            Vertex data file address in vorpy format
        balls_file : str, optional
            Voronota Ball file, triggers Voronota reading of the verts file
        network_file : str, optional
            Network data file address in vorpy format
        index_file : str, optional
            Index file address in GROMACS index format
        frame_files : list, optional
            Files for atom movements
        output_directory : str, optional
            Directory for export files to be output to
        gui : object, optional
            The GUI object (tkinter) associated with loading the system
        root_dir : str, optional
            Root directory for the system
        print_actions : bool, optional
            Tells the system to print or not
        atoms : list, optional
            List holding the atom objects
        residues : list, optional
            List of residues (lists of atoms)
        chains : list, optional
            List of the chains that make up the molecule
        segments : list, optional
            List of segments in the molecule
        groups : list, optional
            List of groups in the system
        ifaces : list, optional
            List of interface objects between groups
        simple : bool, optional
            Indicates the system is simple and is only a shell
        name : str, optional
            Name describing the system

        Methods
        -------
        make_simple(): 
            Create a simple system with specified components.
        set_files(): 
            Set the file paths for the system components.
        load_files(): 
            Load and initialize all system files specified during initialization.
        load_sys(): 
            Load and set the base file for the system.
        set_radii(): 
            Set the atom radii in the spheres dataframe.
        load_verts(): 
            Load vertices from a file into the system.
        load_net(): 
            Load a previously calculated network into the system.
            
        Notes
        -----
        This initializes a new System object with the specified files and components.
        If no files are provided, an empty system is created.
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
        self.round_to = None                # round_to            :   Decimals setting for the whole system
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
        Create a simple system with specified components.

        Notes
        -----
        This method creates a simple system with the specified components.
        All parameters are optional and can be added later.
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
        """
        Set the file paths for the system components.

        Parameters
        ----------
        base_file : str, optional
            Base system file address
        ball_file : str, optional
            Ball file address
        verts_file : str, optional
            Vertex data file address in vorpy format
        net_file : str, optional
            Network data file address in vorpy format
        ndx_file : str, optional
            Index file address in GROMACS index format
        file_dir : str, optional
            Directory for export files to be output to
        frame_files : list, optional
            Files for atom movements
        root_dir : str, optional
            Root directory for the system

        Notes
        -----
        This method updates the system's file paths dictionary with the provided values.
        If a parameter is not provided, its corresponding entry in the dictionary remains unchanged.
        """
        # Set the defaults
        defaults = {'base_file': base_file, 'ball_file': ball_file, 'verts_file': verts_file, 'net_file': net_file,
                    'ndx_file': ndx_file, 'dir': file_dir, 'frame_files': frame_files, 'vpy_dir': os.getcwd()}
        # Set the files if they aren't set yet
        
        # Get the directory two levels up from this file
        if defaults['vpy_dir'] is None or defaults['vpy_dir'][-5:] != 'vorpy':
            current_file_path = os.path.abspath(__file__)
            two_dirs_up = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_file_path))))
            
            # Update the vpy_dir and root_dir in defaults
            defaults['vpy_dir'] = two_dirs_up
            defaults['root_dir'] = two_dirs_up
        # Set the files if they aren't set yet
        if self.files is None:
            self.files = defaults
        # Go through the files and see if they need to be set
        for file in self.files:
            if self.files[file] is None:
                self.files[file] = defaults[file]

    def load_files(self):
        """
        Load and initialize all system files specified during initialization.

        """
        # Load the system
        if self.files['base_file'] is not None:
            self.load_sys()

        # elif self.user_atoms is not None:
        #     self.load_sys_atoms()
        elif self.atoms is not None:
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

        # Set the output directory 
        if self.files['dir'] is None or 'No System Chosen' in self.files['dir']:
            self.set_output_directory()

    def load_sys(self, file=None, simple=False, make_dir=True):
        """Load and set the base file for the system.

        This method loads a molecular structure file and sets up the system attributes.
        It supports various file formats including PDB, GRO, MOL, and CIF.

        Parameters
        ----------
        file : str, optional
            Path to the molecular structure file (.pdb, .gro, .mol, .cif)
        simple :
        make_dir :
        """
        # We first need to determine what type of file is being loaded
        if file is not None and file[-3:] == 'csv':
            self.load_net(file=file)
            return

        # If a file is given read the file and set the system attributes
        if file is not None:
            # Set the file
            self.files['base_file'] = file

        # Set the name of the system
        if self.files['base_file'] is not None:
            self.name = path.basename(self.files['base_file'])[:-4].capitalize()
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
        if not simple and make_dir:
            self.set_output_directory()

        # If the system wants its actions printed
        if self.print_actions and not simple:
            print("{} loaded - {} atoms, {} residues, {} chain{}, "
                  .format(self.name, len(self.atoms) if self.atoms is not None else len(self.balls),
                          len(self.residues), len(self.chains), 's' if len(self.chains) > 1 else ''))

    def set_radii(self, my_element_radii=None, my_special_radii=None):
        """Set the atom radii in the spheres dataframe.

        This function assigns radii to atoms in the spheres dataframe based on either
        element-specific radii or special radii defined for specific residues and atoms.

        Parameters
        ----------
        my_element_radii : dict, optional
            Dictionary mapping element symbols to their radii
        my_special_radii : dict, optional
            Dictionary mapping residue names to dictionaries of atom names and their radii

        Notes
        -----
        If both parameters are None, the function will use the system's default
        element_radii and special_radii attributes.
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
        """Load vertices from a file into the system.

        This function loads vertex data from either a vorpy-specific vertices file or a Voronota-generated file.
        The function can handle both vorpy-generated vertices and Voronota-generated vertices when provided
        with the appropriate ball file.

        Parameters
        ----------
        file : str, optional
            Path to the main vertices file. Can be either vorpy-generated or Voronota-generated.
        vta_ball_file : str, optional
            Path to the Voronota ball file. If provided, triggers Voronota-specific reading of the vertices file.

        Notes
        -----
        If vta_ball_file is None, the function assumes a vorpy-specific vertices file and creates a group
        if none exists. If vta_ball_file is provided, the function uses Voronota-specific reading methods.
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
            read_vta(self, vert_file=file, ball_file=vta_ball_file)

    def load_net(self, file=None, group=None, rebuild_edges=True, rebuild_surfs=True, analyze=True, store_points=True):
        """
        Load a previously calculated logs/network file into a Group network.
        """
        if file is not None:
            self.files["net_file"] = file

        if self.files["net_file"] is None:
            raise ValueError("No network file provided.")

        if self.groups is None or len(self.groups) == 0:
            self.create_group(make_net=True)

        if group is None:
            group = self.groups[0]

        if group.net is None:
            group.make_net()

        read_net(
            group,
            group.net,
            self.files["net_file"],
            rebuild_edges=rebuild_edges,
            rebuild_surfs=rebuild_surfs,
            analyze=analyze,
            store_points=store_points,
        )

        if self.print_actions:
            print(
                f"\r{self.name} network loaded - "
                f"{len(group.net.verts)} verts, {len(group.net.edges)} edges, {len(group.net.surfs)} surfs\n",
                end=""
            )

        return group.net

    def load_group_logs(self, group_name, file, rebuild_edges=True, rebuild_surfs=True, analyze=True,
                        store_points=True):
        """
        Load a logs file into a specific group.

        Supports logs whose atom indices do not directly align with the currently
        loaded system by building a log-index -> system-index map.
        """
        from vorpy.src.inputs.logs import read_logs

        if self.balls is None or len(self.balls) == 0:
            raise ValueError("Load a ball/base file before loading group logs.")

        log_data = read_logs(file, all_=True)
        log_balls = log_data["atoms"]

        log_to_sys_index = build_log_to_system_index_map(self.balls, log_balls)
        index_offset = infer_index_offset(log_to_sys_index)

        print(f"\nDetected log/system index offset: {index_offset}")

        if self.groups is None:
            self.groups = []

        group = None
        for existing_group in self.groups:
            if existing_group.name == group_name:
                group = existing_group
                break

        if group is None:
            self.create_group(make_net=True)
            group = self.groups[-1]
            group.name = group_name

        if group.net is None:
            group.make_net()

        read_net(
            group,
            group.net,
            file,
            rebuild_edges=rebuild_edges,
            rebuild_surfs=rebuild_surfs,
            analyze=analyze,
            store_points=store_points,
            index_map=log_to_sys_index,
            index_offset=index_offset,
        )

        group.name = log_data["group data"].get("Name", group_name)

        loaded_ball_indices = sorted(log_to_sys_index.values())

        group.loaded_log_ball_indices = loaded_ball_indices
        group.ball_ndxs = loaded_ball_indices
        group.log_to_sys_index = log_to_sys_index

        return group

    def load_ndx(self, file=None):
        """Load GROMACS index files into the system.

        Parameters
        ----------
        file : str, optional
            Path to the GROMACS index file to be loaded.

        Notes
        -----
        This method reads GROMACS index files and stores them in the system's ndxs attribute.
        If print_actions is enabled, it will print confirmation of the loaded indices.
        """
        # Read the ndx file
        read_ndx(self, file=file)

        # If the system wants its actions printed
        if self.print_actions:
            print("{} indices loaded - {} indices total".format(self.name, len(self.ndxs)))

    def print_info(self):
        """
        Print system information including atom, residue, and chain counts.

        Notes
        -----
        Output includes:
        - Total number of atoms
        - Total number of residues
        - Chain information with atom and residue counts
        - Solute information if present
        """
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
        Create a new group in the system.

        Parameters
        ----------
        atoms : list, optional
            List of atoms to include in the group
        residues : list, optional
            List of residues to include in the group
        chains : list, optional
            List of chains to include in the group
        make_net : bool, optional
            Whether to create a network for the group

        Returns
        -------
        Group
            The newly created group object

        Notes
        -----
        If no groups exist, a new groups list will be created.
        The group will be added to the system's groups list.
        """
        # Check to see of any groups have been made
        if self.groups is None:
            self.groups = []
        # Create the group
        self.groups.append(Group(sys=self, atoms=atoms, residues=residues, chains=chains, make_net=make_net))

    def compare_networks(self, group1, group2, data_file=None):
        """Compare two networks and analyze their differences.

        Parameters
        ----------
        group1 : Group
            First group containing a network to compare
        group2 : Group
            Second group containing a network to compare
        data_file : str, optional
            Path to file where comparison data should be saved

        Returns
        -------
        None

        Notes
        -----
        This method compares the networks of two groups and analyzes their structural
        differences. The comparison can be saved to a data file if specified.
        """
        compare_networks(self, group1, group2, data_file)

    def make_interfaces(self, interface_pairs):
        """
        Create and build the requested interfaces.

        Parameters
        ----------
        interface_pairs : iterable
            Iterable of ``(group1, group2)`` pairs. ``group2`` may be None when
            constructing an interface against the surrounding system.
        """
        if self.ifaces is None:
            self.ifaces = []

        for group1, group2 in interface_pairs:
            interface = Interface(
                sys=self,
                group1=group1,
                group2=group2,
            )

            interface.build()
            self.ifaces.append(interface)

    def set_output_directory(self, directory=None):
        """
        Set the output directory for system exports.

        Parameters
        ----------
        directory : str, optional
            Path to the directory where system files will be exported

        Notes
        -----
        If no directory is specified, a default directory will be used.
        The directory will be created if it doesn't exist.
        """
        set_sys_dir(self, dir_name=directory)

    def exports(self, all_=False, pdb=False, set_atoms=False, info=False, mol=False, cif=False, xyz=False, txt=False):
        """Export system data in various formats.

        This method prepares the output directory and system for export in the specified formats.
        It ensures consistent output formatting across different export types.

        Parameters
        ----------
        all_ : bool, optional
            If True, export in all available formats
        pdb : bool, optional
            If True, export as PDB file
        set_atoms : bool, optional
            If True, modify atom properties before export
        info : bool, optional
            If True, export system information
        mol : bool, optional
            If True, export as MOL file
        cif : bool, optional
            If True, export as CIF file
        xyz : bool, optional
            If True, export as XYZ file
        txt : bool, optional
            If True, export as text file

        Returns
        -------
        None

        Notes
        -----
        The output directory must be set before calling this method.
        """
        # Export the system (/System/sys_funcs/output)
        export_sys(self, all_=all_, pdb=pdb, alter_atoms_script=set_atoms, info=info, mol=mol, cif=cif, xyz=xyz, txt=txt)
