from vorpy.src.group.layers import get_layers
from vorpy.src.group.sort import add_balls
from vorpy.src.group.sort import get_info
from vorpy.src.group.export import group_exports
from vorpy.src.network import Network
from vorpy.src.inputs import read_verts


class Group:
    """Group class for managing and analyzing collections of atoms and their properties.

    The Group class serves as a container for selections of atoms, molecules, chains, and residues
    from a molecular system. It provides functionality for network analysis, surface area calculations,
    volume measurements, and layer-based analysis of the selected components.

    Key Features:
    - Network construction and analysis
    - Surface area and volume calculations
    - Center of mass and moment of inertia calculations
    - Layer-based analysis of atoms and surfaces
    - Export capabilities for visualization and analysis

    Attributes:
        sys (Network): The parent system network
        name (str): Name of the group
        dir (str): Output directory for group exports
        net (Network): Network object containing the group's topology
        verts: Vertex information
        ball_ndxs (list): Indices of atoms included in the network
        settings (dict): Network construction settings
        atms (list): List of atoms in the group
        mols (list): List of molecules in the group
        chns (list): List of chains in the group
        rsds (list): List of residues in the group
        sa (float): Surface area of outer surfaces
        vol (float): Volume of atoms' cells
        density (float): Ratio of atom volumes to total space
        mass (float): Total mass of atoms
        com (list): Center of mass coordinates
        vdw_vol (float): Van der Waals volume
        vdw_com (list): Van der Waals center of mass
        spatial_moment (list): Spatial moment tensor
        moi (list): Moment of inertia tensor
        layer_atoms (list): Atoms arranged by layers
        layer_verts (list): Vertices arranged by layers
        layer_edges (list): Edges arranged by layers
        layer_surfs (list): Surfaces arranged by layers
        layer_info (list): Layer-specific information

    Methods:
        __init__: Initialize a new Group instance
        get_settings: Configure network and analysis settings
        set_name: Set a default name for the group
        build: Construct the network for the group
        get_layers: Perform layer-based analysis
        get_info: Gather information about the group
        add_balls: Add atoms to the group
        export: Export group data for visualization

    Examples:
        # Example 1: Basic Group Creation and Analysis
        >>> from vorpy.src.group import Group
        >>> # Create a group from a system with default settings
        >>> group = Group(sys=my_system, name='protein_A')
        >>> # This creates a new group named 'protein_A' from the system 'my_system'
        >>> # The group is initialized with default network and analysis settings

        # Example 2: Custom Group Creation with Specific Settings
        >>> group = Group(
        ...     sys=my_system,
        ...     name='active_site',
        ...     surf_res=0.1,  # Higher resolution surface
        ...     box_size=2.0,  # Larger box size for analysis
        ...     build_type='surface'  # Only build surface network
        ... )
        >>> # This creates a group with custom settings for more detailed surface analysis

        # Example 3: Adding Atoms and Building Network
        >>> group = Group(sys=my_system, name='ligand_binding')
        >>> # Add specific atoms by their indices
        >>> group.add_balls(atom_indices=[1,2,3,4,5])
        >>> # Build the network with the added atoms
        >>> group.build()
        >>> # This creates a network representation of the selected atoms

        # Example 4: Analysis and Property Calculation
        >>> group = Group(sys=my_system, name='protein_interface')
        >>> group.add_balls(atom_indices=range(100))
        >>> group.build()
        >>> # Calculate group properties
        >>> group.get_info()
        >>> # Access calculated properties
        >>> print(f"Surface area: {group.sa}")  # Prints the total surface area
        >>> print(f"Volume: {group.vol}")       # Prints the total volume
        >>> print(f"Center of mass: {group.com}")  # Prints the center of mass coordinates
        >>> # This demonstrates how to calculate and access basic properties

        # Example 5: Layer Analysis
        >>> group = Group(sys=my_system, name='membrane_protein')
        >>> group.build()
        >>> # Perform layer-based analysis
        >>> group.get_layers(max_layers=3)
        >>> # Access layer information
        >>> print(f"Layer atoms: {group.layer_atoms}")
        >>> print(f"Layer surfaces: {group.layer_surfs}")
        >>> # This shows how to analyze the structure in layers

        # Example 6: Export and Visualization
        >>> group = Group(sys=my_system, name='protein_complex')
        >>> group.build()
        >>> # Export group data for visualization
        >>> group.export(output_directory='./visualization')
        >>> # This exports the group data for visualization in external tools
    """
    def __init__(self, sys, name=None, atoms=None, molecules=None, chains=None, residues=None,
                 settings=None, build_net=False, surf_res=0.2, box_size=1.5, max_vert=40, build_type='all', net=None,
                 net_type='aw', surf_col='plasma', surf_scheme='mean', num_splits=None, print_metrics=True,
                 scheme_factor='log', make_net=True, verts=None, vert_col='red', edge_col='grey',
                 output_directory=None):
        # System attributes
        self.sys = sys                  # Network            :    Network of the System
        self.name = name                # Name               :    Name of the group
        self.dir = output_directory     # Directory          :    Directory holding the group export info

        # Network objects attributes
        self.net = net                  # Networks           :    List of Network type objects in the group
        self.verts = verts
        self.ball_ndxs = []             # Group indexes      :    List of the indices that are included in the solve
        self.settings = settings        # Settings           :    List of network settings corresponding to the networks

        # System level classifications involved in the group (must be full)
        self.atms = atoms               # Atoms              :    List of Atoms in the group (Basically spheres)
        self.mols = molecules           # Molecule           :    List of Molecules in the group
        self.chns = chains              # Chains             :    List of molecule objects in the group
        self.rsds = residues            # Residues           :    List of residue objects in the group

        # Analysis attributes
        self.sa = 0                     # Surface Area       :    The surface area of the outer surfaces of the body
        self.vol = 0                    # Volume             :    The volume of the group's atom's cells
        self.density = 0                # Atom vol/space     :    The sum of all the atoms volumes / the total space
        self.mass = 0                   # Mass               :    Mass of the atoms in the group, if foam mass=1 for r=1
        self.com = [0, 0, 0]            # Center of Mass     :
        self.vdw_vol = 0
        self.vdw_com = [0, 0, 0]
        self.spatial_moment = [[0]]
        self.moi = [[0]]

        # Layer attributes
        self.layer_atoms = None         # Layer Atoms        :    List of lists of atoms corresponding to layers
        self.layer_verts = None         # Layer Vertices     :    List of lists of vertices arranged by layer
        self.layer_edges = None         # Layer Edges        :    List of lists of edges arranged by layer
        self.layer_surfs = None         # Layer Surfaces     :    List of lists of surfaces corresponding to layers
        self.layer_info = None          # Layer Information  :    List of information (atoms, SA, vol) for each layer

        # Set the output directory
        if self.sys.files['dir'] is None:
            self.sys.set_output_directory()

        # Get the settings
        self.get_settings(surf_res=surf_res, surf_col=surf_col, surf_scheme=surf_scheme, max_vert=max_vert,
                          box_size=box_size, net_type=net_type, build_type=build_type, num_splits=num_splits,
                          scheme_factor=scheme_factor, print_metrics=print_metrics, ball_type=sys.type,
                          sys_dir=sys.files['dir'], foam_box=sys.foam_box, vert_col=vert_col, edge_col=edge_col)

        # Set the name
        if self.name is None:
            self.set_name()

        # Process the inputs
        self.process_inputs()

        # Set the verts
        if verts is not None and type(verts) == str:
            self.verts = read_verts(self, verts)

        # Make the network
        if make_net:
            self.make_net(verts)

        # Build the Networks
        if build_net:
            self.build()

    def get_settings(self, surf_res=0.2, surf_col='plasma', surf_scheme='mean', scheme_factor='log', max_vert=40,
                     box_size=1.5, net_type='aw', build_type='all', num_splits=1, print_metrics=True, ball_type=None,
                     sys_dir=None, foam_box=None, vert_col='red', edge_col='grey'):
        """
        Sets the settings for the network building
        """
        # Set up the default values
        defaults = {'surf_res': surf_res, 'surf_col': surf_col, 'surf_scheme': surf_scheme, 'max_vert': max_vert,
                    'box_size': box_size, 'net_type': net_type, 'build_type': build_type, 'num_splits': num_splits,
                    'print_metrics': print_metrics, 'ball_type': ball_type, 'sys_dir': sys_dir, 'foam_box': foam_box,
                    'atom_rad': None, 'scheme_factor': scheme_factor, 'vert_col': vert_col, 'edge_col': edge_col}
        # Create the settings dictionary
        if self.settings is None:
            self.settings = defaults
        # Set the settings to their default values
        for setting in defaults:
            if setting not in self.settings or self.settings[setting] is None:
                self.settings[setting] = defaults[setting]

    def set_name(self):
        """
        We are looking for a name that adequately describes the group. For lists of elements > 1, they get a list of
        indices rather than their actual names
        """
        # Set up the names list that is going to be combined
        names = []
        # Get the residue names
        if self.rsds is not None and len(self.rsds) <= 2:
            for res in self.rsds:
                names.append(res.name + str(res.seq))
            self.name = '_'.join(names)
        elif self.sys.name is not None:
            self.name = self.sys.name + '_group_0'
        else:
            self.name = 'group_0'
    # Process inputs method. Goes through the atoms, residues and molecules provided in the group
    def process_inputs(self):
        """
        Processes the inputs to the group and interprets them into atoms
        :return: Sets uo the group for interpretation
        """
        # Set the atoms
        self.atms = self.atms if self.atms is not None else []
        self.rsds = self.rsds if self.rsds is not None else []
        self.chns = self.chns if self.chns is not None else []
        self.mols = self.mols if self.mols is not None else []
        # Check for empty groups
        if len(self.atms + self.rsds + self.chns + self.mols) == 0:
            self.rsds = [i for i in range(len(self.sys.residues))]
        # Add the provided atoms to the self.atoms list
        self.add_balls(self.atms)
        for resid in self.rsds:
            if isinstance(resid, int):
                self.add_balls(self.sys.residues[resid].atoms)
            else:
                self.add_balls(resid.atoms)
        for chain in self.chns:
            if isinstance(chain, int):
                self.add_balls(self.sys.chains[chain].atoms)
            else:
                self.add_balls(chain.atoms)
        # Add the residues and chains to the group
        if self.net is not None and 'res' in self.net.atoms:
            for atom in self.atms:
                if self.net.atoms['res'][atom] not in self.rsds:
                    self.rsds.append(self.net.atoms['res'][atom])
        if self.net is not None and 'chn' in self.net.atoms:
            for atom in self.atms:
                if self.sys.net.atoms['chn'][atom] not in self.chns:
                    self.chns.append(self.net.atoms['chn'][atom])
        # Add a Name If none was provided
        if self.name is None:
            # Or if the group is not in the systems list of groups
            if self not in self.sys.groups:
                # Add the group
                self.sys.groups.append(self)
            # Set the name
            self.name = '{}_group_{}'.format(self.sys.name, self.sys.groups.index(self))

    def make_net(self, verts=None):
        """
        Creates the network without an obligation to necessarily make it
        """
        self.net = Network(locs=self.sys.balls['loc'], rads=self.sys.balls['rad'], group=self.ball_ndxs,
                           settings=self.settings, sort_balls=True, masses=self.sys.balls['mass'], verts=verts)

    def build(self, verts=None):
        """
        Allows user to build the network from the system object.
        """
        self.get_settings()
        if self.net is None:
            self.make_net(verts)
        self.net.build()

    def add_balls(self, ball_list):
        """
        Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
        :param ball_list: List of atom objects expected to be added to the group
        :return: The group will have the new atoms integrated
        """
        add_balls(self, ball_list)

    def get_info(self):
        """
        Gets the info for the group to be able to make an output file with said information and also sorts the network
        """
        get_info(self)

    def get_layers(self, max_layers=50, group_resids=True, build_surfs=True):
        """
        Gets the surrounding layers of the group. Requires the whole network be built
        :param max_layers: The number of layers to go out into the SOL
        :param group_resids: Bool determining whether to keep residues together or not
        :param build_surfs: Bool determining whether to build the surfaces in the network
        :return: All layers with vertices less than the maximum number of layers will be integrated
        """
        get_layers(self, max_layers, group_resids, build_surfs)

    def exports(self, all_=False, atoms=False, atom_surfs=False, atom_edges=False, atom_verts=False, surfs=False,
                sep_surfs=False, shell_surfs=False, edges=False, sep_edges=False, shell_edges=False, verts=False,
                sep_verts=False, shell_verts=False, layers=-1, info=False, surr_atoms=False, logs=False,
                ext_atoms=False, concave_colors=False):
        """
        Exports specified export types for the group
        :param all_: All possible exports for the group will be exported to the group directory
        :param atoms: Exports a new pdb file containing only the atoms of the group
        :param shell: Exports the outer surfaces of the group
        :param surfs: Exports all surfaces in the group as one object
        :param sep_surfs: Exports all surfaces in the group as separate files, named by their atoms
        :param layers: Exports all layers surrounding the group, unless num_layers is specified
        :param info: Exports the information for the group
        :param iface: Exports the interface for the group, bff must be specified first
        :param verts: Exports the vertices of the group as an off file
        :param surr_atoms: Exports the atoms directly surrounding the group (residues intact)
        :param ext_atoms: Exports the outermost atoms in the group's set of atoms (must be a part of shell)
        :param edges: Exports all edges for the group
        :param concave_colors: Exports the concave colors for the surfaces. Default is False
        :return: The specified export is placed in the group's directory
        """
        group_exports(self, all_=all_, atoms=atoms, atom_surfs=atom_surfs, atom_edges=atom_edges, atom_verts=atom_verts,
                      surfs=surfs, sep_surfs=sep_surfs, shell_surfs=shell_surfs, edges=edges, sep_edges=sep_edges,
                      shell_edges=shell_edges, verts=verts, sep_verts=sep_verts, shell_verts=shell_verts, layers=layers,
                      info=info, surr_atoms=surr_atoms, logs=logs, ext_atoms=ext_atoms, concave_colors=concave_colors)
