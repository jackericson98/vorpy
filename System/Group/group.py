from System.Group.build import build_surfs
from System.Group.layers import get_layers
from System.Group.sort import get_surfs, get_edges, get_verts, add_spheres
from System.Group.export import group_exports
from System.Network.network import Network
from System.Network.split_net import split_net_slow
from System.sys_funcs.calcs.sorting import ndx_search
from System.sys_funcs.calcs.surf import calc_surf_sa
from System.sys_funcs.output.net import add_metrics
import numpy as np


class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, sys, spheres=None, atoms=None, name=None, molecules=None, chains=None, residues=None, bff=None,
                 settings=None, build_net=False, surf_res=0.2, box_size=1.5, max_vert=40, build_type='all', net=None,
                 net_type='aw', surf_col='plasma', surf_scheme='curv', num_splits=None, print_metrics=True):
        # System attributes
        self.sys = sys                  # Network            :    Network of the System
        self.name = name                # Name               :    Name of the group
        self.dir = None                 # Directory          :    Directory holding the group export info

        # Network objects attributes
        self.net = net                  # Networks           :    List of Network type objects in the group
        self.ball_ndxs = []             # Group indexes      :    List of the indices that are included in the solve
        self.settings = settings        # Settings           :    List of network settings corresponding to the networks

        # System level classifications involved in the group (must be full)
        self.atms = atoms               # Atoms              :    List of Atoms in the group (Basically spheres)
        self.mols = molecules           # Molecule           :    List of Molecules in the group
        self.chns = chains              # Chains             :    List of molecule objects in the group
        self.rsds = residues            # Residues           :    List of residue objects in the group

        # Analysis attributes
        self.sa = None                  # Surface Area       :    The surface area of the outer surfaces of the body
        self.vol = None                 # Volume             :    The volume of the group's atom's cells
        self.density = None             # Atom vol/space     :    The sum of all the atoms volumes / the total space

        # Layer attributes
        self.layer_atoms = None         # Layer Atoms        :    List of lists of atoms corresponding to layers
        self.layer_verts = None         # Layer Vertices     :    List of lists of vertices arranged by layer
        self.layer_edges = None         # Layer Edges        :    List of lists of edges arranged by layer
        self.layer_surfs = None         # Layer Surfaces     :    List of lists of surfaces corresponding to layers
        self.layer_info = None          # Layer Information  :    List of information (atoms, SA, vol) for each layer

        # Interface attributes
        self.bff = bff                 # BFF                 :   Other group used for comparison
        self.iface_surfs = None        # Interface Surfaces  :   Surfaces that make the interface
        self.iface_edges = None        # Interface Edges     :   Interfacial edges list
        self.iface_verts = None        # Interface Vertices  :   Interfacial vertices list
        self.iface_atoms = None        # Interface atoms     :   Atoms in the group in the interface
        self.iface_sa = None           # Surface area        :   Surface area of the interface
        self.iface_curv = None         # Interface Curvature :   Average curvature from the interface

        # Get the settings
        self.get_settings(surf_res=surf_res, surf_col=surf_col, surf_scheme=surf_scheme, max_vert=max_vert,
                          box_size=box_size, net_type=net_type, build_type=build_type, num_splits=num_splits,
                          print_metrics=print_metrics, ball_type=sys.type, sys_dir=sys.files['dir'],
                          foam_box=sys.foam_box)

        # Process the inputs
        self.process_inputs()

        # Make the Networks
        if build_net:
            self.build_network()

    def get_settings(self, surf_res=0.2, surf_col='plasma', surf_scheme='curv', max_vert=40, box_size=1.5, net_type='aw',
                     build_type='all', num_splits=1, print_metrics=True, ball_type=None,
                     sys_dir=None, foam_box=None):
        """
        Sets the settings for the network building
        """
        # Set up the default values
        defaults = {'surf_res': surf_res, 'surf_col': surf_col, 'surf_scheme': surf_scheme, 'max_vert': max_vert,
                    'box_size': box_size, 'net_type': net_type, 'build_type': build_type, 'num_splits': num_splits,
                    'print_metrics': print_metrics, 'ball_type': ball_type, 'sys_dir': sys_dir, 'foam_box': foam_box}
        # Create the settings dictionary
        if self.settings is None:
            self.settings = defaults
        # Set the settings to their default values
        for setting in self.settings:
            if self.settings[setting] is None:
                self.settings[setting] = defaults[setting]

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
        # Add the provided atoms to the self.atoms list
        self.add_atoms(self.atms)
        for resid in self.rsds:
            self.add_atoms(resid.atoms)
        for chain in self.chns:
            self.add_atoms(chain.atoms)
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

    def build_network(self, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None,
                      calc_verts=None, my_group=None, print_actions=None, num_atoms_sub_net=1000, no_split=True,
                      add_net_metrics=True, min_atom_split=1000):
        """
        Allows user to build the network from the system object.
        """
        if self.net is None:
            self.net = Network(locs=self.sys.spheres['loc'], rads=self.sys.spheres['rad'], group=self.ball_ndxs,
                               settings=self.settings)
        # Small networks and no split option
        if len(self.ball_ndxs) < num_atoms_sub_net or no_split:
            # Build the network
            self.net.build()
            # Add the metrics
            self.net.metrics['splits'] = 1
            if add_net_metrics:
                add_metrics(self)
        else:
            split_net_slow(sys=self, num_atoms_sub_net=num_atoms_sub_net, add_net_metrics=add_net_metrics,
                           min_atom_split=min_atom_split)

    def get_surfs(self):
        """
        Finds and sorts all surfaces in the group without calculating them
        :return: The group will have its surfaces sorted and non-redundant
        """
        get_surfs(self)

    def get_edges(self):
        """
        Finds and sorts the edges in group
        :return: The group will have all edge objects associated with it sorted and non-redundant
        """
        get_edges(self)

    def get_verts(self):
        """
        Finds and sorts all the vertices in the group
        :return: The groups vertices are sorted and non-redundant
        """
        get_verts(self)

    def build_surfs(self, resolution=None):
        """
        Checks the surfaces for points and allows for rebuilds surfaces with incorrect resolutions
        :param resolution: If not None all surfs without this resolution will be rebuilt
        :return: All surfaces in the group will be constructed
        """
        build_surfs(self)

    def add_atoms(self, atom_list):
        """
        Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
        :param atom_list: List of atom objects expected to be added to the group
        :return: The group will have the new atoms integrated
        """
        add_spheres(self, atom_list)

    def get_info(self):
        """
        Gathers information about the group and stores it in a dictionary
        :return:
        """
        net = self.net
        # Get the group objects
        self.get_surfs()
        self.get_edges()
        self.get_verts()
        # Reset the group's data attributes
        self.sa, self.vol, self.density = 0, 0, 0
        tot_atom_vol = 0
        # Get the volume of the group
        for i in self.ball_ndxs:
            atom = self.net.balls.iloc[i]
            if not atom['complete']:
                continue
            # Add the volume to that of the group
            self.vol += atom['vol']
            tot_atom_vol += (4/3)*np.pi*atom['rad']**3
        if self.vol > 0:
            self.density = tot_atom_vol/self.vol
        # Check to see if the first layer has been calculated
        if self.layer_surfs is None or len(self.layer_surfs) == 0:
            self.get_layers(max_layers=1)
        if len(self.layer_surfs) > 0:
            for i in self.layer_surfs[0]:
                surf = self.net.surfs.iloc[i]
                # Check that the surface has a surface area
                if surf['sa'] is None or surf['sa'] == 0:
                    # Get the surface area for the surface
                    edge_ndxss = [ndx_search(net.edge_ndxs, _) for _ in surf['sedges']]
                    edges = np.array([net.edges.iloc[_] for _ in edge_ndxss])
                    surf_sa = calc_surf_sa(edges=edges, com=np.array(surf['com']), tris=surf['tris'], points=surf['points'], flat=surf['flat'])
                else:
                    surf_sa = surf['sa']
                # Add the surface area
                self.sa += surf_sa

    def get_layers(self, max_layers=50, group_resids=True, build_surfs=True):
        """
        Gets the surrounding layers of the group. Requires the whole network be built
        :param max_layers: The number of layers to go out into the SOL
        :param group_resids: Bool determining whether to keep residues together or not
        :param build_surfs: Bool determining whether to build the surfaces in the network
        :return: All layers with vertices less than the maximum number of layers will be integrated
        """
        get_layers(self, max_layers, group_resids, build_surfs)

    def exports(self, all_=False, iface=False, atoms=False, surfs=False, sep_surfs=False, edges=False,
                sep_edges=False, verts=False, sep_verts=False, layers=-1, info=False, surr_atoms=False,
                ext_atoms=False, shell=False):
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
        :return: The specified export is placed in the group's directory
        """
        group_exports(grp=self, all_=all_, iface=iface, atoms=atoms, surfs=surfs, sep_surfs=sep_surfs, edges=edges,
                      sep_edges=sep_edges, verts=verts, sep_verts=sep_verts, layers=layers, info=info,
                      surr_atoms=surr_atoms, ext_atoms=ext_atoms,  shell=shell)
