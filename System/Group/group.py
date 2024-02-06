from System.Group.build import build_surfs
from System.Group.layers import get_layers
from System.Group.sort import get_surfs, get_edges, get_verts, add_atoms
from System.Group.export import group_exports
from System.sys_funcs.calcs.sorting import ndx_search
from System.sys_funcs.calcs.surf import calc_surf_sa
import numpy as np


class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, sys, atoms=None, verts=None, edges=None, surfs=None, name=None, chains=None, residues=None,
                 indices=None, bff=None, vert_color='Reds', vert_scheme='shell', edge_color='white', edge_scheme=None,
                 surf_color=None, surf_scheme=None):

        # System attributes
        self.sys = sys                  # Network            :    Network of the System
        self.name = name                # Name               :    Name of the group
        self.dir = None                 # Directory          :    Directory holding the group export info

        # Network objects attributes
        self.atoms = atoms              # Atoms              :    List of Atom type objects in the group
        self.verts = verts              # Vertices           :    Vertex objects involved in the group
        self.edges = edges              # Edges              :    Edge objects involved in the group's network
        self.surfs = surfs              # Surfaces           :    All surfaces associated with the group

        # Network object tracking attributes
        self.atom_ndxs = []             # Atom indices       :    List of atom indices for checking against (sorted)
        self.vert_ndxs = []             # Vertex indices     :    Tracks the vertices in the group (sorted)
        self.edge_ndxs = []             # Edge indices       :    Tracks the edges in a group (sorted)
        self.surf_ndxs = []             # Surface indices    :    Atom indices of the surfaces associated with the group

        # System level classifications involved in the group (must be full)
        self.chains = chains            # Molecules          :    List of molecule objects in the group
        self.residues = residues        # Residues           :    List of residue objects in the group
        self.ndxs = indices             # Indices            :    List of index objects in the group

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

        # Coloring schemes and maps
        self.vert_color = vert_color    # Vertex coloring    :    Color map for vertices in the group
        self.vert_scheme = vert_scheme  # Vertex Scheme      :    Color Scheme for vertex coloring in the group
        self.edge_color = edge_color    # Edge Coloring      :    Color map for the edges in the group
        self.edge_scheme = edge_scheme  # Edge Scheme        :    Color scheme for the edges in the group
        self.surf_color = surf_color    # Surface Coloring   :    Color map for the surfaces in the group
        self.surf_scheme = surf_scheme  # Surface Scheme     :    Color Scheme for the surfaces in the group
        self.surf_res = None            # Surface resolution :    Surface resolution for the surfaces in the group

        # Interface attributes
        self.bff = bff                 # BFF                 :   Other group used for comparison
        self.iface_surfs = None        # Interface Surfaces  :   Surfaces that make the interface
        self.iface_edges = None        # Interface Edges     :   Interfacial edges list
        self.iface_verts = None        # Interface Vertices  :   Interfacial vertices list
        self.iface_atoms = None        # Interface atoms     :   Atoms in the group in the interface
        self.iface_sa = None           # Surface area        :   Surface area of the interface
        self.iface_curv = None         # Interface Curvature :   Average curvature from the interface

        self.process_inputs()

    # Process inputs method. Goes through the atoms, residues and molecules provided in the group
    def process_inputs(self):
        """
        Processes the inputs to the group and interprets them into atoms
        :return: Sets uo the group for interpretation
        """
        # Set the atoms
        atoms = self.atoms if self.atoms is not None else []
        resids = self.residues if self.residues is not None else []
        chains = self.chains if self.chains is not None else []
        # Set up the atoms list if needed
        self.atoms = []
        if self.residues is None:
            self.residues = []
        if self.chains is None:
            self.chains = []
        # Add the provided atoms to the self.atoms list
        self.add_atoms(atoms)
        for resid in resids:
            self.add_atoms(resid.atoms)
        for chain in chains:
            self.add_atoms(chain.atoms)
        # Add the residues and chains to the group
        if self.sys.net is not None and 'res' in self.sys.net.atoms:
            for atom in self.atoms:
                if self.sys.net.atoms['res'][atom] not in self.residues:
                    self.residues.append(self.sys.net.atoms['res'][atom])
        if self.sys.net is not None and 'chn' in self.sys.net.atoms:
            for atom in self.atoms:
                if self.sys.net.atoms['chn'][atom] not in self.chains:
                    self.chains.append(self.sys.net.atoms['chn'][atom])
        # Add a Name If none was provided
        if self.name is None:
            # Or if the group is not in the systems list of groups
            if self not in self.sys.groups:
                # Add the group
                self.sys.groups.append(self)
            # Set the name
            self.name = '{}_group{}'.format(self.sys.name, self.sys.groups.index(self))
        # Get the surfaces
        self.get_surfs()

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
        add_atoms(self, atom_list)

    def get_info(self):
        """
        Gathers information about the group and stores it in a dictionary
        :return:
        """
        net = self.sys.net
        # Get the group objects
        self.get_surfs()
        self.get_edges()
        self.get_verts()
        # Reset the group's data attributes
        self.sa, self.vol = 0, 0
        tot_atom_vol = 0
        # Get the volume of the group
        for i in self.atoms:
            atom = self.sys.net.atoms.iloc[i]
            # Add the volume to that of the group
            self.vol += atom['vol']
            tot_atom_vol += (4/3)*np.pi*atom['rad']**3
        self.density = tot_atom_vol/self.vol
        # Check to see if the first layer has been calculated
        if self.layer_surfs is None or len(self.layer_surfs) == 0:
            self.get_layers(max_layers=1)
        if len(self.layer_surfs) > 0:
            for i in self.layer_surfs[0]:
                surf = self.sys.net.surfs.iloc[i]
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
