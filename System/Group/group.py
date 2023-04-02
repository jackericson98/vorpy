from System.sys_funcs.output import write_surfs
from System.Group.layers import get_layers
from System.Group.sort import get_surfs, get_edges, get_verts, add_atoms
from System.Group.export import group_exports
from System.sys_funcs.calcs import calc_vol
from System.sys_funcs.input import read_surf_file
from System.Network.net_funcs.build_surf import build_surf
import os


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
        for atom in self.atoms:
            if atom.res not in self.residues:
                self.residues.append(atom.res)
            if atom.chn not in self.chains:
                self.chains.append(atom.chn)
        # Add a Name If none was provided
        if self.name is None:
            # If no system is given
            if self.sys is not None:
                # Or if the group is not in the systems list of groups
                if self not in self.sys.groups:
                    # Add the group
                    self.sys.groups.append(self)
                # Set the name
                self.name = '{}_group{}'.format(self.sys.name, self.sys.groups.index(self))
            else:
                self.name = self.atoms[0].name + "_" + self.atoms[-1].name
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
        # Get the resolution
        if resolution is None:
            resolution = self.sys.net.surf_res
            self.surf_res = resolution
        # Set up the build surfaces list
        build_surfs = []
        # Go through the list of build surfaces checking for
        for surf in self.surfs:
            # Check if the resolution is different from the set resolution or the surface has no points
            if surf.residue != resolution:
                build_surfs.append(surf)
            # Check if there is any sign of missing points or triangles
            elif surf.points is None or surf.tris is None or len(surf.points) <= 2 or len(surf.tris) == 0:
                # If it is possible to load the file
                if surf.file is not None and surf.file not in ["", " "]:
                    test = read_surf_file(surf)
                    if test is None:
                        build_surfs.append(surf)
                # Worst case, add the surface to the list of surfaces to be built
                else:
                    build_surfs.append(surf)
        # Create the system's surface's file if needed
        if len(build_surfs) > 0 and not os.path.exists(self.sys.dir + "/surfs"):
            os.mkdir(self.sys.dir + "/surfs")
            os.chdir(self.sys.dir + '/surfs')
        # Build the surfaces
        for i in range(len(build_surfs)):
            # Print the status of the surfaces being built
            print("\rbuilding " + self.name + " surfaces " + " " * (len(str(len(self.surfs) - 1)) - len(str(i + 1))) +
                  str(i + 1) + "/" + str(len(self.surfs)) + "                   ", end="")
            build_surf(build_surfs[i], res=resolution, flat=self.sys.net.flat_Del)
            if build_surfs[i].file is None:
                write_surfs([build_surfs[i]], "_".join([str(_) for _ in build_surfs[i].ndx]))
        # Change back
        os.chdir(self.sys.dir)

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
        # Get the group objects
        self.get_surfs()
        self.get_edges()
        self.get_verts()
        # Reset the group's data attributes
        self.sa, self.vol = 0, 0
        # Get the volume of the group
        for atom in self.atoms:
            # Check to see that the atom's volume is not 0
            if atom.vol is None or atom.vol == 0:
                # Calculate the volume of the atom
                calc_vol(atom)
            # Add the volume to that of the group
            self.vol += atom.vol
        # Check to see if the first layer has been calculated
        if len(self.layer_surfs) == 0:
            return
        for surf in self.layer_surfs[0]:
            # Check that the surface has a surface area
            if surf.sa is None or surf.sa == 0:
                # Get the surface area for the surface
                surf.calc_sa()
            # Add the surface area
            self.sa += surf.sa

    def get_layers(self, max_layers=50, group_resids=True, build_surfs=True):
        """
        Gets the surrounding layers of the group. Requires the whole network be built
        :param max_layers: The number of layers to go out into the SOL
        :param group_resids: Bool determining whether to keep residues together or not
        :param build_surfs: Bool determining whether to build the surfaces in the network
        :return: All layers with vertices less than the maximum number of layers will be integrated
        """
        get_layers(self, max_layers, group_resids, build_surfs)

    def exports(self, all_=False, atoms=False, shell=False, fill=False, surfaces=False, layers=False, num_layers=50,
                info=False, iface=False, verts=False, surr_atoms=False, ext_atoms=False, shell_edges=False,
                shell_verts=False, edges=False):
        """
        Exports specified export types for the group
        :param all_: All possible exports for the group will be exported to the group directory
        :param atoms: Exports a new pdb file containing only the atoms of the group
        :param shell: Exports the outer surfaces of the group
        :param fill: Exports all surfaces in the group as one object
        :param surfaces: Exports all surfaces in the group as separate files, named by their atoms
        :param layers: Exports all layers surrounding the group, unless num_layers is specified
        :param num_layers: Controls the number of exported layers for the group
        :param info: Exports the information for the group
        :param iface: Exports the interface for the group, bff must be specified first
        :param verts: Exports the vertices of the group as an off file
        :param surr_atoms: Exports the atoms directly surrounding the group (residues intact)
        :param ext_atoms: Exports the outermost atoms in the group's set of atoms (must be a part of shell)
        :param shell_edges: Exports only the outermost edges for the group as an OFF file
        :param shell_verts: Exports the outermost vertices for the group
        :param edges: Exports all edges for the group
        :return: The specified export is placed in the group's directory
        """
        group_exports(grp=self, all_=all_, atoms=atoms, shell=shell, fill=fill, surfaces=surfaces, layers=layers,
                      num_layers=num_layers, info=info, iface=iface, verts=verts, surr_atoms=surr_atoms,
                      ext_atoms=ext_atoms, shell_edges=shell_edges, shell_verts=shell_verts, edges=edges)
