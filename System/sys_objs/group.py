from System.sys_funcs.calcs import ndx_search
from System.sys_funcs.output import write_surfs, write_pdb, write_verts, write_edges
import os


class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, sys, atoms=None, verts=None, edges=None, surfs=None, name=None, mols=None, residues=None,
                 indices=None, bff=None, vert_color='Reds', vert_scheme='shell', edge_color='white', edge_scheme=None,
                 surf_color='plasma', surf_scheme='dist'):

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
        self.mols = mols                # Molecules          :    List of molecule objects in the group
        self.resids = residues          # Residues           :    List of residue objects in the group
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
        self.bff = bff                 # BFF                :    Other group used for comparison
        self.iface_surfs = None        # Interface surfaces :    Surfaces that make the interface
        self.iface_atoms = None        # Interface atoms    :    Atoms in the group in the interface
        self.iface_sa = None           # Surface area       :    Surface area of the interface

        self.process_inputs(atoms=atoms)

    # Process inputs method. Goes through the atoms, residues and molecules provided in the group
    def process_inputs(self, atoms=None, mols=None, resids=None):
        """
        Processes the inputs to the group and interprets them into atoms
        :param atoms: List of atom objects to be added to the group
        :param mols: List of molecule objects to be added to the group
        :param resids: List of residue objects to be added to the group
        :return: Sets uo the group for interpretation
        """
        # Set up the atoms list if needed
        if self.atoms is None:
            self.atoms = []
        # Add the provided atoms to the self.atoms list
        if atoms is not None:
            self.add_atoms(self.atoms)
        # If mols were provided and not entered into the group add them
        if mols is not None and (self.mols is None or len(self.mols) < mols):
            self.mols = mols
        elif self.mols is None:
            self.mols = []
        # Add the molecule atoms
        for mol in self.mols:
            self.add_atoms(mol.atoms)
        # If residues were provided and not entered into the group add them
        if resids is not None and (self.resids is None or len(self.resids) < resids):
            self.resids = resids
        elif self.resids is None:
            self.resids = []
        # Add the residue atoms
        for residue in self.resids:
            self.add_atoms(residue.atoms)
        # Get the surfaces
        self.get_surfs()

    def get_surfs(self):
        """
        Finds and sorts all surfaces in the group without calculating them
        :return: The group will have its surfaces sorted and non-redundant
        """
        # Reset the surfaces lists
        self.surfs, self.surf_ndxs = [], []
        # Go through the atoms in the group
        for atom in self.atoms:
            # Go through the surfaces in the atoms list of surfaces
            for surf in atom.surfs:
                # Get the index of the surface
                surf_ndx = ndx_search(self.surf_ndxs, surf.ndx)
                # Check if the surface has been added yet or not
                if surf_ndx >= len(self.surf_ndxs) or self.surf_ndxs[surf_ndx] == surf.ndx:
                    # Insert the index and the surfaces in their correct place
                    self.surfs.insert(surf_ndx, surf)
                    self.surf_ndxs.insert(surf_ndx, surf.ndx)

    def get_edges(self):
        """
        Finds and sorts the edges in group
        :return: The group will have all edge objects associated with it sirted and non-redundant
        """
        # Reset the surfaces lists
        self.edges, self.edge_ndxs = [], []
        # Go through the surfaces in the atoms list of surfaces
        for edge in self.sys.net.edges:
            # Get the index of the edge
            edge_ndx = ndx_search(self.edge_ndxs, edge.ndx)
            # Check if the edge has been added yet or not
            if edge_ndx >= len(self.edge_ndxs) or self.edge_ndxs[edge_ndx] == edge.ndx:
                # Insert the index and the surfaces in their correct place
                self.edges.insert(edge_ndx, edge)
                self.edge_ndxs.insert(edge_ndx, edge.ndx)

    def get_verts(self):
        """
        Finds and sorts all the vertices in the group
        :return: The groups vertices are sorted and non-redundant
        """
        # Reset the surfaces lists
        self.verts, self.vert_ndxs = [], []
        self.atom_ndxs = [_.num for _ in self.atoms]
        self.atom_ndxs.sort()
        # Go through the surfaces in the atoms list of surfaces
        for vert in self.sys.net.verts:
            # Get the index of the edge
            vert_ndx = ndx_search(self.vert_ndxs, vert.ndx)
            # Check if the edge has been added yet or not
            if vert_ndx >= len(self.vert_ndxs) or self.vert_ndxs[vert_ndx] == vert.ndx:
                # Insert the index and the surfaces in their correct place
                self.verts.insert(vert_ndx, vert)
                self.vert_ndxs.insert(vert_ndx, vert.ndx)

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
            if surf.res != resolution:
                build_surfs.append(surf)
            # Check if there is any sign of missing points or triangles
            elif surf.points is None or surf.tris is None or len(surf.points) <= 2 or len(surf.tris) == 0:
                # If it is possible to load the file
                if surf.file is not None and surf.file not in ["", " "]:
                    test = surf.read_file()
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
            build_surfs[i].build(res=resolution, flat=self.sys.net.flat_Del)
            if build_surfs[i].file is None:
                write_surfs([build_surfs[i]], "_".join([str(_) for _ in build_surfs[i].ndx]), )
        # Change back
        os.chdir(self.sys.dir)

    def add_atoms(self, atom_list):
        """
        Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
        :param atom_list: List of atom objects expected to be added to the group
        :return: The group will have the new atoms integrated
        """
        # Check to see if the atoms list has been instantiated
        if self.atoms is None:
            self.atoms = []
        # Go through the atom_list
        for atom in atom_list:
            atom_ndx = ndx_search(self.atom_ndxs, atom.num)
            # Check to see if we have found this atom before
            if atom_ndx >= len(self.atom_ndxs) or self.atoms[atom_ndx] == atom.num:
                self.atoms.insert(atom_ndx, atom)
                self.atom_ndxs.insert(atom_ndx, atom.num)

    def get_info(self, iface_info=True):
        """
        Gathers information about the group and stores it in a dictionary
        :param iface_info:
        :return:
        """
        # Reset the group's data attributes
        self.sa, self.vol = 0, 0
        # Get the volume of the group
        for atom in self.atoms:
            # Check to see that the atom's volume is not 0
            if atom.vol is None or atom.vol == 0:
                # Calculate the volume of the atom
                atom.calc_vol()
            # Add the volume to that of the group
            self.vol += atom.vol
        # Check to see if the first layer has been calculated
        for surf in self.layer_surfs[0]:
            # Check that the surface has a surface area
            if surf.sa is None or surf.sa == 0:
                # Get the surface area for the surface
                surf.calc_sa()
            # Add the surface area
            self.sa += surf.sa

    def get_iface(self, bff=None):
        # Set the bff
        if bff is not None:
            self.bff = bff
        # Reset the interface attributes for the group, and it's bff
        self.iface_atoms, self.bff.iface_atoms, self.iface_surfs, self.bff.iface_surfs = [], [], [], []
        self.iface_sa = 0
        # Go through the atoms in the group
        for atom in self.atoms:
            # Check to see if the atom is in the bff's list of atoms
            if atom.num in self.bff.atom_ndxs:
                continue
            # Go through the surfaces in the atom's list of surfaces
            for surf in atom.surfs:
                # Check for an interface surf
                if (surf.ndx[0] in self.atom_ndxs and surf.ndx[1] in self.bff.atom_ndxs) or \
                   (surf.ndx[1] in self.atom_ndxs and surf.ndx[0] in self.bff.atom_ndxs):
                    # Get the other atom from the surface's atoms
                    other_atom = [_ for _ in surf.atoms if _ != atom][0]
                    # Add the first atom to the group's list of interface atoms
                    self.iface_atoms.append(atom)
                    self.bff.iface_atoms.append(other_atom)
                    # Add the surface to the list of interface surfs and add the surface area of the surface
                    self.iface_surfs.append(surf)
                    self.bff.iface_surfs.append(surf)
                    self.iface_sa += surf.sa
        # Set the bff's surface area
        self.bff.iface_sa = self.iface_sa

    def get_layers(self, max_layers=50, group_resids=True, build_surfs=True):
        """
        Gets the surrounding layers of the group. Requires the whole network be built
        :param max_layers: The number of layers to go out into the SOL
        :param group_resids: Bool determining whether to keep residues together or not
        :param build_surfs: Bool determining whether to build the surfaces in the network
        :return: All layers with vertices less than the maximum number of layers will be bintegrated
        """
        # Make sure that the group has atoms
        if self.atoms is None:
            return
        # Set up the layer surfs and layer atoms list variables
        counter = 0
        self.layer_atoms = [self.atoms, []]
        layer_atoms_ndxs = [[_.num for _ in self.atoms], []]
        self.layer_surfs = [[]]
        self.layer_verts = [[]]
        self.layer_edges = [[]]
        self.layer_info = [[0, 0]]
        # Set up the loop to keep adding layers
        while counter < max_layers:
            # Go through the atoms in the last layer
            for atom in self.layer_atoms[-2]:
                # Go through the surfaces in the atom's list of surfaces
                for surf in atom.surfs:
                    if surf in self.layer_surfs[-1] or (len(self.layer_surfs) >= 2 and surf in self.layer_surfs[-2]):
                        continue
                    elif surf.ndx[0] in layer_atoms_ndxs[-2] and surf.ndx[1] in layer_atoms_ndxs[-2]:
                        continue
                    self.layer_surfs[-1].append(surf)
                    # Add the vertices
                    for vert in surf.verts:
                        if vert not in self.layer_verts[-1]:
                            self.layer_verts[-1].append(vert)
                    # Add the edges
                    for edge in surf.edges:
                        if edge not in self.layer_edges[-1]:
                            self.layer_edges[-1].append(edge)
                    # Get the index of the surface
                    surf_ndx = ndx_search(self.surf_ndxs, surf.ndx)
                    # Check if the surface has been added yet or not
                    if surf_ndx < len(self.surf_ndxs) and self.surf_ndxs[surf_ndx] != surf.ndx:
                        # Insert the index and the surfaces in their correct place
                        self.surfs.insert(surf_ndx, surf)
                        self.surf_ndxs.insert(surf_ndx, surf.ndx)
                    # Sort the surface's atoms inside or out
                    if surf.ndx[0] in layer_atoms_ndxs[-2] and surf.ndx[1] not in layer_atoms_ndxs[-2]:
                        self.layer_atoms[-1].append(self.sys.atoms[surf.ndx[1]])
                        layer_atoms_ndxs[-1].append(surf.ndx[1])
                    if surf.ndx[1] in layer_atoms_ndxs[-2] and surf.ndx[0] not in layer_atoms_ndxs[-2]:
                        self.layer_atoms[-1].append(self.sys.atoms[surf.ndx[0]])
                        layer_atoms_ndxs[-1].append(surf.ndx[0])
            if build_surfs:
                # Check to make sure the surfaces are built in the layer
                self.build_surfs()
            # Check to see if the residues are supposed to stay together
            if group_resids:
                for atom in self.layer_atoms[-1]:
                    if atom.res is not None:
                        # Get the atoms in the residue that are not already in the layer
                        for resid_atom in atom.res.atoms:
                            # Check if the atom is in the layer or not
                            if resid_atom not in self.layer_atoms[-1]:
                                self.layer_atoms[-1].append(resid_atom)
            # Get the surface area and volume for the layer
            for atom in self.layer_atoms[-1]:
                # Add the volume to the current layer's volume
                self.layer_info[-1][0] += atom.vol
            # Get the surface area of the layer
            for surf in self.layer_surfs[-1]:
                # Add the surface area to the current layer's surface area
                self.layer_info[-1][1] += surf.sa
            # If there is nothing to add leave the layers loop
            if len(self.layer_surfs[-1]) == 0:
                self.layer_surfs.pop()
                break
            # Create the new layer lists
            self.layer_surfs.append([])
            self.layer_atoms.append([])
            self.layer_edges.append([])
            self.layer_verts.append([])
            self.layer_info.append([0, 0])
            layer_atoms_ndxs.append([])
            counter += 1

    def export_iface(self, bff=None, info_file=False, interface_atoms=False):
        """
        Exports the information from the given interface as a txt file
        :param bff: Group for the interface with self
        :param info_file: Whether to export a txt file with info on the interface or not
        :param interface_atoms: Whether to export a pdb file with the atoms around the interface or not
        :return:
        """
        # Check to see if there is a bff or not
        if bff is not None:
            self.bff = bff
        # Check to see that the interface has been calculated
        if self.iface_surfs is None or self.iface_atoms is None:
            self.get_iface()
        # Set the interface name
        interface_name = self.name + "_" + self.bff.name + "_interface"
        # Move to the output directory
        os.chdir(self.sys.dir)
        # Create and move to the interface directory
        os.mkdir(os.getcwd() + "/" + interface_name)
        os.chdir(os.getcwd() + "/" + interface_name)
        # Write the surfaces for the interface
        write_surfs(self.iface_surfs, interface_name)
        # Check to see of the user wants to export the interface's atoms
        if interface_atoms:
            # Get the two sets of interface atoms
            write_pdb(self.iface_atoms, interface_name + "_" + self.name + "_atoms", self.sys)
            write_pdb(self.bff.iface_atoms, interface_name + "_" + self.bff.name + "_atoms", self.bff.sys)
        # Check to see if the user wants to export the interface's information
        if info_file:
            info = open(interface_name + "_info.txt", 'w')
            info.write("Interface between " + self.name + " and " + self.bff.name + " : \n")
            info.write("Number of Surfaces: " + str(len(self.iface_surfs)))
            info.write("Surface Area: " + str(self.iface_sa))
            info.close()

    def exports(self, all_=False, atoms=False, shell=False, fill=False, surfaces=False, layers=False, num_layers=50,
                info=False, iface=False, verts=False, surr_atoms=False, ext_atoms=False, shell_edges=False,
                shell_verts=False, edges=False):
        """
        Exports specified export types for the group
        :param all_: All possible exports for the group will be exported to the group directory
        :param atoms: Exports a new pdb file contasining only the atoms of the group
        :param shell: Exports the outer surfaces of the group
        :param fill: Exports all surfaces in the group as one object
        :param surfaces: Exports all surfaces in the group as seperate files, named by their atoms
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
        # Get the surfaces if they haven't been got
        if self.surfs is None or len(self.surfs) == 0:
            self.build_surfs()
        # Get the surface coloring scheme
        scheme = ""
        if self.surf_scheme is not None:
            scheme = "_" + self.surf_scheme
        # Create the output directory inside the system's directory
        if self.dir is None:
            i = 1
            my_dir = self.sys.dir + "/" + self.name
            first = True
            while os.path.exists(my_dir):
                if first:
                    my_dir += "__"
                    first = False
                my_dir = my_dir[:-(1 + len(str(i)))] + '_' + str(i)
                i += 1
            self.dir = my_dir
            os.mkdir(self.dir)
        os.chdir(self.dir)
        # If the user wants to export the atoms for the group
        if atoms or all_:
            write_pdb(atoms=self.atoms, name=self.name + "_atoms", sys=self.sys)
        # If the user wants to export the shell for the group
        if shell or all_:
            if self.layer_surfs is None:
                # Get the first layer
                self.get_layers(max_layers=1)
            # noinspection PyUnresolvedReferences
            if self.layer_surfs is not None and len(self.layer_surfs) > 0:
                write_surfs(surfs=self.layer_surfs[0], file_name=self.name + "_shell" + scheme, directory=self.dir, color_map=self.surf_color, color_scheme=self.surf_scheme)
        # If the user wants a filled shell for the group
        if fill or all_:
            self.build_surfs()
            write_surfs(surfs=self.surfs, file_name=self.name + "_fill" + scheme, directory=self.dir, color_map=self.surf_color, color_scheme=self.surf_scheme)
        # If the user wants separate surfaces for the group
        if surfaces or all_:
            i = 1
            my_dir = self.dir + "/surfaces"
            while os.path.exists(my_dir):
                if my_dir[-1] == 's':
                    my_dir += '__'
                my_dir  = my_dir[:-2] + str(i)
                i += 1
            os.mkdir(my_dir)
            os.chdir(my_dir)
            for surf in self.surfs:
                write_surfs([surf], file_name="_".join([str(_) for _ in surf.ndx]), directory=my_dir, color_map=self.surf_color, color_scheme=self.surf_scheme)
            os.chdir(self.dir)
        # If the user wants layers
        if layers or all_:
            # First check to see if the number of layers is greater than 1
            if self.layer_atoms is None or len(self.layer_atoms) <= 1:
                self.get_layers(max_layers=num_layers)
            # Create the layers directory
            i = 1
            my_dir = os.getcwd() + "/layers"
            while os.path.exists(my_dir):
                if my_dir[-1] == 's':
                    my_dir += '__'
                my_dir = my_dir[:-2] + str(i)
                i += 1
            os.mkdir(my_dir)
            os.chdir(my_dir)
            # Create the layer and atoms files
            for i in range(len(self.layer_surfs)):
                write_pdb(self.layer_atoms[i + 1], name=str(i) + "_atoms", sys=self.sys)
                write_surfs(self.layer_surfs[i], file_name=str(i) + "_surfs", color_map=self.surf_color, color_scheme=self.surf_scheme)
            # If the user wants info and layers create a layers info file
            if info or all_:
                self.get_info()
                # Create the information file
                info = open(self.name + "_layer_info.txt", 'w')
                info.write(self.name + " body: \n")
                # Go through the layers in the group's layers
                for i in range(len(self.layer_surfs)):
                    info.write("Number of atoms: " + str(len(self.layer_atoms[i])) + "\n")
                    info.write("Volume: " + str(self.layer_info[i][0]) + "\n")
                    info.write("Surface Area: " + str(self.layer_info[i][1]) + "\n")
                info.close()
            # Change back to the group directory
            os.chdir(self.dir)
        # If the user wants to export the interface
        if (iface or all_) and self.bff is not None:
            self.get_iface()
            self.export_iface([self, self.bff], info_file=info)
        # If the user wants a full information file on the group
        if info or all_:
            os.chdir(self.dir)
            self.get_info()
            info = open("cell_" + self.name + "_info.txt", 'w')
            info.write(self.name + " body: \n")
            info.write("Number of atoms: " + str(len(self.atoms)) + "\n")
            info.write("Volume: " + str(self.vol) + "\n")
            info.write("Surface Area: " + str(self.sa) + "\n")
            info.close()
        if verts or all_:
            if self.verts is None:
                self.get_verts()
            write_verts(verts=self.verts, file_name=self.name + "_verts", directory=self.dir)
        if surr_atoms or all_:
            if self.layer_surfs is None:
                # Get the first layer
                self.get_layers(max_layers=1)
            # write the surrounding atoms
            write_pdb(atoms=self.layer_atoms[1], name=self.name + "_surr_atoms", directory=self.dir)
        if ext_atoms or all_:
            if self.layer_surfs is None:
                # Get the first layer
                self.get_layers(max_layers=1)
            # write the surrounding atoms
            write_pdb(atoms=self.layer_atoms[0], name=self.name + "_ext_atoms", directory=self.dir)
        if shell_verts or all_:
            if self.layer_verts is None:
                # Get the first layer
                self.get_layers(max_layers=1, build_surfs=False)
            write_verts(self.layer_verts[0], file_name=self.name + "_shell_verts", directory=self.dir)
        if edges or all_:
            if self.edges is None:
                self.get_edges()
            write_edges(edges=self.edges, file_name=self.name + "_edges", directory=self.dir)
        if shell_edges or all_:
            if self.edges is None:
                self.get_edges()
            if self.layer_edges is None:
                self.get_layers(max_layers=1, build_surfs=False)
            write_edges(self.layer_edges[0], file_name=self.name + "_shell_edges", directory=self.dir)
        os.chdir("..")
        # Change back to the system directory
        os.chdir(self.sys.dir)
