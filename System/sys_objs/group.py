from System.sys_funcs.output import *


class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, sys, atoms=None, name=None, mols=None, residues=None, indices=None, bff=None):

        self.sys = sys                 # Network            :    Network of the System
        self.atoms = atoms             # Atoms              :    List of Atom type objects in the group
        self.mols = mols               # Molecules          :    List of molecule objects in the group
        self.resids = residues         # Residues           :    List of residue objects in the group
        self.ndxs = indices            # Indices            :    List of index objects in the group
        self.atom_ndxs = []            # Atom indices       :    List of atom indices for checking against
        self.name = name               # Name               :    Name of the group
        self.dir = None                # Directory          :    Directory holding the group export info

        self.surfs = None              # Surfaces           :    All surfaces associated with the group
        self.surf_ndxs = []            # Surface indices    :    Atom indices of the surfaces associated with the group
        self.sa = None                 # Surface area       :    The surface area of the outer surfaces of the body
        self.vol = None                # Volume             :    The volume of the group's atom's cells
        self.layer_atoms = None        # Layer atoms        :    List of lists of atoms corresponding to layers
        self.layer_surfs = None        # Layer Surfaces     :    List of lists of surfaces corresponding to layers
        self.layer_info = None         # Layer Information  :    List of information (atoms, SA, vol) for each layer

        self.bff = bff                 # BFF                :    Other group used for comparison
        self.iface_surfs = None        # Interface surfaces :    Surfaces that make the interface
        self.iface_atoms = None        # Interface atoms    :    Atoms in the group in the interface
        self.iface_sa = None           # Surface area       :    Surface area of the interface

        self.process_inputs()

    # Get surfaces method. Finds and sorts all surfaces in the group without needing to calculate them
    def get_surfs(self):
        # Reset the surfaces lists
        self.surfs, self.surf_ndxs = [], []
        # Go through the atoms in the group
        for atom in self.atoms:
            # Go through the surfaces in the atoms list of surfaces
            for surf in atom.surfs:
                # Get the index of the surface
                surf_ndx = ndx_search(self.surf_ndxs, surf.ndx)
                # Check if the surface has been added yet or not
                if surf_ndx < len(self.surf_ndxs) and self.surf_ndxs[surf_ndx] != surf.ndx:
                    # Insert the index and the surfaces in their correct place
                    self.surfs.insert(surf_ndx, surf)
                    self.surf_ndxs.insert(surf_ndx, surf.ndx)

    # Build surfaces method. Checks the surfaces for points and allows for rebuilds surfaces with incorrect resolutions
    def build_surfs(self, resolution=None, surfs=None, name=""):
        # Get the list of surfaces
        if surfs is not None:
            group_surfs = surfs
        # Build all surfaces in the group
        else:
            #
            if self.surfs is None:
                self.get_surfs()
            group_surfs = self.surfs
        # Get the resolution
        if resolution is None:
            resolution = self.sys.net.surf_res
        # Set up the build surfaces list
        build_surfs = []
        # Go through the list of build surfaces checking for
        for surf in group_surfs:
            # Check if the resolution is different from the set resolution or the surface has no points
            if surf.res != resolution or surf.points is None:
                build_surfs.append(surf)
        # Build the surfaces
        for i in range(len(build_surfs)):
            print("\rbuilding " + name + " surfaces " + " " * (len(str(len(surfs) - 1)) - len(str(i + 1))) + str(i + 1)
                  + "/" + str(len(surfs)) + "                   ", end="")
            build_surfs[i].build(res=resolution)

    # Add atoms method. Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
    def add_atoms(self, atom_list):
        # Check to see if the atoms list has been instantiated
        if self.atoms is None:
            self.atoms = []
        # Go through the atom_list
        for atom in atom_list:
            atom_ndx = ndx_search(self.atom_ndxs, atom.num)
            # Check to see if we have found this surface before
            if atom_ndx >= len(self.atom_ndxs) or self.atoms[atom_ndx] != atom.num:
                self.atoms.insert(atom_ndx, atom)
                self.atom_ndxs.insert(atom_ndx, atom.num)

    # Process inputs method. Goes through the atoms, residues and molecules provided in the group
    def process_inputs(self, atoms=None, mols=None, resids=None):
        # Set up the atoms list if needed
        if self.atoms is None:
            self.atoms = []
        # Add the provided atoms to the self.atoms list
        if atoms is not None:
            self.add_atoms(atoms)
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

    # Get information method. Gathers information about the group and stores it in a dictionary
    def get_info(self, iface_info=True):
        # Reset the group's data attributes
        self.sa, self.vol = 0, 0
        # Get the volume of the group
        for atom in self.atoms:
            # Check to see that the atom's volume is not 0
            if atom.vol is None or atom.vol == 0:
                atom.calc_vol()
            self.vol += atom.vol
        # Check to see if the first layer has been calculated
        if self.layer_surfs is None or self.layer_surfs == []:
            # Calculate the first layer
            self.get_layers(max_layers=1)
        # Go through the surfaces in the first layer
        for surf in self.layer_surfs[0]:
            # Add the surface area
            self.sa += surf.sa
        # Check to see if there is an interface in play
        if self.bff is not None and iface_info:
            self.get_iface()

    def get_iface(self, bff=None):
        # Set the bff
        if bff is not None:
            self.bff = bff
        # Reset the interface attributes for the group and it's bff
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

    # Get layers method. Gets the surrounding layers of the group
    def get_layers(self, max_layers=50, group_resids=True):
        # Make sure that the group has atoms
        if self.atoms is None:
            return
        # Set up the layer surfs and layer atoms list variables
        counter, self.layer_atoms, layer_atoms_ndxs, self.layer_surfs, self.layer_info = 0, [self.atoms, []], [[self.sys.atoms.index(_) for _ in self.atoms], []], [[]], [[0, 0]]
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
            # Check to make sure the surfaces are built in the layer
            self.build_surfs(surfs=self.layer_surfs[-1], name="layer {}".format(counter))
            # Check to see if the residues are supposed to stay together
            if group_resids:
                for atom in self.layer_atoms[-1]:
                    if atom.resid is not None:
                        # Get the atoms in the residue that are not already in the layer
                        for resid_atom in atom.resid:
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
        # Set the bff
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


    def exports(self, atoms=False, shell=False, fill=False, surfaces=False, layers=False, num_layers=50, info=False, iface=False):
        # Create the output directory inside the system's directory
        if self.dir is None:
            self.dir = self.sys.dir + "/" + self.name
            os.mkdir(self.dir)
        os.chdir(self.dir)
        # If the user wants to export the atoms for the group
        if atoms:
            write_pdb(atoms=self.atoms, name=self.name, sys=self.sys)
        # If the user wants to export the shell for the group
        if shell:
            if self.layer_surfs is None:
                # Get the first layer
                self.get_layers(max_layers=1)
                # noinspection PyUnresolvedReferences
                write_surfs(surfs=self.layer_surfs[0], file_name=self.name)
        # If the user wants a filled shell for the group
        if fill:
            self.build_surfs()
            write_surfs(surfs=self.surfs, file_name=self.name + "_fill")
        # If the user wants separate surfaces for the group
        if surfaces:
            self.build_surfs()
            os.mkdir(self.dir + "/surfaces")
            os.chdir(self.dir + "/surfaces")
            for surf in self.surfs:
                write_surfs([surf], file_name="_".join([str(_) for _ in surf.ndx]))
            os.chdir(self.dir)
        # If the user wants layers
        if layers:
            # First check to see if the number of layers is greater than 1
            if self.layer_atoms is None or len(self.layer_atoms) <= 1:
                self.get_layers(max_layers=num_layers)
            # Create the layers directory
            os.mkdir(os.getcwd() + "/layers")
            os.chdir(os.getcwd() + "/layers")
            # Create the layer and atoms files
            for i in range(len(self.layer_surfs)):
                write_pdb(self.layer_atoms[i + 1], name=str(i) + "_atoms", sys=self.sys)
                write_surfs(self.layer_surfs[i], file_name=str(i) + "_surfs")
            # If the user wants info and layers create a layers info file
            if info:
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
        if iface and self.bff is not None:
            self.get_iface()
            self.export_iface([self, self.bff], info_file=info)
        # If the user wants a full information file on the group
        if info:
            self.get_info()
            info = open("cell_" + self.name + "_info.txt", 'w')
            info.write(self.name + " body: \n")
            info.write("Number of atoms: " + str(len(self.atoms)) + "\n")
            info.write("Volume: " + str(self.vol) + "\n")
            info.write("Surface Area: " + str(self.sa) + "\n")
            info.close()
        os.chdir("..")
        # Change back to the system directory
        os.chdir(self.sys.dir)
