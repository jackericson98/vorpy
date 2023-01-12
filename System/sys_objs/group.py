from System.sys_funcs.output import *


class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, sys, atoms=None, name=None, mols=None, residues=None, my_atoms=None, ndxs=None, bff=None):

        self.sys = sys                 # Network            :    Network of the System
        self.atoms = atoms             # Atoms              :    List of Atom type objects for the edge
        self.selects = []              # Previous Selection :    List of the previously selected atoms
        self.select_strs = []          # Previous String    :    String from the last group of atoms selected
        self.name = name               # Name               :    Name of the group
        self.dir = None                # Directory          :    Directory holding the group export info

        self.surfs = None              # Surfaces           :    All surfaces associated with the group
        self.surf_ndxs = None          # Surface indices    :    Atom indices of the surfaces associated with the group
        self.body_surfs = None         # Body surfaces      :    The surfaces on the outside of the body
        self.body_sa = None            # Body surface area  :    The surface area of the outer surfaces of the body
        self.body_vol = None           # Body volume        :    The volume of the group's atom's cells
        self.outer_body_atoms = None   # Outer body atoms   :    Body atoms that create outer surfaces
        self.surr_body_atoms = None    # Surrounding atoms  :    Atoms not in the group that create surfaces with it
        self.layer_atoms = None        # Layer atoms        :    List of lists of atoms corresponding to layers
        self.layer_surfs = None        # Layer Surfaces     :    List of lists of surfaces corresponding to layers

        self.bff = bff                 # BFF                :    Other group used for comparison
        self.iface_surfs = None        # Interface surfaces :    Surfaces that make the interface
        self.iface_atoms = None        # Interface atoms    :    Atoms in the group in the interface
        self.iface_sa = None           # Surface area       :    Surface area of the interface

        self.mol_names = None          # Molecule Names     :    Names of the molecules added to the group
        self.res_names = None          # Residue Names      :    Names of the residues added to the group
        self.atom_names = None         # Atom Names         :    Names of the atoms added to the group
        self.ndx_names = None          # Index Names        :    Names of the indices added to the group

        # Check if the group gets initialized with string versions of common identifiers
        if mols is not None:
            self.get_mol_atoms(mols)
        if residues is not None:
            self.get_res_atoms(residues)
        if my_atoms is not None:
            self.get_atom_atoms(my_atoms)
        if ndxs is not None:
            self.get_ndx_atoms(ndxs)
        if type(self.atoms[0]) is int:
            self.get_atom_atoms(self.atoms)

    # Get information method. Gathers the information for the group(s) selected
    def get_info(self):
        # Set the name of the group
        if self.name is None:
            self.set_name()
        # Reset the main information variables
        self.surfs, self.surf_ndxs, self.body_surfs, self.outer_body_atoms, self.surr_body_atoms = [], [], [], [], []
        self.body_vol, self.body_sa = 0, 0
        # Go through the atom in the group
        for atom in self.atoms:
            if atom.vol is None or atom.vol == 0:
                atom.calc_vol()
            # Add the volume of the atom to the group's volume and add the atom to the group
            self.body_vol += atom.vol
            # Check the surfaces of each of the atoms to see if they are on the outside or not
            for surf in atom.surfs:
                # Check if the surface has been set up
                if surf.points is None or len(surf.points) == 0 or surf.tris is None or len(surf.tris) == 0:
                    # Check to see if there is a file to choose from
                    if surf.file is not None:
                        # Try to grab the file
                        surf.read_file()
                    else:
                        surf.build()
                # Get the index of the surface
                my_surf_ndx = ndx_search(self.surf_ndxs, surf.ndx)
                # Check to see if we have found this surface before
                if my_surf_ndx < len(self.surf_ndxs) and self.surf_ndxs[my_surf_ndx] == surf.ndx:
                    continue
                else:
                    self.surfs.insert(my_surf_ndx, surf)
                    self.surf_ndxs.insert(my_surf_ndx, surf.ndx)
                # Check if the surface's first atom is in the group's list of atoms and the second atom is not
                if surf.atoms[0] in self.atoms and surf.atoms[1] not in self.atoms:
                    # Add the surface's body atom to the outer body atoms list if it isn't already in it
                    if surf.atoms[0] not in self.outer_body_atoms:
                        self.outer_body_atoms.append(surf.atoms[0])
                    # Add the surface's non-body atom to the outer body atoms list if it isn't already in it
                    if surf.atoms[1] not in self.surr_body_atoms:
                        self.surr_body_atoms.append(surf.atoms[1])
                    # Add the surface to the list of body surfaces
                    self.body_surfs.append(surf)
                    self.body_sa += surf.sa
                # Check if the surface's second atom is in the group's list of atoms and the first atom is not
                elif surf.atoms[0] not in self.atoms and surf.atoms[1] in self.atoms:
                    # Add the surface's body atom to the outer body atoms list if it isn't already in it
                    if surf.atoms[1] not in self.outer_body_atoms:
                        self.outer_body_atoms.append(surf.atoms[1])
                    # Add the surface to the list of body surfaces
                    if surf.atoms[0] not in self.surr_body_atoms:
                        self.surr_body_atoms.append(surf.atoms[0])
                    # Add the surface to the list of body surfaces
                    self.body_surfs.append(surf)
                    self.body_sa += surf.sa
                else:
                    continue
        # If the group has a bff get that information
        if self.bff is not None:
            # Get the surfaces and the
            self.iface_atoms, self.bff.iface_atoms, self.iface_surfs, self.bff.iface_surfs = [], [], [], []
            self.iface_sa = 0
            # Go through all the atoms in the group
            for atom in self.atoms:
                # Go through the surfaces in the atom's list of surfaces
                for surf in atom.surfs:
                    # Get the other atom from the surface's atoms
                    other_atom = [_ for _ in surf.atoms if _ != atom][0]
                    # Check to see if the other atom is in the self.bff's list of atoms
                    if other_atom in self.bff.atoms:
                        # Add the first atom to the group's list of interface atoms
                        self.iface_atoms.append(atom)
                        self.bff.iface_atoms.append(other_atom)
                        # Add the surface to the list of interface surfs and add the surface area of the surface
                        self.iface_surfs.append(surf)
                        self.bff.iface_surfs.append(surf)
                        self.iface_sa += surf.sa
            # Set the bff's surface area
            self.bff.iface_sa = self.iface_sa
        # Get the layers
        self.get_layers()

    def set_name(self):
        self.name = self.sys.name + "_" + "_".join(self.select_strs) + "group"

    # Add selection method. Adds a new selection to the group
    def add_sele(self, new_sele, new_str):
        if self.atoms is None:
            self.atoms = []
        self.selects.append(new_sele)
        self.select_strs.append(new_str)
        self.atoms += new_sele

    # Undo selection method. Undoes the last selection in the group
    def undo_sele(self):
        last_select = self.selects.pop()
        self.select_strs.pop()
        self.atoms = self.atoms[:len(self.atoms) - len(last_select)]
        self.get_info()

    # Get group function. Interprets the strings from below
    def get_mol_atoms(self, mols):
        # Get the molecules from their names
        for mol in mols:
            mol_atoms = self.sys.mols[self.sys.mol_names.index(int(mol))]
            for atom in mol_atoms:
                self.atoms.append(atom)

    # Get group function. Interprets the strings from below
    def get_res_atoms(self, residues):
        # Get the molecules from their names
        for res in residues:
            res_atoms = self.sys.residues[self.sys.res_names.index(int(res))]
            for atom in res_atoms:
                self.atoms.append(atom)

    # Get group function. Interprets the strings from below
    def get_atom_atoms(self, atoms):
        # Get the molecules from their names
        for atom in atoms:
            self.atoms.append(self.sys.atoms[int(atom)])

    # Get group function. Interprets the strings from below
    def get_ndx_atoms(self, ndxs):
        # Get the molecules from their names
        for ndx in ndxs:
            ndx_atoms = self.sys.ndxs[self.sys.ndx_names.index(ndx)]
            for atom in ndx_atoms:
                self.atoms.append(atom)

    def get_layers(self, max_layers=50):
        # Make sure that the group has atoms
        if self.atoms is None:
            return
        # Set up the layer surfs and layer atoms list variables
        counter, self.layer_atoms, layer_atoms_ndxs, self.layer_surfs = 0, [self.atoms, []], [[self.sys.atoms.index(_) for _ in self.atoms], []], [[]]
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
                    if surf.ndx[0] in layer_atoms_ndxs[-2] and surf.ndx[1] not in layer_atoms_ndxs[-2]:
                        self.layer_atoms[-1].append(self.sys.atoms[surf.ndx[1]])
                        layer_atoms_ndxs[-1].append(surf.ndx[1])

                    if surf.ndx[1] in layer_atoms_ndxs[-2] and surf.ndx[0] not in layer_atoms_ndxs[-2]:
                        self.layer_atoms[-1].append(self.sys.atoms[surf.ndx[0]])
                        layer_atoms_ndxs[-1].append(surf.ndx[0])

            # If there is nothing to add leave the layers loop
            if len(self.layer_surfs[-1]) == 0:
                self.layer_surfs.pop()
                break

            # Create the new layer lists
            self.layer_surfs.append([])
            self.layer_atoms.append([])
            layer_atoms_ndxs.append([])
            counter += 1

    def export_iface(self, g2=None, info_file=True, interface_atoms=True):
        # Check for a second group
        if g2 is None:
            if self.bff is not None:
                g2 = self.bff
            else:
                return
        # Export the interface
        export_iface(groups=[self, g2], info_file=info_file, interface_atoms=interface_atoms)



    def exports(self, atoms=True, shell=True, fill=True, surfaces=True, layers=True):
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
            write_surfs(surfs=self.layer_surfs[0], file_name=self.name + "_shell")
        # If the user wants a filled shell for the group
        if fill:
            write_surfs(surfs=self.surfs, file_name=self.name + "_fill")
        # If the user wants separate surfaces for the group
        if surfaces:
            os.mkdir(self.dir + "/surfaces")
            os.chdir(self.dir + "/surfaces")
            for surf in self.surfs:
                write_surfs([surf], file_name="_".join([str(_) for _ in surf.ndx]))
            os.chdir(self.dir)
        # If the user wants layers
        if layers:
            os.mkdir(os.getcwd() + "/layers")
            os.chdir(os.getcwd() + "/layers")
            for i in range(len(self.layer_surfs)):
                write_pdb(self.layer_atoms[i + 1], name=str(i) + "_atoms", sys=self.sys)
                write_surfs(self.layer_surfs[i], file_name=str(i) + "_surfs")
        os.chdir(self.sys.dir)

