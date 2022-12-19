from System.sys_funcs.calcs import *


class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, net, atoms=None, name=None, mols=None, residues=None, my_atoms=None, ndxs=None, bff=None):

        self.net = net                 # Network            :    Network of the System
        self.atoms = atoms             # Atoms              :    List of Atom type objects for the edge
        self.selects = []              # Previous Selection :    List of the previously selected atoms
        self.select_strs = []          # Previous String    :    String from the last group of atoms selected
        self.name = name               # Name               :    Name of the group

        self.surfs = None              # Surfaces           :    All surfaces associated with the group
        self.surf_ndxs = None          # Surface indices    :    Atom indices of the surfaces associated with the group
        self.body_surfs = None         # Body surfaces      :    The surfaces on the outside of the body
        self.body_sa = None            # Body surface area  :    The surface area of the outer surfaces of the body
        self.body_vol = None           # Body volume        :    The volume of the group's atom's cells
        self.outer_body_atoms = None   # Outer body atoms   :    Body atoms that create outer surfaces
        self.surr_body_atoms = None    # Surrounding atoms  :    Atoms not in the group that create surfaces with it

        self.bff = bff                # BFF                :    Other group used for comparison
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
            # Add the volume of the atom to the group's volume and add the atom to the group
            self.body_vol += atom.cell_vol
            # Check the surfaces of each of the atoms to see if they are on the outside or not
            for surf in atom.surfs:
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

    def set_name(self):
        self.name = self.net.sys.name + "_" + "_".join(self.select_strs) + "group"

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
            self.atoms += self.net.sys.mols[self.net.sys.mol_names.index(mol)]

    # Get group function. Interprets the strings from below
    def get_res_atoms(self, residues):
        # Get the molecules from their names
        for res in residues:
            self.atoms += self.net.sys.residues[self.net.sys.res_names.index(res)]

    # Get group function. Interprets the strings from below
    def get_atom_atoms(self, atoms):
        # Get the molecules from their names
        for atom in atoms:
            self.atoms.append(self.net.sys.atoms[int(atom)])

    # Get group function. Interprets the strings from below
    def get_ndx_atoms(self, ndxs):
        # Get the molecules from their names
        for ndx in ndxs:
            self.atoms += self.net.sys.ndxs[self.net.sys.ndx_names.index(ndx)]

    # Get Sol
    def find_sol_layers(self, mol=None):
        # We want a function that creates full solute layers
        # This variable holds the atoms for their respective layers
        layers_atoms = []
        # This variable holds the surfaces associated with the layers
        layers = []
        # Holds the current set of atoms being surrounded
        current_layer, next_layer = None, None
        sol_atoms = self.net.sys.sol.copy()
        if sol_atoms is None:
            return
        i = 0
        # Keep adding layers until the sol atoms are gone
        while len(sol_atoms) > 0 and i < 10:
            # Set the current layer to what was the next layer
            current_layer, next_layer = next_layer, []
            # If the current layer is none, this is the first layer and we are calculating the molecule's sol layer
            if current_layer is None:
                # If no molecule is provided go off of the system's molecule list
                if mol is None:
                    current_layer = []
                    for my_mol in self.net.sys.mols:
                        if self.net.sys.mol_names[self.net.sys.mols.index(my_mol)].lower() == 'sol':
                            continue
                        current_layer += my_mol
                else:
                    current_layer = mol
            # Set up the shell storage for the surfaces
            shell = []
            # Check each surface to see if it is in the current layer
            for surf in self.surfs:
                if surf.atoms[0] in current_layer and surf.atoms[1] in sol_atoms:
                    my_atom = sol_atoms.pop(sol_atoms.index(surf.atoms[1]))
                    next_layer.append(my_atom)
                    shell.append(surf)
            # Add the correct variables
            layers.append(shell)
            layers_atoms.append(current_layer)

            i += 1
        # Add the last layer of atoms
        layers.append(next_layer)
        return layers, layers_atoms
