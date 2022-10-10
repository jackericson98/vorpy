

class Group:
    """Group class. Used to hold selections of atoms and do analysis on it"""
    def __init__(self, net, atoms):

        self.net = net                 # Network            :    Network of the System
        self.atoms = atoms             # Atoms              :    List of Atom type objects for the edge
        self.prev_sele = None          # Previous Selection :    List of the previously selected atoms
        self.name = None               # Name               :    Name of the group

        self.body_surfs = None         # Body surfaces      :    The surfaces on the outside of the body
        self.body_sa = None            # Body surface area  :    The surface area of the outer surfaces of the body
        self.body_vol = None           # Body volume        :    The volume of the group's atom's cells
        self.all_body_atoms = None     # All body atoms     :    All atoms contained in the group's "body"
        self.outer_body_atoms = None   # Outer body atoms   :    Body atoms that create outer surfaces
        self.surr_body_atoms = None    # Surrounding atoms  :    Atoms not in the group that create surfaces with it

        self.bff = None                # BFF                :    Other group used for comparison
        self.interface_surfs = None    # Interface surfaces :    Surfaces that make the interface
        self.interface_atoms = None    # Interface atoms    :    Atoms in the group in the interface
        self.interface_sa = None       # Surface area       :    Surface area of the interface


    # Get information method. Gathers the information for the group(s) selected
    def get_info(self):
        # Reset the main information variables
        self.body_surfs, self.all_body_atoms, self.outer_body_atoms, self.surr_body_atoms = [], [], [], []
        self.body_vol, self.body_sa = 0, 0
        # Go through the atom in the group
        for atom in self.atoms:
            # Add the volume of the atom to the group's volume and add the atom to the group
            self.body_vol += atom.cell_vol
            self.all_body_atoms.append(atom)
            # Check the surfaces of each of the atoms to see if they are on the outside or not
            for surf in atom.surfs:
                # Check if the surface's first atom is in self.atoms and the surfaces second atom is in self.bff.atoms
                if surf.atoms[0] in self.atoms and surf.atoms[1] not in self.atoms:
                    # Add a0 and a1 to the list of the outer atoms and the list of surrounding atoms respectively
                    self.outer_body_atoms.append(surf.atoms[0])
                    self.surr_body_atoms.append(surf.atoms[1])
                # Check if the surface's first atom is in self.bff.atoms and the surfaces second atom is in self.atoms
                elif surf.atoms[0] not in self.atoms and surf.atoms[1] in self.atoms:
                    # Add a1 and a0 to the list of the outer atoms and the list of surrounding atoms respectively
                    self.outer_body_atoms.append(surf.atoms[1])
                    self.surr_body_atoms.append(surf.atoms[0])
                else:
                    continue
                # Add the surface to the list of interface surfs and add the surface area of the surface
                self.body_surfs.append(surf)
                self.body_sa += surf.sa
        # If the group has a bff get that information
        if self.bff:
            # Get the surfaces and the
            interface = []
            self.interface_sa = 0
            # Check the network's surfaces' atoms looking for one in self.atoms and the other in self.bff.atoms
            for surf in self.net.surfs:
                # Check if the surface's first atom is in self.atoms and the surfaces second atom is in self.bff.atoms
                if surf.atoms[0] in self.atoms and surf.atoms[1] in self.bff.atoms:
                    # Add the first atom to the group's list of interface atoms
                    self.interface_atoms.append(surf.atoms[0])
                # Check if the surface's first atom is in self.bff.atoms and the surfaces second atom is in self.atoms
                elif surf.atoms[0] in self.bff.atoms and surf.atoms[1] in self.atoms:
                    # Add the second atom to the list of interface atoms
                    self.interface_atoms.append(surf.atoms[1])
                # If the surface is not in the interface, continue
                else:
                    continue
                # Add the surface to the list of interface surfs and add the surface area of the surface
                interface.append(surf)
                self.interface_sa += surf.sa
