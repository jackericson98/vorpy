from System.sys_funcs.calcs import ndx_search


def add_atoms(grp, atom_list):
    """
    Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
    :param grp:
    :param atom_list: List of atom objects expected to be added to the group
    :return: The group will have the new atoms integrated
    """
    # Check to see if the atoms list has been instantiated
    if grp.atoms is None:
        grp.atoms = []
    # Go through the atom_list
    for atom in atom_list:
        # Get the atom's location
        atom_ndx = ndx_search(grp.atom_ndxs, atom.num)
        # Check to see if we have found this atom before
        if atom_ndx >= len(grp.atom_ndxs) or grp.atoms[atom_ndx].num != atom.num:
            grp.atoms.insert(atom_ndx, atom)
            grp.atom_ndxs.insert(atom_ndx, atom.num)


def get_surfs(grp):
    """
    Finds and sorts all surfaces in the group without calculating them
    :return: The group will have its surfaces sorted and non-redundant
    """
    # Reset the surfaces lists
    grp.surfs, grp.surf_ndxs = [], []
    # Go through the atoms in the group
    for atom in grp.atoms:
        # Go through the surfaces in the atoms list of surfaces
        for surf in atom.surfs:
            # Get the index of the surface
            surf_ndx = ndx_search(grp.surf_ndxs, surf.ndx)
            # Check if the surface has been added yet or not
            if surf_ndx >= len(grp.surf_ndxs) or grp.surf_ndxs[surf_ndx] == surf.ndx:
                # Insert the index and the surfaces in their correct place
                grp.surfs.insert(surf_ndx, surf)
                grp.surf_ndxs.insert(surf_ndx, surf.ndx)


def get_edges(grp):
    """
    Finds and sorts the edges in group
    :return: The group will have all edge objects associated with it sirted and non-redundant
    """
    # Reset the surfaces lists
    grp.edges, grp.edge_ndxs = [], []
    # Go through the surfaces in the atoms list of surfaces
    for edge in grp.sys.net.edges:
        # Get the index of the edge
        edge_ndx = ndx_search(grp.edge_ndxs, edge.ndx)
        # Check if the edge has been added yet or not
        if edge_ndx >= len(grp.edge_ndxs) or grp.edge_ndxs[edge_ndx] == edge.ndx:
            # Insert the index and the surfaces in their correct place
            grp.edges.insert(edge_ndx, edge)
            grp.edge_ndxs.insert(edge_ndx, edge.ndx)


def get_verts(grp):
    """
    Finds and sorts all the vertices in the group
    :return: The groups vertices are sorted and non-redundant
    """
    # Reset the surfaces lists
    grp.verts, grp.vert_ndxs = [], []
    grp.atom_ndxs = [_.num for _ in grp.atoms]
    grp.atom_ndxs.sort()
    # Go through the surfaces in the atoms list of surfaces
    for vert in grp.sys.net.verts:
        # Get the index of the edge
        vert_ndx = ndx_search(grp.vert_ndxs, vert.ndx)
        # Check if the edge has been added yet or not
        if vert_ndx >= len(grp.vert_ndxs) or grp.vert_ndxs[vert_ndx] == vert.ndx:
            # Insert the index and the surfaces in their correct place
            grp.verts.insert(vert_ndx, vert)
            grp.vert_ndxs.insert(vert_ndx, vert.ndx)


def get_iface(grp, bff=None):
    # Set the bff
    if bff is not None:
        grp.bff = bff
    # Reset the interface attributes for the group, and it's bff
    grp.iface_atoms, grp.bff.iface_atoms, grp.iface_surfs, grp.bff.iface_surfs = [], [], [], []
    grp.iface_sa = 0
    # Go through the atoms in the group
    for atom in grp.atoms:
        # Check to see if the atom is in the bff's list of atoms
        if atom.num in grp.bff.atom_ndxs:
            continue
        # Go through the surfaces in the atom's list of surfaces
        for surf in atom.surfs:
            # Check for an interface surf
            if (surf.ndx[0] in grp.atom_ndxs and surf.ndx[1] in grp.bff.atom_ndxs) or \
               (surf.ndx[1] in grp.atom_ndxs and surf.ndx[0] in grp.bff.atom_ndxs):
                # Get the other atom from the surface's atoms
                other_atom = [_ for _ in surf.atoms if _ != atom][0]
                # Add the first atom to the group's list of interface atoms
                grp.iface_atoms.append(atom)
                grp.bff.iface_atoms.append(other_atom)
                # Add the surface to the list of interface surfs and add the surface area of the surface
                grp.iface_surfs.append(surf)
                grp.bff.iface_surfs.append(surf)
                grp.iface_sa += surf.sa
    # Set the bff's surface area
    grp.bff.iface_sa = grp.iface_sa
