from System.sys_funcs.calcs.sorting import ndx_search


def get_layers(grp, max_layers=50, group_resids=True, build_surfs=True):
    """
    Gets the layers surrounding the group if they are calculated
    :param grp: Group to find layers from
    :param max_layers: Number of layers to go out
    :param group_resids: Whether to group together surrounding residues for layers, i.e. keep waters together
    :param build_surfs: To build the needed surfaces or not
    """

    net = grp.net
    # Make sure that the group has atoms
    if grp.atms is None:
        return
    # Set up the layer surfs and layer atoms list variables
    counter = 0
    grp.layer_atoms = [grp.atms[:], []]

    layer_atoms_ndxs = [grp.atms[:], []]
    grp.layer_surfs = [[]]
    grp.layer_verts = [[]]
    grp.layer_edges = [[]]
    grp.layer_info = [[0, 0]]
    # Set up the loop to keep adding layers
    while counter < max_layers:
        # Go through the atoms in the last layer
        for i in grp.layer_atoms[-2]:
            atom = net.balls.iloc[i]
            # Go through the surfaces in the atom's list of surfaces
            for j in atom['asurfs']:
                surf = grp.net.surfs.iloc[j]
                if j in grp.layer_surfs[-1] or (len(grp.layer_surfs) >= 2 and j in grp.layer_surfs[-2]):
                    continue
                elif surf['satoms'][0] in layer_atoms_ndxs[-2] and surf['satoms'][1] in layer_atoms_ndxs[-2]:
                    continue
                grp.layer_surfs[-1].append(j)
                # Add the vertices
                for k in surf['sverts']:
                    if k not in grp.layer_verts[-1]:
                        grp.layer_verts[-1].append(k)
                for edge in surf['sedges']:
                    if edge not in grp.layer_edges[-1]:
                        grp.layer_edges[-1].append(edge)
                # Get the index of the surface
                surf_ndx = ndx_search(grp.surf_ndxs, surf['satoms'])
                # Check if the surface has been added yet or not
                if surf_ndx < len(grp.surf_ndxs) and grp.surf_ndxs[surf_ndx] != surf['satoms']:
                    # Insert the index and the surfaces in their 181L place
                    grp.surfs.insert(surf_ndx, j)
                    grp.surf_ndxs.insert(surf_ndx, surf['satoms'])
                # Sort the surface's atoms inside or out
                if surf['satoms'][0] in layer_atoms_ndxs[-2] and surf['satoms'][1] not in layer_atoms_ndxs[-2]:
                    grp.layer_atoms[-1].append(surf['satoms'][1])
                    layer_atoms_ndxs[-1].append(surf['satoms'][1])
                if surf['satoms'][1] in layer_atoms_ndxs[-2] and surf['satoms'][0] not in layer_atoms_ndxs[-2]:
                    grp.layer_atoms[-1].append(surf['satoms'][0])
                    layer_atoms_ndxs[-1].append(surf['satoms'][0])
        if build_surfs and grp.sys.cmnds['xpt'] != [['logs']]:
            # Check to make sure the surfaces are built in the layer
            grp.build_surfs()
        # Check to see if the residues are supposed to stay together
        if group_resids:
            for my_atom in grp.layer_atoms[-1]:
                atom = net.balls.iloc[my_atom]
                if atom['res'] is not None:
                    # Get the atoms in the residue that are not already in the layer
                    for resid_atom in atom['res'].atoms:
                        # Check if the atom is in the layer or not
                        if resid_atom not in grp.layer_atoms[-1]:
                            grp.layer_atoms[-1].append(resid_atom)
        # Get the surface area and volume for the layer
        for my_atom in grp.layer_atoms[-1]:
            atom = net.balls.iloc[my_atom]
            # Add the volume to the current layer's volume
            grp.layer_info[-1][0] += atom['vol']
        # Get the surface area of the layer
        for surf in grp.layer_surfs[-1]:
            # Add the surface area to the current layer's surface area
            grp.layer_info[-1][1] += net.surfs['sa'][surf]
        # If there is nothing to add leave the layers loop
        if len(grp.layer_surfs[-1]) == 0:
            grp.layer_surfs.pop()
            break
        # Create the new layer lists
        grp.layer_surfs.append([])
        grp.layer_atoms.append([])
        grp.layer_edges.append([])
        grp.layer_verts.append([])
        grp.layer_info.append([0, 0])
        layer_atoms_ndxs.append([])
        counter += 1
