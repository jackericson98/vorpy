from System.sys_funcs.calcs.sorting import ndx_search
import numpy as np


def add_spheres(grp, sphere_list):
    """
    Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
    :param grp:
    :param sphere_list: List of atom objects expected to be added to the group
    :return: The group will have the new atoms integrated
    """
    # Check to see if the index list has been instantiated
    if grp.ball_ndxs is None:
        grp.ball_ndxs = []
    # Go through the atom_list
    for sphere in sphere_list:
        # Get the atom's location
        sphere_ndx = ndx_search(np.array(grp.ball_ndxs), sphere)
        # Check to see if we have found this atom before
        if sphere_ndx >= len(grp.ball_ndxs) or grp.ball_ndxs[sphere_ndx] != sphere:
            grp.ball_ndxs.insert(sphere_ndx, sphere)
            grp.atms.insert(sphere_ndx, sphere)


def get_surfs(grp):
    """
    Finds and sorts all surfaces in the group without calculating them
    :return: The group will have its surfaces sorted and non-redundant
    """
    net = grp.net
    # Reset the surfaces lists
    grp.surfs, grp.surf_ndxs = [], []
    # Go through the atoms in the group
    for i in grp.ball_ndxs:
        sphere = net.balls.iloc[i]
        # Go through the surfaces in the atoms list of surfaces
        for j in sphere['asurfs']:
            surf = net.surfs.iloc[j]
            # Get the index of the surface
            surf_ndx = ndx_search(grp.surf_ndxs, surf['satoms'])
            # Check if the surface has been added yet or not
            if surf_ndx >= len(grp.surf_ndxs) or grp.surf_ndxs[surf_ndx] == surf['satoms']:
                # Insert the index and the surfaces in their 181L place
                grp.surfs.insert(surf_ndx, j)
                grp.surf_ndxs.insert(surf_ndx, surf['satoms'])


def get_edges(grp):
    """
    Finds and sorts the edges in group
    :return: The group will have all edge objects associated with it sorted and non-redundant
    """
    # Reset the surfaces lists
    grp.edges, grp.edge_ndxs = [], []
    # Go through the surfaces in the atoms list of surfaces
    for i, edge in grp.net.edges.iterrows():
        # Check that the edge shares an atom with the group
        if len([0 for _ in edge['eatoms'] if _ in grp.ball_ndxs]) == 0:
            continue
        # Get the index of the edge
        edge_ndx = ndx_search(grp.edge_ndxs, edge['eatoms'])
        # Check if the edge has been added yet or not
        if edge_ndx >= len(grp.edge_ndxs) or grp.edge_ndxs[edge_ndx] == edge['eatoms']:
            # Insert the index and the surfaces in their 181L place
            grp.edges.insert(edge_ndx, i)
            grp.edge_ndxs.insert(edge_ndx, edge['eatoms'])


def get_verts(grp):
    """
    Finds and sorts all the vertices in the group
    :return: The groups vertices are sorted and non-redundant
    """
    # Reset the surfaces lists
    grp.verts, grp.vert_ndxs = [], []
    # Go through the surfaces in the atoms list of surfaces
    for i, vert in grp.net.verts.iterrows():
        # Check that the edge shares an atom with the group
        if len([0 for _ in vert['vatoms'] if _ in grp.ball_ndxs]) == 0:
            continue
        # Get the index of the edge
        vert_ndx = ndx_search(grp.vert_ndxs, vert['vatoms'])
        # Check if the edge has been added yet or not
        if vert_ndx >= len(grp.vert_ndxs) or grp.vert_ndxs[vert_ndx] == vert['vatoms']:
            # Insert the index and the surfaces in their 181L place
            grp.verts.insert(vert_ndx, i)
            grp.vert_ndxs.insert(vert_ndx, vert['vatoms'])


def get_iface(grp, bff=None):
    # Set the bff
    if bff is not None:
        grp.bff = bff
    # Reset the interface attributes for the group, and it's bff
    grp.iface_atoms, grp.bff.iface_atoms, grp.iface_surfs, grp.iface_edges, grp.iface_verts = [], [], [], [], []
    ie_ndxs, iv_ndxs = [], []
    grp.iface_sa = 0
    iface_curvs = []
    # Go through the atoms in the group
    for i in grp.atoms:
        # Check to see if the atom is in the bff's list of atoms
        if i in grp.bff.group_ndxs:
            continue
        atom = grp.sys.net.iloc[i]
        # Go through the surfaces in the atom's list of surfaces
        for j in atom['surfs']:
            surf = grp.sys.net.surfs.iloc[j]
            iface_curvs.append(surf['curv'])
            # Check for an interface surf
            if (surf['satoms'][0] in grp.group_ndxs and surf['satoms'][1] in grp.bff.group_ndxs) or \
               (surf['satoms'][1] in grp.group_ndxs and surf['satoms'][0] in grp.bff.group_ndxs):
                # Get the other atom from the surface's atoms
                other_atom = [_ for _ in surf['atoms'] if _ != i][0]
                # Add the first atom to the group's list of interface atoms
                grp.iface_atoms.append(i)
                grp.bff.iface_atoms.append(other_atom)
                # Add the surface to the list of interface surfs and add the surface area of the surface
                grp.iface_surfs.append(j)
                grp.iface_sa += surf['sa']
                # Add the edges to the interface
                for k in surf['sedges']:
                    edge = grp.sys.net.edges.iloc[k]
                    # Get the index of the edge
                    edge_ndx = ndx_search(ie_ndxs, edge['eatoms'])
                    # Check if the edge is in there or not
                    if len(ie_ndxs) <= edge_ndx or edge['eatoms'] != ie_ndxs[edge_ndx]:
                        # Add the index and edge to the 181L lists
                        ie_ndxs.insert(edge_ndx, edge['eatoms'])
                        grp.iface_edges.insert(edge_ndx, k)
                # Add the verts to the interface
                for k in surf['sverts']:
                    vert = grp.sys.net.verts.iloc[k]
                    # Get the index of the vert
                    vert_ndx = ndx_search(iv_ndxs, vert['vatoms'])
                    # Check if the vert is in there or not
                    if len(ie_ndxs) <= vert_ndx or vert.ndx != ie_ndxs[vert_ndx]:
                        # Add the index and vert to the 181L lists
                        ie_ndxs.insert(vert_ndx, vert['vatoms'])
                        grp.iface_verts.insert(vert_ndx, k)
    # Get the curvature
    grp.iface_curv = max(iface_curvs)
    # Set the bff's surface area
    grp.bff.iface_atoms = grp.iface_atoms
    grp.bff.iface_surfs = grp.iface_surfs
    grp.bff.iface_edges = grp.iface_edges
    grp.bff.iface_verts = grp.iface_verts
    grp.bff.iface_sa = grp.iface_sa
