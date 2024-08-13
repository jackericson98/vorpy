from System.sys_funcs.calcs.sorting import ndx_search
from System.sys_funcs.calcs.surf import calc_surf_sa
import numpy as np


def get_info(group):
    """
        Gathers information about the group and stores it in a dictionary
    """
    net = group.net
    # Get the group objects
    group.get_surfs()
    group.get_edges()
    group.get_verts()
    # Reset the group's data attributes
    group.sa, group.vol, group.density = 0, 0, 0
    tot_atom_vol = 0
    # Get the volume of the group
    for i in group.ball_ndxs:
        atom = group.net.balls.iloc[i]
        if not atom['complete']:
            continue
        # Add the volume to that of the group
        group.vol += atom['vol']
        tot_atom_vol += (4 / 3) * np.pi * atom['rad'] ** 3
    if group.vol > 0:
        group.density = tot_atom_vol / group.vol
    # Check to see if the first layer has been calculated
    if group.layer_surfs is None or len(group.layer_surfs) == 0:
        group.get_layers(max_layers=1)
    if len(group.layer_surfs) > 0:
        for i in group.layer_surfs[0]:
            surf = group.net.surfs.iloc[i]
            # Check that the surface has a surface area
            if surf['sa'] is None or surf['sa'] == 0:
                # Get the surface area for the surface
                edge_ndxss = [ndx_search(net.edge_ndxs, _) for _ in surf['sedges']]
                edges = np.array([net.edges.iloc[_] for _ in edge_ndxss])
                surf_sa = calc_surf_sa(edges=edges, com=np.array(surf['com']), tris=surf['tris'], points=surf['points'],
                                       flat=surf['flat'])
            else:
                surf_sa = surf['sa']
            # Add the surface area
            group.sa += surf_sa


def add_balls(grp, ball_list):
    """
    Adds the atoms from a list (mol.atoms, res.atoms, atoms, etc) to the group checking duplicates
    :param grp:
    :param ball_list: List of atom objects expected to be added to the group
    :return: The group will have the new atoms integrated
    """
    # Check to see if the index list has been instantiated
    if grp.ball_ndxs is None:
        grp.ball_ndxs = []
    # Go through the atom_list
    for sphere in ball_list:
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
        for j in sphere['surfs']:
            surf = net.surfs.iloc[j]
            # Get the index of the surface
            surf_ndx = ndx_search(grp.surf_ndxs, surf['balls'])
            # Check if the surface has been added yet or not
            if surf_ndx >= len(grp.surf_ndxs) or grp.surf_ndxs[surf_ndx] == surf['balls']:
                # Insert the index and the surfaces in their 181L place
                grp.surfs.insert(surf_ndx, j)
                grp.surf_ndxs.insert(surf_ndx, surf['balls'])


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
        if len([0 for _ in edge['balls'] if _ in grp.ball_ndxs]) == 0:
            continue
        # Get the index of the edge
        edge_ndx = ndx_search(grp.edge_ndxs, edge['balls'])
        # Check if the edge has been added yet or not
        if edge_ndx >= len(grp.edge_ndxs) or grp.edge_ndxs[edge_ndx] == edge['balls']:
            # Insert the index and the surfaces in their 181L place
            grp.edges.insert(edge_ndx, i)
            grp.edge_ndxs.insert(edge_ndx, edge['balls'])


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
        if len([0 for _ in vert['balls'] if _ in grp.ball_ndxs]) == 0:
            continue
        # Get the index of the edge
        vert_ndx = ndx_search(grp.vert_ndxs, vert['balls'])
        # Check if the edge has been added yet or not
        if vert_ndx >= len(grp.vert_ndxs) or grp.vert_ndxs[vert_ndx] == vert['balls']:
            # Insert the index and the surfaces in their 181L place
            grp.verts.insert(vert_ndx, i)
            grp.vert_ndxs.insert(vert_ndx, vert['balls'])
