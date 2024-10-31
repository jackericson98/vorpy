from System.sys_funcs.calcs.calcs import combine_inertia_tensors, calc_total_inertia_tensor
from System.sys_funcs.calcs.sorting import ndx_search
from System.sys_funcs.calcs.surf import calc_surf_sa
import numpy as np


def get_info(group):
    """
        Gathers information about the group and stores it in a dictionary
    """
    # Reset the group's data attributes
    group.sa, group.vol, group.vdw_vol, group.density, group.mass = 0, 0, 0, 0, 0
    com, vdw_com = [0, 0, 0], [0, 0, 0]
    # Get the balls in the group
    group_balls = group.net.balls.iloc[group.ball_ndxs].to_dict(orient='records')
    # Get the volume of the group
    for i, ball in enumerate(group_balls):
        # Check for the ball to be complets
        if not ball['complete']:
            continue
        # Add the volume to that of the group
        group.vol += ball['vol']
        group.vdw_vol += ball['vdw_vol']
        group.mass += ball['mass']
        # Add to the coms
        com = [com[j] + ball['com'][j] * ball['vol'] for j in range(3)]
        vdw_com = [vdw_com[j] + ball['loc'][j] * ball['mass'] for j in range(3)]
    if group.vol > 0:
        group.density = group.vdw_vol / group.vol
        group.com = np.array([com[j] / group.vol for j in range(3)])
        group.vdw_com = [vdw_com[j] / group.vdw_vol for j in range(3)]
    if 'moi' in group.net.balls.iloc[group.ball_ndxs[0]]:
        group.spatial_moment = combine_inertia_tensors([_['moi'] for _ in group_balls], [_['com'] for _ in group_balls],
                                                       group.com, [_['vol'] for _ in group_balls])
    if group.vdw_vol > 0:
        group.moi = calc_total_inertia_tensor(group_balls, group.vdw_com)
    # Check to see if the first layer has been calculated
    if group.layer_surfs is None or len(group.layer_surfs) == 0:
        group.get_layers(max_layers=1)
    if len(group.layer_surfs) > 0:
        for i in group.layer_surfs[0]:
            surf = group.net.surfs.iloc[i]
            # Check that the surface has a surface area
            if surf['sa'] is None or surf['sa'] == 0:
                surf_sa = calc_surf_sa(tris=surf['tris'], points=surf['points'])
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
