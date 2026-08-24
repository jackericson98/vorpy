import numpy as np
from vorpy.src.calculations import combine_inertia_tensors
from vorpy.src.calculations import calc_total_inertia_tensor
from vorpy.src.calculations import ndx_search
from vorpy.src.calculations import calc_surf_sa


def get_info(group):
    """
    Gathers and calculates comprehensive information about a molecular group, including geometric, physical, and structural properties.

    This function performs a detailed analysis of a molecular group, calculating:
    - Surface area and volume metrics
    - Center of mass (both geometric and van der Waals)
    - Mass and density properties
    - Moment of inertia tensors
    - Layer-based surface properties

    Parameters:
        group (Group): The Group object to analyze. Must have a valid network (group.net) and ball indices (group.ball_ndxs).

    Returns:
        None. Results are stored in the Group object's attributes:
        - sa (float): Total surface area of the group
        - vol (float): Total volume of the group
        - vdw_vol (float): Van der Waals volume
        - density (float): Ratio of van der Waals volume to total volume
        - mass (float): Total mass of atoms
        - com (numpy.ndarray): Center of mass coordinates
        - vdw_com (list): Van der Waals center of mass
        - spatial_moment (list): Spatial moment tensor
        - moi (list): Moment of inertia tensor

    Examples:
        >>> from vorpy.src.group import Group
        >>> # Create and analyze a group
        >>> group = Group(sys=my_system, name='protein_A')
        >>> group.add_balls(atom_indices=range(100))
        >>> group.build()
        >>> # Calculate group properties
        >>> get_info(group)
        >>> # Access calculated properties
        >>> print(f"Surface area: {group.sa:.2f} Å²")
        >>> print(f"Volume: {group.vol:.2f} Å³")
        >>> print(f"Center of mass: {group.com}")
        >>> print(f"Moment of inertia: {group.moi}")
    """
    # Reset the group's data attributes
    group.sa, group.vol, group.vdw_vol, group.density, group.mass = 0, 0, 0, 0, 0
    # Reset the group's curvature geometry
    group.avg_mean_curv = 0.0
    group.avg_gauss_curv = 0.0

    group.int_mean_curv = 0.0
    group.int_mean_curv_sq = 0.0
    group.int_gauss_curv = 0.0
    # Center of masses
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
        # Add the vdw volume to that of the group
        group.vdw_vol += ball['vdw_vol']
        # Add the mass to that of the group
        group.mass += ball['mass']
        # Add to the coms
        com = [com[j] + ball['com'][j] * ball['vol'] for j in range(3)]
        vdw_com = [vdw_com[j] + ball['loc'][j] * ball['mass'] for j in range(3)]
    # Check to see if the volume is greater than 0
    if group.vol > 0:
        # Calculate the density
        group.density = group.vdw_vol / group.vol
        # Calculate the center of mass
        group.com = np.array([com[j] / group.vol for j in range(3)])
        # Calculate the vdw center of mass
        group.vdw_com = [vdw_com[j] / group.vdw_vol for j in range(3)]
    # Check to see if the moi has been calculated
    if 'moi' in group.net.balls.iloc[group.ball_ndxs[0]]:
        # Calculate the spatial moment
        group.spatial_moment = combine_inertia_tensors([_['moi'] for _ in group_balls], [_['com'] for _ in group_balls],
                                                       group.com, [_['vol'] for _ in group_balls])
    if group.vdw_vol > 0:
        group.moi = calc_total_inertia_tensor(group_balls, group.vdw_com)
    # Check to see if the first layer has been calculated
    if group.layer_surfs is None or len(group.layer_surfs) == 0:
        group.get_layers(max_layers=1)

    # Check to see if there are any layers
    if group.layer_surfs is not None and len(group.layer_surfs) > 0:

        # The first surface layer defines the exposed group boundary.
        for surf_ndx in group.layer_surfs[0]:

            # Get the surface
            surf = group.net.surfs.iloc[surf_ndx]

            # --------------------------------------------------------------
            # Surface area
            # --------------------------------------------------------------

            if surf['sa'] is None or surf['sa'] == 0:
                surf_sa = calc_surf_sa(
                    tris=surf['tris'],
                    points=surf['points']
                )
            else:
                surf_sa = float(surf['sa'])

            group.sa += surf_sa

            # --------------------------------------------------------------
            # Integrated curvature
            # --------------------------------------------------------------

            int_mean_curv = surf.get('int_mean_curv', 0.0)
            int_mean_curv_sq = surf.get('int_mean_curv_sq', 0.0)
            int_gauss_curv = surf.get('int_gauss_curv', 0.0)

            # Protect against missing/None values
            if int_mean_curv is not None:
                group.int_mean_curv += float(int_mean_curv)

            if int_mean_curv_sq is not None:
                group.int_mean_curv_sq += float(int_mean_curv_sq)

            if int_gauss_curv is not None:
                group.int_gauss_curv += float(int_gauss_curv)

    # ----------------------------------------------------------------------
    # Area-weighted group curvature
    # ----------------------------------------------------------------------

    if group.sa > 0.0:
        group.avg_mean_curv = (
                group.int_mean_curv / group.sa
        )

        group.avg_gauss_curv = (
                group.int_gauss_curv / group.sa
        )


def add_balls(grp, ball_list):
    """
    Adds atoms to a group while maintaining sorted order and preventing duplicates.

    This function efficiently integrates new atoms into a group's existing atom list. It uses binary search
    to maintain a sorted order of atom indices and ensures no duplicate atoms are added. The function
    handles various input types including lists of atoms from molecules, residues, or direct atom selections.

    Parameters:
        grp (Group): The Group object to which atoms will be added
        ball_list (list): List of atom indices to be added to the group. Can be from various sources
                         (e.g., molecule.atoms, residue.atoms, or direct atom selections)

    Returns:
        None. The group's ball_ndxs and atms attributes are updated with the new atoms.
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
