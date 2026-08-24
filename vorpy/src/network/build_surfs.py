import time
from vorpy.src.calculations import calc_surf_sa
from vorpy.src.calculations import calc_tetra_vol
from vorpy.src.network.build_surf import build_surf


def build_surfs(net, store_points=True):
    """
    Builds surfaces for all edges in a network and calculates their geometric properties.

    This function processes each surface in the network, constructing the surface mesh and computing
    various geometric properties including:
    - Surface area
    - Volume contributions to adjacent atoms
    - Mean and Gaussian curvatures
    - Surface points and triangulation
    - Center of mass locations

    Parameters
    ----------
    net : Network
        The network object containing surfaces to be built
    store_points : bool, optional
        If True, stores detailed surface point and triangle data. If False, only stores
        essential metrics. Default is True.

    Returns
    -------
    None
        The function modifies the network object in place, adding surface properties to
        the network's surfaces DataFrame.

    Notes
    -----
    - Progress is displayed during surface construction
    - Invalid surfaces are removed from the network
    - Surface properties are stored in lists and later added to the network DataFrame
    - For large networks, setting store_points=False can reduce memory usage
    """
    # Instantiate the lists for storage
    points, tris, mean_tri_curvs, mean_curvs, avg_mean_curvs, gauss_tri_curvs, gauss_curvs, avg_gauss_curvs, funcs, coms, flats, sas, vols, surf_locs = [], [], [], [], [], [], [], [], [], [], [], [], [], []
    # Make each surface
    for i, surf in net.surfs.iterrows():
        # Build the surfaces and print the progress
        percentage = min(100.0, 100.0 * (i + 1) / len(net.surfs))
        net.update_progress("Building surfaces", percentage)
        # Get the radii, locations, and numbers of the balls involved in the surface
        rads = [net.balls['rad'][_] for _ in surf['balls']]
        locs = [net.balls['loc'][_] for _ in surf['balls']]
        nums = [net.balls['num'][_] for _ in surf['balls']]
        # If the first ball has a larger radius than the second, swap them
        if rads[0] > rads[1]:
            rads, locs, nums = [rads[1], rads[0]], [locs[1], locs[0]], [nums[1], nums[0]]
        # Build the surface
        my_surf = build_surf(locs=locs, rads=rads, epnts=[net.edges['points'][_] for _ in surf['edges']],
                             res=net.settings['surf_res'], net_type=net.settings['net_type'])
        # If the surface is not valid, drop it from the network
        if my_surf is None:
            net.surfs.drop(index=i, inplace=True)
            continue
        (surf_points, surf_tris, mean_surf_tri_curvs, mean_surf_curv, avg_mean_surf_curv, gauss_surf_tri_curvs,
         gauss_surf_curv, avg_gauss_surf_curv, surf_func, surf_com, surf_flat, surf_loc) = my_surf
        # full_count = {_: full_count[_] + timer[_] for _ in full_count}
        # Get the surface Volumes
        sv0 = sum([calc_tetra_vol(locs[0], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        sv1 = sum([calc_tetra_vol(locs[1], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        # Calculate the surface area of the surface
        sa = calc_surf_sa(tris=surf_tris, points=surf_points)
        # If we are doing a large export and will need the points later in the process for export and such
        if store_points:
            # Append the surface points and triangles to the lists
            points.append(surf_points)
            tris.append(surf_tris)
            mean_tri_curvs.append(mean_surf_tri_curvs)
            gauss_tri_curvs.append(gauss_surf_tri_curvs)
        else:
            # Append empty lists to the lists
            points.append([])
            tris.append([])
            mean_tri_curvs.append([])
            gauss_tri_curvs.append([])
        # Append the mean and Gaussian curvatures to the lists
        mean_curvs.append(mean_surf_curv)
        avg_mean_curvs.append(avg_mean_surf_curv)
        gauss_curvs.append(gauss_surf_curv)
        avg_gauss_curvs.append(avg_gauss_surf_curv)
        # Append the surface function, center of mass, flatness, surface area, volumes, and location to the lists
        funcs.append(surf_func)
        coms.append(surf_com)
        flats.append(surf_flat)
        sas.append(sa)
        vols.append({nums[0]: sv0, nums[1]: sv1})
        surf_locs.append(surf_loc)
    # Set the dataframe elements
    (net.surfs['points'], net.surfs['tris'], net.surfs['mean_tri_curvs'], net.surfs['mean_curv'],
     net.surfs['avg_mean_curv'], net.surfs['gauss_tri_curvs'], net.surfs['gauss_curv'], net.surfs['avg_gauss_curv'],
     net.surfs['func'], net.surfs['com'], net.surfs['flat'], net.surfs['sa'], net.surfs['vols'], net.surfs['loc']) = \
        (points, tris, mean_tri_curvs, mean_curvs, avg_mean_curvs, gauss_tri_curvs, gauss_curvs, avg_gauss_curvs, funcs,
         coms, flats, sas, vols, surf_locs)

    # Get the curvature in the 95th percentile
    my_surf_curvs = net.surfs['mean_curv'].to_list()
    if net.settings['surf_scheme'] == 'gauss':
        my_surf_curvs = net.surfs['gauss_curv'].to_list()
    # Sort the curvatures
    my_surf_curvs.sort()
    # Get the 95th percentile curvature
    try:
        net.max_curv = my_surf_curvs[min(int(0.99 * len(my_surf_curvs)), len(my_surf_curvs) - 1)]
    except IndexError:
        net.max_curv = 0
    net.update_progress("Building surfaces", 100.0)
    net.metrics['surf'] = time.perf_counter() - net.metrics['start'] - net.metrics['vert'] - net.metrics['con']
