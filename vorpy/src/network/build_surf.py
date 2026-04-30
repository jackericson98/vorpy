import numpy as np
from vorpy.src.calculations import calc_surf_func
from vorpy.src.calculations import calc_surf_tri_curvs
from vorpy.src.network.perimeter import build_perimeter
from vorpy.src.network.fill import calc_surf_point
from vorpy.src.network.fill import calc_surf_point_from_plane
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import project_to_plane
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import unproject_to_3d
from vorpy.src.network.triangulate import triangulate_2D_Surface
from vorpy.src.network.triangulate import is_within


def get_com(locs, rads, perimeter, surf_loc, surf_norm, func, flat, net_type='aw'):
    """
    Calculates the center of mass for a surface's perimeter points, with special handling for different network types.
    
    Parameters
    ----------
    locs : numpy.ndarray
        Array of ball locations
    rads : numpy.ndarray
        Array of ball radii
    perimeter : list
        List of perimeter points defining the surface boundary
    surf_loc : numpy.ndarray
        Location of the surface center
    surf_norm : numpy.ndarray
        Normal vector of the surface plane
    func : callable
        Function used to calculate surface points
    flat : bool
        Flag indicating if the surface is flat
    net_type : str, optional
        Type of network ('aw', 'del', or 'pow'). Default is 'aw'
    
    Returns
    -------
    tuple
        A tuple containing:
        - numpy.ndarray: The calculated center of mass point
        - bool: Flag indicating if the point was found through perimeter point iteration
    
    Notes
    -----
    - For 'del' and 'pow' network types, uses simple center of mass calculation
    - For 'aw' type, attempts multiple strategies to find a valid center point:
        1. Uses surface location if within perimeter
        2. Projects true center of mass onto surface
        3. Uses center of mass of sampled perimeter points
        4. Falls back to closest perimeter point to true center
    """
    # If the surface is flat just get the center of mass
    if net_type in {'del', 'pow'}:
        return calc_com(points=np.array(perimeter)), False
    # Next create the polygon so that we can tell if the center of mass is within the perimeter
    if is_within(perimeter, surf_loc, surf_loc, surf_norm):
        return surf_loc, False
    # First try the center of mass of the 3d points projected onto the surface
    true_com = calc_com(points=np.array(perimeter))
    my_com = calc_surf_point(locs, point=true_com, func=func)
    if my_com is not None:
        if is_within(perimeter, my_com, surf_loc, surf_norm):
            return my_com, False
    # Next try to calculate a center of mass of some of the points
    my_com = calc_surf_point(locs, point=calc_com(points=np.array(perimeter[::5])), func=func)
    if my_com is not None:
        if is_within(perimeter, my_com, surf_loc, surf_norm):
            return my_com, False
    # Loop through the points in the perimeter and choose the point that is the closest to the true center of mass
    min_dist, my_point = np.inf, None
    for point in perimeter:
        dist = calc_dist(point, true_com)
        if dist < min_dist:
            my_point, min_dist = point, dist
    return my_point, True


def project_to_hyperboloid(twoD_points, small_ball_loc, surf_func, plane_normal, plane_location):
    """
    Projects 2D points back onto a hyperboloid surface and selects valid points based on distance criteria.

    Parameters
    ----------
    twoD_points : numpy.ndarray
        Array of 2D points to be projected back to 3D
    small_ball_loc : numpy.ndarray
        Location of the smaller ball used for distance-based point selection
    surf_func : callable
        Function defining the hyperboloid surface
    plane_normal : numpy.ndarray
        Normal vector of the projection plane
    plane_location : numpy.ndarray
        Location point of the projection plane

    Returns
    -------
    list
        List of valid 3D points that lie on the hyperboloid surface

    Notes
    -----
    - Points are first unprojected from 2D to 3D space
    - Each 3D point is then projected onto the hyperboloid surface
    - Only points that successfully project onto the surface are included in the result
    - The process maintains the geometric relationship between the original points
    """
    # First we need to get the 2D points back to 3D
    plane_points = unproject_to_3d(twoD_points, plane_location, plane_normal)

    # Next each point needs to be projected onto the hyperboloid
    new_points = []
    for point in plane_points:
        new_point = calc_surf_point_from_plane(point, plane_normal, surf_func, small_ball_loc)
        if new_point is not None:
            new_points.append(new_point)
    # Return the new points
    return new_points


# Build method. Makes the mesh for the surface and calculates the simplices between them
def build_surf(locs, rads, epnts, res, net_type, sfunc=None, perimeter=None, surf_loc=None, surf_norm=None):
    """
    Main build method for constructing surfaces between two atoms/balls.

    This function handles the complete surface construction process including:
    - Surface function calculation
    - Perimeter construction
    - Center of mass determination
    - Point projection and triangulation
    - Curvature calculations

    Parameters
    ----------
    locs : list of numpy.ndarray
        List containing the 3D coordinates of the two atoms/balls
    rads : list of float
        List containing the radii of the two atoms/balls
    epnts : list of numpy.ndarray
        List of edge points defining the surface boundary
    res : int
        Resolution parameter controlling the density of surface points
    net_type : str
        Type of network being constructed ('prm', 'pow', or other)
    sfunc : callable, optional
        Pre-calculated surface function. If None, will be calculated internally
    perimeter : list, optional
        Pre-calculated perimeter. If None, will be calculated internally
    surf_loc : numpy.ndarray, optional
        Location of the surface plane. If None, will be calculated internally
    surf_norm : numpy.ndarray, optional
        Normal vector of the surface plane. If None, will be calculated internally

    Returns
    -------
    tuple
        A tuple containing:
        - spoints : list of numpy.ndarray
            List of 3D points defining the surface mesh
        - surf_tris : list of tuple
            List of triangles (indices into spoints) defining the surface mesh
        - mean_tri_curvs : list of float
            Mean curvature values for each triangle
        - mean_surf_curv : float
            Average mean curvature of the entire surface
        - gauss_tri_curvs : list of float
            Gaussian curvature values for each triangle
        - gauss_surf_curv : float
            Average Gaussian curvature of the entire surface
        - sfunc : callable
            The surface function used for calculations
        - surf_com : numpy.ndarray
            Center of mass of the surface
        - flat : bool
            Whether the surface is flat (True) or curved (False)
        - surf_loc : numpy.ndarray
            Location point of the surface plane
    """

    # Get the surface function if not already calculated
    if sfunc is None:
        sfunc = calc_surf_func(np.array(locs[0]), rads[0], np.array(locs[1]), rads[1])

    # Check if the surface is flat
    flat = False
    if net_type in {'prm', 'pow'} or rads[0] == rads[1]:
        flat = True

    # Build the perimeter of the surface
    if perimeter is None:
        perimeter, surf_loc, surf_norm = build_perimeter(locs, rads, epnts=epnts, net_type=net_type)
    
    # If the surface location and normal are not provided, calculate them
    if surf_loc is None or surf_norm is None:
        return

    # Get the center of mass for the surface
    surf_com, filter_hard = get_com(locs, rads, perimeter=perimeter, surf_loc=surf_loc, surf_norm=surf_norm, flat=flat,
                                    func=sfunc, net_type=net_type)

    # Calculate the angles to rotate the center point around
    flat_points = project_to_plane(np.array(perimeter), plane_normal=surf_norm, plane_point=surf_loc)

    # Calculate the flat COM
    flat_com, flat_loc = project_to_plane(np.array([surf_com, surf_loc]), plane_normal=surf_norm, plane_point=surf_loc)

    # Filter out the bad triangles
    my_2d_points, surf_tris = triangulate_2D_Surface(flat_points, res=res, center=flat_loc)

    # Calculate the curvature of the triangles and the surface
    if not flat:
        # Project the points onto the surface again
        spoints = project_to_hyperboloid(my_2d_points, locs[0], sfunc, surf_norm, surf_loc)
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv = calc_surf_tri_curvs(sfunc, spoints, surf_tris, curvature_type='mean')
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv = calc_surf_tri_curvs(sfunc, spoints, surf_tris, curvature_type='gauss')
    else:
        spoints = unproject_to_3d(my_2d_points, surf_loc, surf_norm)
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv = [0 for _ in range(len(list(surf_tris)))], 0, 0
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv = [0 for _ in range(len(list(surf_tris)))], 0, 0

    # Return the surface points, triangles, triangle curvatures, total curvature, surface function, com, and flatness
    return (spoints, surf_tris, mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv, gauss_tri_curvs, gauss_surf_curv,
            avg_gauss_surf_curv, sfunc, surf_com, flat, surf_loc)
