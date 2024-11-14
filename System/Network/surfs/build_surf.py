from System.sys_funcs.calcs.surf import calc_surf_func, calc_surf_tri_curvs, calc_surf_point_curv
from System.Network.surfs.perimeter import build_perimeter
from System.Network.surfs.fill import calc_surf_point, calc_surf_point_from_plane
from System.sys_funcs.calcs.calcs import calc_com, project_to_plane, calc_dist, unproject_to_3d
import numpy as np
from shapely import Polygon, Point
from System.Network.surfs.triangulate import triangulate_2D_Surface, is_within


def get_com(locs, rads, perimeter, surf_loc, surf_norm, func, flat, net_type='aw'):
    """
    Finds the center of mass of a surface's perimeter points
    :param net_type:
    :param surf: Surface object holding the perimeter points
    :return: Center of mass
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
    Projects the points back onto the hyperboloid and chooses the correct point based on distance from small loc
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
def build_surf(locs, rads, epnts, res, net_type, sfunc=None):
    """
    Main build method for constructing surfaces
    :param surf:
    :param res: Specifies the resolution the surface is to be constructed with
    :return: The surfaces points and triangles are filled
    """

    # Get the surface function if not already calculated
    if sfunc is None:
        sfunc = calc_surf_func(np.array(locs[0]), rads[0], np.array(locs[1]), rads[1])

    # Check if the surface is flat
    flat = False
    if net_type in {'del', 'pow'} or rads[0] == rads[1]:
        flat = True

    # Build the perimeter of the surface
    perimeter, surf_loc, surf_norm = build_perimeter(locs, rads, epnts=epnts, net_type=net_type)

    # Get the center of mass for the surface
    surf_com, filter_hard = get_com(locs, rads, perimeter=perimeter, surf_loc=surf_loc, surf_norm=surf_norm, flat=flat,
                                    func=sfunc, net_type=net_type)

    # Calculate the angles to rotate the center point around
    flat_points = project_to_plane(np.array(perimeter), plane_normal=surf_norm, plane_point=surf_loc)

    # Calculate the flat COM
    flat_com, flat_loc = project_to_plane(np.array([surf_com, surf_loc]), plane_normal=surf_norm, plane_point=surf_loc)


    perim_poly = Polygon(flat_points)

    # Set the surface curvature to 0
    surf_curv = 0
    # If the network type is voronoi the edges could be curved allowing for triangulations outside the edges
    # Add the normal curvature if possible
    if net_type == 'aw' and not flat and perim_poly.contains(Point(surf_loc)):
        surf_curv = calc_surf_point_curv(sfunc, surf_loc)
    # Filter out the bad triangles
    my_2d_points, surf_tris = triangulate_2D_Surface(flat_points, res=res, center=flat_loc)

    # Calculate the curvature of the triangles and the surface
    if not flat:
        # Project the points onto the surface again
        spoints = project_to_hyperboloid(my_2d_points, locs[0], sfunc, surf_norm, surf_loc)
        tri_curvs, surf_curv = calc_surf_tri_curvs(sfunc, spoints, surf_tris, max_curv=surf_curv)
    else:
        spoints = unproject_to_3d(my_2d_points, surf_loc, surf_norm)
        tri_curvs, surf_curv = [0 for _ in range(len(list(surf_tris)))], 0

    # Return the surface points, triangles, triangle curvatures, total curvature, surface function, com, and flatness
    return spoints, surf_tris, tri_curvs, surf_curv, sfunc, surf_com, flat, surf_loc
