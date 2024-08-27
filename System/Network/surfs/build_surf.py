from System.sys_funcs.calcs.surf import calc_surf_func, calc_surf_tri_curvs, calc_surf_point_curv
from System.Network.surfs.perimeter import build_perimeter
from System.Network.surfs.fill import fill_mesh, calc_surf_point
from System.Network.surfs.triangulate import find_simps, tri_within, filter_tris, tri_within2
from System.sys_funcs.calcs.calcs import calc_com, rotate_points, project_to_plane, map_to_plane
from Visualize.mpl_visualize import plot_balls, plot_surfs, plot_edges
import numpy as np
import triangle as tri
import matplotlib.pyplot as plt
from shapely import Polygon, Point
from System.Network.surfs.ConstrainedDelaunayTriangulation1 import triangulate_2D_Surface
from scipy.spatial import Delaunay


############################################## Triangulate Surface Points  #############################################


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
    flat_perim = project_to_plane(np.array(perimeter), plane_normal=surf_norm, plane_point=surf_loc)
    perim_poly = Polygon(flat_perim)
    # If the surface is flat, the center of mass will not need to be projected
    if perim_poly.contains(Point([0, 0])):
        return surf_loc, False
    # First try the center of mass of the 3d points projected onto the surface
    my_com = calc_surf_point(locs, point=calc_com(points=np.array(perimeter)), func=func)
    if my_com is not None:
        flat_com = project_to_plane(np.array([my_com]), plane_normal=surf_norm, plane_point=surf_loc)
        try:
            if perim_poly.contains(Point(flat_com[0])):
                return my_com, False
        except TypeError:
            pass
    # Next try to calculate a center of mass of some of the points
    my_com = calc_surf_point(locs, point=calc_com(points=np.array(perimeter[::5])), func=func)
    if my_com is not None:
        flat_com = project_to_plane(np.array([my_com]), plane_normal=surf_norm, plane_point=surf_loc)
        try:
            if perim_poly.contains(Point(flat_com[0])):
                return my_com, False
        except TypeError:
            pass

    # # Next choose a point within the flat perimeter
    # point_within = perim_poly.point_on_surface()
    # # Project onto the surface
    # plane_point_mapped = map_to_plane([[point_within.x, point_within.y]], plane_normal=surf_norm, plane_point=surf_loc)
    # # Calculate the surface point
    # on_surface_point_within = calc_surf_point(locs, func=func, point=plane_point_mapped)
    # # Project back onto the plane
    # back_on_plane = project_to_plane([on_surface_point_within], surf_loc, surf_norm)
    # # See if it is inside or not
    # if perim_poly.contains(Point(back_on_plane[0])):
    #     return on_surface_point_within, False
    # If nothing else set the center of mass to the first point in the perimeter
    return perimeter[len(perimeter)//2], True


# Build method. Makes the mesh for the surface and calculates the simplices between them
def build_surf(locs, rads, epnts, res, net_type, sfunc=None, check=False):
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
    if net_type in {'del', 'pow'}:
        flat_points = project_to_plane(np.array(perimeter + [surf_com]), plane_normal=surf_norm, plane_point=surf_loc)
        triangles = Delaunay(flat_points).simplices
        return np.array(perimeter + [surf_com]), triangles, [0 for _ in range(len(triangles))], 0.0, None, surf_com, True
    # Fill the mesh
    spoints = fill_mesh(locs, rads, func=sfunc, surf_loc=surf_loc, surf_norm=surf_norm, perimeter=perimeter,
                        com=surf_com, flat=flat, res=res, check=check)

    # Calculate the angles to rotate the center point around
    flat_points = project_to_plane(np.array(spoints), plane_normal=surf_norm, plane_point=surf_loc)

    # Create the polygon from the flat perimeter
    flat_perim = flat_points[:len(perimeter)]
    perim_poly = Polygon(flat_perim)

    # Set the surface curvature to 0
    surf_curv = 0
    # If the network type is voronoi the edges could be curved allowing for triangulations outside the edges
    # Add the normal curvature if possible
    if net_type == 'aw' and not flat and perim_poly.contains(Point(surf_loc)):
        surf_curv = calc_surf_point_curv(sfunc, surf_loc)
    # Filter out the bad triangles
    surf_points, surf_tris = triangulate_2D_Surface(flat_perim, flat_points, res, surf_loc)
    # plot_surfs([spoints], [surf_tris], True, Show=True)
    # Calculate the curvature of the triangles and the surface
    if not flat:
        tri_curvs, surf_curv = calc_surf_tri_curvs(sfunc, spoints, surf_tris, max_curv=surf_curv)
    else:
        tri_curvs, surf_curv = [0 for _ in range(len(list(surf_tris)))], 0
    # Return the surface points, triangles, triangle curvatures, total curvature, surface function, com, and flatness
    return spoints, surf_tris, tri_curvs, surf_curv, sfunc, surf_com, flat


def build_surf1(locs, rads, epnts, res, net_type, sfunc=None):
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
                                    func=sfunc)
    # Fill the mesh
    # spoints = fill_mesh(locs, rads, func=sfunc, surf_loc=surf_loc, surf_norm=surf_norm, perimeter=perimeter,
    #                     com=surf_com, flat=flat, res=res)

    flat_perims = project_to_plane(np.array(perimeter), np.array(surf_loc), np.array(surf_norm))
    # Triangulate the flat points
    perim_dict = dict(vertices=np.array([[round(_[0], 3), round(_[1], 3)] for _ in flat_perims]),
                      segments=np.array([[i, (i + 1) % len(flat_perims)] for i in range(len(flat_perims))]))

    triangulation = tri.triangulate(perim_dict, 'a' + str(round(np.sqrt(3) / 4 * res ** 2, 3)))
    spoints = triangulation['vertices']
    # Check if we need to filter the points
    new_spoints = spoints
    if filter_hard:
        flat_perim = project_to_plane(np.array(perimeter), np.array(surf_loc), np.array(surf_norm))
        flat_points = project_to_plane(np.array(spoints), np.array(surf_loc), np.array(surf_norm))
        new_spoints = []
        for i, point in enumerate(flat_points):
            if tri_within2(np.array(flat_perim), np.array(point)):
                new_spoints.append(spoints[i])

    # Find the simplices of the surface
    tris, flat_points = find_simps(points=new_spoints, loc=surf_loc, norm=surf_norm)
    # Set the surface curvature to 0
    surf_curv = 0
    # If the network type is voronoi the edges could be curved allowing for triangulations outside the edges
    if net_type == 'aw':
        # Add the normal curvature if possible
        if tri_within(perimeter=perimeter,  loc=surf_loc, norm=surf_norm, flat_points=flat_points, point=surf_loc):
            surf_curv = calc_surf_point_curv(sfunc, surf_loc)
        # Filter out the bad triangles
        surf_tris = filter_tris(tris=tris, flat_points=flat_points, res=res, perimeter=perimeter, loc=surf_loc,
                                norm=surf_norm, filter_hard=filter_hard)
        # Calculate the curvature of the triangles and the surface
        if not flat:
            tri_curvs, surf_curv = calc_surf_tri_curvs(sfunc, new_spoints, surf_tris, max_curv=surf_curv)
        else:
            tri_curvs, surf_curv = [0 for _ in range(len(list(surf_tris)))], 0
    else:
        surf_tris = tris
        tri_curvs, surf_curv = [0 for _ in range(len(surf_tris))], 0

    # Return the surface points, triangles, triangle curvatures, total curvature, surface function, com, and flatness
    return new_spoints, tris, tri_curvs, surf_curv, sfunc, flat


def build_surf1(locs, rads, epnts, res, net_type, sfunc=None):
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
    # Get the 2D Perimeter
    # Copy the perimeter points
    points = perimeter.copy()
    # Move all surf points toward the origin via center point
    # for i in range(len(points)):
    #     points[i] = np.array(points[i]) - np.array(surf_loc)
    # Calculate the angles to rotate the center point around
    flat_points = project_to_plane(np.array(points), np.array(surf_loc), np.array(surf_norm))
    # Get the 2d version of the points and their Delaunay tesselation
    # flat_points = [_[:2] for _ in nps]
    plt.scatter([round(_[0], 3) for _ in flat_points], [round(_[1], 3) for _ in flat_points])
    plt.xlim([-5, 5])
    plt.ylim([-5, 5])
    plt.show()
    # Triangulate the flat points
    perim_dict = dict(vertices=np.array([[round(_[0], 3), round(_[1], 3)] for _ in flat_points]),
                      segments=np.array([[i, (i + 1) % len(flat_points)] for i in range(len(flat_points))]))

    try:
        triangulation = tri.triangulate(perim_dict, 'a' + str(round(np.sqrt(3)/4 * res ** 2, 3)))
    except RuntimeError:
        print("here")
        plt.scatter([round(_[0], 3) for _ in flat_points], [round(_[1], 3) for _ in flat_points])
        plt.show()

    # Rotate the flat points back
    triangulated_points = map_to_plane(triangulation['vertices'], surf_loc, plane_normal=surf_norm)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter([_[0] for _ in triangulated_points], [_[1] for _ in triangulated_points], [_[2] for _ in triangulated_points])
    plt.show()

    triangulated_points = [_ + calc_surf_point(locs, surf_loc, sfunc) for _ in triangulated_points]

    ax.scatter([_[0] for _ in triangulated_points], [_[1] for _ in triangulated_points], [_[2] for _ in triangulated_points])
    ax.scatter([_[0] for _ in perimeter], [_[1] for _ in perimeter], [_[2] for _ in perimeter])
    plot_balls(locs, rads, ['k', 'k'], fig, ax, alpha=0.1)
    spoints = []
    for point in triangulated_points:
        spoints.append(calc_surf_point(locs, point, sfunc))

    plot_surfs([spoints], [triangulation['triangles']], fig=fig, ax=ax, Show=True)
    tri_curvs, surf_curv = calc_surf_tri_curvs(sfunc, spoints, triangulation['triangles'])
    # Return the surface points, triangles, triangle curvatures, total curvature, surface function, com, and flatness
    return np.array(spoints), np.array(triangulation['triangles']), tri_curvs, surf_curv, sfunc, flat
