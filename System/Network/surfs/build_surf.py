import time

from System.sys_funcs.calcs.surf import calc_surf_func, calc_surf_tri_curvs, calc_surf_point_curv
from System.Network.surfs.perimeter import build_perimeter
from System.Network.surfs.fill import fill_mesh, calc_surf_point
from System.sys_funcs.calcs.calcs import calc_com, project_to_plane
import numpy as np
from shapely.plotting import plot_polygon
from shapely import Polygon, Point
from System.Network.surfs.triangulate import triangulate_2D_Surface, test_within
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


############################################## Triangulate Surface Points  #############################################

def plot_points_and_tris(pnts=None, trs=None, pcol=None, tcol=None, plot_points=True, Show=False):

    if trs is not None:
        for tri in trs:
            p0, p1, p2 = [pnts[_] for _ in tri]
            plt.plot([p0[0], p1[0], p2[0], p0[0]], [p0[1], p1[1], p2[1], p0[1]], c=tcol)
    if pnts is not None and plot_points:
        plt.scatter([_[0] for _ in pnts], [_[1] for _ in pnts], c=pcol)
    if Show:
        plt.show()

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
    if test_within(perimeter, surf_loc, surf_loc, surf_norm):
        return surf_loc, False
    # First try the center of mass of the 3d points projected onto the surface
    my_com = calc_surf_point(locs, point=calc_com(points=np.array(perimeter)), func=func)
    if my_com is not None:
        if test_within(perimeter, my_com, surf_loc, surf_norm):
            return my_com, False
    # Next try to calculate a center of mass of some of the points
    my_com = calc_surf_point(locs, point=calc_com(points=np.array(perimeter[::5])), func=func)
    if my_com is not None:
        if test_within(perimeter, my_com, surf_loc, surf_norm):
            return my_com, False

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
def build_surf(locs, rads, epnts, res, net_type, sfunc=None, check=False, timer=False):
    """
    Main build method for constructing surfaces
    :param surf:
    :param res: Specifies the resolution the surface is to be constructed with
    :return: The surfaces points and triangles are filled
    """
    # if timer:
    #     start = time.perf_counter()
    #     clocck = {'calc_func': 0, 'perimeter': 0, 'com': 0, 'fill_mesh': 0, 'spider': 0, 'Delaunay': 0,
    #               'designations': 0, 'reassign': 0}


    # Get the surface function if not already calculated
    if sfunc is None:
        sfunc = calc_surf_func(np.array(locs[0]), rads[0], np.array(locs[1]), rads[1])


    # if timer:
    #     clocck['calc_func'] = time.perf_counter() - start
    #     start = time.perf_counter()


    # Check if the surface is flat
    flat = False
    if net_type in {'del', 'pow'} or rads[0] == rads[1]:
        flat = True

    # Build the perimeter of the surface
    perimeter, surf_loc, surf_norm = build_perimeter(locs, rads, epnts=epnts, net_type=net_type)
    # if check:
    #     print("Flat = {}".format(flat))

    # if timer:
    #     clocck['perimeter'] = time.perf_counter() - start
    #     start = time.perf_counter()
    # Get the center of mass for the surface
    surf_com, filter_hard = get_com(locs, rads, perimeter=perimeter, surf_loc=surf_loc, surf_norm=surf_norm, flat=flat,
                                    func=sfunc, net_type=net_type)
    # if timer:
    #     clocck['com'] = time.perf_counter() - start
    #     start = time.perf_counter()
    if net_type in {'del', 'pow'}:
        flat_points = project_to_plane(np.array(perimeter + [surf_com]), plane_normal=surf_norm, plane_point=surf_loc)
        triangles = Delaunay(flat_points).simplices


        # if timer:
        #     clocck['fill_mesh'] = time.perf_counter() - start
        #     return np.array(perimeter + [surf_com]), triangles, [0 for _ in range(len(triangles))], 0.0, None, surf_com, True, clocck
        return np.array(perimeter + [surf_com]), triangles, [0 for _ in range(len(triangles))], 0.0, None, surf_com, True

    # Fill the mesh
    if not flat:
        spoints = fill_mesh(locs, rads, func=sfunc, surf_loc=surf_loc, surf_norm=surf_norm, perimeter=perimeter,
                            com=surf_com, flat=flat, res=res, check=check)
    else:
        spoints = perimeter + [surf_com]

    # if timer:
    #     clocck['fill_mesh'] = time.perf_counter() - start
    #     start = time.perf_counter()
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
    # plot_polygon(perim_poly)
    # plot_points_and_tris(surf_points, surf_tris, plot_points=False, Show=True)
    # if timer:
    #     for _ in slimer:
    #         clocck[_] = slimer[_]
    # plot_surfs([spoints], [surf_tris], True, Show=True)
    # Calculate the curvature of the triangles and the surface
    if not flat:
        tri_curvs, surf_curv = calc_surf_tri_curvs(sfunc, spoints, surf_tris, max_curv=surf_curv)
    else:
        tri_curvs, surf_curv = [0 for _ in range(len(list(surf_tris)))], 0
    # Return the surface points, triangles, triangle curvatures, total curvature, surface function, com, and flatness
    return spoints, surf_tris, tri_curvs, surf_curv, sfunc, surf_com, flat
