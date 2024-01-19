from System.sys_funcs.calcs.calcs import calc_surf_tri_curvs, calc_surf_point_curv
from System.sys_funcs.calcs.surf import calc_surf_func
from System.Network.surfs.perimeter import build_perimeter
from System.Network.surfs.fill import fill_mesh, calc_surf_point
from System.Network.surfs.triangulate import find_simps, tri_within, filter_tris
from System.sys_funcs.calcs.calcs import calc_com
import numpy as np


############################################## Triangulate Surface Points  #############################################

def get_com(alocs, arads, perimeter, surf_loc, surf_norm, func, flat, net_type='vor'):
    """
    Finds the center of mass of a surface's perimeter points
    :param net_type:
    :param surf: Surface object holding the perimeter points
    :return: Center of mass
    """
    if flat or net_type in {'del', 'pow'}:
        return calc_com(points=np.array(perimeter)), False
        # If the surface is flat, the center of mass will not need to be projected
    if tri_within(perimeter, loc=surf_loc, norm=surf_norm, point=surf_loc):
        return surf_loc, False
    # First try the center of mass of the 3d points projected onto the surface
    my_com = calc_surf_point(alocs, point=calc_com(points=np.array(perimeter[::5])), func=func)
    if my_com is not None and tri_within(perimeter, surf_loc, surf_norm, point=my_com) and arads[0] != arads[1]:
        return my_com, False
    # If nothing else set the center of mass to the first point in the perimeter
    return perimeter[len(perimeter)//2], True


# Build method. Makes the mesh for the surface and calculates the simplices between them
def build_surf(alocs, arads, epnts, res, net_type, sfunc=None):
    """
    Main build method for constructing surfaces
    :param surf:
    :param res: Specifies the resolution the surface is to be constructed with
    :return: The surfaces points and triangles are filled
    """
    # Get the surface function if not already calculated
    if sfunc is None:
        sfunc = calc_surf_func(np.array(alocs[0]), arads[0], np.array(alocs[1]), arads[1])
    # Check if the surface is flat
    flat = False
    if net_type in {'del', 'pow'} or arads[0] == arads[1]:
        flat = True
    # Build the perimeter of the surface
    perimeter, surf_loc, surf_norm = build_perimeter(alocs, arads, epnts=epnts, net_type=net_type)
    # Get the center of mass for the surface
    surf_com, filter_hard = get_com(alocs, arads, perimeter=perimeter, surf_loc=surf_loc, surf_norm=surf_norm, flat=flat, func=sfunc)
    # Fill the mesh
    spoints = fill_mesh(alocs, arads, func=sfunc, surf_loc=surf_loc, surf_norm=surf_norm, perimeter=perimeter, com=surf_com, flat=flat, res=res)
    # Find the simplices of the surface
    tris, flat_points = find_simps(points=spoints, loc=surf_loc, norm=surf_norm)
    # Set the surface curvature to 0
    surf_curv = 0
    # If the network type is voronoi the edges could be curved allowing for triangulations outside the edges
    if net_type == 'vor':
        # Add the normal curvature if possible
        if tri_within(perimeter=perimeter,  loc=surf_loc, norm=surf_norm, flat_points=flat_points, point=surf_loc):
            surf_curv = calc_surf_point_curv(sfunc, surf_loc)
        # Filter out the bad triangles
        surf_tris = filter_tris(tris=tris, flat_points=flat_points, res=res, perimeter=perimeter, loc=surf_loc, norm=surf_norm, filter_hard=filter_hard)
        # Calculate the curvature of the triangles and the surface
        if not flat:
            tri_curvs, surf_curv = calc_surf_tri_curvs(sfunc, spoints, surf_tris, max_curv=surf_curv)
        else:
            tri_curvs, surf_curv = [0 for _ in range(len(list(surf_tris)))], 0
    else:
        surf_tris = tris
        tri_curvs, surf_curv = [0 for _ in range(len(surf_tris))], 0
    return spoints, surf_tris, tri_curvs, surf_curv, sfunc, surf_com, flat
