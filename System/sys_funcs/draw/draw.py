import warnings
import numpy as np
from System.sys_funcs.calcs.surf import calc_surf_func, calc_surf_tri_dists, calc_surf_tri_ins_out, calc_surf_tri_curvs
import matplotlib as mpl
from matplotlib._api.deprecation import MatplotlibDeprecationWarning as MPLDepWarn

warnings.filterwarnings('error')


def color_tris(surf, color_scheme, color_map, max_val=None):
    """
    Colors the triangles in the surface based on the specified coloring scheme and map
    :param surf:
    :param color_scheme: Determines how the colors will be mapped
    :param color_map: Determines the actual colors of triangles
    :return: The triangles in the surface are colored
    """
    # Set up the variable tri_colors for recording the colordesignations for
    tri_colors = None
    # Set up the color map
    try:
        my_cmap = mpl.colormaps.get_cmap(color_map)
    except MPLDepWarn:
        my_cmap = mpl.cm.get_cmap(color_map)
    except AttributeError:
        my_cmap = mpl.cm.get_cmap(color_map)
    # Default is distance based color map
    if color_scheme == 'dist':
        # Check if the tri_dists have been calculated before
        if surf['tri_dists'] is None or len(surf['tri_dists']) == 0 or len(surf['tri_dists']) != len(surf['tris']):
            calc_surf_tri_dists(surf['points'], surf['tris'], surf['loc'])
        tri_colors = [my_cmap(_) for _ in surf['tri_dists']]

    elif color_scheme == 'ins_out':
        # Check if the tri_dists have been calculated before
        if surf['tri_ins_out'] is None or len(surf['tri_ins_out']) == 0 or len(surf['tri_ins_out']) != len(surf['tris']):
            calc_surf_tri_ins_out(surf)
        tri_colors = [my_cmap(_) for _ in surf['tri_ins_out']]

    elif color_scheme == 'curv':
        # Check if the function is None
        if surf['func'] is None:
            a0, a1 = [surf['net'].balls.iloc[_] for _ in surf['balls']]
            func = calc_surf_func(a0['loc'], a0['rad'], a1['loc'], a1['rad'])
        else:
            func = surf['func']
        # Check if the tri_dists have been calculated before
        if surf['tri_curvs'] is None or len(surf['tri_curvs']) == 0 or len(surf['tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], max_curv=surf['curv'])
        else:
            tri_curvs = surf['tri_curvs']
        # First check if the surface is flat
        if surf['flat'] or surf['curv'] == 0:
            my_curvs = [0] * len(surf['tris'])
        else:
            my_curvs = [curv/max_val for curv in tri_curvs]

        # Set the colors
        tri_colors = [my_cmap(_) for _ in my_curvs]

    return tri_colors


def draw_line(points, radius=0.02, color=None, edge_org=None):
    """
    Draws a line from point to points
    :param points: two points
    :param radius: Radius for the line to be drawn
    :param color: Color of the line
    :param edge_org: Vector for the edge orientation
    :return: points and tris for drawing
    """
    if edge_org is None:
        edge_org = [0, 0, 1]
    # Initiate the draw attributes
    draw_points, draw_tris = [], []
    r = None
    # Go through the points
    for i in range(len(points)):
        # If we are at the end of the points list, use the previous point for calibration
        p0 = np.array(points[i])
        if i < len(points) - 1:
            p1 = np.array(points[i + 1])
            r = p1 - p0
        # Find the vector and its normal between the two points
        rn = r / np.linalg.norm(r)
        # In the case that the vector between the points is in the z direction only, move it
        if rn[0] == 0 and rn[1] == 0:
            r = r + np.array([0.001, 0.001, 0])
            rn = r / np.linalg.norm(r)
        # Take the cross product with the +z direction and normalize it
        v0_0x = np.cross(rn, np.array(p0 - edge_org))
        v0_0n = v0_0x / np.linalg.norm(v0_0x)
        # Calculate the location of the first point
        p0_0 = v0_0n * radius + p0
        # Take the cross product of the edge vector and the vector to the first point and normalize it
        v0_1x = np.cross(rn, v0_0n)
        v0_1nx = v0_1x / np.linalg.norm(v0_1x)
        # Find the vectors for the other two points (30/60/90 triangle)
        v0_1 = - 0.5 * radius * v0_0n + 0.5 * np.sqrt(3) * radius * v0_1nx
        v0_2 = - 0.5 * radius * v0_0n - 0.5 * np.sqrt(3) * radius * v0_1nx
        # Get the points and add them to the list of draw points
        p0_1, p0_2 = v0_1 + p0, v0_2 + p0
        draw_points += [p0_0, p0_1, p0_2]
    # Go through the points
    for i in range(len(points) - 1):
        # List the points
        p0_0, p0_1, p0_2, p1_0, p1_1, p1_2 = range(3 * i, 3 * (i + 2))
        # Create the triangles
        draw_tris += [[p0_0, p0_1, p1_0], [p1_0, p1_1, p0_1],
                           [p0_1, p0_2, p1_1], [p1_1, p1_2, p0_2],
                           [p0_2, p0_0, p1_2], [p1_2, p1_0, p0_0]]
    # Return the points and triangles
    return draw_points, draw_tris


# Draw Edge Function. Takes in an edge and updates its attributes draw_points, draw_tris
def draw_edge(edge, radius=0.02, color=None):
    """
    Draws an edge in triangles and points
    :param edge: Edge object for exporting
    :param radius: Radius of the edge to be drawn
    :param color: Color for the edge drawing
    :return: None
    """
    # # Get the edge direction to point away from
    # rads = [_.rad for _ in edge.balls]
    # min_ball = edge['balls'][rads.index(min(rads))]
    # if edge.points is None or len(edge.points) <= 1:
    #     edge.points, edge.vals = build_edge(edge['balls'], edge['verts'], edge.net.surf_res)
    # Calculate the lines
    return draw_line(edge.points, radius, color=color)
