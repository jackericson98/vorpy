import numpy as np
import warnings
import matplotlib as mpl
from matplotlib._api.deprecation import MatplotlibDeprecationWarning as MPLDepWarn
from vorpy.src.calculations.curvature import gaussian_curvature
from vorpy.src.calculations.curvature import mean_curvature
from vorpy.src.calculations.surf import calc_surf_tri_dists
from vorpy.src.calculations.surf import calc_surf_func
from vorpy.src.calculations.calcs import calc_dist
warnings.filterwarnings('error')


def calc_surf_tri_curvs(func, points, tris, curvature_type='gauss'):
    """
    Calculates the curvature values for triangles in a surface defined by an implicit function.

    This function computes either Gaussian or mean curvature values for each triangle in a surface
    by averaging the curvature values at the triangle's vertices. The curvature type can be specified
    to calculate either Gaussian or mean curvature.

    Parameters:
    -----------
    func : list of float
        List of coefficients for the quadratic surface equation:
        [A, B, C, D, E, F, G, H, I, J] representing:
        Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0
    points : list of numpy.ndarray
        List of vertex coordinates [x, y, z] for the surface
    tris : list of list of int
        List of triangles, where each triangle is defined by three indices into the points list
    curvature_type : str, optional
        Type of curvature to calculate: 'gauss' for Gaussian curvature (default) or 'mean' for mean curvature

    Returns:
    --------
    tuple
        A tuple containing:
        - tri_curvs: List of curvature values for each triangle
        - max_curv: Maximum curvature value found across all vertices
    """
    # Set up the curvs variable
    curvs = []
    # Set up the min and max curvature variables
    min_curv, max_curv = np.inf, 0
    # If the surface normal is within the surface,
    # Get the curvature for each point
    for point in points:
        # Get the curvature
        if curvature_type == 'gauss':
            curv = gaussian_curvature(func, point)
        else:
            curv = mean_curvature(func, point)
        # Record the min and max curvatures
        if curv < min_curv:
            min_curv = curv
        elif curv > max_curv:
            max_curv = curv
        # Add the curvature to the list
        curvs.append(curv)
    # Set up the tri_curvs list
    tri_curvs = []
    # Go through the curvature values for each point
    for i in range(len(tris)):
        # Get the triangle
        tri = tris[i]
        # Get the curvatures
        try:
            curv_val = sum([curvs[_] for _ in tri])/3
        except IndexError:
            print(len(points), tri)
        # Add the curve value to the surface's list of curvatures
        tri_curvs.append(curv_val)
    # Return the values
    return tri_curvs, max_curv


def calc_surf_tri_ins_out(b0_loc, b0_rad, surf):
    """
    Determines which triangles in a surface lie within the overlapping region of two spheres.

    This function analyzes each triangle in a surface to determine if it lies within the region
    where two spheres overlap. It does this by checking if all vertices of each triangle are
    within the radius of the first sphere.

    Parameters:
    -----------
    b0_loc : numpy.ndarray
        Center coordinates [x, y, z] of the first sphere
    b0_rad : float
        Radius of the first sphere
    surf : dict
        Dictionary containing surface information including:
        - 'points': List of vertex coordinates [x, y, z]
        - 'tris': List of triangles defined by vertex indices

    Returns:
    --------
    list of float
        List of values (0.25 or 0.75) indicating whether each triangle lies within the
        overlapping region (0.25) or outside it (0.75)
    """
    # Set up a list of tracking
    inside_array = []
    # Go through the points in the surface
    for point in surf['points']:
        # Calculate the distance between the point and the ball
        my_dist = calc_dist(point, b0_loc)
        # Check if the triangle is inside or not
        if my_dist < b0_rad:
            inside_array.append(True)
        else:
            inside_array.append(False)
    # Now add the triangles
    tri_ins_out = []
    # Color the tris based on the inside_array
    for tri in surf.tris:
        # Check if the triangle is inside or not
        if inside_array[tri[0]] and inside_array[tri[1]] and inside_array[tri[2]]:
            # Add the value to the list
            tri_ins_out.append(0.25)
        else:
            # Add the value to the list
            tri_ins_out.append(0.75)
    # Return the list
    return tri_ins_out


def calc_surf_tri_dists(points, tris, loc):
    """
    Calculates the normalized distances between each triangle in a surface and a reference location.

    This function computes the distance between each triangle in a surface and a specified location,
    then normalizes these distances to a range between 0 and 1. The normalization is based on the
    minimum and maximum distances found across all points in the surface.

    Parameters:
    -----------
    points : list of numpy.ndarray
        List of 3D point coordinates [x, y, z] that form the vertices of the triangles
    tris : list of tuples
        List of triangles, where each triangle is represented as a tuple of three indices
        corresponding to points in the points array
    loc : numpy.ndarray
        Reference location coordinates [x, y, z] for distance calculations

    Returns:
    --------
    list of float
        List of normalized distances (0 to 1) corresponding to each triangle in the surface,
        where each distance represents the maximum distance between the triangle's vertices
        and the reference location
    """
    # Set up the distances
    dists = []
    tri_dists = []
    max_dist, min_dist = 0, np.inf
    # Provide value for the points
    for point in points:
        # Calculate the distance
        my_dist = calc_dist(point, loc)
        dists.append(my_dist)
        # Record the minimum and maximum distances
        if my_dist < min_dist:
            min_dist = my_dist
        elif my_dist > max_dist:
            max_dist = my_dist
    # Go through the triangles in the surface
    for i in range(len(tris)):
        # Find the maximum distance point of the triangles
        tri_dists.append(max([dists[_] for _ in tris[i]]))
    # Normalize the tri_dists
    return [(_ - min_dist) / (max_dist - min_dist) for _ in tri_dists]


def color_tris(surf, color_scheme, color_map, color_factor, max_val=None):
    """
    Colors the triangles in a surface based on specified coloring scheme and map.

    This function applies color mapping to surface triangles based on various geometric properties
    such as distance, curvature, or other surface characteristics. The coloring can be modified
    using different color schemes and transformation factors.

    Args:
        surf (dict): Surface data structure containing points, triangles, and other geometric properties
        color_scheme (str): Coloring scheme to apply. Options include:
            - 'dist': Distance-based coloring
            - 'ins_out': Inside/outside coloring
            - 'mean': Mean curvature-based coloring
        color_map (str): Name of the matplotlib colormap to use for coloring
        color_factor (str): Transformation to apply to values before coloring. Options:
            - 'log': Logarithmic transformation
            - 'sqr': Square transformation
            - 'cub': Cubic transformation
            - None: No transformation
        max_val (float, optional): Maximum value for color normalization. If None, uses maximum
            value from the data.

    Returns:
        None: The function modifies the surface data structure in place by adding color information
        to the triangles.

    Notes:
        - The function uses matplotlib colormaps for color generation
        - Color values are stored as RGBA tuples in the range [0,1]
        - The function handles various edge cases including flat surfaces and missing data
        - Color calculations are cached when possible to improve performance
    """
    # Set up the variable tri_colors for recording the color designations for
    tri_colors = None
    # Set up the color map
    try:
        my_cmap = mpl.colormaps.get_cmap(color_map)
    except MPLDepWarn:
        my_cmap = mpl.cm.get_cmap(color_map)
    except AttributeError:
        my_cmap = mpl.cm.get_cmap(color_map)
    except ValueError:
        my_cmap = mpl.cm.get_cmap(color_map.capitalize())
    # Create the surface value multiplier
    if color_factor == 'log':
        def multi(val):
            try:
                return np.log(val + 1)
            except RuntimeWarning:
                return np.log(abs(val + 1))
    elif color_factor == 'sqr':
        def multi(val):
            return val ** 2
    elif color_factor == 'cub':
        def multi(val):
            return val ** 3
    else:
        def multi(val):
            return val
    # Default is distance based color map
    if color_scheme == 'dist':
        # Check if the tri_dists have been calculated before
        if 'tri_dists' not in surf or surf['tri_dists'] is None or len(surf['tri_dists']) == 0 or len(surf['tri_dists']) != len(surf['tris']):
            tri_dists = calc_surf_tri_dists(surf['points'], surf['tris'], surf['loc'])
            tri_colors = [my_cmap(multi(_)) for _ in tri_dists]
        else:
            tri_colors = [my_cmap(multi(_)) for _ in surf['tri_dists']]

    elif color_scheme == 'ins_out':
        # Check if the tri_dists have been calculated before
        tri_colors = [my_cmap(multi(_)) for _ in surf['tris_ins_out']]

    elif color_scheme.lower() == 'mean' or color_scheme.lower() == 'mean_curv' or color_scheme.lower() == 'mean curvature':
        # Check if the function is None
        if surf['func'] is None:
            a0, a1 = [surf['net'].balls.iloc[_] for _ in surf['balls']]
            func = calc_surf_func(a0['loc'], a0['rad'], a1['loc'], a1['rad'])
        else:
            func = surf['func']
        # Check if the tri_dists have been calculated before
        if surf['mean_tri_curvs'] is None or len(surf['mean_tri_curvs']) == 0 or len(surf['mean_tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='mean')
        else:
            tri_curvs = surf['mean_tri_curvs']
        # First check if the surface is flat
        if surf['flat'] or surf['mean_curv'] == 0:
            my_curvs = [0] * len(surf['tris'])
        else:
            my_curvs = [curv/max_val for curv in tri_curvs]

        # Set the colors
        tri_colors = [my_cmap(multi(_)) for _ in my_curvs]
    elif color_scheme.lower() == 'gauss' or color_scheme.lower() == 'gauss_curv' or color_scheme.lower() == 'gaussian curvature' or color_scheme.lower() == 'gaus':
        # Check if the function is None
        if surf['func'] is None:
            a0, a1 = [surf['net'].balls.iloc[_] for _ in surf['balls']]
            func = calc_surf_func(a0['loc'], a0['rad'], a1['loc'], a1['rad'])
        else:
            func = surf['func']
        # Check if the tri_dists have been calculated before
        if surf['gauss_tri_curvs'] is None or len(surf['gauss_tri_curvs']) == 0 or len(surf['gauss_tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='gauss')
        else:
            tri_curvs = surf['gauss_tri_curvs']
        # First check if the surface is flat
        if surf['flat'] or surf['gauss_curv'] == 0:
            my_curvs = [0] * len(surf['tris'])
        elif max_val is None or max_val == 0:
            max_val = max(tri_curvs)
            my_curvs = [curv/max_val for curv in tri_curvs]
        else:
            my_curvs = [curv/max_val for curv in tri_curvs]

        # Set the colors
        tri_colors = [my_cmap(multi(_)) for _ in my_curvs]
    else:
        tri_colors = [np.random.randint(0, 256, size=3) for _ in range(len(surf['tris']))]
    return tri_colors
