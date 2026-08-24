import numpy as np
import warnings
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from vorpy.src.calculations.curvature import gaussian_curvature
from vorpy.src.calculations.curvature import mean_curvature
from vorpy.src.calculations.surf import calc_surf_tri_dists
from vorpy.src.calculations.surf import calc_surf_func
from vorpy.src.calculations.calcs import calc_dist
warnings.filterwarnings('error')


def get_cmap_safe(name: str, fallback: str = "viridis"):
    """
    If user passes: 'solid:red', 'solid #00ff00', 'solid(0.1,0.8,0.3)' etc.
    return a colormap where every entry is that color.
    Otherwise return a normal Matplotlib colormap.
    """
    name = name.strip()

    # --- SOLID COLOR HANDLING ---
    if name.lower().startswith("solid"):
        # extract color after "solid:"
        # examples:
        #   solid:red
        #   solid:#00ff00
        #   solid:(0.1,0.8,0.3)
        if ":" in name:
            color = name.split(":", 1)[1].strip()
        else:
            color = "gray"  # default solid color

        try:
            rgba = mpl.colors.to_rgba(color)
        except ValueError:
            rgba = mpl.colors.to_rgba("gray")

        # single-color colormap
        return ListedColormap([rgba], name=f"solid_{color}")

    # --- NORMAL COLORMAP HANDLING ---
    if hasattr(mpl, "colormaps"):
        getter = mpl.colormaps.get_cmap
    else:
        getter = mpl.cm.get_cmap

    # try raw name and variants
    for n in (name, name.lower(), name.upper(), name.capitalize()):
        try:
            return getter(n)
        except Exception:
            pass

    # fallback
    return getter(fallback)



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


def color_tris(surf, color_scheme, color_map, color_factor, max_val=None, min_val=0, inverse=False, remove_outliers=True):
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
            - 'gauss': Gaussian curvature-based coloring
            - 'avg_mean': Average mean curvature
            - 'avg_gauss': Average Gaussian curvature
            - 'max_mean': Maximum mean curvature
            - 'max_gauss': Maximum Gaussian curvature
        color_map (str): Name of the matplotlib colormap to use for coloring
        color_factor (str): Transformation to apply to values before coloring. Options:
            - 'log': Logarithmic transformation
            - 'sqr': Square transformation
            - 'cub': Cubic transformation
            - 'sqrt': Square root transformation
            - 'lin': Linear transformation
            - None: No transformation
        max_val (float, optional): Maximum value for color normalization. If None, uses maximum
            value from the data.
        min_val (float, optional): Minimum value for color normalization. If None, uses minimum
            value from the data.
        inverse (bool, optional): If True, the color scheme is inverted. Default is False
        remove_outliers (bool, optional): If True, outliers are removed from the color scheme. Default is True
    Returns:
        list: The function modifies the surface data structure in place by adding color information
        to the triangles.

    Notes:
        - The function uses matplotlib colormaps for color generation
        - Color values are stored as RGBA tuples in the range [0,1]
        - The function handles various edge cases including flat surfaces and missing data
        - Color calculations are cached when possible to improve performance
    """
    import matplotlib as mpl
    from matplotlib._api.deprecation import MatplotlibDeprecationWarning as MPLDepWarn
    # Set up the variable tri_colors for recording the color designations for
    tri_colors = None
    # Set up the inverse multiplyer
    inverse_mult = 1
    if inverse:
        inverse_mult = -1
    # Set up the color map
    my_cmap = get_cmap_safe(color_map)
    # Create the surface value multiplier
    if color_factor == 'log':
        def multi(val):
            return np.log(val + 1)
    elif color_factor == 'sqr':
        def multi(val):
            return val ** 2
    elif color_factor == 'cub':
        def multi(val):
            return val ** 3
    elif color_factor == 'sqrt':
        def multi(val):
            return np.sqrt(val)
    elif color_factor == 'lin':
        def multi(val):
            return val
    else:
        def multi(val):
            return val

    def remove_outlrs(tri_curvs):
        arr = np.asarray(tri_curvs, dtype=float)
        if arr.size == 0:
            return arr

        low, high = np.percentile(arr, [10, 90])
        return np.clip(arr, low, high)
    
    # Check if the function is None
    if surf['func'] is None:
        a0, a1 = [surf['net'].balls.iloc[_] for _ in surf['balls']]
        func = calc_surf_func(a0['loc'], a0['rad'], a1['loc'], a1['rad'])
    else:
        func = surf['func']

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
    
    elif color_scheme == 'avg_mean':

        # Check if the tri_dists have been calculated before
        if surf['mean_tri_curvs'] is None or len(surf['mean_tri_curvs']) == 0 or len(surf['mean_tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='mean')
        else:
            tri_curvs = surf['mean_tri_curvs']

        # Check if tri curves are empty
        if len(tri_curvs) == 0:
            tri_colors = [np.random.rand(3) for _ in range(len(surf['tris']))]
            return tri_colors

        # Remove the outliers
        if remove_outliers:
            tri_curvs = remove_outlrs(tri_curvs)

        # Get the average of the curvature values
        avg_curv = sum(tri_curvs) / len(tri_curvs)

        # Put the average color on the scale from 0 to 1
        avg_curv = (avg_curv-min_val)/(max_val-min_val)

        # Map the average curvature to a color
        avg_curv = my_cmap(multi(avg_curv))

        # Set all the colors to the mean of the tri_curvs
        tri_colors = [avg_curv] * len(surf['tris'])

    elif color_scheme == 'avg_gauss':

        # Check if the tri_dists have been calculated before
        if surf['gauss_tri_curvs'] is None or len(surf['gauss_tri_curvs']) == 0 or len(surf['gauss_tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='gauss')
        else:
            tri_curvs = surf['gauss_tri_curvs']

        # Remove the outliers
        if remove_outliers:
            tri_curvs = remove_outlrs(tri_curvs)

        # Get the average of the curvature values
        avg_curv = sum(tri_curvs) / len(tri_curvs)

        # Put the average color on the scale from 0 to 1
        avg_curv = (avg_curv-min_val)/(max_val-min_val)

        # Map the average curvature to a color
        avg_curv = my_cmap(multi(avg_curv))

        # Set all the colors to the mean of the tri_curvs
        tri_colors = [avg_curv] * len(surf['tris'])
    
    elif color_scheme == 'max_mean':

        # Check if the tri_dists have been calculated before
        if surf['mean_tri_curvs'] is None or len(surf['mean_tri_curvs']) == 0 or len(surf['mean_tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='mean')
        else:
            tri_curvs = surf['mean_tri_curvs']
        
        # Remove the outliers
        if remove_outliers:
            tri_curvs = remove_outlrs(tri_curvs)

        # Set the tri curves to a range from 0 to 1
        my_curvs = [(inverse_mult*curv-min_val)/(max_val-min_val) for curv in tri_curvs]

        # Check if we are taking the smallest or the largest
        if inverse:
            max_curve = min(my_curvs)
        else:
            max_curve = max(my_curvs)

        # Map the maximum curvature to a color
        max_curve = my_cmap(multi(max_curve))

        # Set all the colors to the maximum of the tri_curvs
        tri_colors = [max_curve] * len(surf['tris'])

    elif color_scheme == 'max_gauss':

        # Check if the tri_dists have been calculated before
        if surf['gauss_tri_curvs'] is None or len(surf['gauss_tri_curvs']) == 0 or len(surf['gauss_tri_curvs']) != len(surf['tris']):
            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='gauss')
        else:
            tri_curvs = surf['gauss_tri_curvs']

        # Remove the outliers
        if remove_outliers:
            tri_curvs = remove_outlrs(tri_curvs)
        
        # Set the tri curves to a range from 0 to 1
        my_curvs = [(inverse_mult*curv-min_val)/(max_val-min_val) for curv in tri_curvs]

        # Check if we are taking the smallest or the largest
        if inverse:
            max_curve = min(my_curvs)
        else:
            max_curve = max(my_curvs)

        # Map the maximum curvature to a color
        max_curve = my_cmap(multi(max_curve))

        # Set all the colors to the maximum of the tri_curvs
        tri_colors = [max_curve] * len(surf['tris'])


    elif color_scheme.lower() in {'mean', 'mean_curv', 'mean curvature'}:

        if surf['mean_tri_curvs'] is None or len(surf['mean_tri_curvs']) == 0 or len(surf['mean_tri_curvs']) != len(
                surf['tris']):

            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='mean')

        else:

            tri_curvs = surf['mean_tri_curvs']

        if len(tri_curvs) == 0:
            return [np.random.rand(3) for _ in range(len(surf['tris']))]

        tri_curvs = np.asarray(tri_curvs, dtype=float)

        if remove_outliers:
            tri_curvs = remove_outlrs(tri_curvs)

        if surf['flat'] or surf['mean_curv'] == 0:

            my_curvs = np.full(len(surf['tris']), (0 - min_val) / (max_val - min_val))

        else:

            my_curvs = (inverse_mult * tri_curvs - min_val) / (max_val - min_val)

        tri_colors = my_cmap(multi(my_curvs))


    elif color_scheme.lower() in {'gauss', 'gauss_curv', 'gaussian curvature'}:

        if surf['gauss_tri_curvs'] is None or len(surf['gauss_tri_curvs']) == 0 or len(surf['gauss_tri_curvs']) != len(
                surf['tris']):

            tri_curvs, _ = calc_surf_tri_curvs(func, surf['points'], surf['tris'], curvature_type='gauss')

        else:

            tri_curvs = surf['gauss_tri_curvs']

        if len(tri_curvs) == 0:
            return [np.random.rand(3) for _ in range(len(surf['tris']))]

        tri_curvs = np.asarray(tri_curvs, dtype=float)

        if remove_outliers:
            tri_curvs = remove_outlrs(tri_curvs)

        if surf['flat'] or surf['mean_curv'] == 0:

            my_curvs = np.full(len(surf['tris']), (0 - min_val) / (max_val - min_val))

        else:

            my_curvs = (inverse_mult * tri_curvs - min_val) / (max_val - min_val)

        tri_colors = my_cmap(multi(my_curvs))
    # ------------------------------------------------------------------
    # Integrated mean curvature
    # ------------------------------------------------------------------

    elif color_scheme.lower() in {

        'int_mean_curv',

        'int_mean',

        'integrated_mean',

        'integrated mean curvature'

    }:

        value = float(surf['int_mean_curv'])

        norm_value = (

                (inverse_mult * value - min_val)

                / (max_val - min_val)

        )

        norm_value = np.clip(norm_value, 0.0, 1.0)

        surf_color = my_cmap(multi(norm_value))

        tri_colors = [

            surf_color

            for _ in range(len(surf['tris']))

        ]


    # ------------------------------------------------------------------
    # Integrated mean curvature squared
    # ------------------------------------------------------------------

    elif color_scheme.lower() in {
        'int_mean_curv_sq',
        'int_mean_sq',
        'integrated_mean_squared',
        'integrated mean curvature squared'
    }:

        value = float(surf['int_mean_curv_sq'])

        norm_value = (
                (value - min_val)
                / (max_val - min_val)
        )

        norm_value = np.clip(norm_value, 0.0, 1.0)

        surf_color = my_cmap(multi(norm_value))

        tri_colors = [
            surf_color
            for _ in range(len(surf['tris']))
        ]


    # ------------------------------------------------------------------
    # Integrated Gaussian curvature
    # ------------------------------------------------------------------

    elif color_scheme.lower() in {

        'int_gauss_curv',

        'int_gauss',

        'integrated_gauss',

        'integrated gaussian curvature'

    }:
        print("Coloring tris")

        value = float(surf['int_gauss_curv'])

        norm_value = (

                (value - min_val)

                / (max_val - min_val)

        )

        norm_value = np.clip(norm_value, 0.0, 1.0)

        surf_color = my_cmap(multi(norm_value))

        tri_colors = [

            surf_color

            for _ in range(len(surf['tris']))

        ]
    elif color_scheme.lower() == "black":
        tri_colors = [(0, 0, 0) for _ in range(len(surf['tris']))]
    else:
        color = np.random.rand(3)
        tri_colors = [color for _ in range(len(surf['tris']))]

    return tri_colors
