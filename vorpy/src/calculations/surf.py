from numba import jit
import numpy as np
from numba.core.errors import NumbaPendingDeprecationWarning as numba_err
from vorpy.src.calculations.calcs import calc_dist
from vorpy.src.calculations.calcs import calc_tri


@jit(nopython=True)
def calc_tri(points):
    """
    Calculates the area of a triangle formed by three 3D points using the cross product method.

    This function computes the area of a triangle by:
    1. Creating two vectors from the points
    2. Taking their cross product to get a vector perpendicular to the triangle
    3. Taking half the magnitude of this vector

    Parameters:
    -----------
    points : list of array-like
        List containing three vertices of the triangle, each as [x, y, z] coordinates

    Returns:
    --------
    float
        The area of the triangle formed by the three input points

    Notes:
    ------
    - The points should be provided in any order
    - The result is always positive
    - Uses the cross product formula: Area = 0.5 * |AB × AC|
    """
    # Get the two triangles vectors
    ab = [points[0][0] - points[1][0], points[0][1] - points[1][1], points[0][2] - points[1][2]]
    ac = [points[0][0] - points[2][0], points[0][1] - points[2][1], points[0][2] - points[2][2]]

    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm((np.cross(ab, ac)))


def calc_surf_func(l0, r0, l1, r1):
    """
    Calculates the coefficients for the surface between two balls.

    This function computes the mathematical coefficients needed to define the hyperboloid surface
    that represents the boundary between two spheres. The calculation is based on the relative
    positions and radii of the spheres, following the mathematical framework described in Z. Hu's work.

    Parameters:
    -----------
    l0 : numpy.ndarray
        Center coordinates of the first sphere [x, y, z]
    r0 : float
        Radius of the first sphere
    l1 : numpy.ndarray
        Center coordinates of the second sphere [x, y, z]
    r1 : float
        Radius of the second sphere

    Returns:
    --------
    list
        A list containing the coefficients for the hyperboloid equation:
        - ABC: Coefficients for x², y², z² terms
        - DEF: Cross-term coefficients (xy, yz, zx)
        - GHI: Linear coefficients (x, y, z)
        - J: Constant term
        - K: Additional constant term
        - d: Vector between sphere centers
    """
    try:
        vals = calc_surf_func_jit(l0, r0, l1, r1)
    except numba_err:
        vals = calc_surf_func_reg(l0, r0, l1, r1)
    return vals


@jit(nopython=True)
def calc_surf_func_jit(l0, r0, l1, r1):
    """
    Calculates the mathematical coefficients defining the hyperboloid surface between two spheres.

    This function computes the coefficients needed to define the quadratic surface equation
    that represents the boundary between two spheres. The calculation follows the mathematical
    framework described in Z. Hu's work, where the surface is defined by a hyperboloid equation
    of the form Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0.

    Parameters:
    -----------
    l0 : numpy.ndarray
        Center coordinates of the first sphere [x, y, z]
    r0 : float
        Radius of the first sphere
    l1 : numpy.ndarray
        Center coordinates of the second sphere [x, y, z]
    r1 : float
        Radius of the second sphere

    Returns:
    --------
    list
        A list containing the coefficients for the hyperboloid equation:
        - ABC: Coefficients for x², y², z² terms
        - DEF: Cross-term coefficients (xy, yz, zx)
        - GHI: Linear coefficients (x, y, z)
        - J: Constant term
        - K: Additional constant term
        - d: Vector between sphere centers
    """
    # Check the smaller ball is first
    if r1 < r0:
        l0, r0, l1, r1 = l1, r1, l0, r0
    # Grab the centers of the spheres
    x1, y1, z1 = l0
    x2, y2, z2 = l1
    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = r0 - r1
    K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
    d = [x1 - x2, y1 - y2, z1 - z2]
    J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2
    # Instantiate/reset the hyperboloid coefficient vector lists
    ABC, DEF, GHI = [], [], []
    # Calculate hyperboloid coefficients
    for i in range(3):
        ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
        DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
        GHI.append(-8 * R ** 2 * l0[i] - 4 * K * d[i])
    # Return the function coefficients
    return ABC + DEF + GHI + [J] + [K] + d


def calc_surf_func_reg(l0, r0, l1, r1):
    """
    Calculates the mathematical coefficients defining the hyperboloid surface between two spheres.

    This function computes the coefficients needed to define the quadratic surface equation
    that represents the boundary between two spheres. The calculation follows the mathematical
    framework described in Z. Hu's work, where the surface is defined by a hyperboloid equation
    of the form Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0.

    Parameters:
    -----------
    l0 : array-like
        Center coordinates of the first sphere [x, y, z]
    r0 : float
        Radius of the first sphere
    l1 : array-like
        Center coordinates of the second sphere [x, y, z]
    r1 : float
        Radius of the second sphere

    Returns:
    --------
    list
        A list containing the coefficients for the hyperboloid equation:
        - ABC: Coefficients for x², y², z² terms
        - DEF: Cross-term coefficients (xy, yz, zx)
        - GHI: Linear coefficients (x, y, z)
        - J: Constant term
        - K: Additional constant term
        - d: Vector between sphere centers
    """
    # Check the smaller ball is first
    if r1 < r0:
        l0, r0, l1, r1 = l1, r1, l0, r0
    # Grab the centers of the spheres
    x1, y1, z1 = l0
    x2, y2, z2 = l1
    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = r0 - r1
    K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
    d = [x1 - x2, y1 - y2, z1 - z2]
    J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2
    # Instantiate/reset the hyperboloid coefficient vector lists
    ABC, DEF, GHI = [], [], []
    # Calculate hyperboloid coefficients
    for i in range(3):
        ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
        DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
        GHI.append(-8 * R ** 2 * l0[i] - 4 * K * d[i])
    # Return the function coefficients
    return ABC + DEF + GHI + [J] + [K] + d


def calc_2d_surf_sa(tris, points):
    """
    Calculates the surface area of a 2D surface defined by triangles and points.

    This function computes the total surface area of a 2D surface by summing the areas
    of individual triangles that make up the surface. The area calculation uses the
    shoelace formula (also known as the surveyor's formula) for computing the area of
    a triangle given its vertices in 2D space.

    Parameters:
    -----------
    tris : list of tuples
        List of triangles, where each triangle is represented as a tuple of three indices
        corresponding to points in the points array
    points : list of numpy.ndarray
        List of 2D point coordinates [x, y] that form the vertices of the triangles

    Returns:
    --------
    float
        The total surface area of the 2D surface
    """
    # Set up the sa variable
    sa = 0
    # Loop through the triangles
    for tri in tris:
        x1, x2, x3 = [points[_][0] for _ in tri]
        y1, y2, y3 = [points[_][1] for _ in tri]

        sa += 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    return sa


def calc_surf_sa(tris, points):
    """
    Calculates the surface area of a 3D surface defined by triangles and points.

    This function computes the total surface area of a 3D surface by summing the areas
    of individual triangles that make up the surface. The area calculation uses the
    cross product method to compute the area of each triangle in 3D space.

    Parameters:
    -----------
    tris : list of tuples
        List of triangles, where each triangle is represented as a tuple of three indices
        corresponding to points in the points array
    points : list of numpy.ndarray
        List of 3D point coordinates [x, y, z] that form the vertices of the triangles

    Returns:
    --------
    float
        The total surface area of the 3D surface
    """
    # Create the surface area variable
    sa = 0
    # Go through the triangles in the surface
    for tri in tris:
        tri1 = np.array([points[tri[_]] for _ in range(3)])
        sa += calc_tri(tri1)
    # Return the surface area
    return sa


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
        if my_dist > max_dist:
            max_dist = my_dist
        if my_dist < min_dist:
            min_dist = my_dist
    # Normalize the distances
    for dist in dists:
        tri_dists.append((dist - min_dist) / (max_dist - min_dist))
    return tri_dists
