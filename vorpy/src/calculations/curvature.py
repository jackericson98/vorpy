import numpy as np
from vorpy.src.calculations.calcs import calc_dist


def gaussian_curvature(func, point):
    """
    Calculates the Gaussian curvature at a point on a surface.

    This function computes the Gaussian curvature at a given point on a surface defined
    by a quadratic function. The Gaussian curvature is a measure of the intrinsic curvature
    of the surface at that point.

    Parameters
    ----------
    func : list
        List of coefficients defining the quadratic surface equation
    point : numpy.ndarray
        Point coordinates [x, y, z] where the curvature is to be calculated

    Returns
    -------
    float
        The Gaussian curvature at the specified point

    Examples
    --------
    >>> func = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Example coefficients
    >>> point = np.array([0, 0, 0])
    >>> K = gaussian_curvature(func, point)
    """
    # Unpack the function coefficients
    A, B, C, D, E, F, G, H, I, J, K, dx, dy, dz = func
    x, y, z = point

    # Calculate first derivatives
    fx = 2*A*x + D*y + F*z + G
    fy = 2*B*y + D*x + E*z + H
    fz = 2*C*z + F*x + E*y + I

    # Calculate second derivatives
    fxx = 2*A
    fyy = 2*B
    fzz = 2*C
    fxy = D
    fxz = F
    fyz = E

    # Calculate the gradient magnitude
    grad_mag = np.sqrt(fx**2 + fy**2 + fz**2)

    # Calculate the Hessian matrix
    H = np.array([[fxx, fxy, fxz],
                  [fxy, fyy, fyz],
                  [fxz, fyz, fzz]])

    # Calculate the Gaussian curvature
    K = np.linalg.det(H) / (grad_mag**4)

    return K


def mean_curvature(func, point):
    """
    Calculates the mean curvature at a point on a surface.

    This function computes the mean curvature at a given point on a surface defined
    by a quadratic function. The mean curvature is a measure of the extrinsic curvature
    of the surface at that point.

    Parameters
    ----------
    func : list
        List of coefficients defining the quadratic surface equation
    point : numpy.ndarray
        Point coordinates [x, y, z] where the curvature is to be calculated

    Returns
    -------
    float
        The mean curvature at the specified point

    Examples
    --------
    >>> func = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Example coefficients
    >>> point = np.array([0, 0, 0])
    >>> H = mean_curvature(func, point)
    """
    # Unpack the function coefficients
    A, B, C, D, E, F, G, H, I, J, K, dx, dy, dz = func
    x, y, z = point

    # Calculate first derivatives
    fx = 2*A*x + D*y + F*z + G
    fy = 2*B*y + D*x + E*z + H
    fz = 2*C*z + F*x + E*y + I

    # Calculate second derivatives
    fxx = 2*A
    fyy = 2*B
    fzz = 2*C
    fxy = D
    fxz = F
    fyz = E

    # Calculate the gradient magnitude
    grad_mag = np.sqrt(fx**2 + fy**2 + fz**2)

    # Calculate the Hessian matrix
    H = np.array([[fxx, fxy, fxz],
                  [fxy, fyy, fyz],
                  [fxz, fyz, fzz]])

    # Calculate the mean curvature
    H = (np.trace(H) * grad_mag**2 - np.dot(np.array([fx, fy, fz]), np.dot(H, np.array([fx, fy, fz])))) / (2 * grad_mag**3)
    return H


def calc_surf_tri_curvs(func, points, tris, curvature_type='gauss'):
    """
    Calculates the curvature values for each triangle in a surface.

    This function computes either Gaussian or mean curvature values for each triangle
    in a surface by evaluating the curvature at the triangle's centroid.

    Parameters
    ----------
    func : list
        List of coefficients defining the quadratic surface equation
    points : list of numpy.ndarray
        List of 3D point coordinates [x, y, z] that form the vertices of the triangles
    tris : list of tuples
        List of triangles, where each triangle is represented as a tuple of three indices
        corresponding to points in the points array
    curvature_type : {'gauss', 'mean'}, optional
        Type of curvature to calculate. Default is 'gauss'.

    Returns
    -------
    tuple
        A tuple containing:
        - List of curvature values for each triangle
        - Maximum curvature value

    Examples
    --------
    >>> func = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    >>> points = [np.array([0, 0, 0]), np.array([1, 0, 0]), np.array([0, 1, 0])]
    >>> tris = [(0, 1, 2)]
    >>> curvs, max_curv = calc_surf_tri_curvs(func, points, tris)
    """
    # Initialize lists to store curvatures and centroids
    tri_curvs = []
    tri_centroids = []

    # Calculate curvature for each triangle
    for tri in tris:
        # Get the triangle vertices
        v1, v2, v3 = [points[i] for i in tri]
        
        # Calculate the centroid
        centroid = (v1 + v2 + v3) / 3
        tri_centroids.append(centroid)
        
        # Calculate the curvature at the centroid
        if curvature_type == 'gauss':
            curv = gaussian_curvature(func, centroid)
        else:  # mean curvature
            curv = mean_curvature(func, centroid)
            
        tri_curvs.append(curv)
    
    # Calculate the max of the tri curves if it isn't empty
    max_tcs = max(tri_curvs) if tri_curvs else None

    return tri_curvs, max_tcs
