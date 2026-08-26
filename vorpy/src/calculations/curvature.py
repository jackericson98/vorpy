import numpy as np
from vorpy.src.calculations.calcs import calc_tri


def gaussian_curvature(func, point, tol=1e-12):
    """
    Calculate the Gaussian curvature at a point on an implicit quadratic surface.

    For an implicit surface F(x, y, z) = 0,

        K = (grad(F)^T adj(Hess(F)) grad(F)) / |grad(F)|^4

    where adj(Hess(F)) is the adjugate of the Hessian matrix.

    Parameters
    ----------
    func : list
        Coefficients defining the quadratic surface equation.
    point : numpy.ndarray
        Point coordinates [x, y, z] where the curvature is calculated.
    tol : float, optional
        Numerical tolerance used to detect degenerate gradients.

    Returns
    -------
    float
        Gaussian curvature at the specified point.
    """
    A, B, C, D, E, F, G, Hc, Ic, J, Kc, dx, dy, dz = func
    x, y, z = point

    # First derivatives / gradient
    fx = 2 * A * x + D * y + F * z + G
    fy = 2 * B * y + D * x + E * z + Hc
    fz = 2 * C * z + F * x + E * y + Ic

    grad = np.array([fx, fy, fz], dtype=float)
    grad_mag = np.linalg.norm(grad)

    # Guard against degenerate / ill-posed points
    if not np.isfinite(grad_mag) or grad_mag < tol:
        return 0.0

    # Hessian matrix
    hess = np.array([
        [2 * A, D, F],
        [D, 2 * B, E],
        [F, E, 2 * C]
    ], dtype=float)

    # Explicit adjugate. This avoids matrix inversion, which would fail
    # for valid surfaces having a singular Hessian.
    a, b, c = hess[0]
    d, e, f = hess[1]
    g, h, i = hess[2]

    adj_hess = np.array([
        [e * i - f * h, c * h - b * i, b * f - c * e],
        [f * g - d * i, a * i - c * g, c * d - a * f],
        [d * h - e * g, b * g - a * h, a * e - b * d]
    ], dtype=float)

    numerator = grad @ adj_hess @ grad
    denominator = grad_mag ** 4

    if not np.isfinite(denominator) or denominator < tol:
        return 0.0

    K = numerator / denominator

    if not np.isfinite(K):
        return 0.0

    return K


def mean_curvature(func, point, tol=1e-12):
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
    A, B, C, D, E, F, G, Hc, Ic, J, Kc, dx, dy, dz = func
    x, y, z = point

    # First derivatives
    fx = 2 * A * x + D * y + F * z + G
    fy = 2 * B * y + D * x + E * z + Hc
    fz = 2 * C * z + F * x + E * y + Ic

    # Second derivatives
    fxx = 2 * A
    fyy = 2 * B
    fzz = 2 * C
    fxy = D
    fxz = F
    fyz = E

    # Gradient magnitude
    grad_mag = np.sqrt(fx ** 2 + fy ** 2 + fz ** 2)

    # Guard against degenerate / ill-posed points
    if not np.isfinite(grad_mag) or grad_mag < tol:
        return 0.0

    H_mat = np.array([[fxx, fxy, fxz],
                      [fxy, fyy, fyz],
                      [fxz, fyz, fzz]])

    # Mean curvature formula
    grad_vec = np.array([fx, fy, fz])
    num = np.trace(H_mat) * grad_mag ** 2 - grad_vec @ (H_mat @ grad_vec)
    denom = 2.0 * grad_mag ** 3

    # denom is protected by the tol-check above, but keep it explicit
    if not np.isfinite(denom) or abs(denom) < tol:
        return 0.0

    H_val = num / denom
    # Optional: clamp crazy values
    if not np.isfinite(H_val):
        return 0.0

    return H_val


def calc_avg_curv(points, tris, curvs):
    """
    Calculates the average curvature of the surface

    This function computes the average surface curvature based on the sa of each triangle and the total surface area of
    the surface to compute a standardized value.

    Parameters
    ----------
    points: list
        List of points associates with the vertices of the triangles of a surface.
    tris: list
        List of triangles; each triangle is three integer indices into the points.
    curvs: list
        List of points that correspond to the curvature on the surface at the centroid of the triangle

    Returns
    -------
    float
        The average curvature of the entire surface
    """

    sa = 0.0
    tot_curv = 0.0

    if len(tris) != len(curvs):
        raise ValueError("tris and curvs must have the same length.")

    for i, tri in enumerate(tris):
        p0, p1, p2 = [points[_] for _ in tri]

        tri_sa = calc_tri(p0, p1, p2)

        tot_curv += tri_sa * curvs[i]
        sa += tri_sa

    if sa == 0:
        return 0.0  # or np.nan, depending on what makes more sense

    return tot_curv / sa


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
    max_tcs = max(map(abs, tri_curvs)) if tri_curvs else None

    # Calculate the average of the tri curves
    avg_curv = calc_avg_curv(points, tris, tri_curvs)

    return tri_curvs, max_tcs, avg_curv



def calc_surf_tri_curvs_both(func, points, tris, tol=1e-12):
    """
    Calculate mean and Gaussian curvature together for every surface triangle.

    This version preserves the existing formulas while:
    - computing each centroid once,
    - computing the gradient once,
    - using scalar arithmetic instead of tiny NumPy matrix operations,
    - calculating triangle area once,
    - accumulating both area-weighted averages in the same triangle pass.

    Returns
    -------
    tuple
        mean_tri_curvs, mean_max_abs, mean_avg,
        gauss_tri_curvs, gauss_max_abs, gauss_avg
    """
    A, B, C, D, E, F, G, Hc, Ic, J, Kc, dx, dy, dz = func

    # Hessian terms are constant across a quadratic surface.
    h00 = 2.0 * A
    h01 = D
    h02 = F
    h11 = 2.0 * B
    h12 = E
    h22 = 2.0 * C
    hess_trace = h00 + h11 + h22

    # Adjugate of the symmetric Hessian, also constant per surface.
    a00 = h11 * h22 - h12 * h12
    a01 = h02 * h12 - h01 * h22
    a02 = h01 * h12 - h02 * h11
    a11 = h00 * h22 - h02 * h02
    a12 = h01 * h02 - h00 * h12
    a22 = h00 * h11 - h01 * h01

    mean_tri_curvs = []
    gauss_tri_curvs = []

    mean_max = None
    gauss_max = None

    total_area = 0.0
    weighted_mean = 0.0
    weighted_gauss = 0.0

    for tri in tris:
        i0, i1, i2 = tri
        p0 = points[i0]
        p1 = points[i1]
        p2 = points[i2]

        # Triangle centroid.
        x = (p0[0] + p1[0] + p2[0]) / 3.0
        y = (p0[1] + p1[1] + p2[1]) / 3.0
        z = (p0[2] + p1[2] + p2[2]) / 3.0

        # Shared gradient.
        fx = 2.0 * A * x + D * y + F * z + G
        fy = 2.0 * B * y + D * x + E * z + Hc
        fz = 2.0 * C * z + F * x + E * y + Ic

        grad_sq = fx * fx + fy * fy + fz * fz

        if not np.isfinite(grad_sq) or grad_sq < tol * tol:
            mean_val = 0.0
            gauss_val = 0.0
        else:
            grad_mag = np.sqrt(grad_sq)

            # Mean curvature:
            # H = (tr(Hess)*|grad|^2 - grad^T Hess grad) / (2|grad|^3)
            hg0 = h00 * fx + h01 * fy + h02 * fz
            hg1 = h01 * fx + h11 * fy + h12 * fz
            hg2 = h02 * fx + h12 * fy + h22 * fz
            grad_hess_grad = fx * hg0 + fy * hg1 + fz * hg2

            mean_num = hess_trace * grad_sq - grad_hess_grad
            mean_denom = 2.0 * grad_sq * grad_mag

            if not np.isfinite(mean_denom) or abs(mean_denom) < tol:
                mean_val = 0.0
            else:
                mean_val = mean_num / mean_denom
                if not np.isfinite(mean_val):
                    mean_val = 0.0

            # Gaussian curvature:
            # K = grad^T adj(Hess) grad / |grad|^4
            ag0 = a00 * fx + a01 * fy + a02 * fz
            ag1 = a01 * fx + a11 * fy + a12 * fz
            ag2 = a02 * fx + a12 * fy + a22 * fz
            gauss_num = fx * ag0 + fy * ag1 + fz * ag2
            gauss_denom = grad_sq * grad_sq

            if not np.isfinite(gauss_denom) or gauss_denom < tol:
                gauss_val = 0.0
            else:
                gauss_val = gauss_num / gauss_denom
                if not np.isfinite(gauss_val):
                    gauss_val = 0.0

        mean_tri_curvs.append(mean_val)
        gauss_tri_curvs.append(gauss_val)

        abs_mean = abs(mean_val)
        abs_gauss = abs(gauss_val)
        if mean_max is None or abs_mean > mean_max:
            mean_max = abs_mean
        if gauss_max is None or abs_gauss > gauss_max:
            gauss_max = abs_gauss

        # Triangle area, computed once and shared by both averages.
        ax = p1[0] - p0[0]
        ay = p1[1] - p0[1]
        az = p1[2] - p0[2]

        bx = p2[0] - p0[0]
        by = p2[1] - p0[1]
        bz = p2[2] - p0[2]

        cx = ay * bz - az * by
        cy = az * bx - ax * bz
        cz = ax * by - ay * bx

        tri_area = 0.5 * np.sqrt(cx * cx + cy * cy + cz * cz)

        if np.isfinite(tri_area) and tri_area > 0.0:
            total_area += tri_area
            weighted_mean += mean_val * tri_area
            weighted_gauss += gauss_val * tri_area

    if total_area > 0.0:
        mean_avg = weighted_mean / total_area
        gauss_avg = weighted_gauss / total_area
    else:
        mean_avg = 0.0
        gauss_avg = 0.0

    return (
        mean_tri_curvs, mean_max, mean_avg,
        gauss_tri_curvs, gauss_max, gauss_avg,
    )

def calc_avg_surface_curvature(func, s_points, s_tris, curvature_type='gauss'):
    """
    Compute the area-weighted average curvature over a triangulated surface.

    Parameters
    ----------
    func : list
        Coefficients defining the quadratic surface (passed directly to
        gaussian_curvature / mean_curvature).
    s_points : sequence of array_like, shape (3,)
        List/array of 3D points (x, y, z) defining the vertices of the surface.
    s_tris : sequence of array_like or tuple, length 3
        List of triangles; each triangle is three integer indices into s_points.
    curvature_type : {'gauss', 'mean'}, optional
        Which curvature to use at each triangle centroid.

    Returns
    -------
    float
        Area-weighted average curvature over the surface.
    """

    # Ensure we have numpy arrays for the points
    points = [np.asarray(p, dtype=float) for p in s_points]

    total_area = 0.0
    weighted_sum = 0.0

    for tri in s_tris:
        i, j, k = tri
        p1, p2, p3 = points[i], points[j], points[k]

        # Triangle centroid
        centroid = (p1 + p2 + p3) / 3.0

        # Curvature at centroid
        if curvature_type == 'gauss':
            curv = gaussian_curvature(func, centroid)
        elif curvature_type == 'mean':
            curv = mean_curvature(func, centroid)
        else:
            raise ValueError("curvature_type must be 'gauss' or 'mean'")

        # Triangle area (using your existing calc_tri)
        area = calc_tri(p1, p2, p3)

        total_area += area
        weighted_sum += curv * area

    if total_area == 0.0:
        # Degenerate surface; you can also return np.nan if you prefer
        return 0.0

    # Area-weighted average curvature (Riemann sum → surface integral / area)
    return weighted_sum / total_area