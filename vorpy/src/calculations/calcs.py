import numpy as np
import warnings
from numba import jit
import math
from numba.core.errors import TypingError

try:
    import importlib

    _native_calc = importlib.import_module("vorpy._native._calc")
except Exception:
    _native_calc = None

warnings.filterwarnings("error")


def round_func(round_to):
    """
    Creates a configurable rounding function that can handle both single values and iterables.

    This function returns a closure that maintains the specified rounding precision and can be
    reused for consistent rounding across multiple values. The returned function handles both
    single numeric values and iterables of values, applying the same rounding precision to all.

    Parameters
    ----------
    round_to : int
        The number of decimal places to round to. A positive value rounds to that many decimal
        places, while a negative value rounds to the left of the decimal point.

    Returns
    -------
    function
        A closure that takes a value (or iterable) and optionally a new rounding precision,
        returning the rounded value(s) with the specified precision.
    """

    # Define the inner round function
    def round_(val, new_num=None):
        """
        Inner round function operating on outer defined round to value
        :param val: float/iterable - val(s) to be rounded
        :param new_num: New round to value
        :return: float/list - rounded values
        """
        # Set the new round to number if specified
        if new_num is None:
            new_num = round_to
        # Return the values
        try:
            return round(val, new_num)
        except TypeError:
            return [round(_, new_num) for _ in val]

    # Return the function for the outer function
    return round_


def calc_dist_seq_py(l0, l1):
    if len(l0) != len(l1):
        raise ValueError("calc_dist_seq_py expects same-length inputs")

    acc = 0.0

    for a, b in zip(l0, l1):
        d = float(a) - float(b)
        acc += d * d

    return math.sqrt(acc)


def calc_dist_seq(l0, l1):
    """
    Calculate Euclidean distance between two 1D sequences (lists/tuples/arrays).

    This is a safe fallback for when native acceleration isn't available.
    Accepts array-like inputs; coerces to 1D float64 arrays and enforces same length.
    """
    a = np.asarray(l0, dtype=np.float64)
    b = np.asarray(l1, dtype=np.float64)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("calc_dist_seq expects 1D array-like inputs")

    if a.shape[0] != b.shape[0]:
        raise ValueError("calc_dist_seq expects same-length inputs")

    d = a - b

    return float(np.sqrt(np.dot(d, d)))


def calc_dist(l0, l1):
    """Calculate the Euclidean distance between two points in n-dimensional space.

    Parameters
    ----------
    l0 : array-like
        First point coordinates as an n-dimensional array or list
    l1 : array-like
        Second point coordinates as an n-dimensional array or list with same dimensionality as l0

    Returns
    -------
    float
        The Euclidean distance between the two points

    Examples
    """
    use_native = (_native_calc is not None) and hasattr(_native_calc, "calc_dist")

    if use_native:
        return float(_native_calc.calc_dist(l0, l1))

    return calc_dist_seq(l0, l1)


@jit(nopython=True)
def calc_dist_numba(l0, l1):
    """Calculate the Euclidean distance between two points in n-dimensional space.

    This function computes the straight-line distance between two points using the
    Pythagorean theorem generalized to n dimensions. The function is optimized with
    Numba's JIT compilation for improved performance.

    Parameters
    ----------
    l0 : numpy.ndarray
        First point coordinates as an n-dimensional array
    l1 : numpy.ndarray
        Second point coordinates as an n-dimensional array with same dimensionality as l0

    Returns
    -------
    float
        The Euclidean distance between the two points

    Notes
    -----
    - Both input points must have the same dimensionality
    - Uses numpy's square and sqrt functions for efficient computation
    - JIT compiled for performance optimization
    """
    # Pythagorean theorem
    return np.sqrt(sum(np.square(l0 - l1)))


@jit(nopython=True)
def calc_angle_jit(p0, p1, p2=None):
    """Calculates the angle between three points in radians using vector geometry.

    This function computes the angle between vectors formed by the points in two possible ways:
    1. If p2 is not provided: Angle between vectors from origin to p0 and p1
    2. If p2 is provided: Angle between vectors from p0 to p1 and p0 to p2

    Parameters
    ----------
    p0 : array-like
        First point coordinates [x, y, z, ...]
    p1 : array-like
        Second point coordinates [x, y, z, ...]
    p2 : array-like, optional
        Third point coordinates [x, y, z, ...]
        If not provided, the origin (0,0,0) is used as the reference point

    Returns
    -------
    float
        The angle in radians between the vectors formed by the points

    Notes
    -----
    - All points must have the same dimensionality
    - Uses numpy's arccos function for angle calculation
    - Handles edge cases where vectors are parallel or antiparallel
    """
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = p0, p1
    else:
        v0, v1 = p1 - p0, p2 - p0

    # Check for zero-length vectors
    norm_v0 = np.linalg.norm(v0)
    norm_v1 = np.linalg.norm(v1)

    if norm_v0 == 0.0 or norm_v1 == 0.0:
        return 0.0  # Return 0 for degenerate cases

    n0, n1 = v0 / norm_v0, v1 / norm_v1
    # Calculate the angle between the two vectors with catches for 180 and 0
    my_dot = np.dot(n0, n1)
    if my_dot <= -1.0:
        my_dot = -1.0
    elif my_dot >= 1.0:
        my_dot = 1.0
    angle = np.arccos(my_dot)
    return angle


def calc_angle(p0, p1, p2=None):
    """Calculate the angle between three points in radians.

    This function computes the angle between vectors formed by the points in two possible ways:
    1. If p2 is not provided: Angle between vectors from origin to p0 and p1
    2. If p2 is provided: Angle between vectors from p0 to p1 and p0 to p2

    Parameters
    ----------
    p0 : array-like
        First point coordinates [x, y, z, ...]
    p1 : array-like
        Second point coordinates [x, y, z, ...]
    p2 : array-like, optional
        Third point coordinates [x, y, z, ...]
        If not provided, the origin (0,0,0) is used as the reference point

    Returns
    -------
    float
        The angle in radians between the vectors formed by the points

    Examples
    --------
    >>> import numpy as np
    >>> p0 = np.array([1, 0, 0])
    >>> p1 = np.array([0, 1, 0])
    >>> calc_angle(p0, p1)  # Angle between x and y axes
    1.5707963267948966
    """
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = p0, p1
    else:
        v0, v1 = p1 - p0, p2 - p0

    # Check for zero-length vectors
    norm_v0 = np.linalg.norm(v0)
    norm_v1 = np.linalg.norm(v1)

    if norm_v0 == 0.0 or norm_v1 == 0.0:
        return np.nan  # Return NaN for degenerate cases

    n0, n1 = v0 / norm_v0, v1 / norm_v1
    # Calculate the angle between the two vectors with catches for 180 and 0
    my_dot = np.dot(n0, n1)
    if my_dot <= -1.0:
        my_dot = -1.0
    elif my_dot >= 1.0:
        my_dot = 1.0
    angle = np.arccos(my_dot)
    return angle


@jit(nopython=True)
def calc_tetra_vol(p0, p1, p2, p3):
    """Calculate the volume of a tetrahedron defined by four vertices in 3D space.

    This function uses the scalar triple product formula to compute the volume:
    V = (1/6) * |(p3-p0) · ((p1-p0) × (p2-p0))|
    where · denotes the dot product and × denotes the cross product.

    Parameters
    ----------
    p0 : array-like
        First vertex of the tetrahedron [x, y, z]
    p1 : array-like
        Second vertex of the tetrahedron [x, y, z]
    p2 : array-like
        Third vertex of the tetrahedron [x, y, z]
    p3 : array-like
        Fourth vertex of the tetrahedron [x, y, z]

    Returns
    -------
    float
        The volume of the tetrahedron formed by the four vertices

    Examples
    --------
    >>> import numpy as np
    >>> p0 = np.array([0, 0, 0])
    >>> p1 = np.array([1, 0, 0])
    >>> p2 = np.array([0, 1, 0])
    >>> p3 = np.array([0, 0, 1])
    >>> calc_tetra_vol(p0, p1, p2, p3)  # Volume of unit tetrahedron
    0.16666666666666666
    """
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = np.array([p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]])

    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    return (1 / 6) * abs(np.dot(r03, np.cross(r01, r02)))


def calc_tetra_inertia(ps, mass):
    """Calculate the moment of inertia tensor of a tetrahedron about its centroid.

    This function computes the inertia tensor for a tetrahedron with uniform density distribution.
    The calculation is based on the parallel axis theorem and the inertia tensor of a tetrahedron
    about its centroid.

    Parameters
    ----------
    ps : list of array-like
        List containing four vertices of the tetrahedron, each as [x, y, z] coordinates
    mass : float
        Total mass of the tetrahedron

    Returns
    -------
    numpy.ndarray
        A 3x3 inertia tensor matrix where:
        - Diagonal elements represent moments of inertia about x, y, and z axes
        - Off-diagonal elements represent products of inertia

    Examples
    --------
    >>> import numpy as np
    >>> ps = [
    ...     np.array([0, 0, 0]),
    ...     np.array([1, 0, 0]),
    ...     np.array([0, 1, 0]),
    ...     np.array([0, 0, 1])
    ... ]
    >>> mass = 1.0
    >>> calc_tetra_inertia(ps, mass)
    array([[ 0.1, -0.05, -0.05],
           [-0.05,  0.1, -0.05],
           [-0.05, -0.05,  0.1]])
    """
    # Placeholder for inertia tensor calculation.
    # For simplicity, this uses an approximate inertia formula for a solid tetrahedron.
    # More accurate calculations can be done by integrating over the volume.
    inertia_tensor = np.zeros((3, 3))

    # Sum contributions from the vertices
    for i in range(4):
        x, y, z = ps[i]
        inertia_tensor[0, 0] += mass * (y ** 2 + z ** 2) / 10.0
        inertia_tensor[1, 1] += mass * (x ** 2 + z ** 2) / 10.0
        inertia_tensor[2, 2] += mass * (x ** 2 + y ** 2) / 10.0
        inertia_tensor[0, 1] -= mass * x * y / 10.0
        inertia_tensor[0, 2] -= mass * x * z / 10.0
        inertia_tensor[1, 2] -= mass * y * z / 10.0

    # Symmetric tensor: fill in the other values
    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    return inertia_tensor


@jit(nopython=True)
def calc_tri_py(points):
    """Calculate the area of a triangle formed by three 3D points.

    This function computes the area of a triangle by:
    1. Creating two vectors from the points
    2. Taking their cross product to get a vector perpendicular to the triangle
    3. Taking half the magnitude of this vector

    Parameters
    ----------
    points : list of array-like
        List containing three vertices of the triangle, each as [x, y, z] coordinates

    Returns
    -------
    float
        The area of the triangle formed by the three input points

    Notes
    -----
    - The points should be provided in any order
    - The result is always positive
    - Uses the cross product formula: Area = 0.5 * |AB × AC|

    Examples
    --------
    >>> points = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    >>> calc_tri(points)
    0.5
    """
    # Get the two triangles vectors
    ab = [points[0][0] - points[1][0], points[0][1] - points[1][1], points[0][2] - points[1][2]]
    ac = [points[0][0] - points[2][0], points[0][1] - points[2][1], points[0][2] - points[2][2]]

    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm(np.cross(ab, ac))


def calc_tri(p0, p1=None, p2=None):
    # single-arg (3,3) form
    if p1 is None and p2 is None:
        tri = np.asarray(p0, dtype=np.float64)
        if _native_calc is not None and hasattr(_native_calc, "calc_tri"):
            return float(_native_calc.calc_tri(tri))

        # python fallback
        e1 = tri[1] - tri[0]
        e2 = tri[2] - tri[0]
        return 0.5 * float(np.linalg.norm(np.cross(e1, e2)))

    # three-arg form
    if _native_calc is not None and hasattr(_native_calc, "calc_tri"):
        return float(_native_calc.calc_tri(p0, p1, p2))

    p0 = np.asarray(p0, dtype=np.float64)
    p1 = np.asarray(p1, dtype=np.float64)
    p2 = np.asarray(p2, dtype=np.float64)

    e1 = p1 - p0
    e2 = p2 - p0
    return 0.5 * float(np.linalg.norm(np.cross(e1, e2)))


def calc_com(points, masses=None):
    """Calculate the center of mass for a set of points.

    This function computes the center of mass (centroid) for a collection of points in 3D space.
    If masses are provided, the calculation is weighted by the masses. If no masses are provided,
    all points are assumed to have equal mass.

    Parameters
    ----------
    points : numpy.ndarray
        Array of points, where each point is a [x, y, z] coordinate
    masses : numpy.ndarray, optional
        Array of masses corresponding to each point. If None, all points are assumed
        to have equal mass.

    Returns
    -------
    numpy.ndarray
        The center of mass coordinates [x, y, z]

    Examples
    --------
    >>> import numpy as np
    >>> points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    >>> masses = np.array([1, 2, 1])
    >>> calc_com(points, masses)
    array([0.5, 0.25, 0.])
    """
    if masses is None:
        return np.mean(points, axis=0)
    else:
        return np.average(points, weights=masses, axis=0)


@jit(nopython=True)
def calc_length(points):
    """Calculates the total length of a path defined by a sequence of points.

    This function computes the sum of Euclidean distances between consecutive points in the input sequence.
    The points are assumed to be ordered in the sequence they should be connected.

    Parameters
    ----------
    points : list of array-like
        List of point coordinates in n-dimensional space. Each point should be a list or array
        of coordinates [x1, x2, ..., xn]. Points must be ordered in the sequence they should
        be connected.

    Returns
    -------
    float
        The total length of the path formed by connecting consecutive points in the input sequence.

    Notes
    -----
    - Points must be ordered in the sequence they should be connected
    - Uses Euclidean distance between consecutive points
    - Returns 0 if the input list contains fewer than 2 points
    """
    # Reset the length
    length = 0
    # Go through the points in the list
    for m, point in enumerate(points):
        # Make sure not to index error
        if m + 1 < len(points):
            # Add the length to the total
            length += calc_dist_numba(point, points[m + 1])
    return length


def calc_sphericity(volume, surface_area):
    """Calculate the sphericity of a geometric object based on its volume and surface area.

    Parameters
    ----------
    volume : float
        The volume of the object.
    surface_area : float
        The surface area of the object.

    Returns
    -------
    float
        The sphericity of the object.
    """
    if volume <= 0 or surface_area <= 0:
        raise ValueError("Volume and surface area must be positive numbers.")

    # Calculate sphericity using the geometric formula
    sphericity = (np.pi ** (1 / 3) * (6 * volume) ** (2 / 3)) / surface_area
    return sphericity


def calc_isoperimetric_quotient(volume, surface_area):
    """Calculate the isoperimetric quotient

    Parameters
    ----------
    volume : float
        The volume of the object.
    surface_area : float
        The surface area of the object.

    Returns
    -------
    float
        The isoperimetric quotient of the object
    """
    if volume <= 0 or surface_area <= 0:
        raise ValueError("Volume and surface area must be positive numbers.")

    return (36 * np.pi * volume ** 2) / (surface_area ** 3)


@jit(nopython=True, cache=True)
def _cell_point_properties_surface_kernel(ball_loc, points):
    """Return spike-distance extrema and XYZ bounds for one surface."""
    min_d2 = np.inf
    max_d2 = 0.0

    min_x = np.inf
    min_y = np.inf
    min_z = np.inf
    max_x = -np.inf
    max_y = -np.inf
    max_z = -np.inf

    x0 = ball_loc[0]
    y0 = ball_loc[1]
    z0 = ball_loc[2]

    for i in range(points.shape[0]):
        x = points[i, 0]
        y = points[i, 1]
        z = points[i, 2]

        dx = x - x0
        dy = y - y0
        dz = z - z0
        d2 = dx * dx + dy * dy + dz * dz

        if d2 < min_d2:
            min_d2 = d2
        if d2 > max_d2:
            max_d2 = d2

        if x < min_x:
            min_x = x
        if y < min_y:
            min_y = y
        if z < min_z:
            min_z = z
        if x > max_x:
            max_x = x
        if y > max_y:
            max_y = y
        if z > max_z:
            max_z = z

    return min_d2, max_d2, min_x, min_y, min_z, max_x, max_y, max_z


def calc_cell_point_properties(ball_loc, surfs):
    """
    Calculate spike extrema and the axis-aligned cell bounding box in one pass.

    This is mathematically equivalent to calling ``calc_spikes`` and
    ``calc_cell_box`` separately, but each surface point is classified only once
    and the point loop runs inside a cached Numba kernel.
    """
    ball_loc_arr = np.asarray(ball_loc, dtype=np.float64)

    min_d2 = np.inf
    max_d2 = 0.0
    mins = np.array([np.inf, np.inf, np.inf], dtype=np.float64)
    maxs = np.array([-np.inf, -np.inf, -np.inf], dtype=np.float64)

    found_points = False

    for surf in surfs:
        points = np.asarray(surf['points'], dtype=np.float64)
        if points.size == 0:
            continue
        if points.ndim != 2:
            points = points.reshape((-1, 3))

        found_points = True
        vals = _cell_point_properties_surface_kernel(ball_loc_arr, points)

        if vals[0] < min_d2:
            min_d2 = vals[0]
        if vals[1] > max_d2:
            max_d2 = vals[1]

        if vals[2] < mins[0]:
            mins[0] = vals[2]
        if vals[3] < mins[1]:
            mins[1] = vals[3]
        if vals[4] < mins[2]:
            mins[2] = vals[4]
        if vals[5] > maxs[0]:
            maxs[0] = vals[5]
        if vals[6] > maxs[1]:
            maxs[1] = vals[6]
        if vals[7] > maxs[2]:
            maxs[2] = vals[7]

    if not found_points:
        # Historical helpers would fail on a point-free cell; this defensive
        # fallback keeps the return shape predictable for malformed input.
        return 0.0, 0.0, [[np.inf, np.inf, np.inf], [-np.inf, -np.inf, -np.inf]]

    return math.sqrt(min_d2), math.sqrt(max_d2), [mins.tolist(), maxs.tolist()]


def calc_cell_point_properties_cached(ball_loc, surf_ids, surf_points):
    """Fast spike extrema and bounding box using pre-normalized surface point arrays."""
    ball_loc_arr = np.asarray(ball_loc, dtype=np.float64)

    min_d2 = np.inf
    max_d2 = 0.0
    mins = np.array([np.inf, np.inf, np.inf], dtype=np.float64)
    maxs = np.array([-np.inf, -np.inf, -np.inf], dtype=np.float64)
    found_points = False

    for surf_id in surf_ids:
        points = surf_points[int(surf_id)]
        if points.size == 0:
            continue
        found_points = True
        vals = _cell_point_properties_surface_kernel(ball_loc_arr, points)
        if vals[0] < min_d2:
            min_d2 = vals[0]
        if vals[1] > max_d2:
            max_d2 = vals[1]
        if vals[2] < mins[0]:
            mins[0] = vals[2]
        if vals[3] < mins[1]:
            mins[1] = vals[3]
        if vals[4] < mins[2]:
            mins[2] = vals[4]
        if vals[5] > maxs[0]:
            maxs[0] = vals[5]
        if vals[6] > maxs[1]:
            maxs[1] = vals[6]
        if vals[7] > maxs[2]:
            maxs[2] = vals[7]

    if not found_points:
        return 0.0, 0.0, [[np.inf, np.inf, np.inf], [-np.inf, -np.inf, -np.inf]]

    return math.sqrt(min_d2), math.sqrt(max_d2), [mins.tolist(), maxs.tolist()]


def calc_spikes(ball_loc, surfs):
    """Calculate the minimum and maximum distances (spikes) from a ball's center to all surface points.

    This function measures the distances from a ball's center location to all points on its surrounding
    surfaces, which helps characterize the shape and extent of the ball's influence region.

    Parameters
    ----------
    ball_loc : list or numpy.ndarray
        The 3D coordinates of the ball's center location
    surfs : list of dict
        List of surface dictionaries, where each surface contains a 'points' key with
        a list of 3D coordinates representing surface points

    Returns
    -------
    tuple
        A tuple containing:
        - min_spike (float): The minimum distance from the ball center to any surface point
        - max_spike (float): The maximum distance from the ball center to any surface point

    Notes
    -----
    - Uses calc_dist function to compute Euclidean distances
    - Useful for analyzing the shape and extent of a ball's influence region
    """
    spikes = []
    for surf in surfs:
        for point in surf['points']:
            spikes.append(calc_dist(ball_loc, point))

    return min(spikes), max(spikes)


def calc_cell_box(surfs):
    """Calculate the bounding box of a cell defined by its surfaces.

    This function computes the minimum and maximum coordinates in each dimension
    (x, y, z) that fully enclose all points of the cell's surfaces, effectively
    creating a rectangular prism that bounds the cell.

    Parameters
    ----------
    surfs : list of dict
        List of surface dictionaries, where each surface contains a 'points' key with
        a list of 3D coordinates representing surface points

    Returns
    -------
    list
        A list containing two 3D coordinate lists:
        - [0]: Minimum coordinates [x_min, y_min, z_min]
        - [1]: Maximum coordinates [x_max, y_max, z_max]

    Notes
    -----
    - Useful for determining the spatial extent of a cell
    - Can be used for visualization or spatial analysis
    - Returns a bounding box that may not be axis-aligned if the cell is rotated
    """
    # Create the mins and maxs varaibles
    mins, maxs = [np.inf, np.inf, np.inf], [-np.inf, -np.inf, -np.inf]
    # Loop through the surfaces
    for surf in surfs:
        for point in surf['points']:
            for i in range(3):
                if point[i] < mins[i]:
                    mins[i] = point[i]
                if point[i] > maxs[i]:
                    maxs[i] = point[i]
    # Return the bounding box for the cell
    return [mins, maxs]


@jit(nopython=True, cache=True)
def _tetra_volume_scalar(p0, p1, p2, p3):
    """Fast scalar tetrahedron volume used inside compiled hot loops."""
    a0 = p1[0] - p0[0]
    a1 = p1[1] - p0[1]
    a2 = p1[2] - p0[2]

    b0 = p2[0] - p0[0]
    b1 = p2[1] - p0[1]
    b2 = p2[2] - p0[2]

    c0 = p3[0] - p0[0]
    c1 = p3[1] - p0[1]
    c2 = p3[2] - p0[2]

    cx0 = a1 * b2 - a2 * b1
    cx1 = a2 * b0 - a0 * b2
    cx2 = a0 * b1 - a1 * b0

    triple = c0 * cx0 + c1 * cx1 + c2 * cx2
    if triple < 0.0:
        triple = -triple
    return triple / 6.0


@jit(nopython=True, cache=True)
def _triangle_area_scalar(p0, p1, p2):
    """Fast scalar triangle area used inside compiled hot loops."""
    a0 = p1[0] - p0[0]
    a1 = p1[1] - p0[1]
    a2 = p1[2] - p0[2]

    b0 = p2[0] - p0[0]
    b1 = p2[1] - p0[1]
    b2 = p2[2] - p0[2]

    cx0 = a1 * b2 - a2 * b1
    cx1 = a2 * b0 - a0 * b2
    cx2 = a0 * b1 - a1 * b0
    return 0.5 * math.sqrt(cx0 * cx0 + cx1 * cx1 + cx2 * cx2)


@jit(nopython=True, cache=True)
def _surface_mass_properties_kernel(ball_loc, points, tris, density):
    """
    Accumulate COM numerator and the historical VorPy MOI approximation for one
    triangulated surface.

    The formulas intentionally reproduce the existing calc_cell_com /
    calc_cell_moi definitions; this is a performance rewrite, not a change in
    the reported physical quantity.
    """
    mx = 0.0
    my = 0.0
    mz = 0.0

    ixx = 0.0
    iyy = 0.0
    izz = 0.0
    ixy = 0.0
    ixz = 0.0
    iyz = 0.0

    x0 = ball_loc[0]
    y0 = ball_loc[1]
    z0 = ball_loc[2]

    for t in range(tris.shape[0]):
        i0 = tris[t, 0]
        i1 = tris[t, 1]
        i2 = tris[t, 2]

        p1 = points[i0]
        p2 = points[i1]
        p3 = points[i2]

        tet_vol = _tetra_volume_scalar(ball_loc, p1, p2, p3)
        tet_mass = density * tet_vol

        # Historical tetrahedron centroid.
        cx = (x0 + p1[0] + p2[0] + p3[0]) / 4.0
        cy = (y0 + p1[1] + p2[1] + p3[1]) / 4.0
        cz = (z0 + p1[2] + p2[2] + p3[2]) / 4.0

        mx += tet_vol * cx
        my += tet_vol * cy
        mz += tet_vol * cz

        # Historical calc_tetra_inertia approximation. Accumulate the six
        # independent tensor elements directly to avoid per-triangle matrices.
        # Vertex 0: ball_loc
        scale = tet_mass / 10.0

        ixx += scale * (y0 * y0 + z0 * z0)
        iyy += scale * (x0 * x0 + z0 * z0)
        izz += scale * (x0 * x0 + y0 * y0)
        ixy -= scale * x0 * y0
        ixz -= scale * x0 * z0
        iyz -= scale * y0 * z0

        # Vertices 1-3.
        for p in (p1, p2, p3):
            x = p[0]
            y = p[1]
            z = p[2]
            ixx += scale * (y * y + z * z)
            iyy += scale * (x * x + z * z)
            izz += scale * (x * x + y * y)
            ixy -= scale * x * y
            ixz -= scale * x * z
            iyz -= scale * y * z

        # Historical parallel-axis shift from tetrahedron centroid to ball_loc.
        rx = cx - x0
        ry = cy - y0
        rz = cz - z0

        ixx += tet_mass * (ry * ry + rz * rz)
        iyy += tet_mass * (rx * rx + rz * rz)
        izz += tet_mass * (rx * rx + ry * ry)
        ixy -= tet_mass * rx * ry
        ixz -= tet_mass * rx * rz
        iyz -= tet_mass * ry * rz

    return mx, my, mz, ixx, iyy, izz, ixy, ixz, iyz


def calc_cell_mass_properties(ball_loc, surfs, volume, density=1.0):
    """
    Calculate cell center of mass and moment of inertia in one surface traversal.

    This uses the same tetrahedral decomposition and historical VorPy MOI
    approximation as calc_cell_com() and calc_cell_moi(), but moves the
    triangle-level work into a Numba-compiled kernel.
    """
    ball_loc_arr = np.asarray(ball_loc, dtype=np.float64)

    mx = my = mz = 0.0
    ixx = iyy = izz = ixy = ixz = iyz = 0.0

    for surf in surfs:
        points = np.asarray(surf['points'], dtype=np.float64)
        tris = np.asarray(surf['tris'], dtype=np.int64)

        if tris.size == 0:
            continue
        if tris.ndim != 2:
            tris = tris.reshape((-1, 3))

        vals = _surface_mass_properties_kernel(ball_loc_arr, points, tris, float(density))
        mx += vals[0]
        my += vals[1]
        mz += vals[2]
        ixx += vals[3]
        iyy += vals[4]
        izz += vals[5]
        ixy += vals[6]
        ixz += vals[7]
        iyz += vals[8]

    com = np.array([mx / volume, my / volume, mz / volume], dtype=np.float64)
    moi = np.array([
        [ixx, ixy, ixz],
        [ixy, iyy, iyz],
        [ixz, iyz, izz],
    ], dtype=np.float64)

    return com, moi


def calc_cell_mass_properties_cached(ball_loc, surf_ids, surf_points, surf_tris, volume, density=1.0):
    """COM and historical VorPy MOI using pre-normalized surface geometry arrays."""
    ball_loc_arr = np.asarray(ball_loc, dtype=np.float64)
    mx = my = mz = 0.0
    ixx = iyy = izz = ixy = ixz = iyz = 0.0

    for surf_id in surf_ids:
        sid = int(surf_id)
        points = surf_points[sid]
        tris = surf_tris[sid]
        if tris.size == 0:
            continue
        vals = _surface_mass_properties_kernel(ball_loc_arr, points, tris, float(density))
        mx += vals[0];
        my += vals[1];
        mz += vals[2]
        ixx += vals[3];
        iyy += vals[4];
        izz += vals[5]
        ixy += vals[6];
        ixz += vals[7];
        iyz += vals[8]

    com = np.array([mx / volume, my / volume, mz / volume], dtype=np.float64)
    moi = np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]], dtype=np.float64)
    return com, moi


def calc_cell_com(ball_loc, surfs, volume):
    """Backward-compatible fast center-of-mass calculation."""
    return calc_cell_mass_properties(ball_loc, surfs, volume, density=1.0)[0]


def calc_cell_moi(ball_loc, surfs, volume, density=1.0):
    """Backward-compatible fast historical VorPy moment-of-inertia calculation."""
    return calc_cell_mass_properties(ball_loc, surfs, volume, density=density)[1]


def combine_inertia_tensors(inertia_tensors, centroids, common_centroid, masses):
    """Combines multiple inertia tensors into a single inertia tensor about a common reference point.

    This function implements the parallel axis theorem to shift each inertia tensor from its local
    centroid to a common reference point, then sums them to get the total inertia tensor.

    Parameters
    ----------
    inertia_tensors : list of numpy.ndarray
        List of 3x3 inertia tensors for each element, where each tensor is about its local centroid
    centroids : list of numpy.ndarray
        List of 3D centroid coordinates for each element
    common_centroid : numpy.ndarray
        The reference point to which all inertia tensors will be shifted
    masses : list of float
        List of masses (or volumes if uniform density) for each element

    Returns
    -------
    numpy.ndarray
        A 3x3 inertia tensor representing the combined moment of inertia about the common centroid

    Notes
    -----
    - Uses the parallel axis theorem: I_total = I_local + m(d^2*I - d*d^T)
    - All input arrays should be numpy arrays
    - The function assumes consistent units across all inputs
    """
    # Initialize the total inertia tensor as a zero matrix
    I_total = np.zeros((3, 3))

    # Loop over each element
    for I_i, C_i, m_i in zip(inertia_tensors, centroids, masses):
        # Calculate the displacement vector from the element's centroid to the common centroid
        d = C_i - common_centroid
        d_squared = np.dot(d, d)  # Squared magnitude of the displacement vector

        # Compute the parallel axis theorem adjustment term
        shift_tensor = m_i * (d_squared * np.eye(3) - np.outer(d, d))

        # Shift the inertia tensor of the element to the common centroid and add to total
        I_shifted = I_i + shift_tensor
        I_total += I_shifted

    return I_total


def calc_total_inertia_tensor(spheres, common_point):
    """
    Calculates the total moment of inertia tensor for a collection of spheres about a common reference point.

    This function computes the combined moment of inertia tensor by:
    1. Calculating each sphere's local inertia tensor about its center
    2. Using the parallel axis theorem to shift each tensor to the common reference point
    3. Summing all shifted tensors to obtain the total inertia tensor

    Parameters
    ----------
    spheres : list of dict
        List of sphere dictionaries containing:
        - 'mass' : float
            Mass of the sphere
        - 'rad' : float
            Radius of the sphere
        - 'loc' : numpy.ndarray
            3D coordinates of the sphere's center
    common_point : numpy.ndarray
        3D coordinates of the reference point about which the total inertia tensor is calculated

    Returns
    -------
    numpy.ndarray
        3x3 inertia tensor representing the total moment of inertia about the common point

    Notes
    -----
    - Uses the parallel axis theorem for rigid body mechanics
    - Assumes uniform density spheres
    - All input arrays should be numpy arrays
    - Units should be consistent across all inputs
    """
    # Initialize the total inertia tensor as a 3x3 zero matrix
    I_total = np.zeros((3, 3))

    # Iterate through each sphere
    for sphere in spheres:
        m = sphere['mass']
        r = sphere['rad']
        loc = sphere['loc']

        # Moment of inertia tensor of the sphere about its own center (3x3 identity scaled by (2/5) * m * r^2)
        I_center = (2 / 5) * m * r ** 2 * np.eye(3)

        # Calculate the displacement vector from the sphere's center to the common point
        d = loc - common_point
        d_squared = np.dot(d, d)  # Squared magnitude of the displacement vector

        # Calculate the parallel axis shift tensor: m * (d^2 * I3 - d * d^T)
        shift_tensor = m * (d_squared * np.eye(3) - np.outer(d, d))

        # Shift the inertia tensor to the common reference point
        I_shifted = I_center + shift_tensor

        # Add the shifted inertia tensor to the total inertia tensor
        I_total += I_shifted

    return I_total


@jit(nopython=True, cache=True)
def _contact_surface_kernel(loc, rad, points, tris):
    """Compiled contact-area and contribution-volume calculation for one surface."""
    n_points = points.shape[0]
    projected = np.empty((n_points, 3), dtype=np.float64)
    inside = np.empty(n_points, dtype=np.uint8)

    lx = loc[0]
    ly = loc[1]
    lz = loc[2]
    rad2 = rad * rad

    # Classify/project every point once.
    for i in range(n_points):
        dx = points[i, 0] - lx
        dy = points[i, 1] - ly
        dz = points[i, 2] - lz
        dist2 = dx * dx + dy * dy + dz * dz

        # Historical code compares sqrt(dist2) <= rad. For non-negative
        # radii this squared comparison is mathematically equivalent.
        if dist2 <= rad2:
            inside[i] = 1
            projected[i, 0] = points[i, 0]
            projected[i, 1] = points[i, 1]
            projected[i, 2] = points[i, 2]
        else:
            inside[i] = 0
            norm = math.sqrt(dist2)
            if norm > 0.0:
                scale = rad / norm
                projected[i, 0] = lx + scale * dx
                projected[i, 1] = ly + scale * dy
                projected[i, 2] = lz + scale * dz
            else:
                projected[i, 0] = points[i, 0]
                projected[i, 1] = points[i, 1]
                projected[i, 2] = points[i, 2]

    contact_area = 0.0
    contribution_vol = 0.0

    for t in range(tris.shape[0]):
        i0 = tris[t, 0]
        i1 = tris[t, 1]
        i2 = tris[t, 2]

        n_inside = int(inside[i0]) + int(inside[i1]) + int(inside[i2])

        p0 = points[i0]
        p1 = points[i1]
        p2 = points[i2]

        if n_inside == 3:
            contact_area += _triangle_area_scalar(p0, p1, p2)
            contribution_vol += _tetra_volume_scalar(loc, p0, p1, p2)

        elif n_inside == 0:
            contribution_vol += _tetra_volume_scalar(
                loc, projected[i0], projected[i1], projected[i2]
            )

        else:
            # Historical mixed case: preserve inside vertices and use projected
            # coordinates for outside vertices.
            q0 = p0 if inside[i0] else projected[i0]
            q1 = p1 if inside[i1] else projected[i1]
            q2 = p2 if inside[i2] else projected[i2]

            contribution_vol += _tetra_volume_scalar(loc, q0, q1, q2)
            contact_area += _triangle_area_scalar(p0, p1, p2)

    return contact_area, contribution_vol


def calc_contacts(loc, rad, surfs, surf_ndxs):
    """
    Fast backward-compatible contact-area / contribution-volume calculation.

    The public API and historical triangle classification rules are unchanged;
    point projection and triangle processing now execute inside a compiled
    kernel instead of Python loops.
    """
    loc_arr = np.asarray(loc, dtype=np.float64)
    contact_areas = {}
    contribution_vol = 0.0

    for i, surf in enumerate(surfs):
        points = np.asarray(surf['points'], dtype=np.float64)
        tris = np.asarray(surf['tris'], dtype=np.int64)

        if tris.size == 0:
            contact_areas[surf_ndxs[i]] = 0.0
            continue
        if tris.ndim != 2:
            tris = tris.reshape((-1, 3))

        contact_area, surf_vol = _contact_surface_kernel(
            loc_arr, float(rad), points, tris
        )
        contact_areas[surf_ndxs[i]] = contact_area
        contribution_vol += surf_vol

    return contact_areas, contribution_vol


def calc_contacts_cached(loc, rad, surf_ids, surf_points, surf_tris):
    """Contact areas/volume using pre-normalized surface geometry arrays."""
    loc_arr = np.asarray(loc, dtype=np.float64)
    contact_areas = {}
    contribution_vol = 0.0

    for surf_id in surf_ids:
        sid = int(surf_id)
        points = surf_points[sid]
        tris = surf_tris[sid]
        if tris.size == 0:
            contact_areas[sid] = 0.0
            continue
        contact_area, surf_vol = _contact_surface_kernel(loc_arr, float(rad), points, tris)
        contact_areas[sid] = contact_area
        contribution_vol += surf_vol

    return contact_areas, contribution_vol


def rotate_points(vec, points, reverse=False):
    """Rotates a set of points around a given vector using rotation matrices.

    This function performs a 3D rotation of points around a specified vector by:
    1. Calculating the rotation angles (phi and theta) needed to align the vector with the z-axis
    2. Creating rotation matrices for both z-axis and y-axis rotations
    3. Combining the rotations in the correct order
    4. Applying the combined rotation to all input points

    Parameters
    ----------
    vec : numpy.ndarray
        The vector about which to rotate the points
    points : list of numpy.ndarray
        List of 3D points to be rotated
    reverse : bool, optional
        If True, performs the inverse rotation (default: False)

    Returns
    -------
    list of numpy.ndarray
        List of rotated 3D points

    Notes
    -----
    - Uses standard rotation matrices for 3D transformations
    - Handles both forward and reverse rotations
    - Maintains point positions relative to the rotation vector
    """
    if reverse:
        vec = - vec
    vx, vy, vz = vec
    mag = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    phi = np.arctan2(vy, vx)
    theta = np.arccos(vz / mag)
    if reverse:
        theta, phi = -theta, -phi

    # Forward rotations to align with z-axis
    Rz = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])

    # Combine rotations to align vector with +z direction
    if reverse:
        # Correct sequence for inverse rotation
        rotation_matrix = np.dot(Rz, Ry)
    else:
        # Correct sequence for forward rotation
        rotation_matrix = np.dot(Ry, Rz)

    # Apply rotation to all points
    rotated_points = [np.dot(rotation_matrix, p) for p in points]
    return rotated_points


@jit(nopython=True)
def get_time(seconds):
    """Converts a duration in seconds into hours, minutes, and remaining seconds.

    Parameters
    ----------
    seconds : float
        Total duration in seconds to be converted

    Returns
    -------
    tuple
        A tuple containing (hours, minutes, seconds) where:
        - hours: Number of complete hours
        - minutes: Number of complete minutes after hours
        - seconds: Remaining seconds after hours and minutes

    Examples
    --------
    >>> get_time(3661)
    (1, 1, 1)  # 1 hour, 1 minute, 1 second
    """
    # Divide up the values
    hours = seconds // 3600
    minutes = (seconds - (hours * 3600)) // 60
    seconds = seconds - hours * 3600 - minutes * 60
    # Return the values
    return hours, minutes, seconds


def calc_vol(a_loc, surfs_points, surfs_tris):
    """Calculates the volume of a ball by summing the volumes of tetrahedrons formed between the ball's center and the triangular faces of its surfaces.

    Parameters
    ----------
    a_loc : numpy.ndarray
        The 3D coordinates of the ball's center point
    surfs_points : list of numpy.ndarray
        List of arrays containing the 3D coordinates of points for each surface
    surfs_tris : list of list of tuples
        List of lists containing triangle indices for each surface, where each tuple contains
        three indices referencing points in the corresponding surfs_points array

    Returns
    -------
    tuple
        A tuple containing:
        - float: Total volume of the ball
        - list: List of volumes for each individual surface
    """
    # Create the volume variable
    surf_vols = []
    # Go through each surface on the ball
    for i in range(len(surfs_points)):
        # Calculate the volume of the
        surf_vol = 0
        for tri in surfs_tris[i]:
            # Calculate the tetrahedron volume between the balls' location and the surface triangle's points
            surf_vol += calc_tetra_vol(np.array(a_loc), surfs_points[i][tri[0]], surfs_points[i][tri[1]],
                                       surfs_points[i][tri[2]])
        # Add the surface's volume to the list
        surf_vols.append(surf_vol)
    # Get the total volume by summing the surfaces volumes
    vol = sum(surf_vols)
    # Set the volume and return it
    return vol, surf_vols


def calc_curvature(points, normals):
    """Calculate the mean curvature at each point using neighboring points and normals.

    This function computes the mean curvature at each point in a surface by analyzing
    the local geometry defined by neighboring points and their surface normals.

    Parameters
    ----------
    points : numpy.ndarray
        Array of points on the surface, where each point is a [x, y, z] coordinate
    normals : numpy.ndarray
        Array of surface normals at each point, where each normal is a [nx, ny, nz] vector

    Returns
    -------
    numpy.ndarray
        Array of mean curvature values at each point

    Examples
    --------
    >>> import numpy as np
    >>> points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    >>> normals = np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1]])
    >>> curvatures = calc_curvature(points, normals)
    >>> print(f"Curvatures: {curvatures}")
    Curvatures: [0. 0. 0.]
    """
    # Create the curvature variable
    n_points = len(points)
    curvatures = np.zeros(n_points)

    for i in range(n_points):
        # Find neighboring points (excluding self)
        neighbors = [j for j in range(n_points) if j != i]

        if not neighbors:
            continue

        # Calculate curvature based on normal variations
        normal_variations = []
        for j in neighbors:
            # Project the difference vector onto the normal plane
            diff = points[j] - points[i]
            proj_diff = diff - np.dot(diff, normals[i]) * normals[i]

            if np.linalg.norm(proj_diff) > 0:
                # Calculate the angle between normals
                cos_angle = np.dot(normals[i], normals[j])
                angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))

                # Add to variations list
                normal_variations.append(angle / np.linalg.norm(proj_diff))

        if normal_variations:
            # Mean curvature is the average of normal variations
            curvatures[i] = np.mean(normal_variations)

    return curvatures


def calc_aw_center(r1, r2, l1, l2):
    """
    Calculate the distance between two spheres using the AW method.

    This function calculates the distance between two spheres based on their radii and locations.
    It uses the formula:

    Parameters
    ----------
    r1 : float
        The radius of the first sphere
    r2 : float
        The radius of the second sphere
    l1 : numpy.ndarray
        The location of the first sphere
    l2 : numpy.ndarray
        The location of the second sphere

    Returns
    -------
    tuple
        A tuple containing:
        - float: The aw distance between the two spheres
        - numpy.ndarray: The aw center point between the two spheres
    """
    # Calculate the distance between the two spheres
    dist = np.linalg.norm(l1 - l2)
    # Calculate the aw distance
    aw_dist = dist / 2 - (r2 - r1) / 2
    # Calculate the aw center point
    aw_center = l1 + (l2 - l1) * (aw_dist / dist)
    # Return the aw distance and center point
    return aw_dist, aw_center


def calc_pw_center(r1, r2, l1, l2):
    """
    Calculate the distance between two spheres using the PW method.

    This function calculates the distance between two spheres based on their radii and locations.
    It uses the formula:

    Parameters
    ----------
    r1 : float
        The radius of the first sphere
    r2 : float
        The radius of the second sphere
    l1 : numpy.ndarray
        The location of the first sphere
    l2 : numpy.ndarray
        The location of the second sphere

    Returns
    -------
    tuple
        A tuple containing:
        - float: The pw distance between the two spheres
        - numpy.ndarray: The pw center point between the two spheres
    """
    # Calculate the distance between the two spheres
    dist = np.linalg.norm(l1 - l2)
    # Calculate the pw distance
    pw_dist = dist / 2 - (r2 ** 2 - r1 ** 2) / (2 * dist)
    # Calculate the pw center point
    pw_center = l1 + (l2 - l1) * (pw_dist / dist)
    # Return the pw distance and center point
    return pw_dist, pw_center
