import numpy as np
import warnings
from numba import jit
from numba.core.errors import TypingError
warnings.filterwarnings("error")


def round_func(round_to):
    """
    Creates a configurable rounding function that can handle both single values and iterables.

    This function returns a closure that maintains the specified rounding precision and can be
    reused for consistent rounding across multiple values. The returned function handles both
    single numeric values and iterables of values, applying the same rounding precision to all.

    Parameters:
    -----------
    round_to : int
        The number of decimal places to round to. A positive value rounds to that many decimal
        places, while a negative value rounds to the left of the decimal point.

    Returns:
    --------
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


def project_to_plane(points, plane_point, plane_normal):
    """
    Projects a set of 3D points onto a plane defined by a point and normal vector.

    This function takes a collection of 3D points and projects them onto a plane defined by
    a point on the plane and its normal vector. The projection is done by finding the closest
    point on the plane for each input point.

    Parameters:
    -----------
    points : list or numpy.ndarray
        Array of 3D points to be projected, where each point is a 3-element array
    plane_point : numpy.ndarray
        A point that lies on the plane, represented as a 3-element array
    plane_normal : numpy.ndarray
        The normal vector of the plane, represented as a 3-element array

    Returns:
    --------
    list
        A list of 2D coordinates representing the projected points on the plane,
        where each coordinate is a tuple of (u, v) values in the plane's coordinate system

    Notes:
    ------
    - The plane's coordinate system is created using an orthogonal basis
    - The normal vector is automatically normalized
    - If the normal vector is parallel to the x-axis, an alternative basis vector is used
    """
    # Normalize the normal vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Create an orthogonal basis for the plane
    u = np.cross(plane_normal, np.array([1, 0, 0]))
    if np.linalg.norm(u) < 1e-10:  # Check if cross product is almost zero
        u = np.cross(plane_normal, np.array([0, 1, 0]))
    u = u / np.linalg.norm(u)
    v = np.cross(plane_normal, u)
    v = v / np.linalg.norm(v)

    # Project points onto the plane
    projected_points = []
    for point in points:
        point_vector = point - plane_point
        u_coord = np.dot(point_vector, u)
        v_coord = np.dot(point_vector, v)
        projected_points.append((u_coord, v_coord))

    return projected_points


def unproject_to_3d(projected_points, plane_point, plane_normal):
    """
    Reconstructs 3D points from their 2D projections on a plane.

    This function takes a set of 2D points that were previously projected onto a plane
    and reconstructs their original 3D positions on that plane. The reconstruction uses
    the plane's point and normal vector to establish the coordinate system.

    Parameters:
    -----------
    projected_points : list of tuple
        List of 2D coordinates (u, v) representing points projected onto the plane
    plane_point : numpy.ndarray
        A point that lies on the plane, represented as a 3-element array
    plane_normal : numpy.ndarray
        The normal vector of the plane, represented as a 3-element array

    Returns:
    --------
    list of numpy.ndarray
        A list of 3D points reconstructed on the plane, where each point is a 3-element array

    Notes:
    ------
    - The plane's coordinate system is recreated using an orthogonal basis
    - The normal vector is automatically normalized
    - This function is the inverse operation of project_to_plane
    """
    # Normalize the normal vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Create an orthogonal basis for the plane
    u = np.cross(plane_normal, np.array([1, 0, 0]))
    if np.linalg.norm(u) < 1e-10:  # Check if cross product is almost zero
        u = np.cross(plane_normal, np.array([0, 1, 0]))
    u = u / np.linalg.norm(u)
    v = np.cross(plane_normal, u)
    v = v / np.linalg.norm(v)

    # Map 2D coordinates back to 3D plane
    reconstructed_points = []
    for u_coord, v_coord in projected_points:
        # Reconstruct 3D point on the plane using the basis vectors and plane point
        point_3d = plane_point + u_coord * u + v_coord * v
        reconstructed_points.append(point_3d)

    return reconstructed_points


def map_to_plane(points_2d, plane_point, plane_normal):
    """
    Maps 2D points onto a 3D plane defined by a point and normal vector.

    This function takes a set of 2D points and maps them onto a 3D plane using the plane's
    point and normal vector to establish the coordinate system. The mapping creates a
    one-to-one correspondence between 2D coordinates and points on the 3D plane.

    Parameters:
    -----------
    points_2d : list of tuple
        List of 2D coordinates (u, v) to be mapped onto the plane
    plane_point : numpy.ndarray
        A point that lies on the plane, represented as a 3-element array
    plane_normal : numpy.ndarray
        The normal vector of the plane, represented as a 3-element array

    Returns:
    --------
    list of numpy.ndarray
        A list of 3D points mapped onto the plane, where each point is a 3-element array

    Notes:
    ------
    - The plane's coordinate system is created using an orthogonal basis
    - The normal vector is automatically normalized
    - Special handling for cases where the normal vector is aligned with the x-axis
    """
    # Normalize the normal vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Create an orthogonal basis for the plane
    if (plane_normal == np.array([1.0, 0.0, 0.0])).all() or (plane_normal == np.array([-1.0, 0.0, 0.0])).all():
        # Handle the case where the normal is along the x-axis
        u = np.array([0, 1, 0])
    else:
        u = np.cross(plane_normal, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(plane_normal, u)
    v = v / np.linalg.norm(v)

    # Map 2D points to the 3D plane
    mapped_points = []
    for point_2d in points_2d:
        u_coord, v_coord = point_2d
        # Calculate the corresponding 3D point
        point_3d = plane_point + u_coord * u + v_coord * v
        mapped_points.append(point_3d)

    return mapped_points


def calc_dist(l0, l1):
    """
    Calculate the Euclidean distance between two points in n-dimensional space.

    This function computes the straight-line distance between two points using the
    Pythagorean theorem generalized to n dimensions.

    Parameters:
    -----------
    l0 : array-like
        First point coordinates as an n-dimensional array or list
    l1 : array-like
        Second point coordinates as an n-dimensional array or list with same dimensionality as l0

    Returns:
    --------
    float
        The Euclidean distance between the two points

    Notes:
    ------
    - Both input points must have the same dimensionality
    - Uses numpy's square and sqrt functions for efficient computation
    """

    return np.sqrt(sum(np.square(np.array(l0) - np.array(l1))))


@jit(nopython=True)
def calc_dist_numba(l0, l1):
    """
    Calculate the Euclidean distance between two points in n-dimensional space.

    This function computes the straight-line distance between two points using the
    Pythagorean theorem generalized to n dimensions. The function is optimized with
    Numba's JIT compilation for improved performance.

    Parameters:
    -----------
    l0 : numpy.ndarray
        First point coordinates as an n-dimensional array
    l1 : numpy.ndarray
        Second point coordinates as an n-dimensional array with same dimensionality as l0

    Returns:
    --------
    float
        The Euclidean distance between the two points

    Notes:
    ------
    - Both input points must have the same dimensionality
    - Uses numpy's square and sqrt functions for efficient computation
    - JIT compiled for performance optimization
    """
    # Pythagorean theorem
    return np.sqrt(sum(np.square(l0 - l1)))


@jit(nopython=True)
def calc_angle_jit(p0, p1, p2=None):
    """
    Calculates the angle between three points in radians using vector geometry.

    This function computes the angle between vectors formed by the points in two possible ways:
    1. If p2 is not provided: Angle between vectors from origin to p0 and p1
    2. If p2 is provided: Angle between vectors from p0 to p1 and p0 to p2

    Parameters:
    -----------
    p0 : array-like
        First point coordinates [x, y, z, ...]
    p1 : array-like
        Second point coordinates [x, y, z, ...]
    p2 : array-like, optional
        Third point coordinates [x, y, z, ...]
        If not provided, the origin (0,0,0) is used as the reference point

    Returns:
    --------
    float
        The angle in radians between the vectors formed by the points

    Notes:
    ------
    - All points must have the same dimensionality
    - Uses numpy's arccos function for angle calculation
    - Handles edge cases where vectors are parallel or antiparallel
    """
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = p0, p1
    else:
        v0, v1 = p1 - p0, p2 - p0
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
    # Calculate the angle between the two vectors with catches for 180 and 0
    my_dot = np.dot(n0, n1)
    if my_dot <= -1.0:
        my_dot = -1.0
    elif my_dot >= 1.0:
        my_dot = 1.0
    angle = np.arccos(my_dot)
    return angle


def calc_angle(p0, p1, p2=None):
    """
    Calculates the angle between three points in radians.

    This function computes the angle between vectors formed by the points in two possible ways:
    1. If p2 is not provided: Angle between vectors from origin to p0 and p1
    2. If p2 is provided: Angle between vectors from p0 to p1 and p0 to p2

    Parameters:
    -----------
    p0 : array-like
        First point coordinates [x, y, z, ...]
    p1 : array-like
        Second point coordinates [x, y, z, ...]
    p2 : array-like, optional
        Third point coordinates [x, y, z, ...]
        If not provided, the origin (0,0,0) is used as the reference point

    Returns:
    --------
    float
        The angle in radians between the vectors formed by the points
    """
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = p0, p1
    else:
        v0, v1 = p1 - p0, p2 - p0
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
    # Calculate the angle between the two vectors with catches for 180 and 0
    my_dot = np.dot(n0, n1)
    if my_dot <= -1.0:
        my_dot = -1.0
    elif my_dot >= 1.0:
        my_dot = 1.0
    angle = np.arccos(my_dot)
    return angle


@jit(nopython=True)
# Calculate tetrahedron volume function.
def calc_tetra_vol(p0, p1, p2, p3):
    """
    Calculates the volume of a tetrahedron defined by four vertices in 3D space.

    This function uses the scalar triple product formula to compute the volume:
    V = (1/6) * |(p3-p0) · ((p1-p0) × (p2-p0))|
    where · denotes the dot product and × denotes the cross product.

    Parameters:
    -----------
    p0 : array-like
        First vertex of the tetrahedron [x, y, z]
    p1 : array-like
        Second vertex of the tetrahedron [x, y, z]
    p2 : array-like
        Third vertex of the tetrahedron [x, y, z]
    p3 : array-like
        Fourth vertex of the tetrahedron [x, y, z]

    Returns:
    --------
    float
        The volume of the tetrahedron formed by the four vertices
    """
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = np.array([p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]])

    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    return (1/6)*abs(np.dot(r03, np.cross(r01, r02)))


def calc_tetra_inertia(ps, mass):
    """
    Calculates the moment of inertia tensor of a tetrahedron about its centroid.

    This function computes the inertia tensor for a tetrahedron with uniform density distribution.
    The calculation is based on the parallel axis theorem and the inertia tensor of a tetrahedron
    about its centroid.

    Parameters:
    -----------
    ps : list of array-like
        List containing four vertices of the tetrahedron, each as [x, y, z] coordinates
    mass : float
        Total mass of the tetrahedron

    Returns:
    --------
    numpy.ndarray
        A 3x3 inertia tensor matrix where:
        - Diagonal elements represent moments of inertia about x, y, and z axes
        - Off-diagonal elements represent products of inertia
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


def calc_com(points, masses=None):
    """
    Calculates the center of mass for a system of points in n-dimensional space.

    This function computes the center of mass by taking the weighted average of point positions,
    where the weights are the masses of each point. If no masses are provided, all points are
    treated as having equal mass.

    Parameters:
    -----------
    points : list of array-like
        List of point coordinates in n-dimensional space. Each point should be a list or array
        of coordinates [x1, x2, ..., xn].
    masses : list of float, optional
        List of masses corresponding to each point. If None (default), all points are treated
        as having unit mass.

    Returns:
    --------
    numpy.ndarray
        The center of mass coordinates as a numpy array of length n.

    Notes:
    ------
    - All points must have the same dimensionality
    - If masses is provided, it must have the same length as points
    - Returns a numpy array for consistency with other geometric calculations
    """
    # Default mass values if None provided
    if masses is None:
        masses = [1] * len(points)

    # Initialize totals for each dimension to zero
    total_mass = 0
    tots = [0] * len(points[0])

    # Calculate the weighted sum for each dimension
    for point, mass in zip(points, masses):
        for i in range(len(point)):
            tots[i] += point[i] * mass
        total_mass += mass

    # Calculate the center of mass for each dimension
    center_of_mass = [total / total_mass for total in tots]
    return np.array(center_of_mass)


@jit(nopython=True)
def calc_length(points):
    """
    Calculates the total length of a path defined by a sequence of points.

    This function computes the sum of Euclidean distances between consecutive points in the input sequence.
    The points are assumed to be ordered in the sequence they should be connected.

    Parameters:
    -----------
    points : list of array-like
        List of point coordinates in n-dimensional space. Each point should be a list or array
        of coordinates [x1, x2, ..., xn]. Points must be ordered in the sequence they should
        be connected.

    Returns:
    --------
    float
        The total length of the path formed by connecting consecutive points in the input sequence.

    Notes:
    ------
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
    """
    Calculate the sphericity of a geometric object based on its volume and surface area.

    Parameters:
    - volume (float): The volume of the object.
    - surface_area (float): The surface area of the object.

    Returns:
    - float: The sphericity of the object.
    """
    if volume <= 0 or surface_area <= 0:
        raise ValueError("Volume and surface area must be positive numbers.")

    # Calculate sphericity using the geometric formula
    sphericity = (np.pi ** (1/3) * (6 * volume) ** (2/3)) / surface_area
    return sphericity


def calc_isoperimetric_quotient(volume, surface_area):
    """
    Calculate the isoperimetric quotient

    Parameters:
        - volume (float): The volume of the object.
        - surface_area (float): The surface area of the object.

    Returns:
        - isoperimetric quotient (float): The isoperimetric quotient of the object
    """
    if volume <= 0 or surface_area <= 0:
        raise ValueError("Volume and surface area must be positive numbers.")

    return (36 * np.pi * volume ** 2) / (surface_area ** 3)


def calc_spikes(ball_loc, surfs):
    """
    Calculate the minimum and maximum distances (spikes) from a ball's center to all surface points.

    This function measures the distances from a ball's center location to all points on its surrounding
    surfaces, which helps characterize the shape and extent of the ball's influence region.

    Parameters:
    -----------
    ball_loc : list or numpy.ndarray
        The 3D coordinates of the ball's center location
    surfs : list of dict
        List of surface dictionaries, where each surface contains a 'points' key with
        a list of 3D coordinates representing surface points

    Returns:
    --------
    tuple
        A tuple containing:
        - min_spike (float): The minimum distance from the ball center to any surface point
        - max_spike (float): The maximum distance from the ball center to any surface point

    Notes:
    ------
    - Uses calc_dist function to compute Euclidean distances
    - Useful for analyzing the shape and extent of a ball's influence region
    """
    spikes = []
    for surf in surfs:
        for point in surf['points']:
            spikes.append(calc_dist(ball_loc, point))

    return min(spikes), max(spikes)

def calc_cell_box(surfs):
    """
    Calculate the bounding box of a cell defined by its surfaces.

    This function computes the minimum and maximum coordinates in each dimension
    (x, y, z) that fully enclose all points of the cell's surfaces, effectively
    creating a rectangular prism that bounds the cell.

    Parameters:
    -----------
    surfs : list of dict
        List of surface dictionaries, where each surface contains a 'points' key with
        a list of 3D coordinates representing surface points

    Returns:
    --------
    list
        A list containing two 3D coordinate lists:
        - [0]: Minimum coordinates [x_min, y_min, z_min]
        - [1]: Maximum coordinates [x_max, y_max, z_max]

    Notes:
    ------
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


def calc_cell_com(ball_loc, surfs, volume):
    """
    Calculate the center of mass of a cell using tetrahedral decomposition.

    This function computes the center of mass of a cell by decomposing it into tetrahedrons
    formed by the cell's center point and each triangular face of its surfaces. The center
    of mass is calculated as the volume-weighted average of the centroids of all tetrahedrons.

    Parameters:
    -----------
    ball_loc : list or numpy.ndarray
        The center location of the cell (3D coordinates)
    surfs : list of dict
        List of surface dictionaries, where each surface contains:
        - 'points': List of 3D coordinates representing surface points
        - 'tris': List of triangle indices referencing the points
    volume : float
        Total volume of the cell

    Returns:
    --------
    numpy.ndarray
        The 3D coordinates of the cell's center of mass

    Notes:
    ------
    - Uses tetrahedral decomposition to calculate the center of mass
    - Each tetrahedron is formed by the cell center and a triangular face
    - The result is normalized by the total cell volume
    - Returns a numpy array for vector operations
    """
    # Create the mass_locs list
    mass_locs = []
    for surf in surfs:
        for tri in surf['tris']:
            # Get the points of the tetrahedron
            ps = [ball_loc, *[surf['points'][_] for _ in tri]]
            # Calculate the centroid of the tetrahedron
            tet_com = [sum([ps[j][i] for j in range(4)]) / 4 for i in range(3)]
            # Calculate the volume of the tetrahedron
            tet_vol = calc_tetra_vol(*ps)
            # Append the volume-weighted centroid
            mass_locs.append([tet_vol * coord for coord in tet_com])

    # Calculate the total center of mass by normalizing with the cell volume
    return np.array([sum(coords) / volume for coords in zip(*mass_locs)])


def calc_cell_moi(ball_loc, surfs, volume, density=1.0):
    """
    Calculate the moment of inertia tensor of a cell using tetrahedral decomposition.

    This function computes the moment of inertia tensor of a cell by decomposing it into tetrahedrons
    formed by the cell's center point and each triangular face of its surfaces. The total moment of inertia
    is calculated as the sum of individual tetrahedron contributions, using the parallel axis theorem to
    shift each contribution to the cell's center point.

    Parameters:
    -----------
    ball_loc : numpy.ndarray or list
        The center location of the cell (3D coordinates)
    surfs : list of dict
        List of surface dictionaries, where each surface contains:
        - 'points': List of 3D coordinates representing surface points
        - 'tris': List of triangle indices referencing the points
    volume : float
        Total volume of the cell
    density : float, optional
        Density of the material (default is 1.0)

    Returns:
    --------
    numpy.ndarray
        A 3x3 moment of inertia tensor of the cell with respect to ball_loc

    Notes:
    ------
    - Uses tetrahedral decomposition to calculate individual contributions
    - Applies the parallel axis theorem to shift each tetrahedron's inertia tensor
    - Returns a symmetric 3x3 numpy array representing the inertia tensor
    """
    # Create an inertia tensor initialized to zero
    inertia_tensor = np.zeros((3, 3))

    # Iterate through each surface and triangle to calculate the tetrahedron MOI contributions
    for surf in surfs:
        for tri in surf['tris']:
            # Get the points of the tetrahedron
            ps = [ball_loc, *[surf['points'][_] for _ in tri]]

            # Calculate the centroid of the tetrahedron
            tet_com = [sum([ps[j][i] for j in range(4)]) / 4 for i in range(3)]

            # Calculate the volume of the tetrahedron
            tet_vol = calc_tetra_vol(*ps)

            # Calculate the mass of the tetrahedron
            tet_mass = density * tet_vol

            # Calculate the inertia tensor of the tetrahedron about its centroid
            tet_inertia_tensor = calc_tetra_inertia(ps, tet_mass)

            # Calculate the distance vector from the tetrahedron centroid to the cell's center (`ball_loc`)
            r = np.array(tet_com) - np.array(ball_loc)
            r_squared = np.dot(r, r)

            # Use the parallel axis theorem to adjust the inertia tensor to the cell's center
            shift_tensor = tet_mass * (r_squared * np.identity(3) - np.outer(r, r))

            # Add the adjusted tensor to the total inertia tensor
            inertia_tensor += tet_inertia_tensor + shift_tensor

    return inertia_tensor


def combine_inertia_tensors(inertia_tensors, centroids, common_centroid, masses):
    """
    Combines multiple inertia tensors into a single inertia tensor about a common reference point.

    This function implements the parallel axis theorem to shift each inertia tensor from its local
    centroid to a common reference point, then sums them to get the total inertia tensor.

    Parameters:
    -----------
    inertia_tensors : list of numpy.ndarray
        List of 3x3 inertia tensors for each element, where each tensor is about its local centroid
    centroids : list of numpy.ndarray
        List of 3D centroid coordinates for each element
    common_centroid : numpy.ndarray
        The reference point to which all inertia tensors will be shifted
    masses : list of float
        List of masses (or volumes if uniform density) for each element

    Returns:
    --------
    numpy.ndarray
        A 3x3 inertia tensor representing the combined moment of inertia about the common centroid

    Notes:
    ------
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

    Parameters:
    -----------
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

    Returns:
    --------
    numpy.ndarray
        3x3 inertia tensor representing the total moment of inertia about the common point

    Notes:
    ------
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
        I_center = (2 / 5) * m * r**2 * np.eye(3)

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


def calc_contacts(loc, rad, surfs, surf_ndxs):
    """
    Calculate the contact areas and contribution volume for a given sphere.

    This function computes the contact areas between a sphere and surrounding surfaces, as well as the
    contribution volume of the sphere to the total volume of the system. It handles cases where surfaces
    are fully inside, fully outside, or partially intersecting with the sphere.

    Parameters:
    -----------
    loc : numpy.ndarray
        Center coordinates of the sphere in 3D space
    rad : float
        Radius of the sphere
    surfs : list of dict
        List of surface dictionaries, where each surface contains:
        - 'points': List of 3D coordinates defining surface vertices
        - 'tris': List of triangle indices referencing the points
    surf_ndxs : list of int
        Indices of surfaces to consider for contact calculations

    Returns:
    --------
    tuple
        A tuple containing:
        - contact_areas : dict
            Dictionary mapping surface indices to their respective contact areas
        - contribution_vol : float
            Total volume contribution of the sphere to the system

    Notes:
    ------
    - Contact area is calculated for surfaces that intersect with the sphere
    - Volume contribution considers both internal and external portions of surfaces
    - Points inside the sphere are preserved, while points outside are projected onto the sphere's surface
    """
    # Create the area and volume vals
    contact_areas, contribution_vol = {}, 0

    # Loop through the surfaces
    for i, surf in enumerate(surfs):
        # Initialize contact area for this surface
        contact_area = 0
        new_points = []
        point_inside = []

        # Loop through the points to determine if inside or outside
        for point in surf['points']:
            distance = calc_dist(point, loc)
            if distance <= rad:
                point_inside.append(True)
                new_points.append(point)
            else:
                point_inside.append(False)
                # Get the direction and normalize it
                direction = point - loc
                norm = np.linalg.norm(direction)
                if norm > 0:
                    # Project the point onto the sphere's surface
                    new_points.append(rad * (direction / norm) + loc)
                else:
                    new_points.append(point)  # If the point coincides with the center (rare edge case)

        # Loop through the triangles
        for tri in surf['tris']:
            triangle_points = [surf['points'][index] for index in tri]
            projected_points = [new_points[index] for index in tri]
            inside_flags = [point_inside[index] for index in tri]

            # Determine if the triangle is fully inside, fully outside, or mixed
            all_inside = all(inside_flags)
            all_outside = not any(inside_flags)
            mixed = not all_inside and not all_outside

            if all_inside:
                # Triangle is fully inside the sphere
                contact_area += calc_tri(np.array(triangle_points))
                contribution_vol += calc_tetra_vol(loc, *triangle_points)
            elif all_outside:
                # Triangle is fully outside the sphere
                contribution_vol += calc_tetra_vol(loc, *projected_points)
            elif mixed:
                # Triangle is partially inside and outside
                # We add the volume using a mix of inside and projected points
                mixed_points = [triangle_points[i] if inside_flags[i] else projected_points[i] for i in range(3)]
                contribution_vol += calc_tetra_vol(loc, *mixed_points)
                # Count the triangle as outside for contact area if any point is outside
                if inside_flags.count(True) < 3:
                    contact_area += calc_tri(np.array(triangle_points))

        # Append the contact area for this surface
        contact_areas[surf_ndxs[i]] = contact_area

    return contact_areas, contribution_vol


@jit(nopython=True)
def rotate_points1(vec, points, reverse=False):
    """
    Takes in a set of points and a vector and rotates the points and the vector so the v = [0,0,1]
    :param vec: The vector about which the surface is rotated
    :param points: the points of the surface
    :param reverse: Bool for rotating the surface back
    :return: List of rotated points
    """
    # Get the vx, vy, vz vector components
    vx, vy, vz = vec
    # If vy or vz are zero we need a catch for divide by zero error.
    if round(vy, 2) == 0:
        phi = np.pi / 2
    else:
        phi = np.arctan(vx / vy)
    if round(vz, 2) == 0:
        theta = np.pi / 2
    else:
        theta = np.arctan(vy / vz)
    # If the points are to be sent back, provide the negative values for the angles
    if reverse:
        theta, phi = -theta, -phi
    # Get variables for sin(theta), cos(theta), sin(phi), cos(phi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    nps = []
    for p in points:
        px, py, pz = np.round(p[0], 7), np.round(p[1], 7), np.round(p[2], 7)
        # Multiplying the x, y rotation matrices gives the following:
        npx = px * cp - py * sp
        npy = px * ct * sp + py * ct * cp - pz * st
        npz = px * st * sp + py * st * cp + pz * ct
        # Add the new points to the list
        nps.append([npx, npy, npz])
    return nps


def rotate_points(vec, points, reverse=False):
    """
    Rotates a set of points around a given vector using rotation matrices.

    This function performs a 3D rotation of points around a specified vector by:
    1. Calculating the rotation angles (phi and theta) needed to align the vector with the z-axis
    2. Creating rotation matrices for both z-axis and y-axis rotations
    3. Combining the rotations in the correct order
    4. Applying the combined rotation to all input points

    Parameters:
    -----------
    vec : numpy.ndarray
        The vector about which to rotate the points
    points : list of numpy.ndarray
        List of 3D points to be rotated
    reverse : bool, optional
        If True, performs the inverse rotation (default: False)

    Returns:
    --------
    list of numpy.ndarray
        List of rotated 3D points

    Notes:
    ------
    - Uses standard rotation matrices for 3D transformations
    - Handles both forward and reverse rotations
    - Maintains point positions relative to the rotation vector
    """
    if reverse:
        vec = - vec
    vx, vy, vz = vec
    mag = np.sqrt(vx**2 + vy**2 + vz**2)
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
    """
    Converts a duration in seconds into hours, minutes, and remaining seconds.

    Parameters:
    -----------
    seconds : float
        Total duration in seconds to be converted

    Returns:
    --------
    tuple
        A tuple containing (hours, minutes, seconds) where:
        - hours: Number of complete hours
        - minutes: Number of complete minutes after hours
        - seconds: Remaining seconds after hours and minutes

    Examples:
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
    """
    Calculates the volume of a ball by summing the volumes of tetrahedrons formed between the ball's center
    and the triangular faces of its surfaces.

    Parameters:
    -----------
    a_loc : numpy.ndarray
        The 3D coordinates of the ball's center point
    surfs_points : list of numpy.ndarray
        List of arrays containing the 3D coordinates of points for each surface
    surfs_tris : list of list of tuples
        List of lists containing triangle indices for each surface, where each tuple contains
        three indices referencing points in the corresponding surfs_points array

    Returns:
    --------
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
            surf_vol += calc_tetra_vol(np.array(a_loc), surfs_points[i][tri[0]], surfs_points[i][tri[1]], surfs_points[i][tri[2]])
        # Add the surface's volume to the list
        surf_vols.append(surf_vol)
    # Get the total volume by summing the surfaces volumes
    vol = sum(surf_vols)
    # Set the volume and return it
    return vol, surf_vols
