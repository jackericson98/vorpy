import numpy as np


def project_to_plane(points, plane_point, plane_normal):
    """
    Projects a set of 3D points onto a plane defined by a point and normal vector.

    This function takes a collection of 3D points and projects them onto a plane defined by
    a point on the plane and its normal vector. The projection is done by finding the closest
    point on the plane for each input point.

    Parameters
    ----------
    points : list or numpy.ndarray
        Array of 3D points to be projected, where each point is a 3-element array
    plane_point : numpy.ndarray
        A point that lies on the plane, represented as a 3-element array
    plane_normal : numpy.ndarray
        The normal vector of the plane, represented as a 3-element array

    Returns
    -------
    list
        A list of 2D coordinates representing the projected points on the plane,
        where each coordinate is a tuple of (u, v) values in the plane's coordinate system

    Notes
    -----
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

    Parameters
    ----------
    projected_points : list of tuple
        List of 2D coordinates (u, v) representing points projected onto the plane
    plane_point : numpy.ndarray
        A point that lies on the plane, represented as a 3-element array
    plane_normal : numpy.ndarray
        The normal vector of the plane, represented as a 3-element array

    Returns
    -------
    list of numpy.ndarray
        A list of 3D points reconstructed on the plane, where each point is a 3-element array

    Notes
    -----
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

    Parameters
    ----------
    points_2d : list of tuple
        List of 2D coordinates (u, v) to be mapped onto the plane
    plane_point : numpy.ndarray
        A point that lies on the plane, represented as a 3-element array
    plane_normal : numpy.ndarray
        The normal vector of the plane, represented as a 3-element array

    Returns
    -------
    list of numpy.ndarray
        A list of 3D points mapped onto the plane, where each point is a 3-element array

    Notes
    -----
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