import numpy as np
import warnings
from numba import jit
from numba.core.errors import TypingError
warnings.filterwarnings("error")


def round_func(round_to):
    """
    Nested round function for defining round schemes and rounding multiple values
    :param round_to: int - number of decimal places
    :return: Round function set to round to value
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


def calc_dist(l0, l1):
    return np.sqrt(sum(np.square(l0 - l1)))


@jit(nopython=True)
def calc_dist_numba(l0, l1):
    """
    Calculate distance function used to simplify code
    :param l0: Point 0 list, array, n-dimensional must match point 1
    :param l1: Point 1 list, array, n-dimensional must match point 0
    :return: float distance between the two points
    """
    # Pythagorean theorem
    return np.sqrt(sum(np.square(l0 - l1)))


@jit(nopython=True)
def calc_angle_jit(p0, p1, p2=None):
    """
    Finds the angle (in rads) between three points
    :param p0: Point 0 list, array, n-dimensional must match points 1 and 2
    :param p1: Point 1 list, array, n-dimensional must match points 0 and 2
    :param p2: (optional) Point 2 list, array, n-dimensional must match points 0 and 1
    :return: Angle between (p0, O) and (p1, O) or (p0, p1) and (p0, p2)
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
    Finds the angle (in rads) between three points
    :param p0: Point 0 list, array, n-dimensional must match points 1 and 2
    :param p1: Point 1 list, array, n-dimensional must match points 0 and 2
    :param p2: (optional) Point 2 list, array, n-dimensional must match points 0 and 1
    :return: Angle between (p0, O) and (p1, O) or (p0, p1) and (p0, p2)
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
    Calculates the volume of a tetrahedron defined by its vertices
    :param p0: Point 0
    :param p1: Point 1
    :param p2: Point 2
    :param p3: Point 3
    :return: Volume of the tetrahedron made by the points
    """
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = np.array([p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]])

    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    return (1/6)*abs(np.dot(r03, np.cross(r01, r02)))


@jit(nopython=True)
def calc_tri(points):
    """
    Takes in 3 points and returns the area of the triangle created by them
    :param points: 3D points
    :return: Area of the triangle made by the three points
    """
    # Get the two triangles vectors
    ab = [points[0][0] - points[1][0], points[0][1] - points[1][1], points[0][2] - points[1][2]]
    ac = [points[0][0] - points[2][0], points[0][1] - points[2][1], points[0][2] - points[2][2]]

    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm((np.cross(ab, ac)))


def calc_com(points):
    """
    Takes in a set of points and returns the coordinates of the center of mass
    :param points: lists of locations in n-dimensions
    :return: Center of mass of the inputs
    """

    # Set the running sum for the x, y, z values to 0
    tots = [0 for _ in range(len(points[0]))]
    for point in points:
        for i in range(len(points[0])):
            tots[i] += point[i]

    # Return the center of mass of inputs
    return [tots[i]/len(points) for i in range(len(points[0]))]


@jit(nopython=True)
def calc_length(points):
    """
    Calculates the total length of the points assuming they are in order
    :param points: Points for length calculations
    :return: float total length between consecutive points
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


@jit(nopython=True)
def rotate_points(vec, points, reverse=False):
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


@jit(nopython=True)
def get_time(seconds):
    """
    Turns seconds into hours, minutes and seconds
    :param seconds: Number of seconds in the counter
    :return: hours, minutes, seconds
    """
    # Divide up the values
    hours = seconds // 3600
    minutes = (seconds - (hours * 3600)) // 60
    seconds = seconds - hours * 3600 - minutes * 60
    # Return the values
    return hours, minutes, seconds


def calc_vol(aloc, surfs_points, surfs_tris):
    """
    Calculates the volume of an atom using its surfaces
    :param atom: Atom object for volume calculation
    :return: returns the volume for the atom object
    """
    # Create the volume variable
    surf_vols = []
    # Go through each surface on the atom
    for i in range(len(surfs_points)):
        # Calculate the volume of the
        surf_vol = 0
        for tri in surfs_tris[i]:
            # Calculate the tetrahedron volume between the atoms' location and the surface triangle's points
            surf_vol += calc_tetra_vol(np.array(aloc), surfs_points[i][tri[0]], surfs_points[i][tri[1]], surfs_points[i][tri[2]])
        # Add the surface's volume to the list
        surf_vols.append(surf_vol)
    # Get the total volume by summing the surfaces volumes
    vol = sum(surf_vols)
    # Set the volume and return it
    return vol, surf_vols
