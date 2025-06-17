import numpy as np
from numpy import array, dot, isreal, linalg, roots
from numba import jit
from vorpy.src.calculations.calcs import calc_dist, calc_dist_numba


@jit(nopython=True)
def calc_vert_abcfs(locs, rads):
    """
    Calculate and organize coefficients for solving the system of equations that determine additively weighted vertices.

    This function calculates the coefficients necessary for finding vertices of the inscribed sphere from the
    locations and radii of four spheres. It adjusts all sphere locations relative to the first sphere's location
    for simpler calculation and computes coefficients for a system of linear equations derived from geometric
    properties.

    Parameters
    ----------
    locs : numpy.ndarray of arrays
        Coordinates of the centers of the four spheres
    rads : numpy.ndarray of floats
        Radii of the four spheres

    Returns
    -------
    tuple
        Contains arrays of calculated coefficients (fs, abcdfs), an array of radii (rs), and the base location (l0)

    Notes
    -----
    The function adjusts all sphere locations relative to the first sphere's location for simpler calculation.
    It then calculates the coefficients of a system of linear equations derived from the geometric properties
    of the spheres.
    """

    # Unpack the radii of the four spheres
    r0, r1, r2, r3 = rads

    # Calculate the square of the first sphere's radius for use in equations
    r0_2 = r0 ** 2

    # Adjust locations relative to the first sphere's location to simplify the system of equations
    l0, l1, l2, l3 = locs[0], locs[1] - locs[0], locs[2] - locs[0], locs[3] - locs[0]

    # Calculate the coefficients for the system of linear equations
    a1, b1, c1, d1, f1 = 2 * l1[0], 2 * l1[1], 2 * l1[2], 2 * (r1 - r0), r0_2 - r1 ** 2 + l1[0] ** 2 + l1[1] ** 2 + l1[2] ** 2
    a2, b2, c2, d2, f2 = 2 * l2[0], 2 * l2[1], 2 * l2[2], 2 * (r2 - r0), r0_2 - r2 ** 2 + l2[0] ** 2 + l2[1] ** 2 + l2[2] ** 2
    a3, b3, c3, d3, f3 = 2 * l3[0], 2 * l3[1], 2 * l3[2], 2 * (r3 - r0), r0_2 - r3 ** 2 + l3[0] ** 2 + l3[1] ** 2 + l3[2] ** 2

    # Calculate determinant and other coefficients for solving the vertex positions
    F = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1
    F_2 = F ** 2
    F10 = b1 * c2 * f3 - b1 * c3 * f2 - b2 * c1 * f3 + b2 * c3 * f1 + b3 * c1 * f2 - b3 * c2 * f1
    F11 = -b1 * c2 * d3 + b1 * c3 * d2 + b2 * c1 * d3 - b2 * c3 * d1 - b3 * c1 * d2 + b3 * c2 * d1
    F20 = -a1 * c2 * f3 + a1 * c3 * f2 + a2 * c1 * f3 - a2 * c3 * f1 - a3 * c1 * f2 + a3 * c2 * f1
    F21 = a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1
    F30 = a1 * b2 * f3 - a1 * b3 * f2 - a2 * b1 * f3 + a2 * b3 * f1 + a3 * b1 * f2 - a3 * b2 * f1
    F31 = -a1 * b2 * d3 + a1 * b3 * d2 + a2 * b1 * d3 - a2 * b3 * d1 - a3 * b1 * d2 + a3 * b2 * d1

    # Store the calculated coefficients in arrays for easy access
    fs = array([F, F_2, F10, F11, F20, F21, F30, F31])
    abcdfs = array([[a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f3]])
    rs = array([r0, r1, r2, r3])

    # Return the necessary values for vertex calculation
    return fs, abcdfs, rs, l0


def calc_vert_case_1(Fs, l0, r0):
    """
    Calculate vertices for Case 1 in a vertex calculation scenario involving spheres.

    This function solves a quadratic equation to determine possible radii (R values) and their corresponding
    vertex coordinates. It handles the case where the vertex calculation involves solving a quadratic equation
    to determine valid radii and uses these radii to compute vertex coordinates.

    Parameters
    ----------
    Fs : list
        List of polynomial coefficients F, F_2, F10, F11, etc., that define the conditions for vertex calculation
    l0 : array
        The original location of the sphere center used to adjust the calculated vertices back to the actual position
    r0 : float
        The radius component used in the calculation of polynomial coefficients

    Returns
    -------
    list
        A list of vertices, where each vertex is represented as a list containing its x, y, z coordinates and the radius R

    Notes
    -----
    The function solves a quadratic equation to determine valid radii and uses these radii to compute vertex
    coordinates. Only real and positive roots of the quadratic equation are considered for vertex calculation.
    """

    # Unwrap the polynomial coefficients from Fs for convenience
    F, F_2, F10, F11, F20, F21, F30, F31 = Fs

    # # Compute the coefficients of the quadratic equation for the radius R
    # a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F_2) - 1  # Quadratic term
    # b = 2 * (((F10 * F11 + F20 * F21 + F30 * F31) / F_2) - r0)  # Linear term
    # c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F_2) - r0 ** 2  # Constant term

    # Compute the coefficients of the quadratic equation for the radius R
    a = (F11 ** 2 + F21 ** 2 + F31 ** 2) - F_2  # Quadratic term
    b = 2 * (((F10 * F11 + F20 * F21 + F30 * F31)) - r0 * F_2)  # Linear term
    c = (F10 ** 2 + F20 ** 2 + F30 ** 2) - F_2 *r0 ** 2  # Constant term
    # Initialize an empty list to store the vertices
    verts = []

    # Check if the quadratic equation has real solutions (discriminant >= 0)
    if -4 * a * c + b ** 2 >= 0:
        # Solve the quadratic equation and filter real roots
        Rs = [R for R in roots(array([a, b, c])) if isreal(R)]
    else:
        return  # Exit if the discriminant is negative (no real solutions)

    # Ensure there are valid roots to process
    if Rs is not None and len(Rs) > 0:
        # Loop through each valid radius (R) and calculate the corresponding vertex coordinates
        for R in Rs:
            x = F10 / F + R * F11 / F + l0[0]  # x-coordinate of the vertex
            y = F20 / F + R * F21 / F + l0[1]  # y-coordinate of the vertex
            z = F30 / F + R * F31 / F + l0[2]  # z-coordinate of the vertex
            # Append the calculated vertex (including radius R) to the list
            verts.append([x, y, z, R])

    # Return the list of vertices
    return verts


@jit(nopython=True)
def calc_vert_case_2(Fs, r0, l0):
    """
    Calculate vertices for Case 2 in a vertex calculation scenario involving spheres.

    This function computes vertices based on polynomial roots derived from given coefficients, which describe
    the geometric and algebraic conditions for sphere intersections. It handles three subcases based on the
    values of the coefficients F31, F21, and F11.

    Parameters
    ----------
    Fs : list
        List of polynomial coefficients F, F_2, F10, F11, etc., that define the conditions for vertex calculation
    r0 : float
        The radius component used in the calculation of polynomial coefficients
    l0 : array
        The original location of the sphere center used to adjust the calculated vertices back to the actual position

    Returns
    -------
    list
        A list of vertices, each represented as a tuple containing the vertex coordinates and a corresponding radius

    Notes
    -----
    This function handles three subcases within Case 2 based on the values of the coefficients F31, F21, and F11.
    It checks for real roots of the polynomial defined by the coefficients a, b, and c, calculated from the input Fs.
    """

    # Unpack the F values for easier handling
    F, F_2, F10, F11, F20, F21, F30, F31 = Fs

    # Calculate polynomial coefficients for the vertex equation
    a = F_2 + F11 ** 2 + F21 ** 2 - F31 ** 2
    b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * r0)
    c = F10 ** 2 + F20 ** 2 - (F30 + F * r0)**2

    # Initialize list to store vertices
    verts = []

    # Calculate the discriminant to check for real roots
    disc = b**2 - 4 * a * c

    # Proceed only if the discriminant is non-negative, indicating real roots
    if disc < 0:
        return  # Exit if roots would be complex

    # Solve the quadratic equation to find potential z, y, or x values (roots)
    rts = [root for root in roots([a, b, c]) if isreal(root)]

    # Handle different subcases based on the non-zero coefficient
    if F31 != 0:  # Case 2.1
        for z in rts:
            x = F10 / F + z * F11 / F
            y = F20 / F + z * F21 / F
            R = F30 / F + z * F31 / F
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])

    elif F21 != 0:  # Case 2.2
        for y in rts:
            x = F10 / F + y * F11 / F
            R = F20 / F + y * F21 / F
            z = F30 / F + y * F31 / F
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])

    elif F11 != 0:  # Case 2.3
        for x in rts:
            R = F10 / F + x * F11 / F
            y = F20 / F + x * F21 / F
            z = F30 / F + x * F31 / F
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])

    return verts


def filter_vert_locrads(verts, rs):
    """
    Filter and sort vertices based on their radii.

    This function ensures that encapsulating vertices are removed and the smallest vertex is listed first.
    It is typically used in geometric processing where vertices represent possible solutions that need to be
    validated based on physical or geometric constraints.

    Parameters
    ----------
    verts : list
        List of vertices to be filtered and sorted
    rs : array-like
        Array of radii corresponding to the vertices

    Returns
    -------
    list
        Filtered and sorted list of vertices, with encapsulating vertices removed and the smallest vertex first

    Notes
    -----
    The function removes vertices that are completely encapsulated by other vertices and sorts the remaining
    vertices by their radii in ascending order.
    """

    # Initialize return variables for location and radii
    loc, rad, loc2, rad2 = None, None, None, None

    # Handle the case with a single vertex
    if len(verts) == 1:
        loc, rad = verts[0][0], verts[0][1]  # Directly assign the location and radius

    # Handle the case with two vertices
    elif len(verts) == 2:
        max_ball_rad = max(rs)  # Determine the largest radius from the original spheres for comparison

        # Ensure the smaller vertex is listed first based on absolute radius
        if abs(verts[0][1]) > abs(verts[1][1]):
            verts[0], verts[1] = verts[1], verts[0]  # Swap if necessary

        # Extract locations and radii after potential swap
        locs, rads = [verts[0][0], verts[1][0]], [verts[0][1], verts[1][1]]

        # Validate vertices based on their radii
        if rads[0] < 0 or rads[1] < 0:
            # Check first vertex
            if rads[0] > 0 or abs(rads[0]) < max_ball_rad:
                loc, rad = locs[0], rads[0]  # Assign if valid
                # Check second vertex
                if rads[1] > 0 or abs(rads[1]) < max_ball_rad:
                    loc2, rad2 = locs[1], rads[1]  # Assign if also valid
            # If first vertex wasn't valid, check the second
            elif rads[1] > 0 or abs(rads[1]) < max_ball_rad:
                loc, rad = locs[1], rads[1]
        else:
            # If both radii are positive, assign both with the smallest listed first
            loc, loc2, rad, rad2 = locs[0], locs[1], rads[0], rads[1]

    # Return sorted and validated vertex information
    return loc, loc2, rad, rad2


def calc_vert(locs, rads):
    """
    Calculate the geometrically inscribed or additively weighted vertex between four spheres.

    This function calculates the vertex between four spheres based on their locations and radii.
    It handles different geometrical configurations by applying appropriate computational cases.

    Parameters
    ----------
    locs : list of arrays
        A list of coordinates for the centers of the four spheres
    rads : list of floats
        A list of radii for the four spheres

    Returns
    -------
    tuple
        Returns a tuple of vertices locations and their respective radii calculated for the inscribed sphere

    Notes
    -----
    The function uses different computational cases based on the geometric configuration of the spheres.
    It first calculates vertex coefficients using calc_vert_abcfs, then determines the appropriate case
    based on matrix ranks and coefficient conditions.
    """

    # Attempt to calculate vertex coefficients using a JIT-accelerated function
    Fs, abcdfs, rs, l0 = calc_vert_abcfs(array(locs), array(rads))

    # Initialize matrix ranks needed for determining the computational case
    m_rank, f_rank = 3, 3  # Default ranks if F != 0

    # Adjust ranks based on the coefficients if the first F coefficient is zero
    if Fs[0] == 0:
        my_mtx = [abcdfs]  # Construct matrix from coefficients
        m_rank = linalg.matrix_rank(array(my_mtx[:-1]))  # Calculate rank excluding the last element
        if m_rank != 3:
            f_rank = linalg.matrix_rank(array(my_mtx))  # Calculate full matrix rank

    # Initialize a list to store vertices
    verts = []

    # Case 1: Standard case where the first coefficient of F is non-zero
    if Fs[0] != 0:
        verts = calc_vert_case_1(Fs, l0, rs[0])  # Calculate vertices for case 1
        if verts is not None:
            verts = [[vert[:3], vert[3]] for vert in verts]  # Format vertices
        else:
            verts = []  # Reset vertices if None found

    # Case 2: Special case based on matrix ranks and specific coefficient conditions
    elif abcdfs[0][0] * abcdfs[1][1] - abcdfs[0][1] * abcdfs[1][0] != 0 and m_rank == 3 and f_rank == 3 and Fs[0] > 0:
        verts = calc_vert_case_2(Fs, rs[0], l0)  # Calculate vertices for case 2

    # Filter and sort vertices to find the appropriate geometric solution
    loc, loc2, rad, rad2 = filter_vert_locrads(verts, rs)

    # Return the first and second vertex locations and their corresponding radii
    return loc, rad, loc2, rad2


def calc_flat_vert(locs, rads, power=False):
    """
    Calculate the vertex at the intersection of planes bisecting line segments between balls.

    This function calculates the vertex at the intersection of the planes bisecting the line segments
    between the first ball and each of the other three balls. This vertex represents the geometric
    solution where these planes intersect, which can be interpreted as the center of a circumsphere
    in Delaunay triangulation or as a power center in Laguerre (power) diagrams.

    Parameters
    ----------
    locs : list of arrays
        Coordinates of the centers of the four balls
    rads : list of floats
        Radii of the four balls
    power : bool, optional
        If True, calculates using the power diagram method, which accounts for the radii differences;
        otherwise, uses the Delaunay triangulation method

    Returns
    -------
    tuple
        A tuple containing the coordinates of the calculated vertex and its associated radius or power distance

    Notes
    -----
    The function first sorts the balls by their radii to consistently define the plane equations.
    Plane equations are derived from the midpoints of the line segments (or their power equivalents).
    The intersection of these planes is found by solving a linear system derived from the plane equations.
    """
    # Sort the locations and radii in terms of radii and retun a list of loc, rad tuples
    ball_rads = [(x, _) for _, x in sorted(zip(rads, locs), key=lambda pair: pair[0])]
    # Get the plane equations
    coeffs = []
    # Go through the balls to make the planes
    for an in ball_rads[1:]:
        # Get the point between the balls
        r = array(an[0]) - array(ball_rads[0][0])
        norm = linalg.norm(r)
        rn = r / norm
        if power:
            d0 = 0.5 * (norm ** 2 + ball_rads[0][1] ** 2 - an[1] ** 2) / norm
            center = ball_rads[0][0] + d0 * rn
        else:
            center = 0.5 * r + array(ball_rads[0][0])
        coeffs.append(rn.tolist() + [dot(rn, center)])
    # Unpack the coefficients for the planes
    a1, b1, c1, d1 = coeffs[0]
    a2, b2, c2, d2 = coeffs[1]
    a3, b3, c3, d3 = coeffs[2]
    # Find the discriminant?
    disc = c1 * b2 * a3 - b1 * c2 * a3 - c1 * a2 * b3 + a1 * c2 * b3 + b1 * a2 * c3 - a1 * b2 * c3
    # Calculate the intersection numerators
    x_numerator = d1 * c2 * b3 - c1 * d2 * b3 - d1 * b2 * c3 + b1 * d2 * c3 + c1 * b2 * d3 - b1 * c2 * d3
    y_numerator = - d1 * c2 * a3 + c1 * d2 * a3 + d1 * a2 * c3 - a1 * d2 * c3 - c1 * a2 * d3 + a1 * c2 * d3
    z_numerator = d1 * b2 * a3 - b1 * d2 * a3 - d1 * a2 * b3 + a1 * d2 * b3 + b1 * a2 * d3 - a1 * b2 * d3
    # Calculate the location of the intersection of the planes
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    # Get the radius
    if power:
        # Calculate the power distance between the vertex and an arbitrary ball
        rad = calc_dist(array([x, y, z]), array(ball_rads[0][0])) ** 2 - ball_rads[0][1] ** 2
    else:
        # Calculate the distance between the vertex and an arbitrary ball
        rad = calc_dist(array([x, y, z]), array(ball_rads[0][0]))
    return [x, y, z], rad


@jit(nopython=True)
def verify_aw(loc, rad, test_locs, test_rads):
    """
    Verify if a sphere does not encroach within the radius of any other spheres.

    This function determines if a given sphere (defined by its center 'loc' and radius 'rad') does not
    encroach within the radius of any other spheres in a given list, adjusted for their radii. This
    function is tailored for applications in atomic weaving network calculations and is optimized
    with Numba for high performance.

    Parameters
    ----------
    loc : numpy.ndarray
        The center of the sphere to verify
    rad : float
        The radius of the sphere to verify
    test_locs : numpy.ndarray
        An array of centers of other spheres to check against
    test_rads : numpy.ndarray or list
        An array or list of radii corresponding to the centers in 'test_locs'

    Returns
    -------
    bool
        Returns True if the sphere does not encroach within the radii of any other spheres in the list,
        otherwise False

    Notes
    -----
    The function checks for non-encroachment by ensuring the distance between 'loc' and each 'test_loc'
    minus the respective 'test_rad' is greater than 'rad'. This method is suited for verifying spatial
    configurations in models where spheres represent atoms or particles and their interactions or
    separations are critical. The function is optimized with Numba's nopython mode, which ensures it
    is compiled to machine code for faster execution.
    """

    # Iterate through each sphere in the list to check for encroachment
    for i, b_loc in enumerate(test_locs):
        b_rad = test_rads[i]  # Get the radius for the current sphere
        # Calculate if the center 'loc' encroaches within the adjusted radius of any sphere
        if calc_dist_numba(b_loc, loc) - b_rad < rad:
            return False  # Encroachment detected, return False

    return True  # No encroachments found, return True


@jit(nopython=True)
def verify_prm(loc, rad, test_locs):
    """
    Verify if a location does not fall within the power radius of any other locations.

    This function verifies if a given location 'loc' with a specified 'rad' does not fall within the
    power radius of any other locations in 'test_locs'. This function is intended for use in solving
    the power diagram of a system of balls and is optimized with Numba for high performance.

    Parameters
    ----------
    loc : numpy.ndarray
        The center of the location to be verified
    rad : float
        The radius within which no other centers should exist
    test_locs : numpy.ndarray
        An array of centers to check against

    Returns
    -------
    bool
        Returns True if no other centers are within the radius 'rad' from 'loc', otherwise returns False

    Notes
    -----
    This function iterates over each center in 'test_locs' to check if 'loc' is outside the specified
    'rad'. The function is optimized with Numba's nopython mode for faster execution.
    """

    # Iterate through each location in the list to check for proximity
    for i, b_loc in enumerate(test_locs):
        # Check if the distance between 'loc' and the current location 'b_loc' is less than 'rad'
        if calc_dist_numba(b_loc, loc) < rad:
            return False  # If within radius, return False indicating an invalid position

    return True  # If no overlaps are found, return True indicating a valid position


@jit(nopython=True)
def verify_pow(loc, rad, test_locs, test_rads):
    """
    Verify if a sphere does not overlap with any other spheres.

    This function determines if a given sphere (defined by its center 'loc' and 'rad') does not overlap
    with any other spheres in a given list. This function is optimized for use in power diagram
    computations and is compiled with Numba for performance.

    Parameters
    ----------
    loc : numpy.ndarray
        The center of the sphere to verify
    rad : float
        The radius of the sphere to verify
    test_locs : numpy.ndarray
        An array of centers of other spheres to check against
    test_rads : numpy.ndarray or list
        An array or list of radii corresponding to the centers in 'test_locs'

    Returns
    -------
    bool
        Returns True if the sphere does not overlap with any other spheres in the list, otherwise False

    Notes
    -----
    The function iterates over a list of spheres defined by 'test_locs' and 'test_rads'. It checks for
    non-overlapping conditions by comparing the squared distance between sphere centers to the squared
    sum of radii. This function is suitable for high-performance computational needs due to its
    compilation with Numba, which translates Python functions to optimized machine code at runtime.
    """

    # Iterate through each sphere in the list to check for overlaps
    for i, b_loc in enumerate(test_locs):
        b_rad = test_rads[i]  # Get the radius for the current sphere
        # Calculate the squared distance and compare it to the squared sum of radii
        if calc_dist_numba(b_loc, loc) ** 2 - b_rad ** 2 < rad:
            return False  # Overlap detected, return False

    return True  # No overlaps found, return True


def verify_site(loc, rad, test_locs, test_rads, net_type='aw'):
    """
    Check if a site (vertex) overlaps with other sites.

    This function checks if a given site (vertex) specified by its location and radius overlaps with
    other sites. It can adapt to different network types by selecting appropriate verification methods.

    Parameters
    ----------
    loc : array-like or numpy.ndarray
        The location of the vertex as coordinates
    rad : float
        The radius of the vertex
    test_locs : list or numpy.ndarray
        A collection of locations for other sites to test against
    test_rads : list or numpy.ndarray
        Radii corresponding to each location in test_locs
    net_type : str, optional
        Type of network to use for verification. Options include 'aw' for atomic weaving,
        'prm' for probabilistic roadmaps, and 'pow' for power diagrams

    Returns
    -------
    bool
        True if the site is verified (does not overlap or meets criteria specific to the network type),
        False otherwise

    Notes
    -----
    The function first ensures that the 'loc' parameter is a numpy.ndarray. It then delegates the
    actual overlap checking to specific functions based on the network type: 'aw' for atomic weaving
    networks, 'prm' for probabilistic roadmaps, and 'pow' for power diagrams. These specific functions
    check for conditions like overlapping or proximity based on network-specific rules.
    """

    # Ensure the location is in numpy array format for consistency in mathematical operations
    if not isinstance(loc, np.ndarray):
        loc = np.array(loc)

    # Call the appropriate function to verify the site based on the type of network
    if net_type == 'aw':
        return verify_aw(loc, rad, test_locs, test_rads)
    elif net_type == 'prm':
        return verify_prm(loc, rad, test_locs)
    elif net_type == 'pow':
        return verify_pow(loc, rad, test_locs, test_rads)


# def calc_vert1(locs, rads):
#     """
#     Calculates the Voronoi vertex from the centers and radii of four spheres.
#     The Voronoi vertex is the center of a sphere tangential to all four spheres.
#     """
#     # Sorting the locations and radii in terms of radii and return a list of loc, rad tuples
#     ball_rads = sorted(zip(locs, rads), key=lambda pair: pair[1])
#
#     # Define the symbols
#     x, y, z, R = symbols('x y z R')
#
#     # Matrix M and vector f initialization
#     M = []
#     f = []
#
#     # Iterate over the sorted ball radii and locations to populate matrix M and vector f
#     for i, (loc, rad) in enumerate(ball_rads[1:], start=2):
#         xi, yi, zi = symbols(f'x{i} y{i} z{i}')
#         Ri = symbols(f'R{i}')
#
#         # Compute coefficients for the system of linear equations
#         a = 2 * (loc[0] - ball_rads[0][0][0])
#         b = 2 * (loc[1] - ball_rads[0][0][1])
#         c = 2 * (loc[2] - ball_rads[0][0][2])
#         d = 2 * (rad - rads[0])
#
#         # Compute constants for the system
#         fi = rad**2 - rads[0]**2 + sum((loc[j] - ball_rads[0][0][j])**2 for j in range(3))
#
#         # Append to the system matrix and vector
#         M.append([a, b, c, d])
#         f.append(fi)
#
#     # Convert to SymPy Matrix
#     M = Matrix(M)
#     f = Matrix(f)
#
#     # Solve the linear system
#     system = M.row_join(f)
#     loc_ = list(linsolve(system, (x, y, z, R)))
#     print(loc_[0])
#
#     loc = np.array([float(loc_[0][i] + locs[0][i]) for i in range(3)])
#     print(loc, locs[0])
#     R = calc_dist(loc, locs[0]) - rads[0]
#
#     return loc, R, None, None
