from System.sys_funcs.calcs.calcs import calc_dist
from numpy import array, dot, isreal, linalg, roots
from numba import jit

################################## Calc vert functions #################################################################


# Speed Boost decorator
@jit(nopython=True)
def calc_vert_abcfs_jit(locs, rads):
    """
    Calculate the fs, abcfs, rs and l0 for an additively weighted vertex (set) from the locations and radii of 4 spheres
    """
    # Unpack the atom radii
    r0, r1, r2, r3 = rads
    # Calculate the square of the base sphere's radius
    r0_2 = r0 ** 2
    # Move all atoms' locations to the base sphere for simpler calculation
    l0, l1, l2, l3 = locs[0], locs[1] - locs[0], locs[2] - locs[0], locs[3] - locs[0]
    # Calculate our System of linear equations coefficients
    a1, b1, c1, d1, f1 = 2 * l1[0], 2 * l1[1], 2 * l1[2], 2 * (r1 - r0), r0_2 - r1 ** 2 + l1[0] ** 2 + l1[
        1] ** 2 + l1[2] ** 2
    a2, b2, c2, d2, f2 = 2 * l2[0], 2 * l2[1], 2 * l2[2], 2 * (r2 - r0), r0_2 - r2 ** 2 + l2[0] ** 2 + l2[
        1] ** 2 + l2[2] ** 2
    a3, b3, c3, d3, f3 = 2 * l3[0], 2 * l3[1], 2 * l3[2], 2 * (r3 - r0), r0_2 - r3 ** 2 + l3[0] ** 2 + l3[
        1] ** 2 + l3[2] ** 2
    # Calculate the F values
    F = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1
    F_2 = F ** 2
    F10 = b1 * c2 * f3 - b1 * c3 * f2 - b2 * c1 * f3 + b2 * c3 * f1 + b3 * c1 * f2 - b3 * c2 * f1
    F11 = -b1 * c2 * d3 + b1 * c3 * d2 + b2 * c1 * d3 - b2 * c3 * d1 - b3 * c1 * d2 + b3 * c2 * d1
    F20 = -a1 * c2 * f3 + a1 * c3 * f2 + a2 * c1 * f3 - a2 * c3 * f1 - a3 * c1 * f2 + a3 * c2 * f1
    F21 = a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1
    F30 = a1 * b2 * f3 - a1 * b3 * f2 - a2 * b1 * f3 + a2 * b3 * f1 + a3 * b1 * f2 - a3 * b2 * f1
    F31 = -a1 * b2 * d3 + a1 * b3 * d2 + a2 * b1 * d3 - a2 * b3 * d1 - a3 * b1 * d2 + a3 * b2 * d1
    # Place the F values in an array
    fs = array([F, F_2, F10, F11, F20, F21, F30, F31])
    # Place the abcdfs in an array
    abcdfs = array([[a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f3]])
    # Place the radii in an array
    rs = array([r0, r1, r2, r3])
    # Return values needed for further calculation
    return fs, abcdfs, rs, l0


def calc_vert_abcfs(locs, rads):
    """
    Calculate the fs, abcfs, rs and l0 for an additively weighted vertex (set) from the locations and radii of 4 spheres
    """
    # Unpack the atom radii
    r0, r1, r2, r3 = rads
    # Calculate the square of the base sphere's radius
    r0_2 = r0 ** 2
    # Move all atoms' locations to the base sphere for simpler calculation
    l0, l1, l2, l3 = locs[0], locs[1] - locs[0], locs[2] - locs[0], locs[3] - locs[0]
    # Calculate our System of linear equations coefficients
    a1, b1, c1, d1, f1 = 2 * l1[0], 2 * l1[1], 2 * l1[2], 2 * (r1 - r0), r0_2 - r1 ** 2 + l1[0] ** 2 + l1[
        1] ** 2 + l1[2] ** 2
    a2, b2, c2, d2, f2 = 2 * l2[0], 2 * l2[1], 2 * l2[2], 2 * (r2 - r0), r0_2 - r2 ** 2 + l2[0] ** 2 + l2[
        1] ** 2 + l2[2] ** 2
    a3, b3, c3, d3, f3 = 2 * l3[0], 2 * l3[1], 2 * l3[2], 2 * (r3 - r0), r0_2 - r3 ** 2 + l3[0] ** 2 + l3[
        1] ** 2 + l3[2] ** 2
    # Calculate the F values
    F = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1
    F_2 = F ** 2
    F10 = b1 * c2 * f3 - b1 * c3 * f2 - b2 * c1 * f3 + b2 * c3 * f1 + b3 * c1 * f2 - b3 * c2 * f1
    F11 = -b1 * c2 * d3 + b1 * c3 * d2 + b2 * c1 * d3 - b2 * c3 * d1 - b3 * c1 * d2 + b3 * c2 * d1
    F20 = -a1 * c2 * f3 + a1 * c3 * f2 + a2 * c1 * f3 - a2 * c3 * f1 - a3 * c1 * f2 + a3 * c2 * f1
    F21 = a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1
    F30 = a1 * b2 * f3 - a1 * b3 * f2 - a2 * b1 * f3 + a2 * b3 * f1 + a3 * b1 * f2 - a3 * b2 * f1
    F31 = -a1 * b2 * d3 + a1 * b3 * d2 + a2 * b1 * d3 - a2 * b3 * d1 - a3 * b1 * d2 + a3 * b2 * d1
    # Place the F values in an array
    fs = array([F, F_2, F10, F11, F20, F21, F30, F31])
    # Place the abcdfs in an array
    abcdfs = array([[a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f3]])
    # Place the radii in an array
    rs = array([r0, r1, r2, r3])
    # Return values needed for further calculation
    return fs, abcdfs, rs, l0


# @jit(nopython=True)   <---  Throws error in large systems for the roots call. Negative discriminant, though filtered
def calc_vert_case_1(Fs, l0, r0):
    """
    Calculates the case 1 vertices from Fs, r0 and l0
    """
    # Unwrap the F values
    F, F_2, F10, F11, F20, F21, F30, F31 = Fs
    # Calculate the radius polynomial coefficients
    a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F_2) - 1
    b = 2 * (((F10 * F11 + F20 * F21 + F30 * F31) / F_2) - r0)
    c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F_2) - r0 ** 2
    # Set up the list of vertices (0, 1, or 2)
    verts = []
    # If the discriminant is positive, find the real positive roots of the quadratic
    if -4 * a * c + b ** 2 >= 0:
        Rs = [R for R in roots(array([a, b, c])) if isreal(R)]
    else:
        return
    # Ensure that a valid R is available
    if Rs is not None and len(Rs) > 0:
        # Go through each radius and calculate the vertex
        for R in Rs:
            x = F10 / F + R * F11 / F + l0[0]
            y = F20 / F + R * F21 / F + l0[1]
            z = F30 / F + R * F31 / F + l0[2]
            # Move the vertex back to the actual location of the atoms
            verts.append([x, y, z, R])
    # Return all calculated roots
    return verts


@jit(nopython=True)
def calc_vert_case_2(Fs, r0, l0):
    """
    Calculates the case 2 vertices from Fs, r0 and l0
    """
    # Unpack the F values
    F, F_2, F10, F11, F20, F21, F30, F31 = Fs

    # Calculate the _ polynomial coefficients
    a = F_2 + F11 ** 2 + F21 ** 2 - F31 ** 2
    b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * r0)
    c = F10 ** 2 + F20 ** 2 - (F30 + F * r0)
    # Instantiate the roots and verts lists
    verts = []
    # Check the discriminant
    disc = -4 * a * c + b ** 2
    # If the discriminant is negative escape, we don't want complex roots
    if disc <= 0:
        return

    # Get the roots of the abc values
    rts = [root for root in roots([a, b, c]) if isreal(root)]

    # Case 2 subcases:
    # Case 2.1
    if F31 != 0:
        # Go through each radius and calculate the vertex
        for z in rts:
            x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
            # Move the vertex back to the actual location of the atoms
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
    # Case 2.2
    elif F21 != 0:
        # Go through each radius and calculate the vertex
        for y in rts:
            x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
            # Move the vertex back to the actual location of the atoms
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
    # Case 2.3
    elif F11 != 0:
        # Go through each radius and calculate the vertex
        for x in rts:
            R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
            # Move the vertex back to the actual location of the atoms
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
    return verts


def filter_vert_locrads(verts, rs):
    """
    Filters and sorts additively weighted vertex pairs removing encapsulating verts and returning the smaller vert first
    """
    # Returning these four variables either way, so defining them here makes that easier if no vertices result
    loc, rad, loc2, rad2 = None, None, None, None

    # If one root exists return it
    if len(verts) == 1:
        loc, rad = verts[0][0], verts[0][1]

    # If two roots exist, we need to return the locs and rads with most likely (smaller) vert first
    elif len(verts) == 2:

        # Get the largest atom's radius to rule out negative encapsulating vertices
        max_atom_rad = max(rs)
        # Set the locations and radii, so that the smaller vertex is first
        if abs(verts[0][1]) > abs(verts[1][1]):
            verts[0], verts[1] = verts[1], verts[0]
        # Set the locations and radii variables
        locs, rads = [verts[0][0], verts[1][0]], [verts[0][1], verts[1][1]]

        # If either radii are negative we need to check them
        if rads[0] < 0 or rads[1] < 0:
            # Check the first vert for validity. Positive verts are ok and negative verts under max_atom_rad are ok
            if rads[0] > 0 or abs(rads[0]) < max_atom_rad:
                loc, rad = locs[0], rads[0]
                # Check for a possible second root under the same criteria
                if rads[1] > 0 or abs(rads[1]) < max_atom_rad:
                    loc2, rad2 = locs[1], rads[1]
            # If vert 1 is bad, check vert2 for positiveness or negative relative size
            elif rads[1] > 0 or abs(rads[1]) < max_atom_rad:
                loc, rad = locs[1], rads[1]
        # If both radii are positive we have a doublet. Choose the smaller vertex to be the lead vertex and set loc2
        else:
            loc, loc2, rad, rad2 = locs[0], locs[1], rads[0], rads[1]
    # Return the sorted filtered vert values
    return loc, loc2, rad, rad2


# Calculate vertex function. Takes in 4 atoms, calculates the loc and rad of the inscribed sphere and adds the
def calc_vert(locs, rads):
    """
    Calculates the additively weighted vertex between the locations and radii of four spheres
    """

    # Get the first set of major coefficients
    try:
        Fs, abcdfs, rs, l0 = calc_vert_abcfs_jit(array(locs), array(rads))
    except AssertionError:
        Fs, abcdfs, rs, l0 = calc_vert_abcfs(array(locs), array(rads))

    # Calculate the ranks of the coefficient matrices to determine which vert calculation case to use. F != 0 means 3, 3
    m_rank, f_rank = 3, 3

    # Other rank cases
    if Fs[0] == 0:
        my_mtx = [abcdfs]
        m_rank = linalg.matrix_rank(array(my_mtx[:-1]))
        if m_rank != 3:
            f_rank = linalg.matrix_rank(array(my_mtx))

    # Instantiate the vertices list
    verts = []

    # Case 1:
    if Fs[0] != 0:
        verts = calc_vert_case_1(Fs, l0, rs[0])
        if verts is not None:
            verts = [[vert[:3], vert[3]] for vert in verts]
        else:
            verts = []

    # Case 2:
    elif abcdfs[0][0] * abcdfs[1][1] - abcdfs[0][1] * abcdfs[1][0] != 0 and m_rank == 3 and f_rank == 3 and Fs[0] > 0:
        verts = calc_vert_case_2(Fs, rs[0], l0)

    # Filter out vertices that encapsulate and sort the smaller vertex first
    loc, loc2, rad, rad2 = filter_vert_locrads(verts, rs)

    # Return the locations and radii
    return loc, rad, loc2, rad2


def calc_flat_vert(locs, rads, power=False):
    """
    Calculates the flat vertex between 4 atoms by finding the intersection of the mid-point planes between the first
    atom and the others
    :param locs: List of atomic locations
    :param rads: List of atomic radii
    :param power: Whether to calculate the vertex as a power or Delaunay
    :return: Location and radius tuple
    """
    # Sort the locations and radii in terms of radii and retun a list of loc, rad tuples
    atom_rads = [(x, _) for _, x in sorted(zip(rads, locs), key=lambda pair: pair[0])]
    # Get the plane equations
    coeffs = []
    # Go through the atoms to make the planes
    for an in atom_rads[1:]:
        # Get the point between the atoms
        r = array(an[0]) - array(atom_rads[0][0])
        norm = linalg.norm(r)
        rn = r / norm
        if power:
            d0 = 0.5 * (norm ** 2 + atom_rads[0][1] ** 2 - an[1] ** 2) / norm
            center = atom_rads[0][0] + d0 * rn
        else:
            center = 0.5 * r + array(atom_rads[0][0])
        coeffs.append(rn.tolist() + [dot(rn, center)])
    # Unpack the coefficients for the planes
    x1, y1, z1, c1 = coeffs[0]
    x2, y2, z2, c2 = coeffs[1]
    x3, y3, z3, c3 = coeffs[2]
    # Find the discriminant?
    disc = z1 * y2 * x3 - y1 * z2 * x3 - z1 * x2 * y3 + x1 * z2 * y3 + y1 * x2 * z3 - x1 * y2 * z3
    # Calculate the intersection numerators
    x_numerator = c1 * z2 * y3 - z1 * c2 * y3 - c1 * y2 * z3 + y1 * c2 * z3 + z1 * y2 * c3 - y1 * z2 * c3
    y_numerator = - c1 * z2 * x3 + z1 * c2 * x3 + c1 * x2 * z3 - x1 * c2 * z3 - z1 * x2 * c3 + x1 * z2 * c3
    z_numerator = c1 * y2 * x3 - y1 * c2 * x3 - c1 * x2 * y3 + x1 * c2 * y3 + y1 * x2 * c3 - x1 * y2 * c3
    # Calculate the location of the intersection of the planes
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    # Get the radius
    if power:
        # Calculate the power distance between the vertex and an arbitrary atom
        rad = calc_dist(array([x, y, z]), array(atom_rads[0][0])) ** 2 - atom_rads[0][1] ** 2
    else:
        # Calculate the distance between the vertex and an arbitrary atom
        rad = calc_dist(array([x, y, z]), array(atom_rads[0][0]))
    return [x, y, z], rad
