import numpy as np
from numpy import array, dot, linalg
from numba import jit
import warnings
from vorpy.src.calculations.calcs import calc_dist, calc_dist_numba
from contextlib import contextmanager
warnings.simplefilter('error', RuntimeWarning)


@contextmanager
def _numeric_guard():
    """
    Locally escalate dangerous FP conditions and numpy RuntimeWarnings to exceptions.
    - Only active inside the `with` block.
    - Restores the previous NumPy error state on exit.
    """
    old = np.seterr()  # save previous behavior
    try:
        # Raise on divide-by-zero / invalid / overflow; underflow is usually benign
        np.seterr(divide='raise', invalid='raise', over='raise', under='ignore')
        with warnings.catch_warnings():
            warnings.filterwarnings('error', category=RuntimeWarning)
            yield
    finally:
        np.seterr(**old)


def _real_roots_quadratic(a, b, c, tol=1e-12):
    """
    Return the real roots of:

        a*x**2 + b*x + c = 0

    Uses the analytic quadratic formula instead of np.roots(), avoiding
    a full eigenvalue calculation for every candidate Voronoi vertex.
    """

    a = float(a)
    b = float(b)
    c = float(c)

    coefficient_scale = max(abs(a), abs(b), abs(c), 1.0)
    coefficient_tol = tol * coefficient_scale

    # Degenerate quadratic: solve as a linear equation.
    if abs(a) <= coefficient_tol:
        if abs(b) <= coefficient_tol:
            return []

        return [-c / b]

    discriminant = b * b - 4.0 * a * c

    discriminant_scale = max(
        abs(b * b),
        abs(4.0 * a * c),
        1.0,
    )
    discriminant_tol = tol * discriminant_scale

    # Definitely complex roots.
    if discriminant < -discriminant_tol:
        return []

    # Treat a small negative value caused by floating-point error as zero.
    if discriminant < 0.0:
        discriminant = 0.0

    sqrt_discriminant = np.sqrt(discriminant)

    # Stable quadratic formula. This avoids cancellation when b and
    # sqrt(discriminant) have similar magnitudes.
    if b >= 0.0:
        q = -0.5 * (b + sqrt_discriminant)
    else:
        q = -0.5 * (b - sqrt_discriminant)

    if abs(q) <= coefficient_tol:
        roots = [-b / (2.0 * a)]
    else:
        roots = [
            q / a,
            c / q,
        ]

    # Remove duplicate roots from a zero discriminant.
    if len(roots) == 2 and np.isclose(
        roots[0],
        roots[1],
        rtol=tol,
        atol=coefficient_tol,
    ):
        return [roots[0]]

    return roots


def _safe_div(num, den, name="denominator", eps=1e-15):
    """Divide while rejecting zero, near-zero, and non-finite denominators."""
    den = float(den)

    if not np.isfinite(den) or abs(den) <= eps:
        raise ValueError(f"{name} is zero or non-finite (|{den}| <= {eps}).")

    return num / den


@jit(nopython=True, cache=True)
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

    # Build quadratic in a numerically safer form
    a = (F11**2 + F21**2 + F31**2) - F_2
    b = 2.0 * ((F10*F11 + F20*F21 + F30*F31) - r0*F_2)
    c = (F10**2 + F20**2 + F30**2) - (r0**2)*F_2

    Rs = _real_roots_quadratic(a, b, c)
    if not Rs:
        return []

    verts = []
    for R in Rs:
        # use safe division by F (already checked F != 0 in the outer dispatcher)
        x = _safe_div(F10 + R*F11, F, name="F") + l0[0]
        y = _safe_div(F20 + R*F21, F, name="F") + l0[1]
        z = _safe_div(F30 + R*F31, F, name="F") + l0[2]
        verts.append([x, y, z, R])
    return verts


@jit(nopython=True, cache=True)
def calc_vert_case_1_numba(Fs, l0, r0, tol=1e-12):
    """
    Numba-compiled implementation of the standard AW vertex solution.

    This function is numerically equivalent to ``calc_vert_case_1`` but uses
    fixed-size NumPy arrays and an analytic quadratic solution so the hot
    four-sphere calculation can execute efficiently in compiled code.

    Parameters
    ----------
    Fs : numpy.ndarray
        Determinant coefficients returned by ``calc_vert_abcfs``.
    l0 : numpy.ndarray
        Reference sphere location.
    r0 : float
        Reference sphere radius.
    tol : float, optional
        Relative tolerance used for degenerate coefficients and duplicate roots.

    Returns
    -------
    tuple
        ``(verts, n_roots)`` where ``verts`` is a preallocated ``(2, 4)``
        array containing candidate ``x, y, z, radius`` values and ``n_roots``
        gives the number of valid rows.
    """
    F, F_2, F10, F11, F20, F21, F30, F31 = Fs
    # Construct the quadratic in the vertex radius.
    a = (F11**2 + F21**2 + F31**2) - F_2
    b = 2.0 * ((F10*F11 + F20*F21 + F30*F31) - r0*F_2)
    c = (F10**2 + F20**2 + F30**2) - (r0**2)*F_2
    # Handle degenerate linear and standard quadratic cases.
    coefficient_scale = max(abs(a), abs(b), abs(c), 1.0)
    coefficient_tol = tol * coefficient_scale
    # Convert each radius root into its corresponding global vertex location.
    roots_out = np.empty(2, dtype=np.float64)
    n_roots = 0

    if abs(a) <= coefficient_tol:
        if abs(b) > coefficient_tol:
            roots_out[0] = -c / b
            n_roots = 1
    else:
        discriminant = b*b - 4.0*a*c
        discriminant_scale = max(abs(b*b), abs(4.0*a*c), 1.0)
        discriminant_tol = tol * discriminant_scale

        if discriminant >= -discriminant_tol:
            if discriminant < 0.0:
                discriminant = 0.0

            sqrt_discriminant = np.sqrt(discriminant)

            if b >= 0.0:
                q = -0.5 * (b + sqrt_discriminant)
            else:
                q = -0.5 * (b - sqrt_discriminant)

            if abs(q) <= coefficient_tol:
                roots_out[0] = -b / (2.0*a)
                n_roots = 1
            else:
                r1 = q / a
                r2 = c / q

                roots_out[0] = r1
                n_roots = 1

                if not np.isclose(r1, r2, rtol=tol, atol=coefficient_tol):
                    roots_out[1] = r2
                    n_roots = 2

    verts = np.empty((2, 4), dtype=np.float64)

    for i in range(n_roots):
        R = roots_out[i]
        verts[i, 0] = (F10 + R*F11) / F + l0[0]
        verts[i, 1] = (F20 + R*F21) / F + l0[1]
        verts[i, 2] = (F30 + R*F31) / F + l0[2]
        verts[i, 3] = R

    return verts, n_roots


def calc_vert_case_2(Fs, r0, l0, tol=1e-12):
    """
    Calculate the legacy Case 2 AW vertex solutions.

    Case 2 parameterizes the vertex equations by one coordinate (z, y, or x,
    depending on which radius coefficient is non-zero) and solves the resulting
    quadratic analytically. The helper intentionally remains Python-level:
    its legacy return contract is a heterogeneous nested list of
    ``[[x, y, z], radius]`` records, which Numba cannot type in nopython mode.
    The active ``calc_vert`` hot path uses the compiled Case 1 solver.

    Notes
    -----
    This routine requires ``F != 0`` because reconstruction of x, y, z and R
    divides by F.  Therefore it is not a valid fallback for determinant-zero
    systems.  ``calc_vert`` uses the standard Case 1 solver for F != 0 and
    rejects determinant-zero systems that do not have a unique solution under
    the current formulation.
    """
    F, F_2, F10, F11, F20, F21, F30, F31 = Fs

    # Case 2 reconstruction divides by F throughout.
    coefficient_scale = max(
        abs(F), abs(F_2), abs(F10), abs(F11),
        abs(F20), abs(F21), abs(F30), abs(F31), 1.0
    )
    coefficient_tol = tol * coefficient_scale
    if not np.isfinite(F) or abs(F) <= coefficient_tol:
        return []

    # Quadratic in the free coordinate.
    a = F_2 + F11 ** 2 + F21 ** 2 - F31 ** 2
    b = 2.0 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * r0)
    c = F10 ** 2 + F20 ** 2 - (F30 + F * r0) ** 2

    quad_scale = max(abs(a), abs(b), abs(c), 1.0)
    quad_tol = tol * quad_scale

    # Fixed-size root storage keeps the function straightforward for Numba.
    rts = np.empty(2, dtype=np.float64)
    n_roots = 0

    # Degenerate quadratic -> linear equation.
    if abs(a) <= quad_tol:
        if abs(b) <= quad_tol:
            return []
        rts[0] = -c / b
        n_roots = 1
    else:
        disc = b * b - 4.0 * a * c
        disc_scale = max(abs(b * b), abs(4.0 * a * c), 1.0)
        disc_tol = tol * disc_scale

        if disc < -disc_tol:
            return []

        # Small negative discriminants are roundoff at a repeated root.
        if disc < 0.0:
            disc = 0.0

        sqrt_disc = np.sqrt(disc)

        # Numerically stable quadratic formula.
        if b >= 0.0:
            q = -0.5 * (b + sqrt_disc)
        else:
            q = -0.5 * (b - sqrt_disc)

        if abs(q) <= quad_tol:
            rts[0] = -b / (2.0 * a)
            n_roots = 1
        else:
            r1 = q / a
            r2 = c / q
            rts[0] = r1
            n_roots = 1

            if not np.isclose(r1, r2, rtol=tol, atol=quad_tol):
                rts[1] = r2
                n_roots = 2

    verts = []

    # Reconstruct the full vertex from the free coordinate.
    if abs(F31) > coefficient_tol:  # Case 2.1: z is free
        for i in range(n_roots):
            z = rts[i]
            x = (F10 + z * F11) / F
            y = (F20 + z * F21) / F
            R = (F30 + z * F31) / F
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])

    elif abs(F21) > coefficient_tol:  # Case 2.2: y is free
        for i in range(n_roots):
            y = rts[i]
            x = (F10 + y * F11) / F
            R = (F20 + y * F21) / F
            z = (F30 + y * F31) / F
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])

    elif abs(F11) > coefficient_tol:  # Case 2.3: x is free
        for i in range(n_roots):
            x = rts[i]
            R = (F10 + x * F11) / F
            y = (F20 + x * F21) / F
            z = (F30 + x * F31) / F
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])

    return verts


def filter_vert_locrads(verts, rs):
    """
    Filter candidate AW vertices and order the surviving solutions.

    Mathematical solutions with radii smaller than the negative radius of the
    smallest defining ball are physically invalid and are removed. Remaining
    solutions are ordered by absolute radius so the geometrically nearest
    solution is returned first.

    Parameters
    ----------
    verts : list
        Candidate vertices in the form ``[[x, y, z], radius]``.
    rs : array-like
        Radii of the four defining balls.

    Returns
    -------
    tuple
        ``(loc, loc2, rad, rad2)`` containing the primary and optional
        secondary vertex locations and radii. Missing solutions are returned
        as ``None``.
    """
    loc, rad, loc2, rad2 = None, None, None, None

    if verts is None or len(verts) == 0:
        return loc, loc2, rad, rad2

    # A negative vertex cannot contract farther than the smallest defining ball.
    min_allowed_rad = -min(rs)
    verts = [vert for vert in verts if vert[1] >= min_allowed_rad]

    if len(verts) == 0:
        return loc, loc2, rad, rad2

    # Return the solution nearest to zero radius first.
    verts.sort(key=lambda vert: abs(vert[1]))

    loc, rad = verts[0][0], verts[0][1]

    if len(verts) >= 2:
        loc2, rad2 = verts[1][0], verts[1][1]

    return loc, loc2, rad, rad2


def calc_vert(locs, rads):
    """
    Calculate the additively weighted Voronoi vertex defined by four spheres.

    The four sphere centers and radii are converted into the linearized AW
    coefficient system. The standard non-degenerate case is solved by the
    Numba-compiled Case 1 implementation. Determinant-zero configurations do
    not have a unique solution under this formulation and return no candidate.
    Candidate solutions are then filtered according to the physically valid
    radius range.

    Parameters
    ----------
    locs : array-like
        Coordinates of the four defining sphere centers.
    rads : array-like
        Radii of the four defining spheres.

    Returns
    -------
    tuple
        ``(loc, rad, loc2, rad2)`` where ``loc`` and ``rad`` describe the
        primary AW vertex and ``loc2`` and ``rad2`` describe an optional
        secondary solution for doublet configurations.
    """
    # Normalize input data for the compiled coefficient calculation.
    locs_array, rads_array = array(locs), array(rads)

    # Build the determinant coefficients for the four-sphere system.
    Fs, abcdfs, rs, l0 = calc_vert_abcfs(locs_array, rads_array)

    # The standard AW formulation requires a non-zero determinant F.
    # The former Case 2 fallback was unreachable (it was guarded by both
    # F == 0 and F > 0) and, more importantly, Case 2 itself divides by F.
    # Preserve the existing determinant-zero behavior explicitly: no candidate
    # vertex is produced when the defining system has no unique solution.
    verts = []

    if Fs[0] != 0:
        numba_verts, numba_n = calc_vert_case_1_numba(Fs, l0, rs[0])
        verts = [[numba_verts[i, :3].tolist(), numba_verts[i, 3]]
                 for i in range(numba_n)]

    # Remove physically invalid roots and order any remaining solutions.
    loc, loc2, rad, rad2 = filter_vert_locrads(verts, rs)

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
    disc = c1 * b2 * a3 - b1 * c2 * a3 - c1 * a2 * b3 + a1 * c2 * b3 + b1 * a2 * c3 - a1 * b2 * c3

    eps = 1e-15
    if not np.isfinite(disc) or abs(disc) <= eps:
        # Singular/near-singular planes → no unique intersection
        return None, None

    # Calculate the intersection numerators
    x_numerator = d1 * c2 * b3 - c1 * d2 * b3 - d1 * b2 * c3 + b1 * d2 * c3 + c1 * b2 * d3 - b1 * c2 * d3
    y_numerator = - d1 * c2 * a3 + c1 * d2 * a3 + d1 * a2 * c3 - a1 * d2 * c3 - c1 * a2 * d3 + a1 * c2 * d3
    z_numerator = d1 * b2 * a3 - b1 * d2 * a3 - d1 * a2 * b3 + a1 * d2 * b3 + b1 * a2 * d3 - a1 * b2 * d3

    # Calculate the location of the intersection of the planes
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc

    if power:
        rad = calc_dist(np.array([x, y, z]), np.array(ball_rads[0][0])) ** 2 - ball_rads[0][1] ** 2
    else:
        rad = calc_dist(np.array([x, y, z]), np.array(ball_rads[0][0]))
    return [x, y, z], rad


@jit(nopython=True, cache=True)
def verify_aw(loc, rad, test_locs, test_rads, skip_ndx=-1):
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
        if i == skip_ndx:
            continue

        b_rad = test_rads[i]
        if calc_dist_numba(b_loc, loc) - b_rad < rad:
            return False

    return True  # No encroachments found, return True


@jit(nopython=True, cache=True)
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


@jit(nopython=True, cache=True)
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


def verify_site(loc, rad, test_locs, test_rads, net_type='aw', skip_ndx=-1):
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
        return verify_aw(loc, rad, test_locs, test_rads, skip_ndx=skip_ndx)
    elif net_type == 'prm':
        return verify_prm(loc, rad, test_locs)
    elif net_type == 'pow':
        return verify_pow(loc, rad, test_locs, test_rads)
