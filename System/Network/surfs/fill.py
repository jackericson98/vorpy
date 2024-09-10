import numba.core.errors
from System.Network.surfs.triangulate import is_within
from System.sys_funcs.calcs.calcs import calc_angle_jit, calc_angle, calc_dist
import numpy as np
from numba import jit


@jit(nopython=True)
def calc_surf_point_abcs(locs, point, func):
    # Set up the unit vector
    vi = point - locs[0]
    vn = vi / np.linalg.norm(vi)
    # Set the atom's location as the root
    vi = locs[0]

    # Solve the surface function's equation for the vector through the given point from the atom's location:

    # Get the a/b/c values for the point(s) that lies on the surface and along the vector from a0 to the given point
    a = func[0] * vn[0] ** 2 + func[1] * vn[1] ** 2 + func[2] * vn[2] ** 2 + func[3] * vn[0] * vn[1] + func[4] * vn[1] \
        * vn[2] + func[5] * vn[2] * vn[0]
    b = 2 * func[0] * vn[0] * vi[0] + 2 * func[1] * vn[1] * vi[1] + 2 * func[2] * vn[2] * vi[2] + func[3] \
        * (vn[0] * vi[1] + vn[1] * vi[0]) + func[4] * (vn[1] * vi[2] + vn[2] * vi[1]) + func[5] \
        * (vn[2] * vi[0] + vn[0] * vi[2]) + func[6] * vn[0] + func[7] * vn[1] + func[8] * vn[2]
    c = func[0] * vi[0] ** 2 + func[1] * vi[1] ** 2 + func[2] * vi[2] ** 2 + func[3] * vi[0] * vi[1] + func[4] * vi[1] \
        * vi[2] + func[5] * vi[2] * vi[0] + func[6] * vi[0] + func[7] * vi[1] + func[8] * vi[2] + func[9]
    return vi, vn, a, b, c


def calc_surf_point(locs, point, func):
    """
    Projects a vector through the reference point and the smaller surface atom's center onto the surface
    :param func: Implicit function for the hyperboloid surface between the atoms
    :param locs: Smaller atom's location used for projection onto the surface
    :param point: Reference point to be projected through
    :return: The point on the surface
    """
    vi, vn, a, b, c = calc_surf_point_abcs(np.array(locs), np.array(point), np.array(func))

    # Choose the 181L root:

    # Check that the discriminant of the solution to at^2 + bt + c = 0, is positive
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        # Calculate the roots of the factoring equation
        roots = np.roots([a, b, c])
        # If one root exists return it
        if len(roots) == 1:
            return vi + roots[0] * vn
        # If the smallest root is negative (i.e. incorrect) return the other root
        if min(roots) < 0:
            return locs[0] + vn * max(roots)
        # Otherwise, return the smaller of the two
        return locs[0] + min(roots) * vn


def find_next_point(locs, func, pn_1, end, d_theta):
    """
    Finds the next point along the given path by projecting a reference point onto the surface
    :param func: Surface's function coefficients
    :param a0_loc: Surface's smaller atom's location
    :param pn_1: Previous path point
    :param end: End path point being moved towards by a d_theta amount
    :param d_theta: Angular increment to move towards the end point
    :return: The new point on the surface
    """
    # Get the first angle
    a0 = d_theta
    # Get the smaller atom's location
    pa = locs[0]
    # Get the location of point b
    pb = np.array(pn_1)
    # Get the distance between pb and pa
    s2 = np.sqrt(sum(np.square(np.array(pa) - np.array(pb))))
    # Get the angle between pa, pb and pv1
    try:
        a1 = calc_angle_jit(pb, pa, end)
    except numba.core.errors.TypingError:
        a1 = calc_angle(pb, pa, end)
    # Get the last angle
    a2 = np.pi - a0 - a1
    # Find a using the law of sines
    s0 = np.sin(a0) * s2 / np.sin(a2)
    # Find the direction of the vector pointing from the previous point to the end point
    rn = end - pb
    # Normalize this vector. Try to supress warnings
    try:
        rn_hat = rn / np.linalg.norm(rn)
    except RuntimeWarning:
        return
    # Find the next projection point by adding the vector with 'a' magnitude and rn_hat direction
    pc = pb + rn_hat * s0
    # Calculate where the point intercepts the surface and return it
    return calc_surf_point(locs, point=pc, func=func)


def fill_mesh(locs, rads, func, surf_loc, surf_norm, perimeter, com, res, flat, check=False):
    """
    Works inward from a set of perimeter points toward a center point filling in equally spaced points
    :param surf: Surface object being filled
    :return: None
    """
    # For each path toward the center of the surface, set up a path list.
    paths = [[_] for _ in perimeter]
    spoints = perimeter[:]

    # Check to see if the atoms have equal radii
    if rads[0] == rads[1] or flat:
        # Go through the paths
        for i in range(len(paths)):
            # Get the
            r = np.array(com) - np.array(paths[i][0])
            if r.all() == 0:
                continue
            norm = np.linalg.norm(r)
            rn = r / norm
            num_steps = max(int(norm / res), 2)
            step = norm / num_steps
            spoints += [paths[i][0] + rn * j * step for j in range(1, num_steps + 1)]
        return spoints
    # Grab the smallest of the 2 surface atoms' location
    pa = locs[0]
    # Get the angles between the edge points and the end points
    dists = []
    angs = []
    for i in range(len(paths)):
        # Calculate the angle for each path
        try:
            angs.append(calc_angle_jit(pa, paths[i][0], com))
        except numba.core.errors.TypingError:
            angs.append(calc_angle(pa, paths[i][0], com))
        # Get the dists from the com to the path
        dists.append(np.sqrt(sum(np.square(np.array(paths[i][0]) - np.array(com)))))
    # Get the maximum path
    max_path_ndx = angs.index(max(angs))
    max_path = paths[max_path_ndx][0]
    # Decide how many rings based off of the ellipticity and density
    num_rings = max(int(np.sqrt(sum(np.square(np.array(max_path) - np.array(com)))) / res), 2)
    # Get the incremental angle increases
    dthetas = [angs[i] / num_rings for i in range(len(angs))]
    # Set the pn_1 point to infinity
    pn_1 = [np.inf, np.inf, np.inf]
    num_paths = len(paths)
    skips = [[False] for _ in paths]
    # Go through ring by ring
    for j in range(num_rings):
        # Go through each of the remaining paths
        i = 0
        # Keep going through the points until the tracker is out
        while i < num_paths:
            # Find the number of previously skipped points
            num_prev_skips = 1
            while True:
                if skips[i][- num_prev_skips]:
                    num_prev_skips += 1
                else:
                    break
            # Get the next point along the path
            pn = find_next_point(locs, func, paths[i][-1], com, num_prev_skips * dthetas[i])
            # Check for edges that start by going outside
            if j == 0 and pn is not None:
                if not is_within(perimeter, pn, surf_loc, surf_norm):
                    paths.pop(i)
                    dthetas.pop(i)
                    skips.pop(i)
                    num_paths -= 1
                    continue
            # Check to see of the new point is too close to the previous point and the path has to end
            if pn is None or calc_dist(pn, pn_1) < 0.5 * res:
                # Add the path to the surfaces points and remove it from the paths list
                spoints += paths.pop(i)[1:]
                dthetas.pop(i)
                skips.pop(i)
                num_paths -= 1
            else:
                # Set the pn_1 to pn and add it to the path
                if not calc_dist(pn, paths[i][-1]) < 0.5 * res:
                    pn_1 = pn
                    paths[i].append(pn)
                    skips[i].append(False)
                else:
                    skips[i].append(True)
                i += 1
    # Add the remaining paths to the surface excluding the first point in the path (i.e. the edge point)
    for path in paths:
        spoints += path[1:]
    # Return the points
    return spoints
