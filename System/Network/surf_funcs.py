from System.sys_funcs import *


# Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smallest of the surfaces 2 atoms through the point into the surface
def calc_surf_point(surf, point):
    # Grab the function's coefficients
    f = surf.func
    # Get the first atoms in the surfaces list of atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    # Set up the unit vector
    vi = np.array(point) - np.array(a0.loc)
    vn = vi/np.linalg.norm(vi)
    # Find the location on the surface of the atom
    vi = np.array(a0.loc) + vn * a0.rad
    # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
    a = f[0] * vn[0] ** 2 + f[1] * vn[1] ** 2 + f[2] * vn[2] ** 2 + f[3] * vn[0] * vn[1] + f[4] * vn[1] * vn[2] + f[5] \
        * vn[2] * vn[0]
    b = 2 * f[0] * vn[0] * vi[0] + 2 * f[1] * vn[1] * vi[1] + 2 * f[2] * vn[2] * vi[2] + f[3] \
        * (vn[0] * vi[1] + vn[1] * vi[0]) + f[4] * (vn[1] * vi[2] + vn[2] * vi[1]) + f[5] \
        * (vn[2] * vi[0] + vn[0] * vi[2]) + f[6] * vn[0] + f[7] * vn[1] + f[8] * vn[2]
    c = f[0] * vi[0] ** 2 + f[1] * vi[1] ** 2 + f[2] * vi[2] ** 2 + f[3] * vi[0] * vi[1] + f[4] * vi[1] * vi[2] + \
        f[5] * vi[2] * vi[0] + f[6] * vi[0] + f[7] * vi[1] + f[8] * vi[2] + f[9]
    # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
    # and add that point to our surface list of points
    if round(b ** 2 - 4 * a * c, 4) >= 0:
        roots = np.roots([a, b, c])
        # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
        if calc_dist(vi, a1.loc) > a1.rad:
            x = (calc_dist(a0.loc, a1.loc)) / a1.rad
            if calc_angle(a0.loc, a1.loc, vi) - np.pi / 2 > np.pi / 8:
                mag = max(abs(roots))
            else:
                mag = min(abs(roots))
        # If the projection point is within the intersection, the magnitude is negative
        else:
            if calc_angle(a0.loc, a1.loc, vi) > np.pi/4:
                mag = - min(abs(roots))
            else:
                mag = -abs(min(roots))
        return vi + mag * vn


# Find next point function. Finds the next point along the given path by projecting a reference point onto the surface
def find_next_point(pn_1, end, d_theta, surf):
    # Get the A angle
    A = d_theta
    # Get the smaller atom's location
    pa = surf.atoms[0].loc
    # Get the location of point b
    pb = np.array(pn_1)
    # Get the distance between pb and pa
    c = calc_dist(pa, pb)
    # Get the angle between pa, pb and pv1
    B = calc_angle(pb, pa, end)
    # Get the last angle
    C = np.pi - A - B
    # Find a using the law of sines
    a = np.sin(A) * c / np.sin(C)
    # Find the intercept point by adding a to pb
    rn = end - pb
    rn_hat = rn / np.linalg.norm(rn)
    pc = pb + rn_hat * a
    # Calculate where the point intercepts the surface
    sp = calc_surf_point(surf, pc)
    return sp


# Edge trace function. I want to update to put the points in order
def edge_trace1(surf, min_dist=None):
    # Instantiate the edge_points list
    edge_points = []

    return edge_points


