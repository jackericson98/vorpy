"""Calculator functions"""
import matplotlib.tri as mtri
from System.sys_calcs import *


# Calculate surface simplices function.
def find_simps(points, atoms):
    # Get the atoms
    a0, a1 = atoms
    # Find the normal to the surface and the magnitude
    r10 = np.array(a0.loc) - np.array(a1.loc)
    d = np.linalg.norm(r10)
    r10_hat = r10/d
    # Get the distance between the surfaces
    ds = d - (a0.rad + a1.rad)
    # Get the center of the surface
    c = np.array(a1.loc) + (0.5 * ds + a0.rad) * r10_hat
    # Move all surf points toward the origin via center point
    for i in range(len(points)):
        points[i] = points[i] - c
    # Calculate the angles to rotate the center point around
    nps = rotate_points(c, points)
    # Get the 2d version of the points
    nps = np.array(nps)
    nps2d = nps[:, 0], nps[:, 1]
    # Get the Delaunay tesselation
    tri = mtri.Triangulation(nps2d[0], nps2d[1])
    # Filter out any connections between the vertices
    for i in range(len(points)):
        points[i] = points[i] + c
    return tri


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
        # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
        roots = np.roots([a, b, c])
        if calc_dist(vi, a1.loc) > a1.rad:
            mag = min(abs(roots))
        else:
            if calc_dist(a0.loc, a1.loc) > a0.rad:
                mag = - abs(min(roots))
            else:
                mag = - min(abs(roots))
        return vi + mag * vn


# Find next point function. Finds the
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
    return calc_surf_point(surf, pc)
