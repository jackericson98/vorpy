"""Calculator functions"""
import numpy as np
from System.checks import *
from scipy.spatial import Delaunay as dl


# Calculate distance function. Finds the distance between 2 points
def calc_dist(l1, l2):
    d = np.sqrt((l1[0]-l2[0])**2+(l1[1]-l2[1])**2+(l1[2]-l2[2])**2)
    return d


# Calculate center of mass function. Takes in a set of points and returns the coordinates of the com
def calc_com(points):
    # Set the running sum for the x, y, z values to 0
    xtot, ytot, ztot = 0, 0, 0
    for point in points:
        xtot = xtot + point[0]
        ytot = ytot + point[1]
        ztot = ztot + point[2]
    return xtot/len(points), ytot/len(points), ztot/len(points)


# Calculate angle function. Finds the angle (in rads) between three points. The first being the common point
def calc_angle(p0, p1, p2=None):
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = np.array(p0), np.array(p1)
    else:
        v0, v1 = np.array(p1) - np.array(p0), np.array(p2) - np.array(p0)
    # Get the unit vectors
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
    # Calculate the angle between the two vectors with catches for 180 and 0
    angle = np.arccos(np.clip(np.dot(n0, n1), -1.0, 1.0))
    return angle


# Rotate points function. takes in a vector and points and rotates the points toward the origin away from the vector
def rotate_points(vec, points):
    vx, vy, vz = vec
    if vy == 0 == vx:
        return
    # Get the x and y angles
    theta = np.arctan(vx/vz)
    phi = np.arctan(vy/vz)
    # Get variables for sin(theta), cos(theta), sin(phi), cos(phi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    nps = []
    for p in points:
        # Multiplying the x, y rotation matrices gives the following:
        npx = p[0] * cp + p[2] * sp
        npy = p[0] * st *sp + p[0] * ct -p[2] * st * cp
        npz = - p[0] * ct * sp + p[1] * st + p[2] * ct * cp
        nps.append([npx, npy, npz])
    return nps


# Calculate surface simplices function.
def calc_surf_simps(surf):
    # Get the atoms
    a0, a1 = surf.atoms
    # Grab the points
    points = surf.points
    # Find the normal to the surface and the magnitude
    r10 = np.array(surf.atoms[0].loc) - np.array(surf.atoms[1].loc)
    d = np.linalg.norm(r10)
    r10_hat = r10/d
    # Get the distance between the surfaces
    ds = d - (a0.rad + a1.rad)
    # Get the center of the surface
    c = np.array(a0.loc) + (0.5 * ds + a0.rad) * r10_hat
    # Move all surf points toward the origin via center point
    for i in range(len(points)):
        points[i] = points[i] - c
    # Calculate the angles to rotate the center point around
    nps = rotate_points(c, points)
    if nps is None:
        nps = points
    # Get the 2d version of the points
    np2ds = []
    for p in nps:
        np2ds.append([p[0], p[1]])
    # Get the Delaunay tesselation
    tri = dl(np2ds)
    return tri.simplices


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
    pb = pn_1
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
