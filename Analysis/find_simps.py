import numpy as np
from anal_calcs import *


# Rotate points function. Takes in a set of points and a vector and returns a set of rotated points about the origin
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
def find_simps(surf):
    # Get the atoms
    a0, a1 = surf.atoms
    # Grab the points
    points = surf.points + surf.edge_points
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
