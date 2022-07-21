import numpy as np
from Analysis.anal_calcs import *
from Presentation.Visualize.visualize import *
import matplotlib.tri as mtri


# Rotate points function. Takes in a set of points and a vector and returns a set of rotated points about the origin
def rotate_points(vec, points):
    # Get the vx, vy, vz vector components
    vx, vy, vz = vec
    # If x and y are 0 no transform is needed
    if vy == 0 == vx:
        return points
    elif vz == 0 == vy:
        theta = 0
        phi = np.pi/2
    elif vz == 0 == vx:
        phi = 0
        theta = np.pi/2
    elif vz == 0:
        theta = np.pi / 2
        phi = np.arctan(vy/vx)
    else:
        theta = np.arctan(vx / vz)
        phi = np.arctan(vy / vz)
    # Get variables for sin(theta), cos(theta), sin(phi), cos(phi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    nps = []
    for p in points:
        px, py, pz = round(p[0], 7), round(p[1], 7), round(p[2], 7)
        # Multiplying the x, y rotation matrices gives the following:
        npx = px * cp + pz * sp
        npy = px * st * sp + py * ct - pz * st * cp
        npz = - px * ct * sp + py * st + pz * ct * cp
        # Add the new points to the list
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
    return tri
