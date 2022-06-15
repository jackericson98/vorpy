from objects import *
import numpy as np


# Bisector function. Creates a bisector surface between 2 atoms
def calc_surface(a1, a2):

    # Set sphere1 to the smaller of the two atoms
    if a1.rad <= a2.rad:
        a1, a2 = a1, a2
    else:
        a2, a1 = a1, a2

    # Grab the centers of our spheres
    C1 = a1.loc
    C2 = a2.loc

    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = a1.rad - a2.rad
    K = (C2[0] ** 2 - C1[0] ** 2) + (C2[1] ** 2 - C1[1] ** 2) + (C2[2] ** 2 - C1[2] ** 2) - R ** 2
    d = C1[0] - C2[0], C1[1] - C2[1], C1[2] - C2[2]
    J = 4 * R ** 2 * (C1[0] ** 2 + C1[1] ** 2 + C1[2] ** 2) - K ** 2

    # Instantiate/reset the hyperboloid coefficient vector lists
    ABC, DEF, GHI = [], [], []
    # Calculate hyperboloid coefficients
    for i in range(3):
        ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
        DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
        GHI.append(-8 * R ** 2 * C1[i] - 4 * K * d[i])

    mySurf = Surface(ABC + DEF + GHI + [J] + [K] + [d], [a1, a2])

    return mySurf


# Make meshes function. Takes in a base sphere (atom) and a set of other spheres to make meshes.
def make_meshes(atom, neighbors):
    # Reset the mesh list
    meshes = []
    for neighbor in neighbors:
        # Set sphere1 to the smaller of the two atoms
        a1 = atom
        if atom.rad > neighbor.rad:
            a1 = neighbor
        # Set the coefficients
        A, B, C, D, E, F, G, H, I, J, K, d = calc_surface(atom, neighbor)
        # Using the coefficients, find the lengths of all sample rays in sphere1, such that they intersect the bisector
        mesh = []
        # Go through each ray sampled from our sphere
        for ray in a1.rays:
            # Get the normal direction and location of the ray
            nx, ny, nz = ray.dir
            Ox, Oy, Oz = ray.loc
            # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
            a = A * nx ** 2 + B * ny ** 2 + C * nz ** 2 + D * nx * ny + E * ny * nz + F * nz * nx
            b = 2 * A * nx * Ox + 2 * B * ny * Oy + 2 * C * nz * Oz + D * (nx * Oy + ny * Ox) + E * (
                    ny * Oz + nz * Oy) + F * (nz * Ox + nx * Oz) + G * nx + H * ny + I * nz
            c = A * Ox ** 2 + B * Oy ** 2 + C * Oz ** 2 + D * Ox * Oy + E * Oy * Oz + F * Oz * Ox + \
                G * Ox + H * Oy + I * Oz + J
            # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
            # and add that point to our surface list of points
            if b ** 2 - 4 * a * c > 0:
                mag = min(abs(np.roots([a, b, c])))
                check = np.dot(ray.dir, d)
                if check < 0:
                    mesh.append([ray.loc[0] + ray.dir[0] * mag, ray.loc[1] + ray.dir[1] * mag,
                                 ray.loc[2] + ray.dir[2] * mag])

        # Add the mesh to the neighbor and
        meshes.append(mesh)

    return meshes

