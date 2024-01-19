from numba import jit
import numpy as np
from System.sys_funcs.calcs.calcs import calc_tri, calc_dist


@jit(nopython=True)
def calc_surf_func(l0, r0, l1, r1):
    """
    Calculates the coefficients for the surface between the two atoms
    :return: Returns a function for the hyperboloid between the atoms
    """
    # Check the smaller atom is first
    if r1 < r0:
        l0, r0, l1, r1 = l1, r1, l0, r0
    # Grab the centers of the spheres
    x1, y1, z1 = l0
    x2, y2, z2 = l1
    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = r0 - r1
    K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
    d = [x1 - x2, y1 - y2, z1 - z2]
    J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2
    # Instantiate/reset the hyperboloid coefficient vector lists
    ABC, DEF, GHI = [], [], []
    # Calculate hyperboloid coefficients
    for i in range(3):
        ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
        DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
        GHI.append(-8 * R ** 2 * l0[i] - 4 * K * d[i])
    # Return the function coefficients
    return ABC + DEF + GHI + [J] + [K] + d


def calc_surf_sa(edges, com, tris, points, flat):
    """
    Calculates the surface area of the input surface
    :param edges: Edges for the surface if the surface is flat
    :param com: Center of mass point for the surface used for calculations
    :param tris: Triangles for summing
    :param points: Points for the surface
    :param flat: Whether the surface is flat or not
    :return: Surface area of the surface
    """
    # Create the surface area variable
    sa = 0
    if flat:
        for edge in edges:
            for i in range(len(edge) - 1):
                tri = np.array([edge[i], edge[i + 1], com])
                sa += calc_tri(tri)
    # Go through the triangles in the surface
    else:
        for tri in tris:
            tri1 = np.array([points[tri[_]] for _ in range(3)])
            sa += calc_tri(tri1)
    # Return the surface area
    return sa


def calc_surf_tri_dists(points, tris, loc):
    """
    Calculate the distances between each triangle and the provided location
    :param points: points from the surface
    :param tris: triangles from the
    :param loc: location for the distance calculations
    :return: List of distance corresponding to the triangles on a surface
    """
    # Set up the distances
    dists = []
    tri_dists = []
    max_dist, min_dist = 0, np.inf
    # Provide value for the points
    for point in points:
        # Calculate the distance
        my_dist = calc_dist(point, loc)
        dists.append(my_dist)
        # Record the minimum and maximum distances
        if my_dist < min_dist:
            min_dist = my_dist
        elif my_dist > max_dist:
            max_dist = my_dist
    # Go through the triangles in the surface
    for i in range(len(tris)):
        # Find the maximum distance point of the triangles
        tri_dists.append(max([dists[_] for _ in tris[i]]))
    # Normalize the tri_dists
    return [(_ - min_dist) / (max_dist - min_dist) for _ in tri_dists]


def calc_surf_point_curv(func, point):
    # Get the function coefficients
    A, B, C, D, E, F, G, H, I, J = func[:10]
    # Label the points
    x, y, z = point
    # Get the gradient of the surface at the point
    delf = [2 * A * x + D * y + F * z + G, 2 * B * y + D * x + E * z + H, 2 * C * z + E * y + F * x + I]
    # Calculate the norm of the gradient
    denominator = np.linalg.norm(delf) ** 4
    # Calculate the determinant of the hessian matrix and the gradient matrix
    numerator = np.linalg.det([[2 * A, D, F, delf[0]], [D, 2 * B, E, delf[1]], [F, E, 2 * C, delf[2]],
                               delf + [0]])
    # Get the curvature
    return - numerator / denominator


def calc_surf_tri_curvs(func, points, tris, max_curv):
    """
    Calculates the curvature of the triangles
    :param calc_max_curv:
    :param func: Surface function
    :param points: points for the surface
    :param tris: triangles in the surface
    :return: A list of curvature values for the triangles
    """
    curvs = []
    min_curv = np.inf
    # If the surface normal is within the surface,
    # Get the curvature for each point
    for point in points:
        curv = calc_surf_point_curv(func, point)
        if curv < min_curv:
            min_curv = curv
        elif curv > max_curv:
            max_curv = curv

        curvs.append(curv)
    # Set up the tri_curvs list
    tri_curvs = []
    # Go through the curvature values for each point
    for i in range(len(tris)):
        # Get the triangle
        tri = tris[i]
        # Get the curvatures
        curv_val = sum([curvs[_] for _ in tri])/3
        # Add the curve value to the surface's list of curvatures
        tri_curvs.append(curv_val)
    # Return the values
    return tri_curvs, max_curv


def calc_surf_tri_ins_out(surf):
    """
    Calculates whether the triangles in the surface are inside the overlapping atoms or not
    :param surf: Surface for calculations
    :return: List of bools for if the triangles in the surface are inside or out
    """
    # Set up a list of tracking
    inside_array = []
    # Go through the points in the surface
    for point in surf.points:
        # Calculate the distance between the point and the atom
        my_dist = calc_dist(point, surf.atoms[0].loc)
        # Check if the triangle is inside or not
        if my_dist < surf.atoms[0].rad:
            inside_array.append(True)
        else:
            inside_array.append(False)
    # Now add the triangles
    surf.tri_ins_out = []
    # Color the tris
    for tri in surf.tris:
        if inside_array[tri[0]] and inside_array[tri[1]] and inside_array[tri[2]]:
            surf.tri_ins_out.append(0.25)
        else:
            surf.tri_ins_out.append(0.75)
