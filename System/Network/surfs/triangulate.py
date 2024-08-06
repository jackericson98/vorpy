from System.sys_funcs.calcs.calcs import rotate_points, calc_com, calc_angle_jit
from scipy.spatial import Delaunay
import numpy as np


def calc_tri_circ(points):
    """
    Finds the circumference of the circumscribed circle for the triangle
    :return: Circumference of the circle circumscribing the triangle
    """
    # Get the points of the triangle
    pa, pb, pc = points
    # Get their squared differences
    a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
    b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
    c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
    # Calculate the area
    s = (a + b + c) / 2.0
    area = np.sqrt(max(s * (s - a) * (s - b) * (s - c), 0))
    # If the triangle is open or has repeat lines or broken in any way, we don't want it so return inf
    if area == 0:
        return np.inf
    # Else, return the circumference of the circle that the triangle inscribes
    circum_r = a * b * c / (4.0 * area)
    return circum_r


def tri_within(perimeter, loc, norm, flat_points=None, my_tri=None, point=None):
    """
    Checks to see if a triangle lies within the perimeter of a surface
    :return: Bool
    """
    # Get the perimeter of the translated and rotated surface
    perimeter = np.array([perimeter[i] - loc for i in range(len(perimeter))])
    # Rotate the point
    perimeter = [point[:2] for point in rotate_points(norm, points=perimeter)]
    if len(perimeter) == 0:
        return False
    center = calc_com(points=perimeter)
    # If we are given a triangle determine the center of mass and use that point
    if my_tri is not None:
        # Copy the triangle, retrieve its points and calculate the center of mass
        tri = my_tri.copy()
        # Get the triangles points
        points = [flat_points[tri[i]] for i in range(len(tri))]
        # Calculate the triangle's center of mass
        point = np.array(calc_com(points=points))
    else:
        # Move the point
        point = point - loc
        # Rotate the point
        new_point = rotate_points(norm, np.array([point]))[0]
        # Get the 2d version
        point = np.array(new_point[:2])
    # Get the projected point
    proj_vec = np.array(center) - np.array(point)
    proj_point = np.array(point) + np.array(proj_vec)

    # Reset the number of intersections
    xings = 0
    # Go through each line segments around the perimeter
    for i in range(len(perimeter)):
        # Get the line segment's points
        p1 = np.array(perimeter[i])
        p2 = np.array(perimeter[(i + 1) % len(perimeter)])
        # Get the angles
        theta = calc_angle_jit(point, p1, p2)
        theta_n = calc_angle_jit(point, p1, proj_point)
        theta_n1 = calc_angle_jit(point, p2, proj_point)
        # If we have a crossing
        if 0 < theta_n < theta and theta_n1 < theta:
            xings += 1

    # If we have an even number of intersections
    if xings % 2 == 0:
        return False
    else:
        return True


def find_simps(points, loc, norm):
    """
    Transforms and rotates surface points to xy-plane and returns the Delaunay simplices
    :param norm:
    :param loc:
    :param points:
    :return: None
    """
    # Copy the surface points
    points = points.copy()
    # Move all surf points toward the origin via center point
    for i in range(len(points)):
        points[i] = np.array(points[i]) - np.array(loc)
    # Calculate the angles to rotate the center point around
    nps = rotate_points(np.array(norm), np.array(points))
    # Get the 2d version of the points and their Delaunay tesselation
    flat_points = [_[:2] for _ in nps]
    tris = Delaunay(flat_points)
    # Add the flat points to the surface's list of flat points
    tris = tris.simplices.tolist()
    return tris, flat_points


def filter_tris(tris, flat_points, res, perimeter, loc, norm, filter_hard):
    """
    Goes through the triangles on the surface measuring the circumference & testing if inside
    :return:
    """
    # Check to see if the surface is flat or not
    # Set up a list of indices to remove for the triangles
    remove_ndxs = []
    # Go through the triangles in the surface
    for i in range(len(tris)):
        # Grab the triangle and calculate its circumference
        tri = tris[i]
        circ = calc_tri_circ(points=[flat_points[_] for _ in tri])
        # If the circumference of the triangle is less than x times the min_dist check to see if tri is within
        if circ > 3 * res and not tri_within(perimeter=perimeter, loc=loc, norm=norm, flat_points=flat_points, my_tri=tri):
            remove_ndxs.append(i)
        elif filter_hard and not tri_within(perimeter, loc, norm, flat_points=flat_points, my_tri=tri):
            remove_ndxs.append(i)
    # Remove the outer triangles
    remove_ndxs.sort()
    for i in range(len(remove_ndxs)):
        tris.pop(remove_ndxs[-(i + 1)])
    # Return the surface's triangles
    return tris
