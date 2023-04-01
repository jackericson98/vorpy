from System.sys_funcs.calcs import calc_angle, calc_com, rotate_points, calc_surf_func, calc_surf_curv
from System.Network.net_funcs.build_edge import build_edge
from scipy.spatial import Delaunay
import numpy as np


################################################# Find Surface Points  #################################################


def calc_surf_point(point, func, a0_loc):
    """
    Projects a vector through the reference point and the smaller surface atom's center onto the surface
    :param func: Implicit function for the hyperboloid surface betwen the atoms
    :param a0_loc: Smaller atom's location used for projection onto the surface
    :param point: Reference point to be projected through
    :return: The point on the surface
    """
    # Set up the unit vector
    vi = np.array(point) - np.array(a0_loc)
    vn = vi / np.linalg.norm(vi)
    # Set the atom's location as the root
    vi = a0_loc

    # Solve the surface function's equation for the vector through the given point from the atom's location:

    # Get the a/b/c values for the point(s) that lies on the surface and along the vector from a0 to the given point
    a = func[0] * vn[0] ** 2 + func[1] * vn[1] ** 2 + func[2] * vn[2] ** 2 + func[3] * vn[0] * vn[1] + func[4] * vn[1] * vn[2] + func[
        5] \
        * vn[2] * vn[0]
    b = 2 * func[0] * vn[0] * vi[0] + 2 * func[1] * vn[1] * vi[1] + 2 * func[2] * vn[2] * vi[2] + func[3] \
        * (vn[0] * vi[1] + vn[1] * vi[0]) + func[4] * (vn[1] * vi[2] + vn[2] * vi[1]) + func[5] \
        * (vn[2] * vi[0] + vn[0] * vi[2]) + func[6] * vn[0] + func[7] * vn[1] + func[8] * vn[2]
    c = func[0] * vi[0] ** 2 + func[1] * vi[1] ** 2 + func[2] * vi[2] ** 2 + func[3] * vi[0] * vi[1] + func[4] * vi[1] * vi[2] + \
        func[5] * vi[2] * vi[0] + func[6] * vi[0] + func[7] * vi[1] + func[8] * vi[2] + func[9]

    # Choose the correct root:

    # Check that the discriminant of the solution to at^2 + bt + c = 0, is positive
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        # Calculate the roots of the factoring equation
        roots = np.roots([a, b, c])
        # If one root exists return it
        if len(roots) == 1:
            return vi + roots[0] * vn
        # If the smallest root is negative (i.e. incorrect) return the other root
        if min(roots) < 0:
            return a0_loc + vn * max(roots)
        # Otherwise, return the smaller of the two
        return a0_loc + min(roots) * vn


def find_next_point(surf, pn_1, end, d_theta):
    """
    Finds the next point along the given path by projecting a reference point onto the surface
    :param surf: Surface object holding the paths
    :param pn_1: Previous path point
    :param end: End path point being moved towards by a d_theta amount
    :param d_theta: Angular increment to move towards the end point
    :return: The new point on the surface
    """
    # Get the A angle
    A = d_theta
    # Get the smaller atom's location
    pa = surf.atoms[0].loc
    # Get the location of point b
    pb = np.array(pn_1)
    # Get the distance between pb and pa
    c = np.sqrt(sum(np.square(np.array(pa) - np.array(pb))))
    # Get the angle between pa, pb and pv1
    B = calc_angle(pb, pa, end)
    # Get the last angle
    C = np.pi - A - B
    # Find a using the law of sines
    a = np.sin(A) * c / np.sin(C)
    # Find the direction of the vector pointing from the previous point to the end point
    rn = end - pb
    # Normalize this vector. Try to supress warnings
    try:
        rn_hat = rn / np.linalg.norm(rn)
    except RuntimeWarning:
        return
    # Find the next projection point by adding the vector with 'a' magnitude and rn_hat direction
    pc = pb + rn_hat * a
    # Calculate where the point intercepts the surface and return it
    return calc_surf_point(point=pc, func=surf.func, a0_loc=surf.atoms[0].loc)


def build_perimeter(surf, net_type='vor'):
    """
    Sorts the edges of the surface to create a list of points in order around the perimeter
    :param net_type:
    :param surf: Surface object for perimeter building
    :return: None
    """
    # Randomly select the first edge
    e0 = surf.edges[0]
    # Add the first edge's vertex location and set of points to the perimeter points list
    perimeter = e0.points.copy()
    # Check to see if the edges have been built yet
    for edge in surf.edges:
        if edge.points is None:
            build_edge(edge, res=surf.res, straight=surf.flat)
    # Make a copy of the edges to organize excluding the first edge
    edges = surf.edges[1:].copy()

    # Keep looping while we haven't gone through the edges
    while edges:

        # Set the max distance to infinity, the index for the intended edge to None and the reverse bool to False
        d, ndx, reverse = np.inf, None, False

        # Go through each of the remaining edges in the list
        for i in range(len(edges)):
            # Calculate the distance between the most recently recorded point and the first/last points in the edge
            d0 = np.sqrt(sum(np.square(np.array(perimeter[-1]) - np.array(edges[i].points[0]))))
            d1 = np.sqrt(sum(np.square(np.array(perimeter[-1]) - np.array(edges[i].points[-1]))))
            # If the first edge point is closer to the last perimeter point and the last isn't closer add that edge
            if d0 < d and d0 < d1:
                d, ndx, reverse = d0, i, False
            # Otherwise, if the last edge point is the closest add the edge in reverse
            elif d1 < d:
                d, ndx, reverse = d1, i, True
        # Pull the edge from the list of edges
        my_edge = edges.pop(ndx)
        # Add the edge's point in the right order and then add the correct vertex
        if not reverse:  # In order
            perimeter += my_edge.points
        else:  # Reverse order
            perimeter += my_edge.points[::-1]

    # Add the perimeter points to the whole set of points
    surf.points += perimeter
    # Get the perimeter flat points
    surf.pflat_points = perimeter.copy()
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    d = np.sqrt(sum(np.square(np.array(a0.loc) - np.array(a1.loc))))
    # Get the center of the surface
    if surf.norm is None:
        r = np.array(a1.loc) - np.array(a0.loc)
        surf.norm = r / np.linalg.norm(r)
    if net_type == 'vor':
        surf.loc = np.array(a0.loc) + (a0.rad + 0.5 * (d - (a0.rad + a1.rad))) * surf.norm
    elif net_type == 'del':
        surf.loc = np.array(a0.loc) + 0.5 * d * surf.norm
    elif net_type == 'pow':
        surf.loc = np.array(a0.loc) + 0.5 * (surf.norm ** 2 + a0.rad ** 2 - a1.rad ** 2) / surf.norm
    for i in range(len(surf.pflat_points)):
        # Move the points
        surf.pflat_points[i] = surf.pflat_points[i] - surf.loc
    # Rotate the point
    surf.pflat_points = rotate_points(surf.norm, surf.pflat_points)
    # Get the 2d version
    surf.pflat_points = [point[:2] for point in surf.pflat_points]
    return perimeter


def get_com(surf, net_type='vor'):
    """
    Finds the center of mass of a surface's perimeter points
    :param surf: Surface object holding the perimeter points
    :return: Center of mass
    """
    if surf.flat or net_type in {'del', 'pow'}:
        return calc_com(points=surf.perimeter)
        # If the surface is flat, the center of mass will not need to be projected
    if tri_within(surf, point=surf.loc):
        return surf.loc
    # First try the center of mass of the 3d points projected onto the surface
    my_com = calc_surf_point(point=calc_com(points=surf.perimeter[::5]), func=surf.func, a0_loc=surf.atoms[0].loc)
    if my_com is not None and tri_within(surf, point=my_com) and surf.atoms[0].rad != surf.atoms[1].rad:
        return my_com
    # If nothing else set the center of mass to the first point in the perimeter
    surf.filter_hard = True
    return surf.perimeter[len(surf.perimeter)//2]


def fill_mesh(surf):
    """
    Works inward from a set of perimeter points toward a center point filling in equally spaced points
    :param surf: Surface object being filled
    :return: None
    """
    # Check to see that the surface has perimeter points
    if len(surf.perimeter) == 0:
        surf.perimeter = build_perimeter(surf)
    # Get the resolution
    res = surf.res
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    # Get the center of mass
    surf.com = get_com(surf, surf.net.type)
    com = surf.com
    # For each path toward the center of the surface, set up a path list.
    paths = [[surf.perimeter[i]] for i in range(len(surf.perimeter))]
    # Check to see if the atoms have equal radii
    if a0.rad == a1.rad or surf.flat:
        # Make sure the surface is flat
        surf.flat = True
        # Go through the paths
        for i in range(len(paths)):
            # Get the
            r = np.array(com) - np.array(paths[i][0])
            if r.all() == 0:
                continue
            norm = np.linalg.norm(r)
            rn = r / norm
            num_steps = max(int(norm / surf.res), 2)
            step = norm / num_steps
            surf.points += [paths[i][0] + rn * j * step for j in range(1, num_steps + 1)]
        return
    # Grab the smallest of the 2 surface atoms' location
    pa = surf.atoms[0].loc
    # Get the angles between the edge points and the end points
    dists = []
    angs = []
    for i in range(len(paths)):
        # Calculate the angle for each path
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
    # Go through ring by ring
    for j in range(num_rings):
        # Go through each of the remaining paths
        i = 0
        # Keep going through the points until the tracker is out
        while i < num_paths:
            # Get the next point along the path
            pn = find_next_point(surf, paths[i][-1], com, dthetas[i])
            # Check for edges that start by going outside
            if j == 0 and pn is not None and not tri_within(surf, point=pn):
                paths.pop(i)
                dthetas.pop(i)
                num_paths -= 1
                continue
            # Check to see if the point is outside the network's box
            if pn is not None and np.array([surf.net.box[0][i] <= pn[i] <= surf.net.box[1][i] for i in range(3)]).all():
                surf.in_box = False
            # Check to see of the new point is too close to the previous point and the path has to end
            if pn is None or (np.sqrt(sum(np.square(np.array(pn) - np.array(pn_1)))) < 0.5 * res and not
               np.sqrt(sum(np.square(np.array(paths[i - 1][-1]) - np.array(pn)))) > res):
                # Add the path to the surfaces points and remove it from the paths list
                surf.points += paths.pop(i)[1:]
                dthetas.pop(i)
                num_paths -= 1
            else:
                # Set the pn_1 to pn and add it to the path
                pn_1 = pn
                paths[i].append(pn)
                i += 1
    # Add the remaining paths to the surface excluding the first point in the path (i.e. the edge point)
    for path in paths:
        surf.points += path[1:]


############################################## Triangulate Surface Points  #############################################

def tri_within(surf, my_tri=None, point=None):
    """
    Checks to see if a triangle lies within the perimeter of a surface
    :param surf: Surface object to check against
    :param my_tri: Triangle to test insideness
    :param point: point to test insideness
    :return: Bool
    """
    # Get the perimeter of the translated and rotated surface
    perimeter = [surf.perimeter[i] - surf.loc for i in range(len(surf.perimeter))]
    # Rotate the point
    perimeter = [point[:2] for point in rotate_points(surf.norm, perimeter)]
    if len(perimeter) == 0:
        return False
    center = calc_com(points=perimeter)
    # If we are given a triangle determine the center of mass and use that point
    if my_tri is not None:
        # Copy the triangle, retrieve its points and calculate the center of mass
        tri = my_tri.copy()
        # Get the triangles points
        points = [surf.flat_points[tri[i]] for i in range(len(tri))]
        # Calculate the triangle's center of mass
        point = calc_com(points=points)
    else:
        # Move the point
        point = point - surf.loc
        # Rotate the point
        new_point = rotate_points(surf.norm, [point])[0]
        # Get the 2d version
        point = new_point[:2]
    # Get the projected point
    proj_vec = np.array(center) - np.array(point)
    proj_point = np.array(point) + np.array(proj_vec)

    # Reset the number of intersections
    xings = 0
    # Go through each line segments around the perimeter
    for i in range(len(perimeter)):
        # Get the line segment's points
        p1 = perimeter[i]
        p2 = perimeter[(i + 1) % len(perimeter)]
        # Get the angles
        theta = calc_angle(point, p1, p2)
        theta_n = calc_angle(point, p1, proj_point)
        theta_n1 = calc_angle(point, p2, proj_point)
        # If we have a crossing
        if 0 < theta_n < theta and theta_n1 < theta:
            xings += 1

    # If we have an even number of intersections
    if xings % 2 == 0:
        return False
    else:
        return True


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


def find_simps(surf):
    """
    Transforms and rotates surface points to xy-plane and returns the Delaunay simplices
    :param surf: Surface object holding the points
    :return: None
    """
    # Copy the surface points
    points = surf.points.copy()
    # Move all surf points toward the origin via center point
    for i in range(len(points)):
        points[i] = np.array(points[i]) - np.array(surf.loc)

    # Calculate the angles to rotate the center point around
    nps = rotate_points(surf.norm, points)

    # Get the 2d version of the points and their Delaunay tesselation
    surf.flat_points = [_[:2] for _ in nps]

    tris = Delaunay(surf.flat_points)
    # Add the flat points to the surface's list of flat points
    surf.tris = tris.simplices.tolist()


def filter_tris(surf):
    """
    Goes through the triangles on the surface measuring the circumference & testing if inside
    :param surf: Surface object holding the triangles for filtration
    :return:
    """
    # Check to see if the surface is flat or not
    # Set up a list of indices to remove for the triangles
    remove_ndxs = []
    # Go through the triangles in the surface
    for i in range(len(surf.tris)):
        # Grab the triangle and calculate its circumference
        tri = surf.tris[i]
        circ = calc_tri_circ(points=[surf.flat_points[_] for _ in tri])
        # If the circumference of the triangle is less than x times the min_dist check to see if tri is within
        if circ > 3 * surf.res and not tri_within(surf, tri):
            remove_ndxs.append(i)
        elif surf.filter_hard and not tri_within(surf, tri):
            remove_ndxs.append(i)
    # Remove the outer triangles
    remove_ndxs.sort()
    for i in range(len(remove_ndxs)):
        surf.tris.pop(remove_ndxs[-(i + 1)])


# Build method. Makes the mesh for the surface and calculates the simplices between them
def build_surf(surf, res=None, color=False):
    """
    Main build method for constructing surfaces
    :param surf:
    :param res: Specifies the resolution the surface is to be constructed with
    :param color: Bool to color the surf or not
    :return: The surfaces points and triangles are filled
    """
    if surf.net.type in {'pow', 'del'} or surf.atoms[0].element == surf.atoms[1].element:
        surf.flat = True
    # Set the resolution value that the surface is built with
    if res is None:
        res = surf.net.surf_res
    surf.res = res
    # Check to see if the function or curvature have been calculated and calculate them if not
    if surf.curv is None:
        calc_surf_curv(surf)
    # Reset the surface's list of points to empty list and reset the vertex indices list
    surf.points = []
    # Build the perimeter of the surface
    surf.perimeter = build_perimeter(surf, surf.net.type)
    # Fill the mesh
    fill_mesh(surf)
    # Find the simplices of the surface
    find_simps(surf)
    # If the network type is voronoi the edges could be curved allowing for triangulations outside the edges
    if surf.net.type == 'vor':
        # Filter out the bad triangles
        filter_tris(surf)
