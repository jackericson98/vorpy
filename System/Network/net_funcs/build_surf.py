from System.sys_funcs.calcs import *
import matplotlib.tri as mtri
from Visualize.mpl_visualize import *

################################################# Find Surface Points  #################################################


def calc_surf_point(surf, point):
    """
    Projects a vector through the reference point and the smaller surface atom's center onto the surface
    :param surf: Surface object to project the point onto
    :param point: Reference point to be projected through
    :return: The point on the surface
    """
    # Check to see if the surface's function has been calculated
    if surf.func is None:
        surf.calc_func()
    # Grab the surfaces function and atoms
    f, a0, a1 = surf.func, surf.atoms[0], surf.atoms[1]
    # Set up the unit vector
    vi = np.array(point) - np.array(a0.loc)
    vn = vi / np.linalg.norm(vi)
    # Set the atom's location as the root
    vi = a0.loc

    # Solve the surface function's equation for the vector through the given point from the atom's location:

    # Get the a/b/c values for the point(s) that lies on the surface and along the vector from a0 to the given point
    a = f[0] * vn[0] ** 2 + f[1] * vn[1] ** 2 + f[2] * vn[2] ** 2 + f[3] * vn[0] * vn[1] + f[4] * vn[1] * vn[2] + f[
        5] \
        * vn[2] * vn[0]
    b = 2 * f[0] * vn[0] * vi[0] + 2 * f[1] * vn[1] * vi[1] + 2 * f[2] * vn[2] * vi[2] + f[3] \
        * (vn[0] * vi[1] + vn[1] * vi[0]) + f[4] * (vn[1] * vi[2] + vn[2] * vi[1]) + f[5] \
        * (vn[2] * vi[0] + vn[0] * vi[2]) + f[6] * vn[0] + f[7] * vn[1] + f[8] * vn[2]
    c = f[0] * vi[0] ** 2 + f[1] * vi[1] ** 2 + f[2] * vi[2] ** 2 + f[3] * vi[0] * vi[1] + f[4] * vi[1] * vi[2] + \
        f[5] * vi[2] * vi[0] + f[6] * vi[0] + f[7] * vi[1] + f[8] * vi[2] + f[9]

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
            return a0.loc + vn * max(roots)
        # Otherwise, return the smaller of the two
        return a0.loc + min(roots) * vn


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
    return calc_surf_point(surf, pc)


def build_perimeter(surf):
    """
    Sorts the edges of the surface to create a list of points in order around the perimeter
    :param surf: Surface object for perimeter building
    :return: None
    """
    # Reset the surface's perimeter points list
    surf.perimeter = []
    e0 = surf.edges[0]
    # Check to see if the edges have been built yet
    for edge in surf.edges:
        if edge.points is None:
            edge.build(res=surf.res, straight=surf.flat)
    # Set the edge's reference surface and range (overwrite the ref in place, to ensure a lighter storage with an e0)
    e0.ref = [surf.net.surfs.index(surf), 0, len(e0.points) - 1]
    # Add the first edge's vertex location and set of points to the perimeter points list
    surf.perimeter = e0.points.copy()
    # Make a copy of the edges to organize excluding the first edge
    edges = surf.edges[1:].copy()

    # Keep looping while we haven't gone through the edges
    while edges:

        # Set the max distance to infinity, the index for the intended edge to None and the reverse bool to False
        d, ndx, reverse = np.inf, None, False

        # Go through each of the remaining edges in the list
        for i in range(len(edges)):
            # Calculate the distance between the most recently recorded point and the first/last points in the edge
            d0 = np.sqrt(sum(np.square(np.array(surf.perimeter[-1]) - np.array(edges[i].points[0]))))
            d1 = np.sqrt(sum(np.square(np.array(surf.perimeter[-1]) - np.array(edges[i].points[-1]))))
            # If the first edge point is closer to the last perimeter point and the last isn't closer add that edge
            if d0 < d and d0 < d1:
                d, ndx, reverse = d0, i, False
            # Otherwise, if the last edge point is the closest add the edge in reverse
            elif d1 < d:
                d, ndx, reverse = d1, i, True
        # Pull the edge from the list of edges
        myEdge = edges.pop(ndx)
        # Add the edge's point in the right order and then add the correct vertex
        if not reverse:  # In order
            surf.perimeter += myEdge.points
        else:  # Reverse order
            surf.perimeter += myEdge.points[::-1]
        if myEdge.ref is None:
            myEdge.ref = [surf.net.surfs.index(surf), len(surf.perimeter) - len(e0.points), len(surf.perimeter)]

    # Add the perimeter points to the whole set of points
    surf.points += surf.perimeter
    # Get the perimeter flat points
    surf.pflat_points = surf.perimeter.copy()
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    d = np.sqrt(sum(np.square(np.array(a0.loc) - np.array(a1.loc))))
    # Get the center of the surface
    if surf.rn is None:
        r = np.array(a1.loc) - np.array(a0.loc)
        surf.rn = r / np.linalg.norm(r)
    surf.center = np.array(a0.loc) + (a0.rad + 0.5 * (d - (a0.rad + a1.rad))) * surf.rn
    for i in range(len(surf.pflat_points)):
        # Move the points
        surf.pflat_points[i] = surf.pflat_points[i] - surf.center
    # Rotate the point
    surf.pflat_points = rotate_points(surf.rn, surf.pflat_points)
    # Get the 2d version
    surf.pflat_points = [point[:2] for point in surf.pflat_points]


def get_com(surf):
    """
    Finds the center of mass of a surface's perimeter points
    :param surf: Surface object holding the perimeter points
    :return: Center of mass
    """
    # First try the center of mass of the 3d points projected onto the surface
    my_com = calc_com(points=surf.perimeter[::5])
    # If the surface is flat, the center of mass will not need to be projected
    if not surf.flat:
        my_com = calc_surf_point(surf, my_com)
    if my_com is not None and tri_within(surf, point=my_com) and surf.atoms[0].rad != surf.atoms[1].rad:
        return my_com
    # Get the center of the surface
    if tri_within(surf, point=surf.center):
        return surf.center
    # If nothing else set the center of mass to the first point in the perimeter
    return surf.perimeter[len(surf.perimeter)//2]


def fill_mesh(surf):
    """
    Works inward from a set of perimeter points toward a center point filling in equally spaced points
    :param surf: Surface object being filled
    :return: None
    """
    # Check to see that the surface has perimeter points
    if len(surf.perimeter) == 0:
        build_perimeter(surf)
    # Get the resolution
    res = surf.res
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    # Get the center of mass
    com = get_com(surf)
    # Check to see if the atoms have equal radii
    if a0.rad == a1.rad or surf.flat:
        surf.flat = True
        surf.points.append(com)
        return
    # For each path toward the center of the surface, set up a path list.
    paths = [[surf.perimeter[i]] for i in range(len(surf.perimeter))]
    # For each ring toward the center of the surface record a list
    surf.rings = [surf.perimeter]
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
        # Reset the ring variable
        ring = []
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
        # Add the ring to the list of rings
        surf.rings.append(ring)

    # Add the remaining paths to the surface excluding the first point in the path (i.e. the edge point)
    for path in paths:
        surf.points += path[1:]


############################################## Triangulate Surface Points  #############################################


def triangulate_rings(surf):
    """
    Finds the triangulation of the surface points by triangulating point rings
    :param surf: Surface object to triangulate
    :return: None
    """
    # Check to see if the surface is flat or not.
    if surf.flat:
        surf.tris = [[i, (i + 1) % len(surf.perimeter), len(surf.points) - 1] for i in range(len(surf.perimeter))]
        return
    # Check to see if the surface has flat points or not
    if surf.flat_points is None or len(surf.flat_points) < len(surf.points):
        surf.find_flat_points()
    # Go through the rings until the last one is used
    ring_num = 0
    surf.ring_tris, surf.tris = [], []
    # Get the point tracking indices for making the triangles
    outer_ndx, inner_ndx = 0, len(surf.rings[0])
    # Loop through the rings
    while ring_num < len(surf.rings) - 1:
        # Set up the ring triangles list variable
        ring_tris = []
        # Copy the rings
        outer_ring, inner_ring = surf.rings[ring_num].copy(), surf.rings[ring_num + 1].copy()
        # Find 2 points in the rings that are close
        min_dist = np.inf
        inner_start = 0
        for k in range(len(inner_ring)):
            if calc_dist(inner_ring[k], outer_ring[0]) < min_dist:
                inner_start = k
        inner_ring = inner_ring[inner_start:] + inner_ring[:inner_start]
        # Tracker variables for the points in the ring lists
        i, j = 0, 0
        # Go through the points in the two lists and get the set of triangles between
        while i + 1 < len(outer_ring) and j + 1 < len(inner_ring):
            # Calculate the circumference of the two new triangles
            oc = calc_dist(inner_ring[j], outer_ring[i + 1])
            ic = calc_dist(outer_ring[i], inner_ring[j + 1])
            # Check to see which of the triangles has the smaller circumference
            if oc < ic:
                tri = [i + outer_ndx, j + inner_ndx, i + outer_ndx + 1]
                i += 1
            else:
                tri = [i + outer_ndx, j + inner_ndx, j + inner_ndx + 1]
                j += 1
            # Add the triangle to the list of triangles
            ring_tris.append(tri)
        # Get the last couple triangles
        last_tris = []
        if i + 1 == len(outer_ring):
            # Keep looping until the points have all been used
            # Once we move past the end point all other triangles are with the starting point
            p1 = True
            while j + 1 < len(inner_ring) and p1:
                # Find the smallest circle (circle end, circle beginning)
                oc = calc_dist(inner_ring[j], outer_ring[0])
                ic = calc_dist(outer_ring[i], inner_ring[j + 1])
                # If the circle made with the next outer ring
                if oc < ic:
                    last_tris.append([inner_ndx + j, outer_ndx + i, outer_ndx])
                    last_tris += [[outer_ndx, inner_ndx + j + k, inner_ndx + j + k + 1] for k in range(len(inner_ring) - j - 1)]
                    p1 = False
                else:
                    last_tris.append([inner_ndx + j, outer_ndx + i, inner_ndx + j + 1])
                j += 1

        elif j + 1 == len(inner_ring):
            # Keep looping until the points have all been used
            # Once we move past the end point all other triangles are with the starting point
            p1 = True
            while i + 1 < len(outer_ring) and p1:
                # Find the smallest circle (circle end, circle beginning)
                oc = calc_dist(outer_ring[i], inner_ring[0])
                ic = calc_dist(inner_ring[j], outer_ring[i + 1])
                # If the circle made with the next outer ring
                if ic < oc:
                    last_tris.append([inner_ndx + j, outer_ndx + i, inner_ndx])
                    last_tris += [[inner_ndx, outer_ndx + i + k, outer_ndx + i + k + 1] for k in
                                  range(len(outer_ring) - i - 1)]
                    p1 = False
                else:
                    last_tris.append([outer_ndx + i, inner_ndx + j, outer_ndx + i + 1])
                i += 1
        else:
            return
        #
        # # Add the triangles to the ring's list
        # ring_tris += last_tris
        # Add the triangles to the surfaces list of triangles and the surfaces list of triangle rings
        surf.tris += ring_tris
        surf.ring_tris.append(ring_tris)
        # Set the indices for the next ring if needed
        outer_ndx, inner_ndx = outer_ndx + len(outer_ring), inner_ndx + len(inner_ring)
        ring_num += 1

    # Plot the results
    plot_surfs([surf], simps=True, Show=True)


def tri_within(surf, myTri=None, point=None):
    """
    Checks to see if a triangle lies within the perimeter of a surface
    :param surf: Surface object to check against
    :param myTri: Triangle to test insideness
    :param point: point to test insideness
    :return: Bool
    """
    # Get the perimeter of the translated and rotated surface
    perimeter = surf.pflat_points
    if len(perimeter) == 0:
        return False
    center = calc_com(points=perimeter)
    # If we are given a triangle determine the center of mass and use that point
    if myTri is not None:
        # Copy the triangle, retrieve its points and calculate the center of mass
        tri = myTri.copy()
        # Get the triangles points
        points = [surf.flat_points[tri[i]] for i in range(len(tri))]
        # Calculate the triangle's center of mass
        point = calc_com(points=points)
    else:
        # Move the point
        point = point - surf.center
        # Rotate the point
        new_point = rotate_points(surf.rn, [point])[0]
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


def calc_tri_circ(surf, tri):
    """
    Finds the circumference of the circumscribed circle for the triangle
    :param surf: Surface from which the triangle was created
    :param tri: Triangle to test
    :return: Circumference of the circle circumscribing the triangle
    """
    # Get the points of the triangle
    pa, pb, pc = surf.flat_points[tri[0]], surf.flat_points[tri[1]], surf.flat_points[tri[2]]
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
    # Check to see if the surface is flat or not.
    if surf.flat:
        surf.tris = [[i, (i + 1) % len(surf.perimeter), len(surf.points) - 1] for i in range(len(surf.perimeter))]
        return
    # Copy the surface points
    points = surf.points.copy()
    # Move all surf points toward the origin via center point
    for i in range(len(points)):
        points[i] = np.array(points[i]) - np.array(surf.center)

    # Calculate the angles to rotate the center point around
    nps = rotate_points(surf.rn, points)

    # Get the 2d version of the points and their Delaunay tesselation
    nps = np.array(nps)
    tris = mtri.Triangulation(nps[:, 0], nps[:, 1])
    # Add the flat points to the surface's list of flat points
    surf.flat_points = [nps[i, :2] for i in range(len(surf.points))]
    surf.tris = tris.triangles.tolist()


def filter_tris(surf):
    """
    Goes through the triangles on the surface measuring the circumference & testing if inside
    :param surf: Surface object holding the triangles for filtration
    :return:
    """
    # Check to see if the surface is flat or not
    if surf.flat:
        return
    # Set up a list of indices to remove for the triangles
    remove_ndxs = []
    # Go through the triangles in the surface
    for i in range(len(surf.tris)):
        # Grab the triangle and calculate its circumference
        tri = surf.tris[i]
        circ = calc_tri_circ(surf, tri)
        # If the circumference of the triangle is less than x times the min_dist check to see if tri is within
        if circ > 5 * surf.res and not tri_within(surf, tri):
            remove_ndxs.append(surf.tris.index(tri))

    # Remove the outer triangles
    remove_ndxs.sort()
    for i in range(len(remove_ndxs)):
        surf.tris.pop(remove_ndxs[-(i + 1)])
