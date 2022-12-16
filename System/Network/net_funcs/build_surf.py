from System.sys_funcs.calcs import *
import matplotlib.tri as mtri
from Visualize.mpl_visualize import *

################################################# Find Surface Points  #################################################


# Calculate surface point function. Finds the projection a surface's small atom through the given point onto the surface
def calc_surf_point(surf, point):
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


# Find next point method. Finds the next point along the given path by projecting a reference point onto the surface
def find_next_point(surf, pn_1, end, d_theta):

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


# Build perimeter function. Sorts the edges of the surface to create a list of points in order around the perimeter
def build_perimeter(surf):

    # Reset the surface's perimeter points list
    surf.perimeter = []
    e0 = surf.edges[0]
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
            d0, d1 = calc_dist(surf.perimeter[-1], edges[i].points[0]), \
                     calc_dist(surf.perimeter[-1], edges[i].points[-1])
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

    # Add the perimeter points to the whole set of points
    surf.points += surf.perimeter
    # Get the perimeter flat points
    surf.pflat_points = surf.perimeter.copy()
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    d = calc_dist(a0.loc, a1.loc)
    # Get the center of the surface
    surf.center = np.array(a0.loc) + (a0.rad + 0.5 * (d - (a0.rad + a1.rad))) * surf.rn
    for i in range(len(surf.pflat_points)):
        # Move the points
        surf.pflat_points[i] = surf.pflat_points[i] - surf.center
    # Rotate the point
    surf.pflat_points = rotate_points(surf.rn, surf.pflat_points)
    # Get the 2d version
    surf.pflat_points = [point[:2] for point in surf.pflat_points]


# Get center of mass function.Finds the center of mass of a surface's perimeter
def get_com(surf):
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


# Fill mesh function. Works inward from a set of perimeter points toward a center point filling in equally spaced points
def fill_mesh(surf):
    # Check to see that the surface has perimeter points
    if len(surf.perimeter) == 0:
        build_perimeter(surf)
    # Get the resolution
    res = surf.net.surf_res
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    # Get the center of mass
    com = get_com(surf)
    # Check to see if the atoms have equal radii
    if a0.rad == a1.rad:
        surf.flat = True
        surf.points.append(com)
        return
    # For each edge point set up a path list.
    paths = [[surf.perimeter[i]] for i in range(len(surf.perimeter))]
    # Grab the smallest of the 2 surface atoms' location
    pa = surf.atoms[0].loc
    # Get the angles between the edge points and the end points
    dists = []
    angs = []
    for i in range(len(paths)):
        # Calculate the angle for each path
        angs.append(calc_angle(pa, paths[i][0], com))
        # Get the dists from the com to the path
        dists.append(calc_dist(paths[i][0], com))
    # Get the maximum path
    max_path_ndx = angs.index(max(angs))
    max_path = paths[max_path_ndx][0]
    # Decide how many rings based off of the ellipticity and density
    num_rings = max(int(calc_dist(max_path, com) / res), 2)
    # Get the incremental angle increases
    dthetas = [angs[i] / num_rings for i in range(len(angs))]
    # Set the pn_1 point to infinity
    pn_1 = [np.inf, np.inf, np.inf]
    num_paths = len(paths)
    # Go through ring by ring
    for j in range(num_rings):
        # Go through each of the remaining paths
        i = 0
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
            if pn is None or (calc_dist(pn, pn_1) < 0.5 * res and not calc_dist(paths[i - 1][-1], pn) > res):
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


# Triangle within the surface function. Checks to see if a triangle lies within the perimeter of a surface
def tri_within(surf, myTri=None, point=None):

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


# Calculate triangle circumference function. Finds the circumference of the circumscribed circle for the triangle
def calc_tri_circ(surf, tri):

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


# Find simplices function. Transforms and rotates surface points to xy-plane and returns the Delaunay simplices
def find_simps(surf):
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


# Filter triangles function. Goes through the triangles on the surface measuring the circumference & testing if inside
def filter_tris(surf):
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
        if circ > 5 * surf.net.surf_res and not tri_within(surf, tri):
            remove_ndxs.append(surf.tris.index(tri))

    # Remove the outer triangles
    remove_ndxs.sort()
    for i in range(len(remove_ndxs)):
        surf.tris.pop(remove_ndxs[-(i + 1)])


# Make mesh method. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf):

    # Prepare Surface

    # Get the surface's function coefficients
    if surf.func is None:
        surf.calc_func()
    # Reset the surface's list of points to empty list and reset the vertex indices list
    surf.points = []

    # Build Surface:

    # Build the perimeter of the surface
    build_perimeter(surf)
    # Fill the mesh
    fill_mesh(surf)
    # Find the simplices of the surface
    find_simps(surf)
    # Filter out the bad triangles
    filter_tris(surf)
