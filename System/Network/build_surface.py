from System.calcs import *
import matplotlib.tri as mtri
from Visualize.visualize import *


# Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smallest of the surfaces 2 atoms through the point into the surface
def calc_surf_point(surf, point, return_roots=False):
    # Grab the surfaces function and atoms
    f, a0, a1 = surf.func, surf.atoms[0], surf.atoms[1]
    # Set up the unit vector
    vi = np.array(point) - np.array(a0.loc)
    vn = vi / np.linalg.norm(vi)
    # Set the atom's location as the root
    vi = a0.loc
    # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
    a = f[0] * vn[0] ** 2 + f[1] * vn[1] ** 2 + f[2] * vn[2] ** 2 + f[3] * vn[0] * vn[1] + f[4] * vn[1] * vn[2] + f[
        5] \
        * vn[2] * vn[0]
    b = 2 * f[0] * vn[0] * vi[0] + 2 * f[1] * vn[1] * vi[1] + 2 * f[2] * vn[2] * vi[2] + f[3] \
        * (vn[0] * vi[1] + vn[1] * vi[0]) + f[4] * (vn[1] * vi[2] + vn[2] * vi[1]) + f[5] \
        * (vn[2] * vi[0] + vn[0] * vi[2]) + f[6] * vn[0] + f[7] * vn[1] + f[8] * vn[2]
    c = f[0] * vi[0] ** 2 + f[1] * vi[1] ** 2 + f[2] * vi[2] ** 2 + f[3] * vi[0] * vi[1] + f[4] * vi[1] * vi[2] + \
        f[5] * vi[2] * vi[0] + f[6] * vi[0] + f[7] * vi[1] + f[8] * vi[2] + f[9]
    # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
    # and add that point to our surface list of points
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        # Calculate the roots
        roots = np.roots([a, b, c])
        # If the roots are requested return them
        if return_roots:
            return roots
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
    # Normalize this vector
    rn_hat = rn / np.linalg.norm(rn)
    # Find the next projection point by adding the vector with 'a' magnitude and rn_hat direction
    pc = pb + rn_hat * a
    # Calculate where the point intercepts the surface and return it
    return calc_surf_point(surf, pc)


# Build perimeter function. Sorts the edges of the surface to create a list of points in order around the perimeter
def build_perimeter(surf):
    # Reset the surface's perimeter points list
    surf.perimeter = []
    # Go through each edge in the surface's list of edges and build it
    for edge in surf.edges:
        if len(edge.points) == 0:
            edge.build()
    # Add the first edge's vertex location and set of points to the perimeter points list
    surf.perimeter = [surf.edges[0].verts[0].loc] + surf.edges[0].points
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
            surf.perimeter += [myEdge.verts[0].loc] + myEdge.points
        else:  # Reverse order
            surf.perimeter += [myEdge.verts[1].loc] + myEdge.points[::-1]
    # Add the perimeter points to the whole set of points
    surf.points += surf.perimeter


# Fill mesh function. Works inward from a set of perimeter points toward a center point filling in equally spaced points
def fill_mesh(surf):
    # Check to see that the surface has perimeter points
    if len(surf.perimeter) == 0:
        build_perimeter(surf)
    # Get the resolution
    res = surf.net.sys.min_dist
    # Get the atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    # Get the center point for the surface
    center = surf.rn * 0.5 * (calc_dist(a0.loc, a1.loc) - (a0.rad + a1.rad)) + a0.loc
    # Get the center of mass for the edges
    com3d = calc_edges_com(surf.edges)
    com3d_proj = calc_surf_point(surf, com3d)
    com3d_trans = com3d_proj - center
    com2d_trans = rotate_points(center, [com3d_trans])
    # Check to see if the center point is inside the perimeter or not
    if tri_within(surf, point=com2d_trans):
        com = com3d_proj
    else:
        if tri_within(surf, point=center):
            com = center
        else:
            print("Bad surface! Surface: ", surf.ndx)
            return
    # Check to see if the atoms have equal radii
    if a0.rad == a1.rad:
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
    # Add the center of mass point to the set of points
    surf.points.append(com)


# Triangle within the surface function. Checks to see if a triangle lies within the perimeter of a surface
def tri_within(surf, myTri=None, point=None):
    # Get the perimeter of the translated and rotated surface
    perimeter = surf.flat_points[:len(surf.perimeter)]
    if len(perimeter) == 0:
        return False
    center = calc_com(points=perimeter)
    # If we are given a triangle determine the center of mass and use that point
    if myTri:
        # Copy the triangle, retrieve its points and calculate the center of mass
        tri = myTri.copy()
        # Get the triangles points
        points = [surf.flat_points[tri[i]] for i in range(len(tri))]
        # Calculate the triangle's center of mass
        point = calc_com(points=points)
    else:
        # The point needs to be rotated and dropped to 2D
        my_point = rotate_points(surf.rn, [point])[0]
        point = my_point[:2]
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
        if theta_n < theta and theta_n1 < theta:
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
    a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
    b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
    c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
    s = (a + b + c) / 2.0
    area = np.sqrt(max(s * (s - a) * (s - b) * (s - c), 0))
    if area == 0:
        return np.inf
    circum_r = a * b * c / (4.0 * area)
    return circum_r


# Find simplices function. Transforms and rotates surface points to xy-plane and returns the Delaunay simplices
def find_simps(surf):
    # Get the atoms
    a0, a1, d = surf.atoms[0], surf.atoms[1], np.linalg.norm(surf.rn)
    # Get the center of the surface
    c = np.array(a1.loc) - (0.5 * (d - (a0.rad + a1.rad)) + a0.rad) * surf.rn
    # Copy the surface points
    points = surf.points.copy()
    # Move all surf points toward the origin via center point
    for i in range(len(points)):
        points[i] = points[i] - c
    # Calculate the angles to rotate the center point around
    nps = rotate_points(surf.rn, points)
    # Get the 2d version of the points
    nps = np.array(nps)
    # Get the Delaunay tesselation
    tris = mtri.Triangulation(nps[:, 0], nps[:, 1])
    # Find the 2d polygon
    surf.flat_points = [nps[i, :2] for i in range(len(surf.points))]
    surf.tris = tris.triangles.tolist()


# Filter triangles function. Goes through the
def filter_tris(surf):
    # Set up a list of indices to remove for the triangles
    remove_ndxs = []
    # Go through the triangles in the surface
    for i in range(len(surf.tris)):
        # Grab the triangle and calculate its circumference
        tri = surf.tris[i]
        circ = calc_tri_circ(surf, tri)
        # If the circumference of the triangle is less than x times the minimum distance check to see if tri is within
        if circ > 5 * surf.net.sys.min_dist and not tri_within(surf, tri):
            remove_ndxs.append(surf.tris.index(tri))
    # Remove the outer triangles
    remove_ndxs.sort()
    for i in range(len(remove_ndxs)):
        surf.tris.pop(remove_ndxs[-(i + 1)])


# Make mesh method. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf):
    # Get the surface's function coefficients
    if surf.func is None:
        surf.calc_func()
    # Get the normal to the surface
    if surf.rn is None:
        r = np.array(surf.atoms[1].loc) - np.array(surf.atoms[0].loc)
        surf.rn = r / np.linalg.norm(r)
    # Reset the surface's list of points to empty list and reset the vertex indices list
    surf.points = []
    # Build the perimeter of the surface
    build_perimeter(surf)
    # Fill the mesh
    fill_mesh(surf)
    # Find the simplices of the surface
    find_simps(surf)
    # Filter out the bad triangles
    filter_tris(surf)
