from System.sys_funcs import *
import matplotlib.tri as mtri


# Calculate surface simplices function.
def find_simps(points, atoms):
    # Get the atoms
    a0, a1 = atoms
    # Find the normal to the surface and the magnitude
    r10 = np.array(a0.loc) - np.array(a1.loc)
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
    # Filter out any connections between the vertices
    for i in range(len(points)):
        points[i] = points[i] + c
    return tri


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
        roots = np.roots([a, b, c])
        # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
        if calc_dist(vi, a1.loc) > a1.rad:
            mag = min(abs(roots))
        # If the projection point is within the intersection, the magnitude is negative
        else:
            #
            if calc_dist(a0.loc, a1.loc) > a0.rad:
                mag = - abs(min(roots))
            else:
                mag = - min(abs(roots))
        return vi + mag * vn


# Find next point function. Finds the next point along the given path by projecting a reference point onto the surface
def find_next_point(pn_1, end, d_theta, surf):
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
    # Find the intercept point by adding a to pb
    rn = end - pb
    rn_hat = rn / np.linalg.norm(rn)
    pc = pb + rn_hat * a
    # Calculate where the point intercepts the surface
    return calc_surf_point(surf, pc)


# Edge trace function. I want to update to put the points in order
def edge_trace1(surf):
    # Instantiate the edge_points list
    edge_points = []
    # Go through each edge in the surface's list of edges
    for edge in surf.edges:
        edge.build(surf=surf)
        # Add the edge's points to the surface's edge points attribute
        edge_points += edge.points
    return edge_points


# Circular edge trace function. When no edges are made make_mesh defaults to this. Creates a circular edge to build from
def circ_edge_trace(surf, radius, min_dist):
    # Grab the surfaces atoms and make sure the smaller one is a0
    a0, a1 = surf.atoms[0], surf.atoms[1]
    if a0.rad > a1.rad:
        a0, a1, = a1, a0
    r01 = np.array(a1.loc) - np.array(a0.loc)
    # Get the normalized vector for the direction toward the center of the surface
    r01_hat = r01/np.linalg.norm(r01)
    # Get the point on the surface of a0 closest to a1
    dist = calc_dist(a0.loc, a1.loc) - (a0.rad + a1.rad)
    # Get the point on the surface corresponding to vc
    center = a0.loc + r01_hat * (dist + a0.rad)
    # Find a vector perpendicular to r01_hat
    if abs(r01_hat[0]) > abs(r01_hat[1]):
        p = np.array([r01_hat[1], -r01_hat[0], 0])
    elif r01_hat[0] == 0 and r01_hat[1] == 0:
        p = np.array([1, 0, 0])
    else:
        p = np.array([-r01_hat[1], r01_hat[0], 0])
    # Normalize it
    p_hat = p/np.linalg.norm(p)
    c_points = [p_hat*radius + center]
    # Get the circumference of the circle and divide by the minimum distance
    num_points = int(np.pi*radius/(2*min_dist))
    # Get the incremental angle change around the circle
    dtheta = 2*np.pi/num_points
    # Find the amount we project on the circle
    proj_dist = radius * np.tan(dtheta)
    for i in range(num_points):
        # Find the binormal vector to r01_hat and the previous circle point by taking their cross products
        bi = np.cross(r01_hat, c_points[-1] - center)
        # Normalize it
        bi_hat = bi/np.linalg.norm(bi)
        # Get the surface point
        samp = c_points[-1] + proj_dist*bi_hat
        rn = samp - center
        rn_hat = rn/np.linalg.norm(rn)
        c_points.append(center + rn_hat*radius)
    edge_points = []
    for i in range(len(c_points)):
        edge_points.append(calc_surf_point(surf, c_points[i]))
    return edge_points


# Make mesh function. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf, min_dist, radius=None, vta=False):
    # Reset the all surface points to empty lists
    surf.points, surf.vert_points, surf.edge_points, surf.surf_points = [], [], [], []
    # If the surface has vertices, add those points to the vert_points attribute of the surface
    if surf.verts and not vta:
        # Go through each vertex on the surface
        for vert in surf.verts:
            # Add the points to the surface's list of vertex points
            surf.vert_points.append(vert.loc)
        # Add the vert points to the surface's points
        surf.points = surf.vert_points
        # Use the edge tracing function to get edges' points
        surf.edge_points = edge_trace1(surf)
        # Calculate the center of mass point of the edge points and where it maps on the surface
        com = calc_surf_point(surf, calc_edges_com(surf.edges))
    # If no edges exist create a circular edge
    elif not surf.edges:
        # If no radius is specified, create one 5x larger than the size of the center atom
        if radius is None:
            radius = surf.atoms[1].rad * 5
        # Add the circular edge points to the surfaces list of edge points
        surf.edge_points = circ_edge_trace(surf, radius, min_dist)
        # Calculate the center of mass point of the edge points and where it maps on the surface
        com = calc_surf_point(surf, calc_edges_com(points=surf.edge_points))
    # This is for the voronota plot. If the first edge in the list of edges has points add all edges' points so the list
    elif vta:
        # Go through the vertices and add their locations to the surface's points
        for vert in surf.verts:
            surf.vert_points.append(vert.loc)
        # Make sure the edge's points are set to just include the vertices
        for edge in surf.edges:
            edge.points = [edge.verts[0].loc, edge.verts[1].loc]
        surf.points += surf.vert_points
        return
    else:
        return
    # Add the edge points to the surface's points
    surf.points += surf.edge_points
    # For each edge point set up a path list.
    paths = [[surf.edge_points[i]] for i in range(len(surf.edge_points))]
    # Grab the smallest of the 2 surface atoms' location
    pa = surf.atoms[0].loc
    # Set up a list of end points
    ends = [com for i in range(len(paths))]
    # Get the angles between the edge points and the end points
    angs = []
    for i in range(len(paths)):
        # Calculate the angle for each path
        angs.append(calc_angle(pa, paths[i][0], ends[i]))
    # Get the maximum path
    max_path_ndx = angs.index(max(angs))
    max_path = paths[max_path_ndx][0]
    # Decide how many rings based off of the ellipticity and density
    num_rings = max(int(calc_dist(max_path, ends[max_path_ndx]) / min_dist), 10)
    # Get the incremental angle increases
    dthetas = [angs[i]/num_rings for i in range(len(angs))]
    # Set the pn_1 point to infinity
    pn_1 = [np.inf, np.inf, np.inf]
    num_paths = len(paths)
    surf_points = []
    # Go through ring by ring
    for j in range(num_rings):
        # Go through each of the remaining paths
        i = 0
        while i < num_paths:
            # Get the next point along the path
            pn = find_next_point(paths[i][-1], ends[i], dthetas[i], surf)
            # Check to see of the new point is too close to the previous point and the path has to end
            if calc_dist(pn, pn_1) < min_dist:
                # Add the path to the surfaces points and remove it from the paths list
                surf_points += paths.pop(i)[1:]
                ends.pop(i)
                dthetas.pop(i)
                num_paths -= 1
            else:
                # Set the pn_1 to pn and add it to the path
                pn_1 = pn
                paths[i].append(pn)
                # Increment i
                i += 1
    # Add the remaining paths to the surface excluding the first point in the path (i.e. the edge point)
    for path in paths:
        surf_points += path[1:]
    surf.surf_points = np.array(surf_points).tolist()
    # Add the surface points to the general list of points
    surf.points += surf.surf_points
