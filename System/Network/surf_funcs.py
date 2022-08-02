from System.sys_funcs import *
import matplotlib.tri as mtri


# Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smallest of the surfaces 2 atoms through the point into the surface
def calc_surf_point(surf, point):
    # Grab the function's coefficients
    f = surf.func
    # Get the first atoms in the surfaces list of atoms
    a0, a1 = surf.atoms[0], surf.atoms[1]
    # Get the vector between the atoms
    r01 = np.array(a1.loc) - np.array(a0.loc)
    r01n = r01/np.linalg.norm(r01)
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
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        roots = np.roots([a, b, c])
        # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
        if calc_dist(vi, a1.loc) > a1.rad:
            if calc_angle(a0.loc, a1.loc, vi) > 5 * np.pi / 8:
                mag = max(abs(roots))
            else:
                mag = min(abs(roots))
        # If the projection point is within the intersection, the magnitude is negative
        else:
            if calc_angle(a0.loc, a1.loc, vi) > np.pi/4:
                mag = - min(abs(roots))
            else:
                mag = -abs(min(roots))
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

