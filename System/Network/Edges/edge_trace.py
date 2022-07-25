from System.Network.Surfaces.surf_calcs import *


# Calculate edge points function. Takes in an edge and a surface and updates the edge's points.
def calc_edge_points(edge, surf, min_dist):
    # Get the location of the base atom
    pa = edge.atoms[0].loc
    # Get the locations of the vertices
    pv0 = edge.verts[0].loc
    pv1 = edge.verts[1].loc
    # Find the angle made between the edges vertices and the atom
    max_ang = calc_angle(pa, pv0, pv1)
    num_points = max(int(calc_dist(pv0, pv1) / min_dist), 10)
    # Set angle A to be the incremental angle decided by num points
    A = max_ang / num_points
    points = []
    # Calculate each point along the way
    for i in range(1, num_points):
        # If the edge points are empty set pb to the start vertex. Else get the previous point in the path
        if not points:
            pb = pv0
        else:
            pb = points[-1]
        point1 = find_next_point(pb, pv1, A, surf)
        points.append([point1[0], point1[1], point1[2]])
    return points


# Edge trace function.
def edge_trace1(surf, min_dist):
    # Instantiate the edge_points list
    edge_points = []
    # Go through each edge in the surface's list of edges
    for edge in surf.edges:
        # If the edge points exist already add them to the surfaces edge points and continue to the next edge
        if not edge.points:
            edge.points = calc_edge_points(edge, surf, min_dist)
        # Add the edge's points to the surface's edge points attribute
        edge_points += edge.points
    return edge_points


# Edge trace function
def edge_trace(surf, min_dist):
    # Go through each edge on the surface
    for edge in surf.edges:
        # If the edge points exist already continue to the next edge
        if not edge.points:
            calc_edge_points(edge, surf, min_dist)
    # Call the recursive edge tracing function
    edges = surf.edges
    # Get the first edge and vert
    v0, vert = edges[0].verts
    surf.edge_points = edges[0].points
    # Find the edge that is closest to vert and add the points accordingly
    v0_found = False
    while not v0_found:
        # Check the edges for similar vertices
        for edge in edges[1:]:
            # If the first vertex in the edge's vertex list is equal to vert, add the points
            if edge.verts[0] == vert:
                surf.edge_points += edge.points
                vert = edge.verts[1]
            elif edge.verts[1] == vert:
                surf.edge_points = edge.points[::-1]
                vert = edge.verts[0]
        if vert.loc == v0.loc:
            v0_found = True


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
    num_points = int(2*np.pi*radius/min_dist)
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
