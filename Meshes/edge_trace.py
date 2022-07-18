from Meshes.mesh_calcs import *


# Calculate edge points function. Takes in an edge and a surface and updates the edge's points.
def calc_edge_points(edge, surf, min_dist):
    # Get the location of the base atom
    pa = edge.atoms[0].loc
    # Get the locations of the vertices
    pv0 = np.array(edge.verts[0].loc)
    pv1 = np.array(edge.verts[1].loc)
    # Find the angle made between the edges vertices and the atom
    max_ang = calc_angle(pa, pv0, pv1)
    num_points = max(int(calc_dist(pv0, pv1) / min_dist), 10)
    # Set angle A to be the incremental angle decided by num points
    A = max_ang / num_points
    # Calculate each point along the way
    for i in range(1, num_points):
        # If the edge points are empty set pb to the start vertex. Else get the previous point in the path
        if not edge.points:
            pb = pv0
        else:
            pb = edge.points[-1]
        edge.points.append(find_next_point(pb, pv1, A, surf))


def edge_trace1(surf, min_dist):
    for edge in surf.edges:
        # If the edge points exist already add them to the surfaces edge points and continue to the next edge
        if edge.points:
            surf.edge_points += edge.points
            continue
        calc_edge_points(edge, surf, min_dist)
        # Add the edge's points to the surface's edge points attribute
        surf.edge_points += edge.points


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