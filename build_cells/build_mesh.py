from build_cells.calculators import *


########################################################################################################################


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


# Make mesh function. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf, min_dist):
    # Set the atoms in the surface to make the smaller one listed first
    if surf.atoms[0].rad > surf.atoms[1].rad:
        surf.atoms = surf.atoms[1], surf.atoms[0]
    # Check to see if the edges' points have been recorded yet
    edge_trace1(surf, min_dist)
    # For each edge point set up a path list.
    paths = [[surf.edge_points[i]] for i in range(len(surf.edge_points))]
    # Grab the smallest of the 2 surface atoms' location
    pa = surf.atoms[0].loc
    # Calculate the center of mass point of the edge points and where it maps on the surface
    com = calc_surf_point(surf, calc_com(surf.edge_points))
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
                surf.points += paths.pop(i)[1:]
                ends.pop(i)
                dthetas.pop(i)
                num_paths -= 1
            else:
                # Set the pn_1 to pn and add it to the path
                pn_1 = pn
                paths[i].append(pn)
            # Increment i
            i += 1
    # Add the remaining paths to the surface
    for path in paths:
        surf.points += path[1:]


# Build meshes function. Runs make_mesh on all surfaces in the network
def build_meshes(sys, min_dist=None):
    # Set the minimum distance
    if min_dist is None:
        min_dist = calc_dist(sys.net.edges[0].verts[0].loc, sys.net.edges[0].verts[1].loc) / 30
    # Make each surface
    for surf in sys.net.surfs:
        make_mesh(surf, min_dist)
