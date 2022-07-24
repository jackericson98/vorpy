from System.Network.Edges.edge_trace import *


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
        surf.edge_points = edge_trace1(surf, min_dist)
    # If no edges exist create a circular edge
    elif not surf.edges:
        # If no radius is specified, create one 5x larger than the size of the center atom
        if radius is None:
            radius = surf.atoms[1].rad * 5
        # Add the circular edge points to the surfaces list of edge points
        surf.edge_points = circ_edge_trace(surf, radius, min_dist)
    # This is for the voronota plot. If the first edge in the list of edges has points add all edges' points so the list
    elif vta:
        # Go through the
        for edge in surf.edges:
            surf.edge_points += edge.points
        surf.points += surf.edge_points
        return
    else:
        return
    # Add the edge points to the surface's points
    surf.points += surf.edge_points
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
