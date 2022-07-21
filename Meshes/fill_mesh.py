from Meshes.edge_trace import *
from System.system import Atom, Surface
from Presentation.Visualize.visualize import plot_surfs, plot_atoms
import matplotlib.pyplot as plt


# Make mesh function. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
def make_mesh(surf, min_dist, circ_mesh=False, radius=None, edges_made=False):
    # Set the atoms in the surface to make the smaller one listed first
    if surf.atoms[0].rad > surf.atoms[1].rad:
        surf.atoms = surf.atoms[1], surf.atoms[0]
    # Check to see if a circular mesh has been requested
    if circ_mesh:
        if radius is None:
            radius = surf.atoms[1].rad * 5
        surf.edge_points = make_bless_mesh(surf, radius, min_dist)
    elif edges_made:
        pass
    else:
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


# Make boundless mesh.
def make_bless_mesh(surf, radius, min_dist):
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

