from System.Network.surf_funcs import *
import matplotlib.tri as mtri


class Surface:
    """Surface object. Holds the mesh data. Used to analyze."""
    def __init__(self, atoms, edges=None, verts=None):
        self.func = None
        self.atoms = atoms  # List of Atom type objects
        self.edges = edges  # List of Edge type objects
        self.verts = verts
        if verts:
            self.vert_points = [verts[i].loc for i in range(len(verts))]
        self.edge_points = []
        self.surf_points = []  # List of points on the surface
        self.points = []
        self.tris = None
        self.sa = None
        self.calc_func()

    # Bisector function. Creates a bisector surface between 2 atoms
    def calc_func(self):
        # Make sure that a0 is the atom with the smaller radius
        if self.atoms[0].rad > self.atoms[1].rad:
            self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
        a0, a1 = self.atoms
        # Grab the centers of the spheres
        x1, y1, z1 = a0.loc
        x2, y2, z2 = a1.loc
        # Calculate the major coefficients (pg. 574 Z. Hu)
        R = a0.rad - a1.rad
        K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
        d = x1 - x2, y1 - y2, z1 - z2
        J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2
        # Instantiate/reset the hyperboloid coefficient vector lists
        ABC, DEF, GHI = [], [], []
        # Calculate hyperboloid coefficients
        for i in range(3):
            ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
            DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
            GHI.append(-8 * R ** 2 * a0.loc[i] - 4 * K * d[i])
        # Set the function attribute
        self.func = ABC + DEF + GHI + [J] + [K] + [d]

    # Make mesh method. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
    def make_mesh(self, min_dist):
        # Get the atoms
        a0, a1 = self.atoms[0], self.atoms[1]
        # Reset the all surface points to empty lists
        self.points, self.vert_points, self.edge_points, self.surf_points = [], [], [], []
        # Go through each vertex on the surface
        for vert in self.verts:
            # Add the points to the surface's list of vertex points
            self.vert_points.append(vert.loc)
        # Add the vert points to the surface's points
        self.points += self.vert_points
        # Go through each edge in the surface's list of edges
        for edge in self.edges:
            edge.build(surf=self, min_dist=min_dist)
            # Add the edge's points to the surface's edge points attribute
            self.edge_points += edge.points
        # Check to see if the atoms have equal radii
        if a0.rad == a1.rad:
            self.points += self.edge_points
            return
        # Get the center of mass of the edges of the surface
        com = calc_edges_com(self.edges)
        # Get the distance between the center of mass and a1
        d_coma1 = calc_dist(com, self.atoms[1].loc)
        # Check to see if the center of mass is
        if d_coma1 > calc_dist(a0.loc, a1.loc):
            theta = calc_angle(a0.loc, com, a1.loc) - np.pi / 2
            mag = 2 * calc_dist(a0.loc, com) * np.sin(theta)
            r = (np.array(a1.loc) - np.array(a0.loc))
            rn = r / np.linalg.norm(r)
            com = com + mag * rn
        # Calculate the center of mass point of the edge points and where it maps on the surface
        com = calc_surf_point(self, com)
        # Add the edge points to the surface's points
        self.points += self.edge_points
        # For each edge point set up a path list.
        paths = [[self.edge_points[i]] for i in range(len(self.edge_points))]
        # Grab the smallest of the 2 surface atoms' location
        pa = self.atoms[0].loc
        # Get the angles between the edge points and the end points
        angs = []
        for i in range(len(paths)):
            # Calculate the angle for each path
            angs.append(calc_angle(pa, paths[i][0], com))
        # Get the maximum path
        max_path_ndx = angs.index(max(angs))
        max_path = paths[max_path_ndx][0]
        # Decide how many rings based off of the ellipticity and density
        num_rings = max(int(calc_dist(max_path, com) / min_dist), 10)
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
                pn = find_next_point(paths[i][-1], com, dthetas[i], self)
                # Check to see of the new point is too close to the previous point and the path has to end
                if calc_dist(pn, pn_1) < min_dist:
                    # Add the path to the surfaces points and remove it from the paths list
                    self.surf_points += paths.pop(i)[1:]
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
            self.surf_points += path[1:]
        # Add the center of mass point to the mesh
        self.surf_points.append(com)
        # Add the surface points to the general list of points
        self.points += self.surf_points

    # Calculate surface simplices function.
    def find_simps(self):
        # Get the atoms
        a0, a1 = self.atoms
        # Find the normal to the surface and the magnitude
        r10 = np.array(a0.loc) - np.array(a1.loc)
        d = np.linalg.norm(r10)
        r10_hat = r10 / d
        # Get the distance between the surfaces
        ds = d - (a0.rad + a1.rad)
        # Get the center of the surface
        c = np.array(a1.loc) + (0.5 * ds + a0.rad) * r10_hat
        # Move all surf points toward the origin via center point
        for i in range(len(self.points)):
            self.points[i] = self.points[i] - c
        # Calculate the angles to rotate the center point around
        nps = rotate_points(c, self.points)
        # Get the 2d version of the points
        nps = np.array(nps)
        nps2d = nps[:, 0], nps[:, 1]
        # Get the Delaunay tesselation
        tri = mtri.Triangulation(nps2d[0], nps2d[1])
        # Move the points back to their original location
        for i in range(len(self.points)):
            self.points[i] = self.points[i] + c
        # Filter out any connections between the vertices or the edges
        self.tris = tri.triangles.tolist()
        # If the surface's atoms have equal radii, were done
        if a0.rad == a1.rad:
            return
        remove_ndxs = []
        # Go through each triangle on the surface
        for i in range(len(self.tris)):
            # Set the counter to 0
            counter = 0
            # Go through each point on the triangle checking to see if it is an edge point
            for j in range(3):
                # If the triangles jth point index is less than the number of vertex & edge points increment the counter
                if self.tris[i][j] < len(self.vert_points) + len(self.edge_points):
                    counter += 1
            # If all three of the points are on an edge we need to check it
            if counter == 3:
                remove_ndxs.append(i)
        # Remove the outer triangles
        remove_ndxs.sort()
        for tri_ndx in remove_ndxs[::-1]:
            if tri_ndx:
                self.tris.pop(tri_ndx)

    # Build method. Makes the mesh for the surface and calculates the simplices between them
    def build(self, min_dist=0.1, simps=True):
        # Build the mesh
        self.make_mesh(min_dist=min_dist)
        # Calculate the simplices
        if simps:
            self.find_simps()

    # Build vta surface function
    def build_vta(self):
        # Add the vertex points to the surface's list of points
        for vert in self.verts:
            self.points.append(vert.loc)
        # Calculate the simplices
        self.find_simps()
