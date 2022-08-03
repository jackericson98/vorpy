from System.sys_funcs import *
import matplotlib.tri as mtri


class Surface:
    """Surface object. Holds the mesh data. Used to analyze."""
    def __init__(self, atoms, edges=None, verts=None, min_dist=0.1):
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
        self.min_dist = min_dist
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

    # Calculate surface point function. Takes in a surface and a point and returns the intersection point of the vector
    # from the center of the smallest of the surfaces 2 atoms through the point into the surface
    def calc_surf_point(self, point):
        # Grab the function's coefficients
        f = self.func
        # Get the first atoms in the surfaces list of atoms
        a0, a1 = self.atoms[0], self.atoms[1]
        # Set up the unit vector
        vi = np.array(point) - np.array(a0.loc)
        vn = vi / np.linalg.norm(vi)
        # Find the location on the surface of the atom
        vi = np.array(a0.loc) + vn * a0.rad
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
        if round(b ** 2 - 4 * a * c, 4) >= 0:
            roots = np.roots([a, b, c])
            # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
            if calc_dist(vi, a1.loc) > a1.rad:
                x = 1 - (calc_dist(a0.loc, a1.loc) * a0.rad) / a1.rad
                if calc_angle(a0.loc, a1.loc, vi) - np.pi / 2 > x:
                    mag = max(abs(roots))
                else:
                    mag = min(abs(roots))
            # If the projection point is within the intersection, the magnitude is negative
            else:
                if calc_angle(a0.loc, a1.loc, vi) > np.pi / 4:
                    mag = - min(abs(roots))
                else:
                    mag = -abs(min(roots))
            return vi + mag * vn

    # Find next point method. Finds the next point along the given path by projecting a reference point onto the surface
    def find_next_point(self, pn_1, end, d_theta):
        # Get the A angle
        A = d_theta
        # Get the smaller atom's location
        pa = self.atoms[0].loc
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
        sp = self.calc_surf_point(pc)
        return sp

    # Make mesh method. Goes in shrinking concentric circles inside the edges of the surface toward the com of the edges
    def make_mesh(self):
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
            edge.build(surf=self)
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
        # Check to see if the center of mass is going to trigger a
        if d_coma1 > calc_dist(a0.loc, a1.loc):
            theta = calc_angle(a0.loc, com, a1.loc) - np.pi / 2
            mag = 2 * calc_dist(a0.loc, com) * np.sin(theta)
            r = (np.array(a1.loc) - np.array(a0.loc))
            rn = r / np.linalg.norm(r)
            com = com + mag * rn
        # Calculate the center of mass point of the edge points and where it maps on the surface
        com = self.calc_surf_point(com)
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
        num_rings = max(int(calc_dist(max_path, com) / self.min_dist), 10)
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
                pn = self.find_next_point(paths[i][-1], com, dthetas[i])
                # Check to see of the new point is too close to the previous point and the path has to end
                if calc_dist(pn, pn_1) < self.min_dist:
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
        self.make_mesh()
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
