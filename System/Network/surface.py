from System.Network.surf_funcs import *


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
        self.simps = None
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
        self.simps = tri

    # Build method. Makes the mesh for the surface and calculates the simplices between them
    def build(self, min_dist=0.1, simps=True):
        # Build the mesh
        make_mesh(self, min_dist=min_dist)
        # Calculate the simplices
        if simps:
            self.find_simps(self.points)

    # Build vta surface function
    def build_vta(self):
        # Add the vertex points to the surface's list of points
        for vert in self.verts:
            self.points.append(vert.loc)
        # Calculate the simplices
        self.find_simps(self.points)
