from System.Network.net_funcs.build_surf import *


class Surface:
    """Surface object. Holds the mesh data. Used to analyze interfaces between atoms."""
    def __init__(self, atoms=None, net=None, edges=None, verts=None, doublet=False, points=None, tris=None, perimeter=None,
                 rn=None, sa=0, function=None, load_ndxs=None):

        # If no network was given have a catch
        if net is not None and net.atoms is not None:
            ndx = [net.atoms.index(atom) for atom in atoms]
            ndx.sort()
            self.ndx = ndx          # Index            : Indices of the atoms of the surface
            self.net = net          # Network          : Network of the System
        self.atoms = atoms          # Atoms            : Atoms of the surface
        self.verts = verts          # Vertices         : Vertices of the surface
        self.edges = edges          # Edges            : Edges of the surface
        self.load_ndxs = load_ndxs  # Load indices     : List of object load indices

        self.func = function        # Surface function : Holds the coefficients of the function describing the surf
        self.perimeter = perimeter  # Perimeter        : The points around the edges of the surface (IN ORDER)
        self.points = points        # Points           : The points that make up the surface
        self.rings = None           # Rings            : Lists of points in rings toward the center
        self.ring_tris = None       # Ring Triangles   : Triangulation of each set of concentric ring points
        self.flat_points = []       # Flattened points : Points projected into 2d based off of the surface normal
        self.pflat_points = []      # Flat perimeter   : Flattened points around the perimeter
        self.tris = tris            # Triangles        : A list of connections between the points
        self.sa = sa                # Surface Area     : The surface area of the
        self.rn = rn                # Surface Normal   : Normal to the center of the surface
        self.center = None          # Center           : Center point of the hyperboloid the surface is made from
        self.com = None             # Center of mass   : The point toward which all building paths travel
        self.doublet = doublet      # Doublet          : Indicates whether a surface is a part of a doublet or not
        self.flat = False           # Flat             : Whether the surface is flat or not

    # Bisector function. Creates a bisector surface between 2 atoms
    def calc_func(self):
        # Make sure that a0 is the atom with the smaller radius
        if self.atoms[0].rad > self.atoms[1].rad:
            self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
        # Create a0, a1 variables
        a0, a1 = self.atoms
        # Set the rn vector for the surface since the atoms are sorted
        r = np.array(a1.loc) - np.array(a0.loc)
        self.rn = r / np.linalg.norm(r)
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

    # Build method. Makes the mesh for the surface and calculates the simplices between them
    def build(self):
        # Build the mesh
        make_mesh(self)

    # Build vta surface function
    def build_vta(self):
        # Instantiate the points attribute
        self.points = []
        # Add the vertex points to the surface's list of points
        for vert in self.verts:
            self.points.append(vert.loc)
        self.perimeter = self.points
        # Calculate the center of the surface
        self.center = calc_com(points=self.points)
        self.points.append(self.center)
        # Calculate the vector between the atoms
        self.rn = np.array(self.atoms[1].loc) - np.array(self.atoms[0].loc)
        # Calculate the simplices
        find_simps(self)

    def find_flat_points(self):
        # Check to see if the surface is flat or not.
        if self.flat:
            self.tris = [[i, (i + 1) % len(self.perimeter), len(self.points) - 1] for i in range(len(self.perimeter))]
            return
        # Copy the surface points
        points = self.points.copy()
        # Move all surf points toward the origin via center point
        for i in range(len(points)):
            points[i] = np.array(points[i]) - np.array(self.center)
        # Calculate the angles to rotate the center point around
        nps = rotate_points(self.rn, points)
        # Get the 2d version of the points and their Delaunay tesselation
        nps = np.array(nps)
        # Add the flat points to the surface's list of flat points
        self.flat_points = [nps[i, :2] for i in range(len(self.points))]
