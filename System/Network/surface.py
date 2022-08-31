from System.Network.build_surface import *


class Surface:
    """Surface object. Holds the mesh data. Used to analyze interfaces between atoms."""
    def __init__(self, atoms, net=None, edges=None, verts=None, min_dist=0.1):
        self.func = None
        self.atoms = atoms  # List of Atom type objects
        self.edges = edges  # List of Edge type objects
        self.verts = verts
        if net is not None:
            self.net = net
            self.ndx = [net.atoms.index(atom) for atom in self.atoms]
        self.perimeter = []
        self.vert_ndxs = []
        self.points = []
        self.tris = None
        self.sa = None
        self.min_dist = min_dist
        self.crit_ang = np.pi
        self.rn = None
        self.calc_func()
        self.flat_points = []

    # Bisector function. Creates a bisector surface between 2 atoms
    def calc_func(self):
        # Make sure that a0 is the atom with the smaller radius
        if self.atoms[0].rad > self.atoms[1].rad:
            self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
        # Create a0, a1 variables
        a0, a1 = self.atoms
        # Set the rn vector for the surface since the atoms are sorted
        self.rn = np.array(a1.loc) - np.array(a0.loc)
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
    def build(self, min_dist=None, simps=True):
        # Set the minimum distance
        if min_dist is None:
            min_dist = self.net.min_dist
        # Build the mesh
        make_mesh(self, min_dist)
        # Calculate the simplices
        if simps:
            find_simps(self)

    # Build vta surface function
    def build_vta(self):
        # Add the vertex points to the surface's list of points
        for vert in self.verts:
            self.points.append(vert.loc)
        # Calculate the simplices
        find_simps(self)
