from System.Network.surf_funcs import *


class Surface:
    """Surface object. Holds the mesh data. Used to analyze."""
    def __init__(self, atoms, edges=None, verts=None):
        self.func = None
        self.atoms = atoms  # List of Atom type objects
        self.edges = edges  # List of Edge type objects
        self.verts = verts
        self.vert_points = []
        self.edge_points = []
        self.surf_points = []  # List of points on the surface
        self.points = []
        self.simps = None
        self.sa = None
        self.calc_surf()

    # Bisector function. Creates a bisector surface between 2 atoms
    def calc_surf(self):
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

        self.func = ABC + DEF + GHI + [J] + [K] + [d]

    # Build method. Makes the mesh for the surface
    def build(self):
        make_mesh(self, min_dist=0.1)
