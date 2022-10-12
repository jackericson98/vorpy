from System.calcs import *


class Vertex:
    """Vertex object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms=None, net=None, location=None, radius=None, loc2=None, rad2=None, doublet=False, ndx=None):

        self.net = net           # Network       :   Network object for the vertex to refer back to
        self.atoms = atoms       # Atoms         :   List of atoms used to construct the vertex
        self.edges = []          # Edges         :   List of Edge type objects connected to the vertex in the network
        self.surfs = []          # Surfaces      :   List of Surface type objects that the vertex is a part of

        self.ndx = ndx           # Index         :   Indices of the atoms in the vertex
        self.loc = location      # Location      :   Where the vertex is located in 3D
        self.rad = radius        # Radius        :   Radius of the vertex's tangential sphere
        self.load_ndxs = []      # Load indices  :   List of object load indices

        self.doublet = doublet   # Doublet       :   Whether the vertex is a doublet
        self.d_type = None       # Doublet type  :   Doublet type it is: 3 edges, 3 surfaces or two edges, one surface
        self.loc2 = loc2         # Location 2    :   Location of the doublet site
        self.rad2 = rad2         # Radius 2      :   Radius of the doublet site's tangential sphere

    # Calculate vertex function. Takes in 4 atoms, calculates the loc and rad of the inscribed sphere and adds the
    def calc_vert(self):
        # If the vertex is mature enough to be calculated, create and sort its indices
        self.ndx = [self.net.atoms.index(atom) for atom in self.atoms]
        self.ndx.sort()
        # The real location and radius of the base sphere
        locs = np.array(self.atoms[0].loc), np.array(self.atoms[1].loc), np.array(self.atoms[2].loc), \
               np.array(self.atoms[3].loc)
        R1, R2, R3, R4 = self.atoms[0].rad, self.atoms[1].rad, self.atoms[2].rad, self.atoms[3].rad
        # Find the recalculated location of the atoms
        l0, l1, l2, l3 = locs[0], locs[1] - locs[0], locs[2] - locs[0], locs[3] - locs[0]
        # Calculate our System of linear equations coefficients
        a1, b1, c1, d1, f1 = 2*l1[0], 2*l1[1], 2*l1[2], 2*(R2 - R1), R1**2 - R2**2 + l1[0]**2 + l1[1]**2 + l1[2]**2
        a2, b2, c2, d2, f2 = 2*l2[0], 2*l2[1], 2*l2[2], 2*(R3 - R1), R1**2 - R3**2 + l2[0]**2 + l2[1]**2 + l2[2]**2
        a3, b3, c3, d3, f3 = 2*l3[0], 2*l3[1], 2*l3[2], 2*(R4 - R1), R1**2 - R4**2 + l3[0]**2 + l3[1]**2 + l3[2]**2
        A, B, C, d, f = [a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f2]
        # Calculate the ranks of the matrices
        ABC_rank = np.linalg.matrix_rank(np.array([A, B, C]))
        m_rank = np.linalg.matrix_rank(np.array([A, B, C, d]))
        f_rank = np.linalg.matrix_rank(np.array([A, B, C, d, f]))
        # Calculate the F values
        F = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1
        F10 = b1 * c2 * f3 - b1 * c3 * f2 - b2 * c1 * f3 + b2 * c3 * f1 + b3 * c1 * f2 - b3 * c2 * f1
        F11 = -b1 * c2 * d3 + b1 * c3 * d2 + b2 * c1 * d3 - b2 * c3 * d1 - b3 * c1 * d2 + b3 * c2 * d1
        F20 = -a1 * c2 * f3 + a1 * c3 * f2 + a2 * c1 * f3 - a2 * c3 * f1 - a3 * c1 * f2 + a3 * c2 * f1
        F21 = a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1
        F30 = a1 * b2 * f3 - a1 * b3 * f2 - a2 * b1 * f3 + a2 * b3 * f1 + a3 * b1 * f2 - a3 * b2 * f1
        F31 = -a1 * b2 * d3 + a1 * b3 * d2 + a2 * b1 * d3 - a2 * b3 * d1 - a3 * b1 * d2 + a3 * b2 * d1
        verts = []
        xs, ys, zs, Rs = [], [], [], []

        # Case 1:
        if ABC_rank == 3 and m_rank == 3 and f_rank == 3:
            # Calculate the radius polynomial coefficients
            a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F ** 2) - 1
            b = (2 * (F10 * F11 + F20 * F21 + F30 * F31) / F ** 2) - 2 * R1
            c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F ** 2) - R1 ** 2
            # If the discriminant is positive, find the real positive roots of the quadratic
            if round(-4 * a * c + b ** 2, 10) >= 0:
                Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for R in Rs:
                x, y, z = F10 / F + R * F11 / F, F20 / F + R * F21 / F, F30 / F + R * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
        # Case 2:
        elif ABC_rank == 2 and m_rank == 3 and f_rank == 3 and F > 0:
            # Calculate the _ polynomial coefficients
            a = F ** 2 + F11 ** 2 + F21 ** 2 - F31 ** 2
            b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F20 ** 2 - (F30 + F * R1)
            roots, verts = [], []
            # Check the discriminant
            disc = -4 * a * c + b ** 2
            if disc > 0:
                roots = [root for root in np.roots([a, b, c]) if np.isreal(root)]
            # Case 2 subcases:
            # Case 2.1
            if np.linalg.matrix_rank([A, B, d]) == 3:
                # Go through each radius and calculate the vertex
                for z in roots:
                    x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
                    # Move the vertex back to the actual location of the atoms
                    verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
            # Case 2.2
            elif np.linalg.matrix_rank([A, d, C]) == 3:
                # Go through each radius and calculate the vertex
                for y in roots:
                    x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
                    # Move the vertex back to the actual location of the atoms
                    verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
            # Case 2.3
            elif np.linalg.matrix_rank([d, B, C]):
                # Go through each radius and calculate the vertex
                for x in roots:
                    R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
                    # Move the vertex back to the actual location of the atoms
                    verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
        # If one root exists return it
        if len(verts) == 1:
            self.loc, self.rad = verts[0][0], verts[0][1]
        # If two roots exist:
        elif len(verts) == 2:
            # If both radii are negative (I'm not sure if this is possible, but let's catch it anyway)
            if verts[0][1] < 0 and verts[1][1] < 0:
                return
            # Check to see if the second vertex is negative (i.e. not a doublet)
            elif verts[0][1] < 0 < verts[1][1]:
                self.loc, self.rad = verts[1][0], verts[1][1]
            elif verts[1][1] < 0 < verts[0][1]:
                self.loc, self.rad = verts[0][0], verts[0][1]
            # If both radii are positive we have a doublet. Choose the smaller vertex to be the lead vertex and set loc2
            elif verts[0][1] < verts[1][1]:
                self.loc, self.rad = verts[0][0], verts[0][1]
                self.loc2, self.rad2 = verts[1][0], verts[1][1]
            elif verts[1][1] < verts[0][1]:
                self.loc, self.rad = verts[1][0], verts[1][1]
                self.loc2, self.rad2 = verts[0][0], verts[0][1]
