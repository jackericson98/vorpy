from System.sys_funcs.calcs import *


class Vertex:
    """Vertex object. Used to build the network and calculate the surfaces"""

    def __init__(self, atoms=None, net=None, location=None, radius=None, loc2=None, rad2=None, doublet=None, ndx=None,
                 load_ndxs=None, flat_faces=False, edges=None, surfaces=None, distance=None):

        self.net = net                # Network       :   Network object for the vertex to refer back to
        self.atoms = atoms            # Atoms         :   List of atoms used to construct the vertex
        self.edges = edges            # Edges         :   List of Edge type objects connected to the vertex
        self.surfs = surfaces         # Surfaces      :   List of Surface type objects that the vertex is a part of

        self.ndx = ndx                # Index         :   Indices of the atoms in the vertex
        self.loc = location           # Location      :   Where the vertex is located in 3D
        self.rad = radius             # Radius        :   Radius of the vertex's tangential sphere
        self.loc_points = None
        self.loc_tris = None
        self.sphere_points = None
        self.sphere_tris = None
        self.box = None               # Box index     :   Sub box index for sorting
        self.dist = distance
        self.load_ndxs = load_ndxs    # Load indices  :   List of object load indices

        self.doublet = doublet        # Doublet       :   Whether the vertex is a doublet
        self.d_type = None            # Doublet type  :   Doublet type it is: 3 edges, 3 surfs or two edges, one surf
        self.loc2 = loc2              # Location 2    :   Location of the doublet site
        self.rad2 = rad2              # Radius 2      :   Radius of the doublet site's tangential sphere

        self.flat_faced = flat_faces  # Flat Faced?   :   Was this vertex constructed with flat faces in mind?
        self.num_edges = 0

        # If the vertex is mature enough to be calculated, create and sort its indices
        if self.net is not None and self.atoms is not None:
            self.ndx = [atom.num for atom in self.atoms]
            self.ndx.sort()

        # Set up the edge and surfaces lists
        if self.edges is None:
            self.edges = []
        if self.surfs is None:
            self.surfs = []

    # Calculate vertex function. Takes in 4 atoms, calculates the loc and rad of the inscribed sphere and adds the
    def calc_vert(self):
        # The real location and radius of the base sphere
        locs = np.array(self.atoms[0].loc), np.array(self.atoms[1].loc), np.array(self.atoms[2].loc), \
               np.array(self.atoms[3].loc)
        r0, r1, r2, r3 = self.atoms[0].rad, self.atoms[1].rad, self.atoms[2].rad, self.atoms[3].rad
        r0_2 = r0 ** 2
        # Find the recalculated location of the atoms
        l0, l1, l2, l3 = locs[0], locs[1] - locs[0], locs[2] - locs[0], locs[3] - locs[0]
        # Calculate our System of linear equations coefficients
        a1, b1, c1, d1, f1 = 2 * l1[0], 2 * l1[1], 2 * l1[2], 2 * (r1 - r0), r0_2 - r1 ** 2 + l1[0] ** 2 + l1[
            1] ** 2 + l1[2] ** 2
        a2, b2, c2, d2, f2 = 2 * l2[0], 2 * l2[1], 2 * l2[2], 2 * (r2 - r0), r0_2 - r2 ** 2 + l2[0] ** 2 + l2[
            1] ** 2 + l2[2] ** 2
        a3, b3, c3, d3, f3 = 2 * l3[0], 2 * l3[1], 2 * l3[2], 2 * (r3 - r0), r0_2 - r3 ** 2 + l3[0] ** 2 + l3[
            1] ** 2 + l3[2] ** 2
        # Calculate the F values
        F = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1
        F_2 = F ** 2
        F10 = b1 * c2 * f3 - b1 * c3 * f2 - b2 * c1 * f3 + b2 * c3 * f1 + b3 * c1 * f2 - b3 * c2 * f1
        F11 = -b1 * c2 * d3 + b1 * c3 * d2 + b2 * c1 * d3 - b2 * c3 * d1 - b3 * c1 * d2 + b3 * c2 * d1
        F20 = -a1 * c2 * f3 + a1 * c3 * f2 + a2 * c1 * f3 - a2 * c3 * f1 - a3 * c1 * f2 + a3 * c2 * f1
        F21 = a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1
        F30 = a1 * b2 * f3 - a1 * b3 * f2 - a2 * b1 * f3 + a2 * b3 * f1 + a3 * b1 * f2 - a3 * b2 * f1
        F31 = -a1 * b2 * d3 + a1 * b3 * d2 + a2 * b1 * d3 - a2 * b3 * d1 - a3 * b1 * d2 + a3 * b2 * d1
        # Calculate the ranks of the matrices
        m_rank, f_rank = 3, 3
        if F == 0:
            my_mtx = [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f3]]
            m_rank = np.linalg.matrix_rank(np.array(my_mtx[:-1]))
            if m_rank != 3:
                f_rank = np.linalg.matrix_rank(np.array(my_mtx))
        verts = []
        xs, ys, zs, Rs = [], [], [], []
        # Case 1:
        if F != 0:
            # Calculate the radius polynomial coefficients
            a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F_2) - 1
            b = 2 * (((F10 * F11 + F20 * F21 + F30 * F31) / F_2) - r0)
            c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F_2) - r0 ** 2
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 >= 0:
                Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for R in Rs:
                x, y, z = F10 / F + R * F11 / F, F20 / F + R * F21 / F, F30 / F + R * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
        # Case 2:
        elif a1 * b2 - a2 * b1 != 0 and m_rank == 3 and f_rank == 3 and F > 0:
            # Calculate the _ polynomial coefficients
            a = F_2 + F11 ** 2 + F21 ** 2 - F31 ** 2
            b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * r0)
            c = F10 ** 2 + F20 ** 2 - (F30 + F * r0)
            roots, verts = [], []
            # Check the discriminant
            disc = -4 * a * c + b ** 2
            if disc > 0:
                roots = [root for root in np.roots([a, b, c]) if np.isreal(root)]
            # Case 2 subcases:
            # Case 2.1
            if F31 != 0:
                # Go through each radius and calculate the vertex
                for z in roots:
                    x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
                    # Move the vertex back to the actual location of the atoms
                    verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
            # Case 2.2
            elif F21 != 0:
                # Go through each radius and calculate the vertex
                for y in roots:
                    x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
                    # Move the vertex back to the actual location of the atoms
                    verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
            # Case 2.3
            elif F11 != 0:
                # Go through each radius and calculate the vertex
                for x in roots:
                    R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
                    # Move the vertex back to the actual location of the atoms
                    verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
        loc, rad, loc2, rad2 = None, None, None, None
        # If one root exists return it
        if len(verts) == 1:
            loc, rad = verts[0][0], verts[0][1]
        # If two roots exist:
        elif len(verts) == 2:
            # Get the largest atom's radius
            max_atom_rad = max([r0, r1, r2, r3])
            # Set the locations and radii, so that the smaller vertex is first
            if abs(verts[0][1]) > abs(verts[1][1]):
                verts[0], verts[1] = verts[1], verts[0]
            # Set the locations and radii variables
            locs, rads = [verts[0][0], verts[1][0]], [verts[0][1], verts[1][1]]
            # If either radii are negative (I'm not sure if this is possible, but let's catch it anyway)
            if rads[0] < 0 or rads[1] < 0:
                if rads[0] > 0 or abs(rads[0]) < max_atom_rad:
                    loc, rad = locs[0], rads[0]
                    if rads[1] > 0 or abs(rads[1]) < max_atom_rad:
                        loc2, rad2 = locs[1], rads[1]
                elif rads[1] > 0 or abs(rads[1]) < max_atom_rad:
                    loc, rad = locs[1], rads[1]
            # If both radii are positive we have a doublet. Choose the smaller vertex to be the lead vertex and set loc2
            else:
                loc, loc2, rad, rad2 = locs[0], locs[1], rads[0], rads[1]
        self.loc, self.rad, self.loc2, self.rad2 = loc, rad, loc2, rad2

    def calc_flat_vert(self, power=False):
        """
        Calculates the flat vertex between 4 atoms by finding the intersection of the mid-point planes between the first
        atom and the others
        :param power:
        :return:
        """
        rads = [_.rad for _ in self.atoms]
        atom_rads = [x for _, x in sorted(zip(rads, self.atoms), key=lambda pair: pair[0])]
        # Get the plane equations
        coeffs = []
        # Go through the atoms to make the planes
        for an in atom_rads[1:]:
            # Get the point between the atoms
            r = np.array(an.loc) - np.array(atom_rads[0].loc)
            norm = np.linalg.norm(r)
            rn = r / norm
            if power:
                d0 = 0.5 * (norm ** 2 + atom_rads[0].rad ** 2 - an.rad ** 2) / norm
                center = atom_rads[0].loc + d0 * rn
            else:
                center = 0.5 * r + np.array(atom_rads[0].loc)
            coeffs.append(rn.tolist() + [np.dot(rn, center)])

        a, b, c, d = coeffs[0]
        e, f, g, h = coeffs[1]
        i, j, k, m = coeffs[2]

        disc = c * f * i - b * g * i - c * e * j + a * g * j + b * e * k - a * f * k
        x_numerator = d * g * j - c * h * j - d * f * k + b * h * k + c * f * m - b * g * m
        y_numerator = - d * g * i + c * h * i + d * e * k - a * h * k - c * e * m + a * g * m
        z_numerator = d * f * i - b * h * i - d * e * j + a * h * j + b * e * m - a * f * m
        x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
        # Get the radius
        if power:
            rad = calc_dist([x, y, z], atom_rads[0].loc) ** 2 - atom_rads[0].rad ** 2
        else:
            rad = calc_dist([x, y, z], atom_rads[0].loc)
        self.loc, self.rad = [x, y, z], rad

    def add_vert_info(self, vert):
        pass
