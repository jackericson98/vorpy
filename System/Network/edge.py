from System.calcs import *
from System.Network.surface import Surface


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms, verts, net, calc_points=True):
        self.atoms = atoms  # List of Atom type objects
        self.verts = verts  # List of Vertex type objects
        self.surfs = []  # List of 2 surfaces attached to the edge
        self.min_dist = net.min_dist
        if net is not None:
            self.net = net
            self.ndx = [net.atoms.index(atom) for atom in self.atoms]
        self.loc = None  # Location of the center of the 3 atoms that make up the edge
        self.rad = None  # Radius of the inscribed circle of the three atoms
        self.dir = None  # Direction along the edge from the edge's first vertex
        self.points = []  # List of points on the edge. These points do not include the vertex points
        if calc_points:  # Build the edge if calc_points is true
            self.build()

    # Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
    def build(self, min_dist=0.1):
        # Grab the vertex points and make numpy arrays out of them
        pv0, pv1 = np.array(self.verts[0].loc), np.array(self.verts[1].loc)
        # If the edge is completely straight add 1 point in the middle and return
        if self.atoms[0].rad == self.atoms[1].rad and self.atoms[1].rad == self.atoms[2].rad:
            r = pv1 - pv0
            self.points.append(pv0 + r/2)
            return
        # If no surface is given, choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
        if round(self.atoms[0].rad, 10) == round(self.atoms[1].rad, 10):
            surf = Surface(self.atoms[1:], self.net)
        else:
            surf = Surface(self.atoms[:2], self.net)
        # Find the point in between the two vertex points
        r01 = pv1 - pv0  # Vector between vertices
        r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
        rn01 = r01 / r_mag  # Normal to the vector between the vertices
        pc01 = pv0 + 0.5 * rn01 * r_mag  # Center point
        # Get the center point of the edge
        circ = calc_circ(self.atoms)
        c, bn = circ[0], circ[1]
        # Determine if the theoretical center of the edge is inside the vertices or not
        dr = 1
        if calc_dist(c, pv0) < r_mag or calc_dist(c, pv1) < r_mag:
            dr = -1
        # Find the vector normal to the projection plane
        P_norm = dr * np.cross(np.array(c) - np.array(pc01), np.array(pv1) - np.array(pc01))
        # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
        rpcr = - np.cross(P_norm, rn01)
        rnpcr = rpcr / np.linalg.norm(rpcr)
        # Calculate the reference point
        pa = pc01 + 0.5 * r_mag * rnpcr
        # Find the number of points
        n = max(int(r_mag / min_dist), 2)
        # Calculate the angle between the vertices and the reference point
        theta = calc_angle(pa, pv0, pv1)
        A = theta / n
        # Add the first vertex to the list of points
        self.points = [pv0.tolist()]
        # Find the edges points. Don't count the vertex
        for i in range(n-1):
            # Set pb to the previous point
            pb = self.points[-1]
            # Get the distance between pb and pa for c
            c = calc_dist(pb, pa)
            # Get the angle between pb, pa and pb + rno1
            B = calc_angle(pb, pb + rn01, pa)
            # Get the last angle
            C = np.pi - B - A
            # Get the distance to our projection point or 'a' on our triangle
            a = np.sin(A) * c / np.sin(C)
            # Use that distance to project rn01 from pb to find our projection point or pc
            pc = pb + a * rn01
            # Get the vector from pa to pc
            rac = np.array(pc) - np.array(pa)
            rnac = rac / np.linalg.norm(rac)
            # Project the vector onto the surface
            surf_point = self.project(rnac, pa, surf)
            if surf_point is None:
                break
            self.points.append(surf_point)

    @staticmethod
    def project(rn, pa, surf):
        # Get the function values
        f, a0, a1 = surf.func, surf.atoms[0], surf.atoms[1]
        # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
        a = f[0] * rn[0] ** 2 + f[1] * rn[1] ** 2 + f[2] * rn[2] ** 2 + f[3] * rn[0] * rn[1] + f[4] * rn[
            1] * rn[
                2] + f[5] \
            * rn[2] * rn[0]
        b = 2 * f[0] * rn[0] * pa[0] + 2 * f[1] * rn[1] * pa[1] + 2 * f[2] * rn[2] * pa[2] + f[3] \
            * (rn[0] * pa[1] + rn[1] * pa[0]) + f[4] * (rn[1] * pa[2] + rn[2] * pa[1]) + f[5] \
            * (rn[2] * pa[0] + rn[0] * pa[2]) + f[6] * rn[0] + f[7] * rn[1] + f[8] * rn[2]
        c = f[0] * pa[0] ** 2 + f[1] * pa[1] ** 2 + f[2] * pa[2] ** 2 + f[3] * pa[0] * pa[1] + f[4] * pa[1] * pa[
            2] + \
            f[5] * pa[2] * pa[0] + f[6] * pa[0] + f[7] * pa[1] + f[8] * pa[2] + f[9]
        # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
        # and add that point to our surface list of points
        if b ** 2 - 4 * a * c > 0:
            # Calculate the roots
            roots = np.roots([a, b, c])
            # If no roots exist return
            if len(roots) == 0:
                return
            # If one root exists return it
            elif len(roots) == 1:
                return pa + roots[0] * rn
            # If the smallest root is negative (i.e. incorrect) return the other root
            if min(roots) < 0:
                return pa + rn * max(roots)
            # Otherwise, return the smaller of the two
            return pa + min(roots) * rn


# Calculate circle function. Takes in 3 atoms, calculates the center and radius of inscribed circle and returns them
def calc_circ(atoms):
    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad
    # Get the relevant variables
    R2, R3 = atoms[1].rad, atoms[2].rad
    x2, y2, z2 = atoms[1].loc[0] - l1[0], atoms[1].loc[1] - l1[1], atoms[1].loc[2] - l1[2]
    x3, y3, z3 = atoms[2].loc[0] - l1[0], atoms[2].loc[1] - l1[1], atoms[2].loc[2] - l1[2]
    # Calculate coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (R1 - R2), R1 ** 2 - R2 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (R1 - R3), R1 ** 2 - R3 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3 = y2*z3 - z2*y3, z2*x3 - x2*z3, x2*y3 - y2*x3
    # More coefficients
    F = a3*b2*c1 - a2*b3*c1 - a3*b1*c2 + a1*b3*c2 + a2*b1*c3 - a1*b2*c3
    Fx0 = b3*c2*f1 - b2*c3*f1 - b3*c1*f2 + b1*c3*f2
    Fx1 = b3*c2*d1 - b2*c3*d1 - b3*c1*d2 + b1*c3*d2
    Fy0 = - a3*c2*f1 + a2*c3*f1 + a3*c1*f2 - a1*c3*f2
    Fy1 = - a3*c2*d1 + a2*c3*d1 + a3*c1*d2 - a1*c3*d2
    Fz0 = a3*b2*f1 - a2*b3*f1 - a3*b1*f2 + a1*b3*f2
    Fz1 = a3*b2*d1 - a2*b3*d1 - a3*b1*d2 + a1*b3*d2
    # Catch for F=0 (i.e. no circle exists)
    if F == 0:
        return
    # Find the radius of the tangential circle using the quadratic formula
    a = (Fx1 ** 2 + Fy1 ** 2 + Fz1 ** 2) / F ** 2 - 1
    b = 2 * (Fx0 * Fx1 + Fy0 * Fy1 + Fz0 * Fz1) / F ** 2 - 2 * R1
    c = (Fx0 ** 2 + Fy0 ** 2 + Fz0 ** 2) / F ** 2 - R1 ** 2
    # Calculate the discriminant.
    disc = b ** 2 - 4 * a * c
    # If the discriminant is negative then the tangential circle does not exist.
    if disc > 0:
        # Grab the two roots
        Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
        if len(Rs) > 1 and abs(Rs[1]) < abs(Rs[0]):
            Rs[0], Rs[1] = Rs[1], Rs[0]
        R = Rs[0]
        # Go through each circle and gather its points
        # Calculate the vertex based off of our coefficient values and the sphere's radius
        x = Fx0 / F + R * Fx1 / F + l1[0]
        y = Fy0 / F + R * Fy1 / F + l1[1]
        z = Fz0 / F + R * Fz1 / F + l1[2]
        return [[x, y, z], R]
