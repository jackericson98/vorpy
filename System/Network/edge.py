from System.sys_calcs import *


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms, verts, surfs=None):
        if surfs is None:
            surfs = []
        self.atoms = atoms  # List of Atom type objects
        self.verts = verts  # List of Vertex type objects
        self.surfs = surfs
        self.loc = None
        self.rad = None
        self.dir = None
        self.points = []  # List of points on the edge. These points do not include the vertex points

    # Calculate points function. Find points long the edge, given it has atoms and vertices
    def calc_points(self, surf=None, min_dist=None):
        # Give the edge a minimum distance
        if min_dist is None:
            min_dist = 0.1
        # I want to be able to calculate a surface here, but don't want to cross pollinate imports on this level
        if surf is None:
            return
        # Grab the function's coefficients
        f = surf.func
        # Grab the vertex points
        pv0, pv1 = np.array(self.verts[0].loc), np.array(self.verts[1].loc)
        # Find the point in between the two vertex points
        r01 = pv1 - pv0
        d = np.linalg.norm(r01)
        rn01 = r01 / d
        # Get the center point of the vertices
        pc01 = pv0 + 0.5 * rn01 * d
        # Get the center point of the edge
        circ = calc_circ(self.atoms)[0]
        pce, rad = circ[0], circ[1]
        # Determine if the center if the edge is inside the vertices or not
        dr = 1
        if calc_dist(pce, pv0) < d or calc_dist(pce, pv1) < d:
            dr = -1
        # Find the vector perpendicular to the vectors between the center of the edge, pv1 and pc01
        P_norm = dr * np.cross(np.array(pce) - np.array(pc01), np.array(pv1) - np.array(pc01))
        # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
        rpcr = np.cross(P_norm, rn01)
        rnpcr = rpcr / np.linalg.norm(rpcr)
        # Calculate the reference point
        pa = pc01 + 0.5 * d * rnpcr
        # Find the number of points
        n = max(int(d / min_dist), 10)
        # Calculate the angle between the vertices and the reference point
        theta = calc_angle(pv0, pv1, pa)
        A = theta / n
        print(A)
        # Set pb to None
        pb = None
        # Find the edges points
        for i in range(n):
            # If this is the first loop set pb equal to the first vertex
            if pb is None:
                pb = pv0
            # Else grab the last point made
            else:
                pb = self.points[-1]
            # Get the distance between pb and pa for c
            c = calc_dist(pb, pa)
            # Get the angle between pb, pa and pb + rno1
            B = calc_angle(pa, pb + rn01, pb)
            # Get the last angle
            C = 180 - B - A
            # Get the distance to our projection point or 'a' on our triangle
            a = np.sin(A) * c / np.sin(C)
            # Use that distance to project rn01 from pb to find our projection point or pc
            pc = pb + a * rn01
            # Get the vector from pa to pc
            rac = np.array(pc) - np.array(pa)
            rnac = rac / np.linalg.norm(rac)
            # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
            a = f[0] * rnac[0] ** 2 + f[1] * rnac[1] ** 2 + f[2] * rnac[2] ** 2 + f[3] * rnac[0] * rnac[1] + f[4] * rnac[1] * rnac[
                2] + f[5] \
                * rnac[2] * rnac[0]
            b = 2 * f[0] * rnac[0] * pa[0] + 2 * f[1] * rnac[1] * pa[1] + 2 * f[2] * rnac[2] * pa[2] + f[3] \
                * (rnac[0] * pa[1] + rnac[1] * pa[0]) + f[4] * (rnac[1] * pa[2] + rnac[2] * pa[1]) + f[5] \
                * (rnac[2] * pa[0] + rnac[0] * pa[2]) + f[6] * rnac[0] + f[7] * rnac[1] + f[8] * rnac[2]
            c = f[0] * pa[0] ** 2 + f[1] * pa[1] ** 2 + f[2] * pa[2] ** 2 + f[3] * pa[0] * pa[1] + f[4] * pa[1] * pa[
                2] + \
                f[5] * pa[2] * pa[0] + f[6] * pa[0] + f[7] * pa[1] + f[8] * pa[2] + f[9]
            # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
            # and add that point to our surface list of points

            if b ** 2 - 4 * a * c >= 0:
                # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
                roots = np.roots([a, b, c])
                if roots == []:
                    print(roots)
                    return
                mag = min([abs(roots[i]) for i in range(len(roots))])
                self.points.append(pa + mag * rnac)
