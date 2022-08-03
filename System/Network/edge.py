from System.sys_funcs import *
from System.Network.surface import Surface


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms, verts, surfs=None, calc_points=True):
        if surfs is None:
            surfs = []
        self.atoms = atoms  # List of Atom type objects
        self.verts = verts  # List of Vertex type objects
        self.surfs = surfs
        self.loc = None
        self.rad = None
        self.dir = None
        self.points = []  # List of points on the edge. These points do not include the vertex points
        if calc_points:
            self.build()

    # Calculate points function. Find points long the edge, given it has atoms and vertices
    def build(self, surf=None, min_dist=None):
        # Get the atoms
        a0, a1, a2 = self.atoms[0], self.atoms[1], self.atoms[2]
        # Give the edge a minimum distance
        if min_dist is None:
            min_dist = 0.1
        # Grab the vertex points
        pv0, pv1 = np.array(self.verts[0].loc), np.array(self.verts[1].loc)
        # Find the point in between the two vertex points
        r01 = pv1 - pv0
        d = calc_dist(pv0, pv1)
        rn01 = r01 / d
        # Determine the number of points to make based on the minimum distance
        num_pts = max(int(d / min_dist), 10)
        # Find the incremental change
        dr = d / num_pts
        proj_pts = [pv0]
        for i in range(num_pts):
            proj_pts.append(proj_pts[-1] + dr * rn01)
        # If no surface is given, choose one
        if surf is None:
            if round(a0.rad, 10) == round(a1.rad, 10):
                surf = Surface([a1, a2])
            elif round(a1.rad, 10) == round(a2.rad, 10):
                surf = Surface([a0, a1])
            else:
                surf = Surface([a0, a2])
        # Grab the function's coefficients
        f = surf.func
        # Get the center point of the vertices
        pc01 = pv0 + 0.5 * rn01 * d
        # Get the center point of the edge
        circ = calc_circ(self.atoms)
        c, bn = circ[0], circ[1]
        # Find the vector perpendicular to the vectors between the center of the edge, pv1 and pc01
        P_norm = dr * np.cross(np.array(c) - np.array(pc01), np.array(pv1) - np.array(pc01))
        # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
        rpcr = - np.cross(P_norm, rn01)
        rnpcr = rpcr / np.linalg.norm(rpcr)
        # Check to see if the unit vector is pointing toward or away from the edge
        if np.dot(rnpcr, np.array(surf.atoms[0].loc) - np.array(surf.atoms[1].loc)) < 0:
            d = -d
        # Calculate the reference point
        pa = pc01 + 0.5 * d * rnpcr
        # Add the first vertex to the list of points
        self.points = [pv0.tolist()]
        # Find the edges points
        for i in range(num_pts):
            # Get the vector between the base point and the projection point
            r_proj = proj_pts[i] - pa
            rn_proj = r_proj/np.linalg.norm(r_proj)
            # Project the vector onto the surface
            surf_point = self.project(rn_proj, pa, f)
            if surf_point is None:
                break
            self.points.append(surf_point)
        self.points.append(pv1.tolist())

    @staticmethod
    def project(rn, pa, f):
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
        if b ** 2 - 4 * a * c >= 0:
            # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
            roots = np.roots([a, b, c])
            if len(roots) == 0:
                return
            mag = min([abs(roots[i]) for i in range(len(roots))])
            return pa + mag * rn
