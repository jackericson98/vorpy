import numpy as np

from System.calcs import *
from System.Network.surface import Surface


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms, verts, net):

        self.atoms = atoms  # List of Atom type objects
        self.verts = verts  # List of Vertex type objects
        self.surfs = []  # List of 2 surfaces attached to the edge
        if net is not None:
            self.net = net
            self.ndx = [net.atoms.index(atom) for atom in self.atoms]
        self.loc = None  # Location of the center of the 3 atoms that make up the edge
        self.rad = None  # Radius of the inscribed circle of the three atoms
        self.points = []  # List of points on the edge. These points do not include the vertex points

    # Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
    def build(self):
        # Reset the edges points
        self.points = []
        # Get the network's minimum distance
        min_dist = self.net.sys.min_dist
        pv0, pv1 = np.array(self.verts[0].loc), np.array(self.verts[1].loc)
        # Doublet catch
        if self.verts[0].doublet or self.verts[1].doublet:
            # If both vertices are doublets we have to find the closest two vertex locations
            if self.verts[0].doublet and self.verts[1].doublet:
                # Get the backup locations for the vertices
                pv0_, pv1_ = np.array(self.verts[0].loc2), np.array(self.verts[1].loc2)
                # Find the minimum distance
                ds = [calc_dist(pv0, pv1), calc_dist(pv0, pv1_), calc_dist(pv0_, pv1), calc_dist(pv0_, pv1_)]
                ndx = ds.index(min(ds))
                # If the minimum distance comes from the second half, replace the vertex location for pv0
                if ndx > 1:
                    pv0 = pv0_
                # If the ndx is odd, replace the vertex location for pv1
                if ndx % 2 == 1:
                    pv1 = pv1_
            # If only v0 is a doublet, find the closest vertex location to v0
            elif self.verts[0].doublet:
                if calc_dist(pv0, pv1) > calc_dist(self.verts[0].loc2, pv1):
                    pv0 = np.array(self.verts[0].loc2)
            # If only v1 is a doublet, find the closest vertex location to v0
            elif self.verts[1].doublet:
                if calc_dist(pv0, pv1) > calc_dist(pv0, self.verts[1].loc2):
                    pv1 = np.array(self.verts[1].loc2)


        # If the edge is completely straight add 1 point in the middle and return
        if self.atoms[0].rad == self.atoms[1].rad and self.atoms[1].rad == self.atoms[2].rad:
            r = pv1 - pv0
            num_points = int(np.linalg.norm(r) // min_dist) + 2
            for i in range(num_points):
                self.points.append(pv0 + r * (i / num_points))
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
        # Get the center point of the edge and the bottleneck
        self.loc, self.rad = calc_circ(self.atoms)
        # Determine if the theoretical center of the edge is inside the vertices or not
        dr = 1
        if calc_dist(self.loc, pv0) < r_mag or calc_dist(self.loc, pv1) < r_mag:
            dr = -1
        # Find the vector normal to the projection plane
        P_norm = dr * np.cross(np.array(self.loc) - np.array(pc01), np.array(pv1) - np.array(pc01))
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

    # Project method. Projects a point onto the surface using a reference point
    def project(self, rn, pa, surf):
        # Check to see if the surface has its function values
        if surf.func is None:
            surf.calc_func()
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
        if round(b ** 2 - 4 * a * c, 10) >= 0:
            # Calculate the roots
            roots = np.roots([a, b, c])
            # If one root exists return it
            if len(roots) == 1:
                return pa + roots[0] * rn
            else:
                p1 = pa + min(roots) * rn
                p2 = pa + max(roots) * rn
            # If the point we are calculating is the first in the edge choose the one closest to the vertex
            if len(self.points) == 1:
                point = p1
                if calc_dist(p2, self.points[0]) <= calc_dist(p1, self.points[0]):
                    point = p2
            # If we have 2 points to choose from, choose the one that makes the angle closer to 180
            else:
                point = p1
                if calc_angle(self.points[-1], self.points[-2], p2) >= calc_angle(self.points[-1], self.points[-2], p1):
                    point = p2
            # Return the point we choose
            return point
