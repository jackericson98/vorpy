from System.sys_funcs.calcs import calc_circ, calc_angle
from System.Network.net_objs.surface import Surface
import numpy as np


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms=None, net=None, verts=None, surfs=None, doublet=False, points=None, loc=None, rad=None,
                 rn=None, dist=None, pv0=None, pv1=None, ndx=None, load_ndxs=None, point_refs=None, straight=False):

        # If no network was given have a catch
        if net is not None and net.atoms is not None and atoms is not None:
            ndx = [atom.num for atom in atoms]
            ndx.sort()
        self.ndx = ndx                   # Index         :   Indices of the atoms of the surface
        self.net = net                   # Network       :   Network of the System
        self.atoms = atoms               # Atoms         :   List of Atom type objects for the edge
        self.verts = verts               # Vertices      :   List of Vertex type objects
        self.surfs = surfs               # Surfaces      :   List of 2 surfaces attached to the edge
        self.load_ndxs = load_ndxs       # Load indices  :   List of object load indices
        self.point_refs = point_refs     # Point refs    :   List of the surface and indices of the points for the edge

        self.loc = loc                   # Location      :   Location of the center of the 3 atoms that make up the edge
        self.rad = rad                   # Radius        :   Radius of the inscribed circle of the three atoms
        self.rn = rn
        self.dist = dist
        self.points = points             # Points        :   List of points along the
        self.pv0 = pv0                   # Vertex pt 0   :   The points on the ends of the edges
        self.pv1 = pv1                   # Vertex pt 1   :   The points on the ends of the edges
        self.pa = None                   # Projection pt :   The projection point from which the edge is built
        self.doublet = doublet           # Doublet       :   Boolean for if the edge is part of a doublet or not
        self.loc2 = None                 # Loc2          :   Allows edges to be checked like vertices

        self.ref = None                  # Reference     :   Tuple holding a surface and a range for efficient storage
        self.straight = straight         # Straight edge :   Straight edge or not

    # Get location method. Calculates the circle made between the atoms
    def get_loc(self):
        # Get the center point of the edge and the bottleneck
        circ = calc_circ(self.atoms)
        if circ is not None:
            self.loc, self.rad = circ

    # Find projection values. Calculates the correct end and projection points for the edge
    def find_pvals(self):

        # Typical case, no doublets
        self.pv0, self.pv1 = np.array(self.verts[0].loc), np.array(self.verts[1].loc)

        if self.straight:
            return

        # Get the projection point
        # Find the point in between the two vertex points
        r01 = self.pv1 - self.pv0  # Vector between vertices
        r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
        rn01 = r01 / r_mag  # Normal to the vector between the vertices
        pc01 = self.pv0 + 0.5 * rn01 * r_mag  # Center point

        # Determine if the theoretical center of the edge is inside the vertices or not
        dr = 1
        if np.sqrt(sum(np.square(np.array(self.loc) - np.array(self.pv0)))) < r_mag or \
                np.sqrt(sum(np.square(np.array(self.loc) - np.array(self.pv1)))) < r_mag:
            dr = -1

        # Find the vector normal to the projection plane
        P_norm = dr * np.cross(np.array(self.loc) - np.array(pc01), np.array(self.pv1) - np.array(pc01))
        # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
        rpcr = - np.cross(P_norm, rn01)
        rnpcr = rpcr / np.linalg.norm(rpcr)
        # Calculate the reference point
        self.pa = pc01 + 2 * r_mag * rnpcr

    # Project method. Projects a point onto the surface using a reference point
    def project(self, rn, pa, surf):
        # Check to see if the surface has its function values
        if surf.func is None:
            surf.calc_func()
        # Get the function values
        f, a0, a1 = surf.func, surf.atoms[0], surf.atoms[1]
        # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
        a = f[0] * rn[0] ** 2 + f[1] * rn[1] ** 2 + f[2] * rn[2] ** 2 + f[3] * rn[0] * rn[1] + f[4] * rn[
            1] * rn[2] + f[5] * rn[2] * rn[0]
        b = 2 * f[0] * rn[0] * pa[0] + 2 * f[1] * rn[1] * pa[1] + 2 * f[2] * rn[2] * pa[2] + f[3] \
            * (rn[0] * pa[1] + rn[1] * pa[0]) + f[4] * (rn[1] * pa[2] + rn[2] * pa[1]) + f[5] \
            * (rn[2] * pa[0] + rn[0] * pa[2]) + f[6] * rn[0] + f[7] * rn[1] + f[8] * rn[2]
        c = f[0] * pa[0] ** 2 + f[1] * pa[1] ** 2 + f[2] * pa[2] ** 2 + f[3] * pa[0] * pa[1] + f[4] * pa[1] * pa[
            2] + f[5] * pa[2] * pa[0] + f[6] * pa[0] + f[7] * pa[1] + f[8] * pa[2] + f[9]
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
                if np.sqrt(sum(np.square(np.array(p2) - np.array(self.points[0])))) <= \
                        np.sqrt(sum(np.square(np.array(p1) - np.array(self.points[0])))):
                    point = p2
            # If we have 2 points to choose from, choose the one that makes the angle closer to 180
            else:
                point = p1
                if calc_angle(self.points[-1], self.points[-2], p2) >= calc_angle(self.points[-1], self.points[-2], p1):
                    point = p2
            # Return the point we choose
            return point

    # Build edge function. Find points along the edge from its first vertex to its second. Has at least 10 points.
    def build(self, surf=None, res=None, straight=False):
        # Get the location and radius of the circle inscribed between the edge atoms
        self.get_loc()
        # Get the pvals
        self.find_pvals()
        if straight:
            self.points = [self.pv0, self.pv1]
            return
        # Reset the edges points
        self.points = []
        # Check to see if a minimum distance has been provided
        if res is None:
            # Get the network's minimum distance
            res = self.net.surf_res
        # Check to see if a surface has been provided
        if surf is None:
            # Choose a curved one to project onto. If the edge isn't straight 2 surfs are curved.
            if round(self.atoms[0].rad, 10) == round(self.atoms[1].rad, 10):
                surf = Surface(self.atoms[1:], self.net)
            else:
                surf = Surface(self.atoms[:2], self.net)

        ################################################# Fill Edge ####################################################

        # If the edge is completely straight add points in a line from pv0 to pv0 and return
        if straight or (self.atoms[0].rad == self.atoms[1].rad and self.atoms[1].rad == self.atoms[2].rad):
            # Get the vector between the two vectors and the number of point in the edge
            r = self.pv1 - self.pv0
            num_points = 5
            # Add the points
            for i in range(num_points + 1):
                self.points.append(self.pv0 + r * (i / num_points))
            return

        # Calculate the points

        # Find the point in between the two vertex points
        r01 = self.pv1 - self.pv0  # Vector between vertices
        r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
        rn01 = r01 / r_mag  # Normal to the vector between the vertices
        # Find the number of points
        n = max(int(r_mag / res), 4)
        # Calculate the angle between the vertices and the reference point
        theta = calc_angle(self.pa, self.pv0, self.pv1)
        # Add the first vertex to the list of points
        self.points = [self.pv0.tolist()]
        # Find the edges points. Don't count the vertex
        for i in range(n + 1):
            if i == 0:
                A = 0.01 * theta / n
            elif i == 1:
                A = 0.99 * theta / n
            else:
                A = theta / n
            # Set pb to the previous point
            pb = self.points[-1]
            # Get the distance between pb and pa for c
            c = np.sqrt(sum(np.square(np.array(pb) - np.array(self.pa))))
            # Get the angle between pb, pa and pb + rno1
            B = calc_angle(pb, pb + rn01, self.pa)
            # Get the last angle
            C = np.pi - B - A
            # Get the distance to our projection point or 'a' on our triangle
            a = np.sin(A) * c / np.sin(C)
            # Use that distance to project rn01 from pb to find our projection point or pc
            pc = pb + a * rn01
            # Get the vector from pa to pc
            rac = np.array(pc) - np.array(self.pa)
            rnac = rac / np.linalg.norm(rac)
            # Project the vector onto the surface
            surf_point = self.project(rnac, self.pa, surf)
            if surf_point is None:
                break
            self.points.append(surf_point)
        # Add the end point
        self.points.append(self.pv1)
