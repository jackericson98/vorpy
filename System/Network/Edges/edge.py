from System.Network.Edges.edge_calcs import *


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

    # Calc
    def calc_points(self, surf=None):
        # I want to be able to calculate a surface here, but don't want to cross pollinate imports on this level
        if surf is None:
            return
        # Grab the function's coefficients
        f = surf.func
        # Grab the vertex points
        pv0, pv1 = np.array(self.verts[0].loc), np.array(self.verts[1].loc)
        # Find the point in between the two vertex points
        r01 = pv1 - pv0
        r01_mag = np.linalg.norm(r01)
        rn01 = r01 / r01_mag
        # Get the center point of the vertices
        cp = pv0 + 0.5 * rn01 * r01_mag

