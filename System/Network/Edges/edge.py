

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
