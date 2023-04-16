class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms=None, net=None, verts=None, surfs=None, points=None, ndx=None, ref=None):

        # If no network was given have a catch
        if net is not None and net.atoms is not None and atoms is not None:
            ndx = [atom.num for atom in atoms]
            ndx.sort()
        self.ndx = ndx                   # Index         :   Indices of the atoms of the surface

        self.net = net                   # Network       :   Network of the System
        self.atoms = atoms               # Atoms         :   List of Atom type objects for the edge
        self.verts = verts               # Vertices      :   List of Vertex type objects
        self.surfs = surfs               # Surfaces      :   List of 2 surfaces attached to the edge

        self.vals = {}                   # Edge values   :   Dictionary of values (loc, rad, length)
        self.ref = ref                    # Reference     :   Holds the surface and the corresponding indices of points
        self.points = points             # Points        :   List of points along the edge

        self.draw_points = None
        self.draw_tris = None
        self.length = None

        if self.ref is None:
            self.ref = {}
