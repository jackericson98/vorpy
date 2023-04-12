class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms=None, net=None, verts=None, surfs=None, points=None, center=None, rad=None, normal=None,
                 pv0=None, pv1=None, ndx=None, straight=False):

        # If no network was given have a catch
        if net is not None and net.atoms is not None and atoms is not None:
            ndx = [atom.num for atom in atoms]
            ndx.sort()
        self.ndx = ndx                   # Index         :   Indices of the atoms of the surface
        self.net = net                   # Network       :   Network of the System
        self.atoms = atoms               # Atoms         :   List of Atom type objects for the edge
        self.verts = verts               # Vertices      :   List of Vertex type objects
        self.surfs = surfs               # Surfaces      :   List of 2 surfaces attached to the edge

        self.loc = center                # Location      :   Location of the center of the 3 atoms that make up the edge
        self.rad = rad                   # Radius        :   Radius of the inscribed circle of the three atoms
        self.norm = normal
        self.points = points             # Points        :   List of points along the
        self.pv0 = pv0                   # Vertex pt 0   :   The points on the ends of the edges
        self.pv1 = pv1                   # Vertex pt 1   :   The points on the ends of the edges
        self.pa = None                   # Projection pt :   The projection point from which the edge is built
        self.loc2 = None                 # Loc2          :   Allows edges to be checked like vertices
        self.draw_points = None
        self.draw_tris = None
        self.straight = straight         # Straight edge :   Straight edge or not
