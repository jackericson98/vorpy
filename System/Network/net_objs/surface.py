class Surface:
    """Surface object. Holds the mesh data. Used to analyze interfaces between atoms."""
    def __init__(self, atoms=None, net=None, edges=None, verts=None, doublet=False, points=None, tris=None,
                 perimeter=None, normal=None, sa=0, curvature=None, function=None, file=None, resolution=None,
                 center=None, tri_colors=None, color_scheme=None, color_map=None, ndx=None):

        # Network objects
        self.net = net                  # Network         : Network of the System
        self.atoms = atoms              # Atoms           : List Atom objects for the surface
        self.verts = verts              # Vertices        : List of Vertex objects for surface
        self.edges = edges              # Edges           : List of Edge objects for the surface

        # Main descriptor attributes
        self.ndx = ndx                  # Index           : Indices of the atoms of the surface
        self.file = file                # File            : File address for the reference file holding points and tris
        self.func = function            # Function        : Holds the coefficients of the function describing the surf
        self.res = resolution           # Resolution      : The resolution with which to build the surface

        # Points and tris
        self.points = points            # Points          : The points that make up the surface
        self.flat_points = []           # Flat points     : Points projected into 2d based off of the surface normal
        self.perimeter = perimeter      # Perimeter       : The points around the edges of the surface (IN ORDER)
        self.pflat_points = []          # Flat perimeter  : Flattened points around the perimeter
        self.tris = tris                # Triangles       : A list of connections between the points
        self.filter_hard = False        # Filter hard     : Whether to filter the triangles extra extra extra hard UWU

        # Coloring values
        self.tri_dists = None           # Distances       : List corresponding to the distance from the center
        self.tri_ins_out = None         # Inside Outside  : List corresponding to inside or outside the atoms
        self.tri_curvs = None           # Curvature       : List corresponding to curvature value for each point
        self.tri_colors = tri_colors    # Tri colors      : Holds the color mapped color for each triangle
        self.scheme = color_scheme      # Color Scheme    : Holds the method by which the color map is mapped
        self.color_map = color_map      # Color Map       : Holds the map applied to the triangles on the surface

        # Calculation attributes
        self.sa = sa                    # Surface Area    : The surface area of the
        self.curv = curvature           # Curvature       : The curvature of the surface between the
        self.vols = [0, 0]              # Volumes         : The volume contributions from each of the atom
        self.norm = normal              # Surface Normal  : Normal to the center of the surface
        self.loc = center               # Location        : Center point of the hyperboloid the surface is made from
        self.com = None                 # Center of mass  : The point toward which all building paths travel
        self.doublet = doublet          # Doublet         : Indicates whether a surface is a part of a doublet or not
        self.flat = False               # Is Flat?        : Whether the surface is flat or not

        # Make sure that a0 is the atom with the smaller radius
        if self.atoms is not None:
            if self.atoms[0].rad > self.atoms[1].rad:
                self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]

        # Create and sort the atom's indices
        if atoms is not None:
            self.ndx = [atom.num for atom in atoms]
            self.ndx.sort()

