from System.sys_funcs.calcs import *


class Vertex:
    """Vertex object. Used to build the network and calculate the surfaces"""

    def __init__(self, atoms=None, net=None, location=None, radius=None, loc2=None, rad2=None, doublet=None, ndx=None,
                 edges=None, surfaces=None):

        self.net = net                # Network       :   Network object for the vertex to refer back to
        self.atoms = atoms            # Atoms         :   List of atoms used to construct the vertex
        self.edges = edges            # Edges         :   List of Edge type objects connected to the vertex
        self.surfs = surfaces         # Surfaces      :   List of Surface type objects that the vertex is a part of

        self.ndx = ndx                # Index         :   Indices of the atoms in the vertex
        self.loc = location           # Location      :   Where the vertex is located in 3D
        self.rad = radius             # Radius        :   Radius of the vertex's tangential sphere
        self.box = None               # Box index     :   Sub box index for sorting

        self.doublet = doublet        # Doublet       :   Whether the vertex is a doublet
        self.loc2 = loc2              # Location 2    :   Location of the doublet site
        self.rad2 = rad2              # Radius 2      :   Radius of the doublet site's tangential sphere

        # If the vertex is mature enough to be calculated, create and sort its indices
        if self.net is not None and self.atoms is not None:
            self.ndx = [atom.num for atom in self.atoms]
            self.ndx.sort()

        # Set up the edge and surfaces lists
        if self.edges is None:
            self.edges = []
        if self.surfs is None:
            self.surfs = []
