from System.Network.edge import Edge
from System.Network.surface import Surface
from System.Network.net_funcs import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, atoms):
        self.atoms = atoms  # List of Atom type objects
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.rad = 50  # Ballpark range for radius needed for the entire network.

    def connect(self):
        # Create edges and add connections between verts and edges
        # Go through each vertex and find its edges
        for vert1 in self.verts:
            # Check every combination of vert atoms as an edge
            for i in range(4):
                # Grab the atoms
                atoms = {vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]}
                verts = []
                # Find the possible verts
                for vert2 in self.verts:
                    if atoms.issubset(vert2.atoms):
                        verts.append(vert2)
                # Find which edge, if any, go nowhere
                if len(verts) == 1:
                    continue
                # Check to see if the edge has been found
                my_edge = check_edge(atoms, self.edges)
                if my_edge is None:
                    # Create the edge
                    my_edge = Edge(list(atoms), verts)
                    # Add the edge to the System
                    self.edges.append(my_edge)
                    # Add the edge to the verts
                    verts[0].edges.append(my_edge)
                    verts[1].edges.append(my_edge)

        # Create surfaces and add connections for edges and verts
        for vert1 in self.verts:
            # Go through each combination of sets atom in the vertices' atom list
            t_ndxs = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
            for ndxs in t_ndxs:
                # Grab the atoms
                t_atoms = {vert1.atoms[ndxs[0]], vert1.atoms[ndxs[1]]}
                # Check to see if we have recorded this surface before
                if check_surf(t_atoms, self.surfs):
                    continue
                # Put together a list of edges that have our atoms
                edges = []
                for edge in self.edges:
                    if t_atoms.issubset(edge.atoms):
                        edges.append(edge)
                # Put together a list of verts that have our atoms
                verts = []
                for vert2 in self.verts:
                    if t_atoms.issubset(vert2.atoms):
                        verts.append(vert2)
                # In order to be a true surface the number of edges need to be equal to the number of verts
                if len(verts) == len(edges):
                    my_surf = Surface(list(t_atoms), verts=verts, edges=edges)
                    self.surfs.append(my_surf)
                    list(t_atoms)[0].surfs.append(my_surf)
                    list(t_atoms)[1].surfs.append(my_surf)
                    list(t_atoms)[0].edges += edges
                    list(t_atoms)[1].edges += edges
                    list(t_atoms)[0].verts += verts
                    list(t_atoms)[1].verts += verts

        # Add the surfaces to the edges
        for edge in self.edges:
            edge.surfs = []
            for surf in self.surfs:
                if set(surf.atoms).issubset(edge.atoms):
                    edge.surfs.append(surf)
        # Add the surfaces to the vertices
        for vert in self.verts:
            vert.surfs = []
            for surf in self.surfs:
                if set(surf.atoms).issubset(vert.atoms):
                    vert.surfs.append(surf)

    def find_vertices(self):
        find_vertices(self)

    # Build network function. Takes in a system and returns a fully connected network
    def build_net(self):
        # Find the vertices of the system
        find_vertices(self)
        # Connect the network of vertices
        self.connect()

    def build_meshes(self, min_dist=None):
        # Set the minimum distance
        if min_dist is None:
            min_dist = 0.5
        num_surfs = len(self.surfs)
        # Make each surface
        for i in range(num_surfs):
            # Calculate and print the running percentage for mesh calculations
            percentage = int(i / num_surfs * 100) + 1
            print("\rBuilding Surfaces: ", percentage, "%", end='')
            self.surfs[i].build(min_dist=min_dist)
