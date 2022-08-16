from System.Network.edge import Edge
from System.Network.surface import Surface
from System.Network.vertex import Vertex
from System.calcs import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, atoms):
        self.atoms = atoms  # List of Atom type objects
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.rad = 50  # Ballpark range for radius needed for the entire network.
        self.vta = False

    # Find v0 function. Finds the first vertex in the network
    def find_v0(self):
        # Find the center of mass of the atoms
        com = calc_atoms_com(self.atoms)
        # First choose an appropriate initial atom based of com proximity
        min_dist = np.inf
        a0 = None
        # Go through each atom determining if it is closer to the com
        for atom in self.atoms:
            # Set the new com distance
            com_dist = calc_dist(atom.loc, com)
            # If is less than the current closest atom's distance to the center of mass update the variables
            if com_dist < min_dist:
                min_dist = com_dist
                a0 = atom
        # Find the set of atoms with the minimum distance between surfaces
        min_dist = np.inf
        a1 = None
        # Go through each atom determining the atom with the minimum distance between it and a0's surfaces
        for atom in self.atoms:
            # Skip a0
            if atom == a0:
                continue
            # Set the new atom distances
            a_dist = calc_dist(a0.loc, atom.loc) - (a0.rad + atom.rad)
            # If the new atom distance is less than the previous minimum distance update the variables
            if a_dist < min_dist:
                min_dist = a_dist
                a1 = atom
        # Find the set of atoms with the minimum inscribed circle
        min_rad = np.inf
        a2 = None
        # Go through each other atom to determine the smallest circle that can be made with our 2 atoms and a third
        for atom in self.atoms:
            # Skip a0, a1
            if atom == a0 or atom == a1:
                continue
            # Calculate the circle made with the 3 atoms
            circ = calc_circ([a0, a1, atom])
            # If the radius of the inscribed circle is smaller than the previous smallest found circle's radius replace
            if circ and abs(circ[1]) < min_rad:
                min_rad = abs(circ[1])
                a2 = atom
        # Find the set of atoms with the minimum inscribed sphere
        min_rad = np.inf
        myVert = None
        # Go through each other atom to determine the smallest possible inscribed sphere
        for atom in self.atoms:
            # Skip a0, a1, a2
            if atom == a0 or atom == a1 or atom == a2:
                continue
            # Get the vertex made from the atoms
            vert = Vertex(atoms=[a0, a1, a2] + [atom])
            # If the radius of the inscribed
            if vert.loc and vert.rad < min_rad:
                min_rad = vert.rad
                myVert = vert
        # Return the vertex
        return myVert

    # Find site function. Takes in an edge and finds the only other vertex that does not overlap with other atoms
    def find_site(self, edge_atoms, vn_1):
        # Instantiate the vertex
        myVert = None
        # Loop through the atoms to see if they create a vertex that doesn't overlap with any other atoms
        for atom in self.atoms:
            # This filters out any of the atoms in the edge or the remaining atom from the previous vertex
            if {atom}.issubset(vn_1.atoms):
                continue
            # Calculate the vertex with atom
            vert = Vertex(atoms=edge_atoms + [atom])
            if vert.loc is None:
                continue
            # Check if the vertex overlaps with any of the networks atoms
            overlap = False
            for a_test in self.atoms:
                if {a_test}.issubset(edge_atoms + [atom]):
                    continue
                if round(calc_dist(a_test.loc, vert.loc) - (a_test.rad + vert.rad), 7) < 0:
                    overlap = True
                    break
            if not overlap:
                myVert = vert
                break
        return myVert

    # Find network function. Keeps searching the network until all verts are found
    def find_vertices(self):
        # Find the first vertex in the System
        v0 = self.find_v0()
        # Add v0 to the System
        self.verts.append(v0)
        # Set up the vertex stack
        vert_stack = [v0]
        # While the verts stack is not empty
        while vert_stack:
            # Running print statement giving an estimate for percentage of the network that has been created
            tot_verts = max(len(self.verts) + int(3 * len(vert_stack) / 4), 6 * len(self.atoms))
            percentage = int(len(self.verts) / tot_verts * 100)
            pertentage = percentage // 10
            print("\rBuilding Network: ", '#' * pertentage + ' ' * (10 - pertentage), percentage, "%", end='')
            # Get the vertex from the top of the stack
            vert = vert_stack.pop()
            # Set up the edge stack
            e_stack = [[[vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]], vert] for i in range(4)]
            # While the edge stack is not empty
            while e_stack:
                # Get the edge from the top of the stack
                edge = e_stack.pop()
                # Find the next site in the network
                myVert = self.find_site(edge[0], edge[1])
                # If the vertex is none continue
                if myVert is None:
                    continue
                # If the vertex exists in the network add the vertex to the edge and move on to the next edge in the stack
                found_vert = check_vert(set(myVert.atoms), self.verts)
                if not found_vert:
                    vert_stack.append(myVert)
                    self.verts.append(myVert)
        print("\rBuilding Network:  ########## 100 %")
        print("\rNetwork Built")

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
                    my_edge = Edge(list(atoms), verts, calc_points=not self.vta)
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

    # Build network function. Takes in a system and returns a fully connected network
    def build(self, min_dist=None, surfs=True):
        # Find the vertices of the system if it is not a voronota system
        if not self.vta and not self.verts:
            self.find_vertices()
        # Connect the network of vertices
        self.connect()
        num_surfs = len(self.surfs)
        if surfs:
            # Make each surface
            for i in range(num_surfs):
                # Calculate and print the running percentage for mesh calculations
                percentage = int((i + 1) / num_surfs * 100)
                pertentage = percentage // 10
                print("\rBuilding Surfaces:", '#' * pertentage + ' ' * (10 - pertentage), percentage,  "%", end='')
                # If the network is a voronota network, use build_vta method
                if self.vta:
                    self.surfs[i].build_vta()
                # Otherwise, proceed with the regular build method
                else:
                    self.surfs[i].build(simps=True, min_dist=min_dist)
            print("\r")
            print("\rSurfaces Built")

    # Analyze system function. Finds the surfaces and volumes of the system
    def analyze(self):
        # Go through each surface in the system and find the simplices and the surface area
        for surf in self.surfs:
            # Get the surfaces simplices
            surf.simps = surf.find_simps()
            # Get the surface area of the surface
            surf.sa = calc_sa(surf)

        # Go through each atom in the system and find the volume
        for atom in self.atoms:
            atom.cell_vol = calc_vol(atom)

