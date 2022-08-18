from System.Network.edge import *
from System.Network.surface import Surface
from System.Network.vertex import Vertex


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, atoms):
        self.atoms = atoms  # List of Atoms
        self.verts = []  # List of Vertices
        self.surfs = []  # List of Surfaces
        self.edges = []  # List of Edges
        self.vta = False  # Indicator for Voronota network or not

    # Find v0 method. Finds the first vertex in the network
    def find_v0(self):
        # Find the center of mass of the atoms
        com = calc_atoms_com(self.atoms)
        # First choose an appropriate initial atom based of com proximity
        min_dist = np.inf
        a0 = None
        # Go through each atom determining if they are closer to the com
        for atom in self.atoms:
            # Set the new com distance
            com_dist = calc_dist(atom.loc, com)
            # If is less than the current closest atom's distance to the center of mass update the variables
            if com_dist < min_dist:
                min_dist = com_dist
                a0 = atom
        # Find the closest atom to a0
        min_dist = np.inf
        a1 = None
        # Go through each atom skipping a0
        for atom in self.atoms:
            if atom == a0:
                continue
            # Find the distance between the surfaces of the atoms
            a_dist = calc_dist(a0.loc, atom.loc) - (a0.rad + atom.rad)
            # If the new atom distance is less than the previous minimum distance update the variables
            if a_dist < min_dist:
                min_dist = a_dist
                a1 = atom
        # Find a2 : a0, a1, a2 have the smallest possible inscribed circle
        min_rad = np.inf
        a2 = None
        # Go through the atoms creating a circle from a0, a1 and third atom
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
            # If the radius of the inscribed sphere is smaller than the previously recorded smallest replace the vars
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
            print("\rFinding Vertices:   ", '#' * pertentage + ' ' * (10 - pertentage), percentage, "%", end='')
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
        print("\rFinding Vertices:   ########## 100 %")
        print("\rVertices Found")

    # Connect network method.
    def connect(self):
        # Create the edges
        for vert1 in self.verts:
            # Check every combination of vert atoms as an edge
            for i in range(4):
                # Grab the atoms
                atoms = {vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]}
                # If the edge has been found before, continue
                if check_edge(atoms, self.edges):
                    continue
                verts = []
                # Find the possible verts (the original vert and the new vert)
                for vert2 in self.verts:
                    if atoms.issubset(vert2.atoms):
                        verts.append(vert2)
                # If the number of valid vertices for the edge is 1
                if len(verts) == 1:
                    continue
                # Create the edge
                my_edge = Edge(list(atoms), verts, calc_points=not self.vta)
                # Add the edge to the System
                self.edges.append(my_edge)

        # Create the surfaces
        for edge in self.edges:
            # Go through the edge's atoms combinations
            for i in range(3):
                atoms = {edge.atoms[i], edge.atoms[(i+1) % 3]}
                # If the surface has been found before continue
                if check_surf(atoms, self.surfs):
                    continue
                # Put together a list of edges that have our atoms
                edges = []
                for edge in self.edges:
                    if atoms.issubset(edge.atoms):
                        edges.append(edge)
                # Put together a list of verts that have our atoms
                verts = []
                for vert2 in self.verts:
                    if atoms.issubset(vert2.atoms):
                        verts.append(vert2)
                # In order to be a true surface the number of edges need to be equal to the number of verts
                if len(verts) == len(edges):
                    my_surf = Surface(list(atoms), verts=verts, edges=edges)
                    self.surfs.append(my_surf)

        # Add the vertices, edges and surfs to the atoms
        for atom in self.atoms:
            # Reset the atom's vert list
            atom.verts = []
            # Go through the verts in the network
            for vert in self.verts:
                # If the atom is in the vertices atoms add the vertex to the atom's list of vertices
                if {atom}.issubset(vert.atoms):
                    atom.verts.append(vert)
            # Reset the atom's edge list
            atom.edges = []
            # Go through the edges in the network
            for edge in self.edges:
                # If the atom is in the edge's list of atoms add the edge to the atoms list of edges
                if {atom}.issubset(edge.atoms):
                    atom.edges.append(edge)
            # Reset the atom's surf list
            atom.surfs = []
            # Go through the surfs in the network
            for surf in self.surfs:
                # If the atom is in the surfs list of atoms add the surf to the atoms list of surfs
                if {atom}.issubset(surf.atoms):
                    atom.surfs.append(surf)


        # Add the edges and surfs to the vertices
        for vert in self.verts:
            # Reset the vertexes edge list
            vert.edges = []
            # Go through the edges in the network
            for edge in self.edges:
                # If the edges atoms are in the vertices atoms add it to the vertex
                if set(edge.atoms).issubset(vert.atoms):
                    vert.edges.append(edge)
            # Reset the vertexes surf list
            vert.surfs = []
            # Go through the surfaces in the network
            for surf in self.surfs:
                # If the surfaces atoms are in the vertexes atoms add it to the vertex
                if set(surf.atoms).issubset(vert.atoms):
                    vert.surfs.append(surf)

        # Add the surfs to the edges
        for edge in self.edges:
            # Reset the edges surf list
            edge.surfs = []
            # Go through the surfaces in the network
            for surf in self.surfs:
                # If the surfaces atoms are in the edges atoms add it to the edge
                if set(surf.atoms).issubset(edge.atoms):
                    edge.surfs.append(surf)

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
                print("\rBuilding Surfaces: ", '#' * pertentage + ' ' * (10 - pertentage), percentage,  "%", end='')
                # If the network is a voronota network, use build_vta method
                if self.vta:
                    self.surfs[i].build_vta()
                # Otherwise, proceed with the regular build method
                else:
                    self.surfs[i].build(simps=True, min_dist=min_dist)
            print("\rBuilding Surfaces:  ########## 100 %")
            print("\rSurfaces Built")

    # Analyze system function. Finds the surfaces and volumes of the system
    def analyze(self):
        # Get the percentage total number
        tot_num = len(self.surfs) + len(self.atoms)
        # Go through each surface in the system and find the simplices and the surface area
        i = 0
        for i in range(len(self.surfs)):
            percentage = int((i + 1) / tot_num * 100)
            pertentage = percentage // 10
            print("\rAnalyzing System:  ", '#' * pertentage + ' ' * (10 - pertentage), percentage, "%", end='')
            # Get the surfaces simplices
            self.surfs[i].simps = self.surfs[i].find_simps()
            # Get the surface area of the surface
            self.surfs[i].sa = calc_sa(self.surfs[i])

        # Go through each atom in the system and find the volume
        for j in range(len(self.atoms)):
            percentage = int((i + j + 1) / tot_num * 100)
            pertentage = percentage // 10
            print("\rAnalyzing System:  ", '#' * pertentage + ' ' * (10 - pertentage), percentage, "%", end='')
            self.atoms[j].vol = calc_vol(self.atoms[j])

        print("\rAnalyzing System:   ########## 100 %")
        print("\rSystem Analyzed")
