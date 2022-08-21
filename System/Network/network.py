from System.Network.edge import *
from System.Network.surface import Surface
from System.Network.vertex import Vertex
from System.calcs import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, sys, atoms, box_size=1.5):
        self.sub_boxes = None
        self.sub_box_size = []
        self.sys = sys
        self.atoms = atoms  # List of Atom type objects
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.rad = 50  # Ballpark range for radius needed for the entire network.
        self.vta = False
        self.box = None
        self.sort_atoms(box_size)

    # Calculate box function. Takes in a System and returns the dimensions of a box x times the size of the atoms
    def calc_box(self, x):
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = np.array([np.inf, np.inf, np.inf])
        max_vert = np.array([-np.inf, -np.inf, -np.inf])
        # Check each atom in the System
        for atom in self.atoms:
            # Go through x, y, z
            for i in range(3):
                # If we find that the x, y, z value is less replace the value in the mins list
                if atom.loc[i] < min_vert[i]:
                    min_vert[i] = atom.loc[i]
                # If we find that the x, y, z value is less replace the value in the mins list
                elif atom.loc[i] > max_vert[i]:
                    max_vert[i] = atom.loc[i]
        # Get the vector between the minimum and maximum vertices for the defining box
        r_box = max_vert - min_vert
        # If the atoms are in the same plane
        for i in range(3):
            if r_box[i] == 0:
                r_box[i] = 4 * self.atoms[0].rad
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = min_vert - r_box * x, max_vert + r_box * x
        # Return the list of array turned list vertices
        return [min_vert.tolist(), max_vert.tolist()]

    # Sort atoms method. Puts the atoms in the network in their respective grid sections
    def sort_atoms(self, box_size=2, num_boxes=1000):
        # First get the box for the
        self.box = self.calc_box(box_size)
        # Number of cells per row/column/aisle
        n = int(np.cbrt(num_boxes))
        # Divide the box into sub_boxes
        self.sub_boxes = [[[[] for i in range(n)] for j in range(n)] for k in range(n)]
        # Get the cell size
        self.sub_box_size = [(self.box[1][0] - self.box[0][0]) / n, (self.box[1][1] - self.box[0][1]) / n,
                     (self.box[1][2] - self.box[0][2]) / n]
        # Sort the atoms
        for atom in self.atoms:
            # Find the box they belong to
            ai = int((atom.loc[0] - self.box[0][0]) / self.sub_box_size[0])
            aj = int((atom.loc[1] - self.box[0][1]) / self.sub_box_size[1])
            ak = int((atom.loc[2] - self.box[0][2]) / self.sub_box_size[2])
            # Add the atom to the box
            self.sub_boxes[ai][aj][ak].append(atom)
            # Add the box to the atom
            atom.box = [ai, aj, ak]

    # Get atoms method. Takes in the cells and the number of additional cells to search and returns an atom list
    def get_atoms(self, cells, rnge, exclusive=False):
        # Get the min and max of the cells
        ndx_min = [np.inf, np.inf, np.inf]
        ndx_max = [-np.inf, -np.inf, -np.inf]
        # Go through the cells and set the minimum and
        for cell in cells:
            for i in range(3):
                if cell[i] < ndx_min[i]:
                    ndx_min[i] = cell[i]
                if cell[i] > ndx_max[i]:
                    ndx_max[i] = cell[i]

        # Set the initial search parameters to the given cells
        xs, ys, zs = [x for x in range(-rnge + ndx_min[0] + 1, rnge + ndx_max[0])], \
                     [y for y in range(-rnge + ndx_min[1] + 1, rnge + ndx_max[1])], \
                     [z for z in range(-rnge + ndx_min[2] + 1, rnge + ndx_max[2])]
        atoms = []
        # First go through the atoms in a0's box
        if exclusive:
            for i in xs:
                for j in ys:
                    for k in zs:
                        if abs(i) != rnge and abs(j) != rnge and abs(k) != rnge:
                            continue
                        try:
                            atoms += self.sub_boxes[i][j][k]
                        except IndexError:
                            continue
        else:
            for i in xs:
                for j in ys:
                    for k in zs:
                        try:
                            atoms += self.sub_boxes[i][j][k]
                        except IndexError:
                            continue
        return atoms

    # Find v0 function. Finds the first vertex in the network
    def find_v0(self):
        # Find the middle sub_box of the set of boxes and
        mid = len(self.sub_boxes) // 2
        atoms = []
        inc = 0
        while not atoms:
            atoms = self.get_atoms([[mid, mid, mid]], inc)
            inc += 1
        a0 = atoms[0]

        # Find the set of atoms with the minimum distance between surfaces
        min_dist = np.inf
        a1 = None
        inc = 0
        while not a1:
            atoms = self.get_atoms([a0.box], inc)
            # Go through each atom determining the atom with the minimum distance between it and a0's surfaces
            for atom in atoms:
                # Skip a0
                if atom == a0:
                    continue
                # Set the new atom distances
                a_dist = calc_dist(a0.loc, atom.loc) - (a0.rad + atom.rad)
                # If the new atom distance is less than the previous minimum distance update the variables
                if a_dist < min_dist:
                    min_dist = a_dist
                    a1 = atom
            inc += 1
        # Find the set of atoms with the minimum inscribed circle
        min_rad = np.inf
        a2 = None
        inc = 0
        while not a2:
            atoms = self.get_atoms([a0.box, a1.box], inc + 1)
            # Go through each other atom to determine the smallest circle that can be made with our 2 atoms and a third
            for atom in atoms:
                # Skip a0, a1
                if atom == a0 or atom == a1:
                    continue
                # Calculate the circle made with the 3 atoms
                circ = calc_circ([a0, a1, atom])
                # If the radius of the inscribed circle is smaller than the previous smallest found circle's radius replace
                if circ and abs(circ[1]) < min_rad:
                    min_rad = abs(circ[1])
                    a2 = atom
            inc += 1
        # Find the set of atoms with the minimum inscribed sphere
        min_rad = np.inf
        myVert = None
        # Go through each other atom to determine the smallest possible inscribed sphere
        for atom in self.atoms:
            # Skip a0, a1, a2
            if atom == a0 or atom == a1 or atom == a2:
                continue
            # Get the vertex made from the atoms
            vert = Vertex(atoms=[a0, a1, a2] + [atom], net=self)
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
        inc = 0
        min_vert = np.inf
        # Loop through the atoms to see if they create a vertex that doesn't overlap with other atoms or whole system
        while inc <= len(self.sub_boxes) + 1:
            # If no vert has been found yet, keep expanding
            if myVert is None or myVert.loc is None:
                # Keep adding boxes around the atoms to check
                atoms = self.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], inc)
                # Go through the atoms in the surrounding boxes and find the smallest vertex that can be created
                for atom in atoms:
                    # This filters out any of the atoms in the edge or the remaining atom from the previous vertex
                    if {atom}.issubset(vn_1.atoms):
                        continue
                    # Calculate the vertex with atom and pass if the vertex location is None
                    vert = Vertex(atoms=edge_atoms + [atom], net=self)
                    # I need to fix Vertex to get to a point where I dont need this
                    if vert is None or vert.loc is None:
                        continue
                    # Sniff out the smallest vertex that can be made in the box and store it
                    if vert.rad < min_vert:
                        # Replace the variables
                        min_vert = vert.rad
                        myVert = vert
            # If no vertex can be made restart the search with a larger set of boxes
            if myVert is None or myVert.loc is None:
                inc += 1
                continue
            # Find the box that the vertex would be in
            vi = int((myVert.loc[0] - self.box[0][0]) / self.sub_box_size[0])
            vj = int((myVert.loc[1] - self.box[0][1]) / self.sub_box_size[1])
            vk = int((myVert.loc[2] - self.box[0][2]) / self.sub_box_size[2])
            # Any atom that can overlap with this vertex is within the 'vertex's radius over the smallest box length'
            # plus the 'maximum radius of an atom over the smallest box length' # of boxes away from the vert.loc box
            atoms = self.get_atoms([[vi, vj, vk]],
                                   int(myVert.rad/min(self.sub_box_size)) + int(5/min(self.sub_box_size)) + 2)
            # Set up an overlap variable
            overlap = False
            # Test the atoms in the new atom list to see if they overlap with the vertex
            for atom in atoms:
                # If the atom is one of the vert atoms move on
                if {atom}.issubset(myVert.atoms):
                    continue
                # If the distance between the vertex and the atom is less than their radii create a new vertex and reset
                if calc_dist(atom.loc, myVert.loc) - (atom.rad + myVert.rad) < 0:
                    myVert = Vertex(edge_atoms + [atom], net=self)
                    overlap = True
                    break
            # Check to see if the vertex had no overlaps --> We found the vertex!
            if not overlap:
                return myVert
            # If we have made it this far increment the counter and keep expanding
            inc += 1

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
            print("\rBuilding Network:  ", '#' * pertentage + ' ' * (10 - pertentage), percentage, "%", end='')
            # Get the vertex from the top of the stack
            vert = vert_stack.pop()
            # Set up the edge stack
            e_stack = [[[vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]], vert] for i in range(4)]
            # While the edge stack is not empty
            while e_stack:
                # Get the edge from the top of the stack
                edge, vert = e_stack.pop()
                # Find the next site in the network
                myVert = self.find_site(edge, vert)
                # If the vertex is none continue
                if myVert is None:
                    continue
                # If the vertex exists in the network add the vertex to the edge and move on to the next edge in the stack
                found_vert = check_vert(set(myVert.atoms), self.verts)
                if not found_vert:
                    vert_stack.append(myVert)
                    self.verts.append(myVert)
        print("\rBuilding Network:   ########## 100 %")
        print("\rNetwork Built")

    # Connect network method.
    def connect(self):
        # Check to see if the voronota data has been loaded
        if self.vta and not self.verts:
            print("Load Voronota data for {} to continue".format(self.sys.name))
            exit()
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
                    my_edge = Edge(list(atoms), verts, calc_points=not self.vta, net=self)
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
                    my_surf = Surface(list(t_atoms), verts=verts, edges=edges, net=self)
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
