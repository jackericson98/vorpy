from System.Network.edge import Edge
from System.Network.surface import Surface
from System.Network.find_vertices import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""

    def __init__(self, sys, atoms, min_dist=0.1, box_size=1.5):
        self.sub_boxes = None
        self.sub_box_size = []
        self.box = None
        self.box_size = box_size
        self.atoms_range = []
        self.sys = sys
        self.min_dist = min_dist
        self.atoms = atoms  # List of Atom type objects
        self.verts = []  # List of Vertex type objects
        self.surfs = []  # List of Surface type objects
        self.edges = []  # List of Edge type objects
        self.vert_ndxs = []
        self.vta = False

    # Calculate box function. Takes in a System and returns the dimensions of a box x times the size of the atoms
    def calc_box(self):
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = np.array([np.inf, np.inf, np.inf])
        max_vert = np.array([-np.inf, -np.inf, -np.inf])
        # Check each atom in the System
        for atom in self.atoms:
            # Go through x, y, z
            for i in range(3):
                # If we find that the x, y, z value is less replace the value in the mins list
                if atom.loc[i] <= min_vert[i]:
                    min_vert[i] = atom.loc[i]
                # If we find that the x, y, z value is less replace the value in the mins list
                elif atom.loc[i] >= max_vert[i]:
                    max_vert[i] = atom.loc[i]
        # Get the vector between the minimum and maximum vertices for the defining box
        r_box = max_vert - min_vert
        # If the atoms are in the same plane
        for i in range(3):
            if r_box[i] == 0 or abs(r_box[i]) == np.inf:
                r_box[i] = 40 * self.atoms[0].rad
        self.atoms_range = [min_vert, max_vert]
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert - r_box * self.box_size, min_vert + r_box * self.box_size
        # Return the list of array turned list vertices
        self.box = [min_vert.tolist(), max_vert.tolist()]

    # Sort atoms method. Puts the atoms in the network in their respective grid sections
    def sort_atoms(self, num_boxes=8000):
        # First get the box for the
        self.calc_box()
        # Number of cells per row/column/aisle
        n = int(np.cbrt(num_boxes))
        # Divide the box into sub_boxes
        self.sub_boxes = [[[[] for _ in range(n)] for _ in range(n)] for _ in range(n)]
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
    def get_atoms(self, cells, reach, exclusive=False):
        # Get the min and max of the cells
        ndx_min = [np.inf, np.inf, np.inf]
        ndx_max = [-np.inf, -np.inf, -np.inf]
        # Go through the cells and set the minimum and maximum indexes for xyz for a rectangle containing the atoms
        for cell in cells:
            # Check each xyz index to see if they are larger or smaller than the max or min
            for i in range(3):
                if cell[i] < ndx_min[i]:
                    ndx_min[i] = cell[i]
                if cell[i] > ndx_max[i]:
                    ndx_max[i] = cell[i]
        # Set the initial search parameters to the given cells
        xs, ys, zs = [x for x in range(-reach + ndx_min[0] + 1, reach + ndx_max[0])], \
                     [y for y in range(-reach + ndx_min[1] + 1, reach + ndx_max[1])], \
                     [z for z in range(-reach + ndx_min[2] + 1, reach + ndx_max[2])]
        atoms = []
        # Go through each box in the range given and add the atoms
        for i in xs:
            for j in ys:
                for k in zs:
                    # If the exclusive parameter was set we only want the outer shell, skip none of the indices are max
                    if exclusive and abs(i) != reach and abs(j) != reach and abs(k) != reach:
                        continue
                    # Easy way around hitting the edge of the box
                    try:
                        atoms += self.sub_boxes[i][j][k]
                    except IndexError:
                        continue
        return atoms

    # Filter vertices function. Filters out any repeat vertices
    def filter_verts(self):
        # Re-sort the atoms
        self.sort_atoms()
        # Set up a list of vertex ndxs and vertices
        vert_ndxs = []
        verts = []
        # Go through the vertices
        for i in range(len(self.verts)):
            # Sort the indices of the vertices
            self.verts[i].ndx.sort()
            vert_inside = all([self.box[0][j] < self.verts[i].loc[j] < self.box[1][j] for j in range(3)])
            if self.verts[i].ndx not in vert_ndxs:
                vert_ndxs.append(self.verts[i].ndx)
                verts.append(self.verts[i])
        # Set the networks vertices
        self.verts = verts

    # Connect network method.
    def connect(self):
        # Filter the vertices
        self.filter_verts()
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
                my_edge = Edge(list(atoms), verts, self, calc_points=not self.vta)
                # Add the edge to the System
                self.edges.append(my_edge)

        # Create the surfaces
        self.surfs = []
        for edge1 in self.edges:
            if edge1.touches_box:
                continue
            # Go through the edge's atoms combinations
            for i in range(3):
                atoms = {edge1.atoms[i], edge1.atoms[(i + 1) % 3]}
                # If the surface has been found before continue
                if check_surf(atoms, self.surfs):
                    continue
                # Put together a list of edges that have our atoms
                edges = []
                for edge2 in self.edges:
                    if atoms.issubset(edge2.atoms):
                        edges.append(edge2)
                # Put together a list of verts that have our atoms
                verts = []
                for vert2 in self.verts:
                    if atoms.issubset(vert2.atoms):
                        verts.append(vert2)
                # In order to be a true surface the number of edges need to be equal to the number of verts
                if len(verts) == len(edges):
                    my_surf = Surface(list(atoms), verts=verts, net=self, edges=edges)
                    self.surfs.append(my_surf)
                else:
                    pass

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

    def find_verts(self):
        # Put the atoms in their place
        self.sort_atoms()
        # Find the vertices of the system if it is not a voronota system or we haven't indicated not to find them
        if not self.vta:
            # Go through each atom in the system
            for i in range(len(self.atoms)):
                # Update the running print statement
                tot_verts = len(self.verts) + (len(self.atoms) - i)
                percentage = int(len(self.verts) / tot_verts * 10000) / 100
                print("\rBuilding Network:  ", '#' * (int(percentage) // 10) + ' ' * (10 - (int(percentage) // 10)),
                      percentage, "%", end='')
                # If the atom has no vertices run the vertex finder on it
                if len(self.atoms[i].verts) == 0:
                    v0 = find_v0(self, self.atoms[i])
                    if v0 is not None:
                        find_vertices(self, v0, i=i)
        print("\rBuilding Network:   ########## 100 %")
        print("\rNetwork Built")

    # Build network function. Takes in a system and returns a fully connected network
    def build(self, get_verts=True, get_surfs=True):
        if get_verts:
            self.find_verts()
        # Connect the network of vertices
        self.connect()
        num_surfs = len(self.surfs)
        if get_surfs:
            # Make each surface
            for i in range(num_surfs):
                # Calculate and print the running percentage for mesh calculations
                percentage = int((i + 1) / num_surfs * 100)
                print("\rBuilding Surfaces: ",
                      '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
                # If the network is a voronota network, use build_vta method
                if self.vta:
                    self.surfs[i].build_vta()
                # Otherwise, proceed with the regular build method
                else:
                    self.surfs[i].build(simps=True, min_dist=self.min_dist)
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
            print("\rAnalyzing System:  ",
                  '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
            # Get the surface area of the surface
            self.surfs[i].sa = calc_sa(self.surfs[i])

        # Go through each atom in the system and find the volume
        for j in range(len(self.atoms)):
            percentage = int((i + j + 1) / tot_num * 100)
            print("\rAnalyzing System:  ",
                  '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
            self.atoms[j].vol = calc_vol(self.atoms[j])

        print("\rAnalyzing System:   ########## 100 %")
        print("\rSystem Analyzed")
