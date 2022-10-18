from System.Network.find_verts import *
from System.Network.build_net import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, sys, index=0, atoms=None, verts=None, edges=None, surfs=None, groups=None,
                 min_dist=0.1, box_size=1.5, beta_val=2, sol_verts=True):
        # Network graph objects
        self.sys = sys              # System         :  Route back to outer system for system attribute access
        self.index = index          # Index          :  Holds the index of the network in the system object
        self.atoms = atoms          # Atoms          :  Atoms of the network. Should be identical to self.sys.atoms
        self.verts = verts          # Vertices       :  Vertices of the network
        self.edges = edges          # Edges          :  Edges of the network
        self.surfs = surfs          # Surfaces       :  Surfaces of the network
        self.groups = groups        # Groups         :  Groups objects for analysis of selected surfaces
        self.name = None            # Name           :  Name of the network. Used to name subnetworks recursively
        # Tools for splitting up the atoms
        self.box = None             # Box            :  Holds a max and min vertex for the retaining box
        self.sub_boxes = None       # Sub boxes      :  Holds atoms in their different relative locations in the grid
        self.vert_sub_boxes = None  # Vert sub-boxes :  Holds the vertex boxes to more easily check for vertices
        self.sub_box_size = None    # Sub box size   :  Holds the size of each sub box
        self.atoms_box = []         # Atoms box      :  Holds the min and max verts for the box containing the atoms
        self.max_atom_rad = 0       # Max atom rad   :  Holds the largest radius of the system for reference
        self.vert_ndxs = []         # Vert indices   :  Holds the indices of the atoms of the vertices in the network
        # Settings
        self.min_dist = min_dist    # Resolution     :  How small the triangles in the surfaces are
        self.beta_val = beta_val    # Beta value     :  The maximum vertex radius for the network
        self.box_size = box_size    # Box size       :  Holds the box multiplier for the system box from the atoms box
        self.parallelize = False    # Parallelize    :  Split the calculations between cores?
        self.sol_verts = sol_verts  # Sol Vertices   :  Solve the solution's vertices?
        self.curved_faces = True    # Curved Faces   :  Create curved faces for surfaces?
        self.flat_faces = False     # Flat Faces     :  Create flat faces for surfaces?
        self.verts_loaded = False   # Verts Loaded   :  Use loaded verts?
        # Run diagnostics
        self.cpu_time = None        # CPU time       :  CPU time taken to calculate the network
        self.my_time = None         # My time        :  Time taken to calculate the network

    # Calculate box function. Takes in a System and returns the dimensions of a box x times the size of the atoms
    def calc_box(self):
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = np.array([np.inf, np.inf, np.inf])
        max_vert = np.array([-np.inf, -np.inf, -np.inf])
        # Check each atom in the network
        for atom in self.atoms:
            # Go through x, y, z
            for i in range(3):
                # If we find that the x, y, z value is less replace the value in the mins list
                if atom.loc[i] <= min_vert[i]:
                    min_vert[i] = atom.loc[i]
                # If we find that the x, y, z value is less replace the value in the mins list
                if atom.loc[i] >= max_vert[i]:
                    max_vert[i] = atom.loc[i]
        # Get the vector between the minimum and maximum vertices for the defining box
        r_box = max_vert - min_vert
        # If the atoms are in the same plane
        for i in range(3):
            if r_box[i] == 0 or abs(r_box[i]) == np.inf:
                r_box[i] = 4 * self.atoms[0].rad
                min_vert[i], max_vert[i] = self.atoms[0].loc[i], self.atoms[0].loc[i]
        self.atoms_box = [min_vert, max_vert]
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert - r_box * self.box_size, min_vert + r_box * self.box_size
        # Return the list of array turned list vertices
        self.box = [min_vert.tolist(), max_vert.tolist()]

    # Sort atoms method. Puts the atoms in the network in their respective grid sections
    def sort_atoms(self, num_boxes=None):
        # Check the length of the atoms list
        if len(self.atoms) < 1:
            return
        # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
        if num_boxes is None:
            n = int(np.sqrt(len(self.atoms))) + 1
        else:
            n = int(np.cbrt(num_boxes)) + 1
        # First get the box for the atoms to be sorted into
        self.calc_box()
        # Instantiate the grid structure of lists is locations representing a grid
        self.sub_boxes = [[[[] for _ in range(n)] for _ in range(n)] for _ in range(n)]
        # Get the cell size
        self.sub_box_size = [(self.box[1][0] - self.box[0][0]) / n, (self.box[1][1] - self.box[0][1]) / n,
                             (self.box[1][2] - self.box[0][2]) / n]
        # Sort the atoms
        for atom in self.atoms:
            # Adjust the maximum radius
            if atom.rad > self.max_atom_rad:
                self.max_atom_rad = atom.rad
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
        xs, ys, zs = [x for x in range(max(0, -reach + ndx_min[0] + 1), reach + ndx_max[0])], \
                     [y for y in range(max(0, -reach + ndx_min[1] + 1), reach + ndx_max[1])], \
                     [z for z in range(max(0, -reach + ndx_min[2] + 1), reach + ndx_max[2])]
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
                        # Add the atoms
                        atoms += self.sub_boxes[i][j][k]
                        # Add a little catch to not go forever
                        if len(atoms) == len(self.atoms):
                            return atoms
                    except IndexError:
                        continue
        return atoms

    # Connect method. Connects the network using the functions in the connect_network.py file
    def connect(self):
        build(self)

    # Find subnetworks method. Divides the network up into smaller networks, finds their vertices and combines them
    def find_subnets(self):
        # Instantiate the vertices
        self.verts = []
        # This gets us to average about 60 atoms per medium box
        n = max(int(np.cbrt(len(self.atoms) // 30)), 2)
        # The range to search for new vertices
        rnge = len(self.sub_boxes) // n
        # Get the atoms in the
        med_boxes = [[[[] for i in range(n)] for j in range(n)] for k in range(n)]
        # Set up the network list and the list of the atom indices from the main network for reference later
        test_nets, test_nets_real_ndxs = [], []
        # Calculate the overlap needed to prevent missing vertices
        overlap = int((self.beta_val + self.max_atom_rad) / min(self.sub_box_size)) + 1
        # Go through each of the boxes in the medium_boxes matrix to create networks
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    # Starting from the middle of the range search out the atoms within 0.75 + the increment
                    med_boxes[i][j][k] = self.get_atoms([[int((i + .5) * rnge), int((j + .5) * rnge),
                                                          int((k + .5) * rnge)]], int(rnge // 2) + overlap)
                    # Skip the net if there are no atoms in the quadrant
                    if len(med_boxes) == 0:
                        continue
                    # Get the increment to spread each net out by
                    inc = overlap
                    # Make sure there are enough atoms to do a valid search
                    while len(med_boxes[i][j][k]) < 30:
                        # Starting from the middle of the range search out the atoms within 0.75 + the increment
                        med_boxes[i][j][k] = self.get_atoms([[int((i + .5) * rnge), int((j + .5) * rnge),
                                                              int((k + .5) * rnge)]], int(rnge // 2) + inc)
                        inc += 1
                    # Set up the indices' tracker list
                    new_atoms, old_ndxs = [], []
                    # Create a copy of each atom in the network
                    for atom in med_boxes[i][j][k]:
                        new_atoms.append(Atom(atom.loc, atom.rad))
                        old_ndxs.append(self.atoms.index(atom))
                    # Create the network and add it to the list of subnetworks
                    test_nets.append(Network(sys=self.sys, atoms=new_atoms, box_size=1.1, min_dist=self.min_dist,
                                             beta_val=self.beta_val, sol_verts=self.sol_verts))
                    # Store the atoms' real indices
                    test_nets_real_ndxs.append(old_ndxs)
        # Find the vertices for each of the networks
        for i in range(len(test_nets)):
            # Get the subnetwork
            net2 = test_nets[i]
            # Sort the atoms in the network
            net2.sort_atoms()
            # Set the name of the subnetwork and print the progress
            net2.name = "Subnetwork " + str(i + 1) + "/" + str(len(test_nets)) + " of " + self.name
            print("\rFinding " + net2.name + " - " + str(len(net2.atoms)) + " atoms", end="")
            # Find the vertices
            net2.find_verts(n)
        # Go through each of the vertices found in each of the subnetworks filtering out repeat/large verts + doublets
        for i in range(len(test_nets)):
            # Grab the subnetwork
            net2 = test_nets[i]
            # Sort through the vertices in each network
            for j in range(len(net2.verts)):
                # Grab the vertex
                vert = net2.verts[j]
                # Set the vertex's actual atoms and use them to set the vertex's index, then sort it
                vert.atoms = [self.atoms[test_nets_real_ndxs[i][ndx]] for ndx in vert.ndx]
                vert.ndx = [self.atoms.index(atom) for atom in vert.atoms]
                vert.ndx.sort()
                # Look for the vertex in the network
                v_ndx = search_verts(self.vert_ndxs, vert.ndx)
                # If found, verify it is not a doublet
                if 0 == len(self.vert_ndxs) or len(self.vert_ndxs) <= v_ndx:
                    vert.myNet = self
                    self.vert_ndxs.append(vert.ndx)
                    self.verts.append(vert)
                # If the indices match we need to check if it is a doublet
                elif vert.ndx == self.vert_ndxs[v_ndx]:
                    # If the location is different, and we haven't indicated the vertex as a doublet already add it
                    if [round(vert.loc[k], 7) for k in range(3)] != \
                            [round(self.verts[v_ndx].loc[k], 7) for k in range(3)] \
                            and not self.verts[v_ndx].doublet:
                        # Add the doublet
                        self.verts[v_ndx].doublet = True
                        self.verts[v_ndx].loc2, self.verts[v_ndx].rad2 = vert.loc, vert.rad
                # If this is a novel vertex, we need to verify it
                else:
                    if verify_site(vert, self):
                        # Add the vertex to the system
                        vert.myNet = self
                        self.vert_ndxs.insert(v_ndx, vert.ndx)
                        self.verts.insert(v_ndx, vert)

    # Find vertices method. Using the functions in find_vertices.py finds the vertices in the network
    def find_verts(self, n=None):
        # For small systems (<= 200) run the normal algorithm
        if len(self.atoms) <= 300 or n == 2:
            find_vertices(self)
        # For large systems, split the atoms into separate smaller networks top search for vertices in
        else:
            self.find_subnets()

    # Build edges function. Builds the edges in the network for use in the surfaces
    def build_edges(self):
        # Go through the edges in the network
        for edge in self.edges:
            edge.build()

    # Build network function. Takes in a system and returns a fully connected network
    def build_surfs(self):
        # Make each surface
        for i in range(len(self.surfs)):
            # If the network is a voronota network, use build_vta method
            if self.flat_faces:
                self.surfs[i].build_vta()
            # Otherwise, proceed with the regular build method
            else:
                print("\rBuilding surface " + str(i + 1) + "/" + str(len(self.surfs)), end="")
                self.surfs[i].build()

    # Analyze system function. Finds the surfaces and volumes of the system
    def analyze(self):
        # Get the percentage total number
        tot_num = len(self.surfs) + len(self.atoms)
        # Go through each surface in the system and find the simplices and the surface area
        i = 0
        for i in range(len(self.surfs)):
            percentage = int((i + 1) / tot_num * 100)
            print("\rAnalyzing: {} %".format(percentage), end="")
            # Get the surface area of the surface
            self.surfs[i].interface_sa = calc_sa(self.surfs[i])
        # Go through each atom in the system and find the volume
        for j in range(len(self.atoms)):
            percentage = int((i + j + 1) / tot_num * 100)
            print("\rAnalyzing: {} %".format(percentage), end="")
            self.atoms[j].cell_vol = calc_vol(self.atoms[j])
