import time
from itertools import chain as chain

from System.Network.net_funcs.find_verts import find_verts
from System.Network.net_funcs.build_net import build, get_time
from System.Network.net_funcs.ffind_verts import ffind_verts
from Visualize.mpl_visualize import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, sys, atoms=None, verts=None, edges=None, surfs=None, surf_res=0.3, box_size=1.25, max_vert=7,
                 build_surfs=False, flat_surfs=False):
        # Network graph objects
        self.sys = sys                 # System          : Route back to outer system for system attribute access
        self.atoms = atoms             # Atoms           : Atoms of the network. Should be identical to self.sys.atoms
        self.verts = verts             # Vertices        : Vertices of the network
        self.edges = edges             # Edges           : Edges of the network
        self.surfs = surfs             # Surfaces        : Surfaces of the network
        # Tools for splitting up the atoms
        self.box = None                 # Box            : Holds a max and min vertex for the retaining box
        self.sub_boxes = None           # Sub boxes      : Holds atoms in their different relative locations in the grid
        self.vert_sub_boxes = None      # Vert Sub Boxes : Holds the vertices of the network by their location
        self.sub_box_size = None        # Sub box size   : Holds the size of each sub box
        self.box_max = None             # Box maxes      : Number of x, y, z boxes or rows, columns, aisles
        self.atoms_box = []             # Atoms box      : Holds the min and max verts for the box containing the atoms
        self.max_atom_rad = 0           # Max atom rad   : Holds the largest radius of the system for reference
        self.max_vert_rad = 0           # Max Vertex Rad : Holds the maximum real vertex recorded
        self.vert_ndxs = []             # Vert indices   : Holds the sorted indices of the atoms of the network's verts
        self.bad_verts = []             # Bad Vertices   : Holds vertex indices for the vertices that arent real
        self.edge_ndxs = []             # Edge indices   : Holds the sorted indices of the atoms of the network's edges
        self.surf_ndxs = []             # Surf indices   : Holds the sorted indices of the atoms of the network's surfs
        self.atom_ndxs = []             # Atom indices   : Used to track atoms that have been used in a vertex
        # Settings
        self.surf_res = surf_res        # Resolution     : How small the triangles in the surfaces are
        self.max_vert = max_vert        # Max vert rad   : The maximum vertex radius for the network
        self.box_size = box_size        # Box size       : Holds the box multiplier for the sys box from the atoms box
        self.build_surfs = build_surfs  # Calc Surfs     : Calculate the network's surfaces? Bool
        self.flat_surfs = flat_surfs    # Flat Faces     : Create flat faces for surfaces. Bool
        self.calc_verts = True          # Calc Verts     : Use loaded verts. Bool
        self.connect_net = True         # Connect net    : Used to differentiate between loaded net and loaded verts
        # Run diagnostics
        self.cpu_time = None            # CPU time       : CPU time taken to calculate the network
        self.my_time = None             # My time        : Time taken to calculate the network

    def calc_box(self):
        """
        Takes in a System and returns the dimensions of a box x times the size of the atoms
        :return:
        """
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

    def sort_atoms(self, num_boxes=None):
        """
        Puts the atoms in the network in their respective grid sections
        :param num_boxes: The number of sub boxes the network is divided into
        :return:
        """
        # Check the length of the atoms list is big enough to make a vertex
        if len(self.atoms) < 4:
            return
        # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
        elif num_boxes is None:
            n = int(np.sqrt(len(self.atoms))) + 1
        else:
            n = int(np.cbrt(num_boxes)) + 1
        # First get the box for the atoms to be sorted into
        self.calc_box()
        # Instantiate the grid structure of lists is locations representing a grid
        self.sub_boxes = [[[[] for _ in range(n)] for _ in range(n)] for _ in range(n)]
        # Get the cell size
        self.sub_box_size = [(self.box[1][i] - self.box[0][i]) / n for i in range(3)]
        # Sort the atoms
        for atom in self.atoms:
            # Adjust the maximum radius
            if atom.rad > self.max_atom_rad:
                self.max_atom_rad = atom.rad
            # Find the box they belong to
            box_ndxs = [int((atom.loc[i] - self.box[0][i]) / self.sub_box_size[i]) for i in range(3)]
            # Add the atom to the box
            self.sub_boxes[box_ndxs[0]][box_ndxs[1]][box_ndxs[2]].append(atom)
            # Add the box to the atom
            atom.box = box_ndxs
        self.box_max = len(self.sub_boxes) - 1, len(self.sub_boxes[0]) - 1, len(self.sub_boxes[0][0]) - 1

    def sort_verts(self):
        """
        Puts the verts in the network in their respective grid sections
        :return:
        """
        # Instantiate the grid structure of lists is locations representing a grid
        self.vert_sub_boxes = [[[[] for _ in range(self.box_max[2] + 1)] for _ in range(self.box_max[1] + 1)] for _ in range(self.box_max[0] + 1)]
        # Sort the atoms
        for vert in self.verts:
            # Adjust the maximum radius
            if vert.rad > self.max_vert_rad:
                self.max_vert_rad = vert.rad
            # Find the box they belong to
            box_ndxs = [int((vert.loc[i] - self.box[0][i]) / self.sub_box_size[i]) for i in range(3)]
            # Add the atom to the box
            self.vert_sub_boxes[box_ndxs[0]][box_ndxs[1]][box_ndxs[2]].append(vert)
            # Add the box to the atom
            vert.box = box_ndxs
        self.box_max = len(self.sub_boxes) - 1, len(self.sub_boxes[0]) - 1, len(self.sub_boxes[0][0]) - 1

    def get_atoms(self, cells, reach=0):
        """
        Takes in the cells and the number of additional cells to search and returns an atom list
        :param cells: The initial boxes in the network to stem from
        :param reach: The number of cells out from the initial set of cells to search
        :return:
        """
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
        xs = [x for x in range(max(0, -reach + ndx_min[0] + 1), reach + ndx_max[0])]
        ys = [y for y in range(max(0, -reach + ndx_min[1] + 1), reach + ndx_max[1])]
        zs = [z for z in range(max(0, -reach + ndx_min[2] + 1), reach + ndx_max[2])]
        atoms = [self.sub_boxes[i][j][k] for k in zs for j in ys for i in xs if 0 <= k <= self.box_max[2] and 0 <= j <= self.box_max[1] and 0 <= i <= self.box_max[0]]
        atoms = list(chain.from_iterable(atoms))
        return atoms

    def get_verts(self, cells, reach=0):
        """
        Takes in the cells and the number of additional cells to search and returns an atom list
        :param cells: The initial boxes in the network to stem from
        :param reach: The number of cells out from the initial set of cells to search
        :return:
        """
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
        xs = [x for x in range(max(0, -reach + ndx_min[0] + 1), reach + ndx_max[0])]
        ys = [y for y in range(max(0, -reach + ndx_min[1] + 1), reach + ndx_max[1])]
        zs = [z for z in range(max(0, -reach + ndx_min[2] + 1), reach + ndx_max[2])]
        verts = [self.vert_sub_boxes[i][j][k] for k in zs for j in ys for i in xs if 0 < k < self.box_max[2] and 0 < j < self.box_max[1] and 0 < i < self.box_max[0]]
        verts = list(chain.from_iterable(verts))

        return verts

    def connect(self, get_edges=True, get_surfs=True):
        """
        Connects the network using the functions in the build_net.py file
        :return:
        """
        print("\rconnecting network", end="")
        build(self, get_edges=get_edges, get_surfs=get_surfs)

    def find_verts(self, time_start=None, process_time_start=None):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        :return:
        """
        # Get the indices of the atoms in the network to keep track of the atoms that haven't been visited
        self.atom_ndxs = [i for i in range(len(self.atoms))]
        # Find curved surfaces verts
        if not self.flat_surfs:
            # Do an initial sweep
            find_verts(self)
            # Check for disconnects in the network
            while len(self.atom_ndxs) > 0:
                find_verts(self, a0=self.atoms[self.atom_ndxs.pop()])

        # Find flat surfaces vertices
        else:
            # Do an initial sweep
            ffind_verts(self)
            # Check for disconnects in the network
            while len(self.atom_ndxs) > 0:
                ffind_verts(self, a0=self.atoms[self.atom_ndxs.pop()])
        # Clear the print statement
        print("\r                                        ", end="")
        # Bit of code for timing the vertex building process
        if time_start is not None:
            self.my_time = time.time() - time_start
            process_time = time.process_time() - process_time_start
            h, m, s = get_time(self.my_time)
            print("\rvertex process ({} verts) = {}:{}:{:.2f} s, cpu time = {}".format(len(self.verts), int(h), int(m), s, process_time))

    def build_edges(self):
        """
        Builds the edges in the network for use in the surfaces
        :return:
        """
        # Go through the edges in the network
        for edge in self.edges:
            edge.build(straight=self.flat_surfs)

    def build_surfaces(self):
        """
        Takes in a system and returns a fully connected network
        :return:
        """
        # Make each surface
        for i in range(len(self.surfs)):
            # Build the surfaces and print the progress
            print("\rbuilding surfaces " + " " * (len(str(len(self.surfs) - 1)) - len(str(i + 1))) + str(i + 1) + "/" +
                  str(len(self.surfs)) + "                   ", end="")
            self.surfs[i].build(flat=self.flat_surfs)

    def analyze(self):
        """
        Analyzes the output surfaces, cells and solute vertices for the network for later reference
        :return:
        """
        # Get the percentage total number
        tot_num = len(self.surfs) + len(self.atoms)
        # Go through each surface in the system and find the simplices and the surface area
        i = 0
        for i in range(len(self.surfs)):
            percentage = int((i + 1) / tot_num * 100)
            print("\ranalyzing: {} %            ".format(percentage), end="")
            # Get the surface area of the surface
            self.surfs[i].calc_sa()
        # Go through each atom in the system and find the volume
        for j in range(len(self.atoms)):
            percentage = int((i + j + 2) / tot_num * 100)
            print("\ranalyzing: {} %          ".format(percentage), end="")
            self.atoms[j].calc_vol()
        # # Get the solute layers
        # self.sol_layers, self.sol_layer_atoms = find_sol_layers(self)

    def build(self, output=True, surf_res=None, max_vert=None, box_size=None, build_surfs=None, flat_surfs=None,
              calc_verts=None):
        """
        Build network function used to calculate the voronoi
        :param output:
        :param surf_res:
        :param max_vert:
        :param box_size:
        :param build_surfs:
        :param flat_surfs:
        :param calc_verts:
        :return:
        """
        # If the system has no name, one needs top be set
        if self.sys.name is None:
            self.sys.name = "User_Atoms"
        # Check for input values for the network build
        if surf_res is not None:
            self.surf_res = surf_res
        if max_vert is not None:
            self.max_vert = max_vert
        if box_size is not None:
            self.box_size = box_size
        if build_surfs is not None:
            self.build_surfs = build_surfs
        if flat_surfs is not None:
            self.flat_surfs = flat_surfs
        if calc_verts is not None:
            self.calc_verts = calc_verts
        # Instantiate the timer variables
        self.my_time, self.cpu_time = 0, 0
        # Start the timer
        start = time.time()
        process_start = time.process_time()
        # Sort the atoms in the network
        self.sort_atoms()
        # Check to see if there are vertices loaded
        if self.calc_verts:
            # Find the vertices
            self.find_verts(start, process_time_start=process_start)
            # Check to see if there are vertices
            if self.verts is None or len(self.verts) == 0:
                return
        # Connect the network
        self.connect()
        # Build the edges in the network
        self.build_edges()
        if self.build_surfs:
            # Build the network
            self.build_surfaces()
            # Analyze the network
            self.analyze()
        else:
            for surf in self.surfs:
                surf.calc_func()
        # Stop the timer and measure the time
        self.my_time = time.time() - start
        self.cpu_time = time.process_time() - process_start
        # Export the network
        if output:
            self.sys.exports(network=True, pdb=True, no_sol_network_object=True, info=self.build_surfs)
        h, m, s = get_time(self.my_time)
        print("\rnetwork built - {} verts, {} surfs - {}:{}:{:.2f} s, cpu time = {}\n".format(len(self.verts), len(self.surfs),
                                                                                int(h), int(m), s, self.cpu_time), end="")

    def rebuild_net(self, resolution=None, flat_faces=None, max_vert=None, box_size=None):
        pass
