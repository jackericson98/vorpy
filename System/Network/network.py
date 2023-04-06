import time
from itertools import chain as chain
from System.Network.net_funcs.find_verts import find_verts
from System.Network.net_funcs.build_net import build, get_time
from System.Network.net_funcs.build_edge import build_edge
from System.Network.net_funcs.build_surf import build_surf
from System.sys_funcs.calcs import calc_vol, calc_surf_func, calc_surf_sa
from Visualize.mpl_visualize import *


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, sys, atoms=None, verts=None, edges=None, surfs=None, surf_res=0.2, box_size=1.25, max_vert=9,
                 calc_verts=True, connect_net=True, build_surfs=True, net_type='vor', surf_col='plasma', surf_scheme='curv'):

        # Main network defining objects
        self.sys = sys                    # System            : Route back to outer system for system attribute access
        self.type = net_type              # Network Type      : String indicating network build type

        # Network element lists
        self.atoms = atoms                # Atoms             :    List of atom objects
        self.verts = verts                # Vertices          :    List of vertex objects
        self.edges = edges                # Edges             :    List of edge objects
        self.surfs = surfs                # Surfaces          :    List of surface objects

        # Index tracking for network elements
        self.atom_ndxs = []                # Atom indices     :    Atom visitation ledger for network building
        self.vert_ndxs = []                # Vert indices     :    Sorted atom indices defining all net verts
        self.edge_ndxs = []                # Edge indices     :    Sorted atom indices defining all net edges
        self.surf_ndxs = []                # Surf indices     :    Sorted atom indices defining all net surfs

        # Tools for splitting up the atoms
        self.box = None                    # Box              :    Holds a max and min vertex for the retaining box
        self.sub_boxes = None              # Sub boxes        :    3D array holding atoms relative locations
        self.vert_sub_boxes = None         # Vert Sub Boxes   :    Holds the vertices of the network by their location
        self.sub_box_size = None           # Sub box size     :    Holds the size of each sub box
        self.box_max = None                # Box maxes        :    number of x, y, z boxes or rows, columns, aisles
        self.atoms_box = []                # Atoms box        :    min/max vals for the box containing the atoms

        # Diagnostic variables
        self.my_time = None                # My time          :    Time taken to calculate the network
        self.max_vert_rad = 0              # Max Vertex Rad   :    Maximum real vertex recorded

        # Build settings
        self.surf_res = surf_res           # Resolution       :    How small the triangles in the surfaces are
        self.surf_col = surf_col           # Color map        :    How the surfaces are colored
        self.surf_scm = surf_scheme        # Coloring scheme  :    How the surfaces will be colored
        self.max_vert = max_vert           # Max vert rad     :    The maximum vertex radius for the network
        self.box_size = box_size           # Box size         :    Retaining box multiplier
        self.calc_verts = calc_verts       # Calc Verts       :    Calculate the vertices
        self.connect_net = connect_net     # Connect net      :    Connect the network's objects
        self.build_surfs = build_surfs     # Calc Surfs       :    Calculate the network's surfaces

        self.sort_atoms()

    def calc_box(self):
        """
        Determines the dimensions of a box x times the size of the atoms
        :return: Sets the self.box attribute with the correct values as well as self.atoms_box
        """
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = np.array([np.inf, np.inf, np.inf])
        max_vert = np.array([-np.inf, -np.inf, -np.inf])
        # Loop through each atom in the network
        for atom in self.atoms:
            # Loop through x, y, z
            for i in range(3):
                # If x, y, z values are less replace the value in the mins list
                if atom.loc[i] < min_vert[i]:
                    min_vert[i] = atom.loc[i]
                # If x, y, z values are greater replace the value in the maxes list
                if atom.loc[i] > max_vert[i]:
                    max_vert[i] = atom.loc[i]
        # Get the vector between the minimum and maximum vertices for the defining box
        r_box = max_vert - min_vert
        # If the atoms are in the same plane adjust the atoms
        for i in range(3):
            if r_box[i] == 0 or abs(r_box[i]) == np.inf:
                r_box[i], min_vert[i], max_vert[i] = 4 * self.atoms[0].rad, self.atoms[0].loc[i], self.atoms[0].loc[i]
        # Set the atoms box value
        self.atoms_box = [min_vert.tolist(), max_vert.tolist()]
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert - r_box * self.box_size, min_vert + r_box * self.box_size
        # Return the list of array turned list vertices
        self.box = [[round(_, 3) for _ in min_vert], [round(_, 3) for _ in max_vert]]

    def sort_atoms(self, num_boxes=None):
        """
        Puts the atoms in the network in their respective grid sections
        :param num_boxes: The number of sub boxes the network is divided into
        :return: Sets the values for self.sub_boxes with the atom objects in their correct locations. Also sets the
        sub-box locations for the atoms themselves
        """
        # Check that the length of the atoms list is big enough to make a vertex
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
        self.sub_box_size = [round((self.box[1][i] - self.box[0][i]) / n, 3) for i in range(3)]
        # Sort the atoms
        for atom in self.atoms:
            # Adjust the maximum radius
            if atom.rad > self.sys.max_atom_rad:
                self.sys.max_atom_rad = atom.rad
            # Find the box they belong to
            box_ndxs = [int((atom.loc[i] - self.box[0][i]) / self.sub_box_size[i]) for i in range(3)]
            # Add the atom to the box
            self.sub_boxes[box_ndxs[0]][box_ndxs[1]][box_ndxs[2]].append(atom)
            # Add the box to the atom
            atom.box = box_ndxs
        # Get the number of rows columns and aisles
        self.box_max = len(self.sub_boxes) - 1, len(self.sub_boxes[0]) - 1, len(self.sub_boxes[0][0]) - 1

    def sort_verts(self):
        """
        Puts the vertices in the network in their respective grid sections
        :return: Places the vertices into their correct sub_boxes
        """
        # Instantiate the grid structure of lists is locations representing a grid
        self.vert_sub_boxes = [[[[] for _ in range(self.box_max[2] + 1)]
                                for _ in range(self.box_max[1] + 1)] for _ in range(self.box_max[0] + 1)]
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
        # Get atoms
        atoms = [self.sub_boxes[i][j][k] for k in zs for j in ys for i in xs
                 if 0 <= k <= self.box_max[2] and 0 <= j <= self.box_max[1] and 0 <= i <= self.box_max[0]]
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
        verts = [self.vert_sub_boxes[i][j][k] for k in zs for j in ys for i in xs if 0 < k < self.box_max[2]
                 and 0 < j < self.box_max[1] and 0 < i < self.box_max[0]]
        verts = list(chain.from_iterable(verts))

        return verts

    def connect(self):
        """
        Connects the network using the functions in the build_net.py file
        :return:
        """
        build(self)

    def find_verts(self, time_start=None, process_time_start=None, my_group=None):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        :return:
        """
        # Check to see if a group has been provided
        if my_group is not None:
            atom_nums = my_group.atom_ndxs
        else:
            atom_nums = [i for i in range(len(self.atoms))]
        # Get the indices of the atoms in the network to keep track of the atoms that haven't been visited
        self.atom_ndxs = [_ for _ in atom_nums]
        # Do an initial sweep
        find_verts(self, my_group=my_group)
        # Check for disconnects in the network
        while len(self.atom_ndxs) > 0:
            find_verts(self, a0=self.atoms[self.atom_ndxs.pop()], my_group=my_group)
        # Clear the print statement
        print("\r                                                                  ", end="")
        # # Bit of code for timing the vertex building process
        # if time_start is not None:
        #     self.my_time = time.time() - time_start
        #     process_time = time.process_time() - process_time_start
        #     h, m, s = get_time(self.my_time)
        #     if self.sys.print_actions:
        #         print("\rvertex process ({} verts) = {}:{}:{:.2f} s, cpu time = {}"
        #               .format(len(self.verts), int(h), int(m), s, process_time))

    def build_edges(self):
        """
        Builds the edges in the network for use in the surfaces
        :return:
        """
        # Go through the edges in the network
        for edge in self.edges:
            # Build the edge depending on if it is straight or not
            build_edge(edge, straight=True if self.type in ['pow', 'flat', 'del'] else False)

    def build_surfaces(self):
        """
        Takes in a system and returns a fully connected network
        :return:
        """
        # Make each surface
        for i in range(len(self.surfs)):
            # Build the surfaces and print the progress
            my_time = time.perf_counter() - self.my_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {:2}:{:2}:{:.2f} - Process: building surfaces {:.2f} %                                 "
                  .format(int(h), int(m), round(s, 2), min(100.0, 100 * round(i/len(self.surfs), 2))), end="")
            build_surf(self.surfs[i])
        print("\r                                                                                             ", end='')

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
            # If the surface area is None calculate it
            if self.surfs[i].sa is None or self.surfs[i].sa == 0:
                # Get the surface area of the surface
                surf = self.surfs[i]
                surf.sa = calc_surf_sa(edges=surf.edges, com=surf.com, tris=surf.tris, points=surf.points, flat=surf.flat)
            if self.sys.print_actions:
                my_time = time.perf_counter() - self.my_time
                h, m, s = get_time(my_time)
                print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                  "
                      .format(int(h), int(m), round(s, 2), percentage), end="")
        # Go through each atom in the system and find the volume
        for j in range(len(self.atoms)):
            percentage = int((i + j + 2) / tot_num * 100)
            if self.atoms[j].vol is None or self.atoms[j].vol == 0:
                calc_vol(self.atoms[j])
            if self.sys.print_actions:
                my_time = time.perf_counter() - self.my_time
                h, m, s = get_time(my_time)
                print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 ".format(int(h), int(m), round(s, 2), percentage), end="")

    def build(self, output=True, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None,
              calc_verts=None, my_group=None, print_actions=None):
        """
        Build network function used to calculate the voronoi
        :param print_actions: Print the network building actions
        :param net_type: Describes the network construction type ('curv', 'del', 'pow')
        :param my_group: Describes the group of atoms (mols, resids, etc) for network construction
        :param output: Output information for the group? If yes all gets outputted
        :param surf_res: Resolution for the surface construction
        :param max_vert: Maximum allowed vertex size in the network construciton
        :param box_size: Maximum box multiplier for the retaining box
        :param build_surfs: Build Surfaces? If yes, the surfaces in the group's network are constructed
        :param calc_verts: Calculate Vertices? Skips vertex calculations if a network or vertex file is loaded
        :return: Builds the network based on the above specifications
        """
        if self.sys.ball_file is None:
            # Reset the network variables in case of rebuild
            self.verts, self.vert_ndxs, self.edges, self.edge_ndxs, self.surfs, self.surf_ndxs, self.atom_ndxs = [], [], [], [], [], [], []
        # If the system has no name, one needs top be set
        if self.sys.name is None:
            self.sys.name = "User_Atoms"
        if print_actions is not None:
            self.sys.print_actions = print_actions
        # Check for input values for the network build
        if surf_res is not None:
            self.surf_res = surf_res
        if max_vert is not None:
            self.max_vert = max_vert
        if box_size is not None:
            self.box_size = box_size
        if build_surfs is not None:
            self.build_surfs = build_surfs
        if net_type is not None:
            self.type = net_type
        if calc_verts is not None:
            self.calc_verts = calc_verts
        # Instantiate the timer variables
        self.my_time = 0
        # Start the timer
        self.my_time = time.perf_counter()
        # Sort the atoms in the network
        self.sort_atoms()
        # Check to see if there are vertices loaded
        if self.calc_verts and self.sys.ball_file is None:
            # Find the vertices
            self.find_verts(my_group=my_group)
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
                surf.func = calc_surf_func(surf.atoms[0].loc, surf.atoms[0].rad, surf.atoms[1].loc, surf.atoms[1].rad)
        # Load the elements to the group
        my_group.get_info()
        # Stop the timer and measure the time
        self.my_time = time.process_time() - self.my_time
        # Export the network
        if output:
            self.sys.exports(network=True, pdb=True, info=self.build_surfs, set_atoms=True)
        h, m, s = get_time(self.my_time)
        print("\rnetwork built - {} verts, {} surfs - {}:{}:{:.2f} s\n"
              .format(len(self.verts), len(self.surfs), int(h), int(m), s), end="")
