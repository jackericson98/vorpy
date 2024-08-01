import time
from _datetime import datetime
import pandas as pd
import csv
import os
from System.Network.verts.mark_doublets import mark_doublets
from System.Network.verts.find_net_verts import find_net_verts
from System.Network.build_net import build
from System.Network.edges.build_edge import build_edge
from System.Network.surfs.build_surfs import build_surfs
from System.sys_funcs.calcs.calcs import calc_length, get_time
from System.sys_funcs.calcs.surf import calc_surf_func
from numpy import array, inf, cbrt, sqrt, pi


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, sys, atoms=None, verts=None, edges=None, surfs=None, surf_res=0.2, box_size=1.5, max_vert=40,
                 calc_verts=True, connect_net=True, build_surfs=True, net_type='vor', surf_col='plasma',
                 surf_scheme='curv', sub_net=False, box=None, vta_verts=None, settings=None):

        # Main network defining objects
        self.num_splits = None
        self.sys = sys                    # System            :   Route back to outer system
        self.id = 0                       # Network id #      :
        self.sub_net = sub_net
        self.settings = settings          # Settings          :    surf_res, surf_col, surf_schm, max_vert, net_type

        # Network element lists
        self.atoms = atoms                # Atoms             :    List of atom objects
        self.verts = verts                # Vertices          :    List of vertex objects
        self.vta_verts = vta_verts        # Voronota Vertices :    List of Voronota vertices
        self.edges = edges                # Edges             :    List of edge objects
        self.surfs = surfs                # Surfaces          :    List of surface objects

        # Index tracking for network elements
        self.atom_ndxs = []                # Atom indices     :    Atom visitation ledger for network building
        self.vert_ndxs = []                # Vert indices     :    Sorted atom indices defining all net verts
        self.edge_ndxs = []                # Edge indices     :    Sorted atom indices defining all net edges
        self.surf_ndxs = []                # Surf indices     :    Sorted atom indices defining all net surfs

        # Tools for splitting up the atoms
        self.box = box                     # Box              :    Holds a max and min vertex for the retaining box
        self.sub_boxes = None              # Sub boxes        :    3D array holding atoms relative locations
        self.vert_sub_boxes = None         # Vert Sub Boxes   :    Holds the vertices of the network by their location
        self.sub_box_size = None           # Sub box size     :    Holds the size of each sub box
        self.box_max = None                # Box maxes        :    number of x, y, z boxes or rows, columns, aisles
        self.atoms_box = []                # Atoms box        :    min/max vals for the box containing the atoms

        # Diagnostic variables
        self.start_time = time.perf_counter()
        self.metrics = {}                  # Build Metrics    :    Holds the time measurements for the build
        self.max_vert_rad = 0              # Max Vertex Rad   :    Maximum real vertex recorded
        self.max_curv = 0

        # Build settings
        self.calc_verts = calc_verts       # Calc Verts       :    Calculate the vertices
        self.connect_net = connect_net     # Connect net      :    Connect the network's objects
        self.build_surfs = build_surfs     # Calc Surfs       :    Calculate the network's surfaces

        # Get the settings
        self.get_settings(surf_res=surf_res, surf_col=surf_col, surf_scheme=surf_scheme, max_vert=max_vert,
                          box_size=box_size, net_type=net_type)

    def get_settings(self, surf_res=None, surf_col=None, surf_scheme=None, max_vert=None, box_size=None, net_type=None):
        """
        Sets the settings for the network building
        """
        # Set up the default values
        defaults = {'surf_res': 0.2, 'surf_col': 'plasma', 'surf_scheme': 'curv', 'max_vert': 40, 'box_size': 1.5,
                    'net_type': 'aw'}
        # Create the settings dictionary
        if self.settings is None:
            self.settings = {'surf_res': None, 'surf_col': None, 'surf_scheme': None, 'max_vert': None, 'box_size': None,
                             'net_type': None}
        # Set the settings to their default values
        for setting in self.settings:
            if self.settings[setting] is None:
                self.settings[setting] = defaults[setting]

    def calc_box(self, locs, rads, return_val=False, box_size=None):
        """
        Determines the dimensions of a box x times the size of the atoms
        :return: Sets the box attribute with the 181L values as well as atoms_box
        """
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = array([inf, inf, inf])
        max_vert = array([-inf, -inf, -inf])
        if box_size is None:
            box_size = self.settings['box_size']
        # Loop through each atom in the network
        for loc in locs:
            # Loop through x, y, z
            for i in range(3):
                # If x, y, z values are less replace the value in the mins list
                if loc[i] < min_vert[i]:
                    min_vert[i] = loc[i]
                # If x, y, z values are greater replace the value in the maxes list
                if loc[i] > max_vert[i]:
                    max_vert[i] = loc[i]
        # Get the vector between the minimum and maximum vertices for the defining box
        r_box = max_vert - min_vert
        # If the atoms are in the same plane adjust the atoms
        for i in range(3):
            if r_box[i] == 0 or abs(r_box[i]) == inf:
                r_box[i], min_vert[i], max_vert[i] = 4 * rads[0], locs[0][i], locs[0][i]
        # Set the atoms box value
        atoms_box = [min_vert.tolist(), max_vert.tolist()]
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert - r_box * box_size, min_vert + r_box * box_size
        # Return the list of array turned list vertices
        box = [[round(_, 3) for _ in min_vert], [round(_, 3) for _ in max_vert]]
        # If the values are to be returned
        if return_val:
            return box
        self.atoms_box, self.box = atoms_box, box

    def sort_atoms(self, num_boxes=None):
        """
        Puts the atoms in the network in their respective grid sections
        :param num_boxes: The number of sub boxes the network is divided into
        :return: Sets the values for self.sub_boxes with the atom objects in their 181L locations. Also sets the
        sub-box locations for the atoms themselves
        """
        # Check that the length of the atoms list is big enough to make a vertex
        if len(self.atoms) < 4:
            return
        # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
        elif num_boxes is None:
            n = int(0.5 * sqrt(len(self.atoms))) + 1
        else:
            n = int(cbrt(num_boxes)) + 1
        self.num_splits = n
        locs, rads = self.atoms['loc'], self.atoms['rad']
        # First get the box for the atoms to be sorted into
        self.calc_box(locs, rads)
        # Instantiate the grid structure of lists is locations representing a grid
        self.sub_boxes = {(-1, -1, -1): [n]}
        # Get the cell size
        self.sub_box_size = [round((self.box[1][i] - self.box[0][i]) / n, 3) for i in range(3)]
        my_boxes = []
        # Sort the atoms
        for i, loc in enumerate(locs):
            # Get the radius
            rad = rads[i]
            # Adjust the maximum radius
            if rad > self.sys.max_atom_rad:
                self.sys.max_atom_rad = rad
            # Find the box they belong to
            box_ndxs = [int((loc[j] - self.box[0][j]) / self.sub_box_size[j]) for j in range(3)]

            # Add the atom to the box
            try:
                self.sub_boxes[box_ndxs[0], box_ndxs[1], box_ndxs[2]].append(i)
            except KeyError:
                self.sub_boxes[box_ndxs[0], box_ndxs[1], box_ndxs[2]] = [i]
            # Add the box to the atom
            my_boxes.append(box_ndxs)
        # Get the number of rows columns and aisles
        self.box_max = n - 1, n - 1, n - 1
        # set the box data
        self.atoms['box'] = my_boxes

    def find_verts(self, my_group=None, print_metrics=False):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        """
        find_net_verts(self, my_group=my_group, print_metrics=print_metrics)

    def connect(self):
        """
        Connects the network using the functions in the build_net.py file
        """
        my_lists = build(self.verts['vatoms'], self.verts['vloc'], self.verts['vdub'], len(self.atoms), self.start_time)
        atom_lists, vert_lists, edge_lists, surf_lists = my_lists
        self.atoms['averts'], self.atoms['aedges'], self.atoms['asurfs'] = atom_lists['averts'], atom_lists['aedges'], \
            atom_lists['asurfs']
        self.verts['vedges'], self.verts['vsurfs'] = vert_lists['vedges'], vert_lists['vsurfs']
        self.edges = pd.DataFrame(edge_lists)
        self.surfs = pd.DataFrame(surf_lists)
        self.metrics['con'] = time.perf_counter() - self.start_time - self.metrics['vert']

    def get_real_verts(self):
        my_name = os.getcwd() + '/Data/user_data/' + self.sys.name + '_Correct/sys/' + self.sys.name + '_logs.csv'
        if not os.path.exists(my_name):
            return
        with open(my_name) as csvfile:
            my_logs = csv.reader(csvfile, delimiter=',')
            at_verts = False
            vert_ndxs = []
            my_i = 0
            for i, line in enumerate(my_logs):
                if line[0] == 'Vertices':
                    at_verts = True
                    my_i = i
                    continue
                if at_verts and i > my_i + 1:
                    vert_ndxs.append([int(_) for _ in line[1:5]])
        return vert_ndxs

    def build_edges(self):
        """
        Builds the edges in the network for use in the surfaces
        """
        # Set the edge points and vals lists
        edges_points, edges_vals, edges_lengths = [], [], []
        # Go through the edges in the network
        for i, edge in self.edges.iterrows():
            # Build the edge depending on if it is straight or not
            straight = True if self.settings['net_type'] in ['pow', 'flat', 'del'] else False
            edge_points, edge_vals = build_edge(alocs=[array(self.atoms['loc'][_]) for _ in edge['eatoms']],
                                                arads=[self.atoms['rad'][_] for _ in edge['eatoms']],
                                                vlocs=[array(self.verts['vloc'][_]) for _ in edge['everts']],
                                                res=self.settings['surf_res'], straight=straight)
            # Add them to the lists
            edges_lengths.append(calc_length(array(edge_points)))
            edges_points.append(edge_points)
            edges_vals.append(edge_vals)
        # Set the dataframe values
        self.edges['points'], self.edges['vals'], self.edges['length'] = edges_points, edges_vals, edges_lengths

    def build_surfaces(self, store_points=True):
        """
        Takes in a system and returns a fully connected network
        """
        build_surfs(self, store_points=store_points)

    def analyze(self):
        """
        Analyzes the output surfaces, cells and solute vertices for the network for later reference
        """
        # Set up the atoms' volumes surface areas, curvatures vars
        avols, asas, acurvs, acell = [], [], [], []
        # Go through each atom in the system and find the volume
        for k, atom in self.atoms.iterrows():
            # Get the percentage for printing
            percentage = int(k / len(self.atoms['loc']) * 100)
            if len(atom['asurfs']) == 0:
                avols.append(0)
                asas.append(0)
                acurvs.append(0)
                acell.append(False)
                continue
            # Calculate the surface area of the atom by summing the surface areas of all it's surfaces
            asas.append(sum([self.surfs['sa'][_] for _ in atom['asurfs']]))
            # Go through the atom's surfaces
            acurvs.append(max([self.surfs['curv'][_] for _ in atom['asurfs']]))
            # Calculate the volume of the atom by the previouslty stored volume data
            avol = sum([self.surfs['vols'][_][atom['num']] for _ in atom['asurfs']])
            # Exclude atoms that have super large volumes (weird edge error)
            bad_atom = False
            if avol > 15 * 4/3 * atom['rad'] ** 3 * pi:
                bad_atom = True
            avols.append(avol)
            # Check for complete cells in the atoms
            complete = True
            # Go through each of the vertices in the in the atom
            for vert in atom['averts']:
                # Check the number of edges from the vertex that hold
                if len([_ for _ in [self.edges['eatoms'][_] for _ in self.verts['vedges'][vert]] if k in _]) != 3 and not bad_atom:
                    complete = False
            # Additional catch for any atom that doesn't have the 181L number of network elements associated with it
            if len(atom['averts']) < 3 or len(atom['aedges']) < 4 or len(atom['asurfs']) < 3:
                complete = False
            # Add the complete designation for the cell
            acell.append(complete)
            # Print the actions
            my_time = time.perf_counter() - self.start_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 "
                  .format(int(h), int(m), round(s, 2), percentage), end="")
        self.atoms['vol'], self.atoms['sa'], self.atoms['curv'], self.atoms['complete'] = avols, asas, acurvs, acell
        self.metrics['anal'] = time.perf_counter() - self.start_time - self.metrics['surf'] - self.metrics['con'] - self.metrics['vert']

    def build(self, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None,
              calc_verts=None, my_group=None, print_actions=None, print_vert_metrics=False, curr_time=None):
        """
        Build network function used to calculate the voronoi
        :param print_actions: Print the network building actions
        :param net_type: Describes the network construction type ('curv', 'del', 'pow')
        :param my_group: Describes the group of atoms (mols, resids, etc) for network construction
        :param surf_res: Resolution for the surface construction
        :param max_vert: Maximum allowed vertex size in the network construciton
        :param box_size: Maximum box multiplier for the retaining box
        :param build_surfs: Build Surfaces? If yes, the surfaces in the group's network are constructed
        :param calc_verts: Calculate Vertices? Skips vertex calculations if a network or vertex file is loaded
        """
        self.start_time = time.perf_counter()
        # If the system has no name, one needs top be set
        if self.sys.name is None:
            self.sys.name = "User_Atoms"
        if print_actions is not None:
            self.sys.print_actions = print_actions
        # Check for input values for the network build
        if surf_res is not None:
            self.settings['surf_res'] = surf_res
        if max_vert is not None:
            self.settings['max_vert'] = max_vert
        if box_size is not None:
            self.settings['box_size'] = box_size
        if build_surfs is not None:
            self.build_surfs = build_surfs
        if net_type is not None:
            self.settings['net_type'] = net_type
        if calc_verts is not None:
            self.calc_verts = calc_verts
        # Check to see if the only output for the exports is logs
        limit_mem = False
        if self.sys.cmnds['xpt'] == [['logs']]:
            limit_mem = True
        # Sort the atoms in the network
        self.sort_atoms()
        # Check to see if there are vertices loaded
        if self.calc_verts and self.sys.ball_file is None:
            # Find the vertices
            self.find_verts(my_group=my_group, print_metrics=print_vert_metrics)
            # Check to see if there are vertices
            if self.verts is None or len(self.verts) == 0:
                return
        elif self.sys.ball_file != 'deez nuts':
            self.metrics['vert'] = 0
            # Filter out the vertices that don't pertain to the group in question
            verts = []
            for i, vert in self.vta_verts.iterrows():
                if any([True if _ in my_group.atoms else False for _ in vert['vatoms']]):
                    verts.append(vert)
            self.verts = pd.DataFrame(verts)
        elif 'vdub' not in self.verts:
            self.metrics['vert'] = 0
            self.verts['vdub'] = mark_doublets(self.verts)
        else:
            self.metrics['vert'] = 0

        # Connect the network
        self.connect()
        # Build the edges in the network
        self.build_edges()
        if self.build_surfs:
            # Build the network
            self.build_surfaces(not limit_mem)
            # Analyze the network
            self.analyze()
        else:
            for surf in self.surfs:
                a0, a1 = self.atoms.iloc[surf['atoms'][0]], self.atoms.iloc[surf['atoms'][1]]
                surf_atoms_vals = a0['loc'], a0['rad'], a1['loc'], a1['rad']
                surf['func'] = calc_surf_func(*surf_atoms_vals)
            self.metrics['surf'], self.metrics['anal'] = 0, 0
        # Load the elements to the group
        if my_group is not None:
            my_group.get_info()
        # Stop the timer and measure the time
        self.metrics['tot'] = time.perf_counter() - self.start_time
        h, m, s = get_time(self.metrics['tot'])
        print("\rnetwork built - {} atoms, {} verts, {} surfs - {}:{}:{:.2f} s - finished at {}\n"
              .format(len([_ for _ in self.atoms['complete'] if _]), len(self.verts), len(self.surfs), int(h), int(m), s, datetime.now()), end="")
