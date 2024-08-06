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
    def __init__(self, group, settings, balls=None, verts=None, edges=None, surfs=None, sub_net=False, box=None,
                 vta_verts=None):

        # Main network defining objects
        self.group = group                # Group             :    Group that made this, references back to the system
        self.id = 0                       # Network id #      :
        self.settings = settings          # Settings          :    surf_res, surf_col, surf_schm, max_vert, net_type

        # Network element lists
        self.balls = balls                # Spheres           :    List of atom objects
        self.verts = verts                # Vertices          :    List of vertex objects
        self.vta_verts = vta_verts        # Voronota Vertices :    List of Voronota vertices
        self.edges = edges                # Edges             :    List of edge objects
        self.surfs = surfs                # Surfaces          :    List of surface objects

        # Index tracking for network elements
        self.ball_ndxs = []                # Ball indices     :    Atom visitation ledger for network building
        self.vert_ndxs = []                # Vert indices     :    Sorted atom indices defining all net verts
        self.edge_ndxs = []                # Edge indices     :    Sorted atom indices defining all net edges
        self.surf_ndxs = []                # Surf indices     :    Sorted atom indices defining all net surfs

        # Tools for splitting up the atoms
        self.box = box                     # Box              :    Holds a max and min vertex for the retaining box
        self.sub_boxes = None              # Sub boxes        :    3D array holding atoms relative locations
        self.vert_sub_boxes = None         # Vert Sub Boxes   :    Holds the vertices of the network by their location
        self.sub_box_size = None           # Sub box size     :    Holds the size of each sub box
        self.box_max = None                # Box maxes        :    number of x, y, z boxes or rows, columns, aisles
        self.ball_box = []                 # Ball box         :    min/max vals for the box containing the atoms

        # Diagnostic variables
        self.start_time = time.perf_counter()
        self.metrics = {}                  # Build Metrics    :    Holds the time measurements for the build
        self.max_vert_rad = 0              # Max Vertex Rad   :    Maximum real vertex recorded
        self.max_curv = 0

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
        ball_box = [min_vert.tolist(), max_vert.tolist()]
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert - r_box * box_size, min_vert + r_box * box_size
        # Return the list of array turned list vertices
        box = [[round(_, 3) for _ in min_vert], [round(_, 3) for _ in max_vert]]
        # If the values are to be returned
        if return_val:
            return box
        self.ball_box, self.box = ball_box, box

    def sort_balls(self, num_boxes=None):
        """
        Puts the atoms in the network in their respective grid sections
        :param num_boxes: The number of sub boxes the network is divided into
        :return: Sets the values for self.sub_boxes with the atom objects in their 181L locations. Also sets the
        sub-box locations for the atoms themselves
        """
        # Check that the length of the spheres list is big enough to make a vertex
        if len(self.balls) < 4:
            return
        # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
        elif num_boxes is None:
            n = int(0.5 * sqrt(len(self.balls))) + 1
        else:
            n = int(cbrt(num_boxes)) + 1
        self.settings['num_splits'] = n
        locs, rads = self.balls['loc'], self.balls['rad']
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
            if rad > self.group.sys.max_atom_rad:
                self.group.sys.max_atom_rad = rad
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
        self.balls['box'] = my_boxes

    def find_verts(self):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        """
        find_net_verts(self)

    def connect(self):
        """
        Connects the network using the functions in the build_net.py file
        """
        my_lists = build(self.verts['vatoms'], self.verts['vloc'], self.verts['vdub'], len(self.balls), self.start_time)
        ball_lists, vert_lists, edge_lists, surf_lists = my_lists
        self.balls['averts'], self.balls['aedges'], self.balls['asurfs'] = ball_lists['averts'], ball_lists['aedges'], \
            ball_lists['asurfs']
        self.verts['vedges'], self.verts['vsurfs'] = vert_lists['vedges'], vert_lists['vsurfs']
        self.edges = pd.DataFrame(edge_lists)
        self.surfs = pd.DataFrame(surf_lists)
        self.metrics['con'] = time.perf_counter() - self.start_time - self.metrics['vert']

    def get_real_verts(self):
        my_name = os.getcwd() + '/Data/user_data/' + self.group.sys.name + '_Correct/sys/' + self.group.sys.name + '_logs.csv'
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
            edge_points, edge_vals = build_edge(locs=[array(self.balls['loc'][_]) for _ in edge['eatoms']],
                                                rads=[self.balls['rad'][_] for _ in edge['eatoms']],
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
        b_vols, b_sas, b_curvs, b_cell = [], [], [], []
        # Go through each atom in the system and find the volume
        for k, ball in self.balls.iterrows():
            # Get the percentage for printing
            percentage = int(k / len(self.balls['loc']) * 100)
            if len(ball['asurfs']) == 0:
                b_vols.append(0)
                b_sas.append(0)
                b_curvs.append(0)
                b_cell.append(False)
                continue
            # Calculate the surface area of the atom by summing the surface areas of all it's surfaces
            b_sas.append(sum([self.surfs['sa'][_] for _ in ball['asurfs']]))
            # Go through the atom's surfaces
            b_curvs.append(max([self.surfs['curv'][_] for _ in ball['asurfs']]))
            # Calculate the volume of the atom by the previouslty stored volume data
            b_vol = sum([self.surfs['vols'][_][ball['num']] for _ in ball['asurfs']])
            # Exclude atoms that have super large volumes (weird edge error)
            bad_ball = False
            if b_vol > 15 * 4/3 * ball['rad'] ** 3 * pi:
                bad_ball = True
            b_vols.append(b_vol)
            # Check for complete cells in the atoms
            complete = True
            # Go through each of the vertices in the in the atom
            for vert in ball['averts']:
                # Check the number of edges from the vertex that hold
                if len([_ for _ in [self.edges['eatoms'][_] for _ in self.verts['vedges'][vert]] if k in _]) != 3 and not bad_ball:
                    complete = False
            # Additional catch for any atom that doesn't have the 181L number of network elements associated with it
            if len(ball['averts']) < 3 or len(ball['aedges']) < 4 or len(ball['asurfs']) < 3:
                complete = False
            # Add the complete designation for the cell
            b_cell.append(complete)
            # Print the actions
            my_time = time.perf_counter() - self.start_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 "
                  .format(int(h), int(m), round(s, 2), percentage), end="")
        self.balls['vol'], self.balls['sa'], self.balls['curv'], self.balls['complete'] = b_vols, b_sas, b_curvs, b_cell
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
        # Check to see if the only output for the exports is logs
        limit_mem = False
        if self.settings['build_type'] == 'logs':
            limit_mem = True
        # Sort the atoms in the network
        self.sort_balls()
        # Check to see if there are vertices loaded
        if self.group.sys.files['ball_file'] is None:
            # Find the vertices
            self.find_verts()
            # Check to see if there are vertices
            if self.verts is None or len(self.verts) == 0:
                return
        elif self.group.sys.files['ball_file'] != 'deez nuts':
            self.metrics['vert'] = 0
            # Filter out the vertices that don't pertain to the group in question
            verts = []
            for i, vert in self.vta_verts.iterrows():
                if any([True if _ in my_group.group_ndxs else False for _ in vert['vatoms']]):
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
        # Build the network
        self.build_surfaces(not limit_mem)
        # Analyze the network
        self.analyze()
        # else:
            # for surf in self.surfs:
            #     a0, a1 = self.atoms.iloc[surf['atoms'][0]], self.atoms.iloc[surf['atoms'][1]]
            #     surf_atoms_vals = a0['loc'], a0['rad'], a1['loc'], a1['rad']
            #     surf['func'] = calc_surf_func(*surf_atoms_vals)
            # self.metrics['surf'], self.metrics['anal'] = 0, 0
        # Load the elements to the group
        if my_group is not None:
            my_group.get_info()
        # Stop the timer and measure the time
        self.metrics['tot'] = time.perf_counter() - self.start_time
        h, m, s = get_time(self.metrics['tot'])
        num_complete = len([_ for _ in self.balls['complete'] if _])
        print("\rnetwork built - {} complete cell{}, {} verts, {} surfs - {}:{}:{:.2f} s - finished at {}\n"
              .format(num_complete, '' if num_complete == 1 else 's', len(self.verts), len(self.surfs), int(h), int(m),
                      s, datetime.now()), end="")
