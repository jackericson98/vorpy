import os
import csv
import time
import pandas as pd
from datetime import datetime
from time import perf_counter as now
from numpy import array, inf, cbrt, sqrt
from vorpy.src.calculations import get_time
from vorpy.src.calculations import calc_length
from vorpy.src.calculations import global_vars
from vorpy.src.calculations import calc_dist
from vorpy.src.network.analyze import analyze
from vorpy.src.network.build_net import build
from vorpy.src.network.build_edge1 import build_edge
from vorpy.src.network.build_surfs import build_surfs
from vorpy.src.network.mark_doublets import mark_doublets
from vorpy.src.network.find_net_verts import find_net_verts
from vorpy.src.network.net_logs_connect import net_logs_connect


class Network:
    """
    A class representing a network of balls, vertices, edges, and surfaces.
    This class handles the construction and analysis of a network structure, including:
    - Ball sorting and organization
    - Vertex finding and connection
    - Edge building
    - Surface construction
    - Network analysis
    """
    def __init__(self, locs, rads, names=None, group=None, iface_grps=None, group_name=None, settings=None, balls=None,
                 verts=None, edges=None, surfs=None, box=None, sort_balls=False, build_net=False,
                 masses=None):
        """
        Initialize a Network object with the given parameters.
        
        Args:
            locs: List of ball locations
            rads: List of ball radii
            names: Optional list of ball names
            group: Optional group of loc and rad indices for calculation
            settings: Optional settings dictionary
            balls: Optional DataFrame of balls
            verts: Optional DataFrame of vertices
            edges: Optional DataFrame of edges
            surfs: Optional DataFrame of surfaces
            box: Optional box dictionary
            sort_balls: Whether to sort balls on initialization
            build_net: Whether to build network on initialization
            masses: Optional list of ball masses
        """
        # Main network defining objects
        self.group = group               # Group          : List of loc and rad indices for calculation
        self.iface_grps = iface_grps     # Interface Grps : Tuple of lists of ball indices from groups in an interface
        self.group_name = group_name     # Group Name     : Name of the group that the network comes from
        self.settings = settings         # Settings       : surf_res, surf_col, surf_schm, max_vert, net_type
        self.metrics = {'start': now()}  # Metrics        : Holds the time measurements for the build
        self.progress_window = None      # Prog. Window   : Progress window for GUI updates
        self.loaded_from_logs = False    # Logs Flag      : Determines how the net was created

        # Network element lists
        self.balls = balls               # Balls          : Ball DF    - (loc, rad, verts, edges, surfs, vol)
        self.verts = verts               # Vertices       : Vertex DF  - (loc, rad, balls, edges, surfs)
        self.edges = edges               # Edges          : Edge DF    - (center, points, balls, verts, surfs, length)
        self.surfs = surfs               # Surfaces       : Surface DF - (center, points, tris, balls, verts, edges, sa)

        # Tool for splitting up the balls
        self.box = box                    # Box           : Dictionary: Ball box, sub_boxes, and

        if iface_grps is not None:
            if len(iface_grps) != 2:
                raise ValueError(
                    "Interface networks require exactly two interface groups."
                )

            self.iface_grps = tuple(
                set(group_indices)
                for group_indices in iface_grps
            )

            overlap = self.iface_grps[0].intersection(self.iface_grps[1])

            if overlap:
                raise ValueError(
                    f"Interface groups overlap by {len(overlap)} balls: "
                    f"{sorted(overlap)}"
                )

        self.network_mode = (
            "interface"
            if self.iface_grps is not None
            else "complete"
        )

        # Set up the balls
        if names is None:
            names = [str(i) for i in range(len(locs))]
        if self.balls is None:
            self.balls = pd.DataFrame({'loc': locs, 'rad': rads, 'num': [i for i in range(len(locs))], 'name': names,
                                       'mass': masses})
        # Sort the balls if need be
        if sort_balls:
            self.sort_balls()
        # If the user wants to build the net
        if build_net:
            if self.settings is None:
                self.default_settings()

    def default_settings(self, surf_res=0.2, box_size=1.5, max_vert=40, build_type='all', net=None,
                         net_type='aw', surf_col='rainbow', surf_scheme='mean', num_splits=None, print_metrics=False,
                         scheme_factor='log', make_net=False, verts=None):
        """
        Set default settings for the network.
        
        Args:
            surf_res: Surface resolution
            box_size: Box size multiplier
            max_vert: Maximum vertex size
            build_type: Type of build
            net: Network type
            net_type: Network construction type
            surf_col: Surface color scheme
            surf_scheme: Surface scheme
            num_splits: Number of box splits
            print_metrics: Whether to print metrics
            scheme_factor: Scheme factor
            make_net: Whether to make network
            verts: Optional vertices
        """
        self.settings = {'surf_res': surf_res, 'surf_col': surf_col, 'surf_scheme': surf_scheme, 'max_vert': max_vert,
                         'box_size': box_size, 'net_type': net_type, 'build_type': build_type, 'num_splits': num_splits,
                         'print_metrics': print_metrics, 'atom_rad': None, 'scheme_factor': scheme_factor,
                         'foam_box': None, 'sys_dir': os.getcwd()}

    def calc_box(self, locs, rads, return_val=False, box_size=None):
        """
        Calculate the dimensions of a box x times the size of the balls.
        
        Args:
            locs: List of ball locations
            rads: List of ball radii
            return_val: Whether to return the box value
            box_size: Optional box size multiplier
            
        Returns:
            If return_val is True, returns the box vertices
        """
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = array([inf, inf, inf])
        max_vert = array([-inf, -inf, -inf])
        if box_size is None:
            box_size = self.settings['box_size']
        # Loop through each ball in the network
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
        # If the balls are in the same plane adjust the balls
        for i in range(3):
            if r_box[i] == 0 or abs(r_box[i]) == inf:
                r_box[i], min_vert[i], max_vert[i] = 4 * rads[0], locs[0][i], locs[0][i]
        # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
        min_vert, max_vert = max_vert - r_box * box_size, min_vert + r_box * box_size
        # Return the list of array turned list vertices
        box = [[round(_, 3) for _ in min_vert], [round(_, 3) for _ in max_vert]]
        # If the values are to be returned
        if return_val:
            return box
        if self.box is None:
            self.box = {}
        self.box['verts'] = box
        return box

    def set_progress_window(self, progress_window):
        """
        Set the progress window for GUI updates.
        
        Args:
            progress_window: Progress window object
        """
        self.progress_window = progress_window

    def update_progress(self, step, progress):
        """
        Update progress in GUI if available, otherwise print to console.
        
        Args:
            step: Current step name
            progress: Progress percentage
        """
        if self.progress_window:
            self.progress_window.update_progress(step, progress)
        else:
            my_time = now() - self.metrics['start']
            h, m, s = get_time(my_time)
            print(f"\rRun Time = {int(h)}:{int(m):02d}:{s:2.2f} - Process: {step} - {progress:.2f}%", end="")

    def sort_balls(self, num_boxes=None):
        """
        Sort balls into their respective grid sections.
        
        Args:
            num_boxes: Optional number of sub-boxes to divide the network into
        """
        # Print the sorting balls prompt
        print("\rRun Time = 0:00:00.00 - Sorting Balls 0.0 %", end="")
        # Check that the length of the spheres list is big enough to make a vertex
        if len(self.balls) < 4:
            return
        # Set the number of boxes to roughly 5x the number of balls must be a cube for the of cells per row/column/aisle
        elif num_boxes is None:
            n = int(0.5 * sqrt(len(self.balls))) + 1
        else:
            n = int(cbrt(num_boxes)) + 1
        self.settings['num_splits'] = n
        locs, rads = self.balls['loc'], self.balls['rad']
        # First get the box for the balls to be sorted into
        self.calc_box(locs, rads)
        # Instantiate the grid structure of lists is locations representing a grid
        self.box['sub_boxes'] = {(-1, -1, -1): [n]}
        # Get the cell size
        self.box['sub_size'] = [round((self.box['verts'][1][i] - self.box['verts'][0][i]) / n, 3) for i in range(3)]
        my_boxes = []
        # Sort the balls
        for i, loc in enumerate(locs):
            # Print the sorting balls prompt
            percentage = min((i + 1) / len(locs) * 100, 100)
            my_time = now() - self.metrics['start']
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{:02d}:{:2.2f} - Process: Sorting Balls - {:.2f} %"
                  .format(int(h), int(m), round(s, 2), percentage), end="")
            # Find the box they belong to
            box_ndxs = [int((loc[j] - self.box['verts'][0][j]) / self.box['sub_size'][j]) for j in range(3)]

            # Add the ball to the box
            try:
                self.box['sub_boxes'][box_ndxs[0], box_ndxs[1], box_ndxs[2]].append(i)
            except KeyError:
                self.box['sub_boxes'][box_ndxs[0], box_ndxs[1], box_ndxs[2]] = [i]
            # Add the box to the ball
            my_boxes.append(box_ndxs)
        # set the box data
        self.balls['box'] = my_boxes
        # Set the global variables
        global_vars(self.box['sub_boxes'], self.box['verts'], self.settings['num_splits'], max(self.balls['rad']),
                    self.box['sub_size'])

    def find_verts(self):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        """
        find_net_verts(self)

    def connect(self):
        """
        Connects the network using the functions in the build_net.py file
        """
        my_lists = build(self.verts['balls'], self.verts['loc'], self.verts['dub'], len(self.balls), self.metrics['start'],
                         group=self.group, interface=self.iface_grps is not None, iface_grps=self.iface_grps)
        ball_lists, vert_lists, edge_lists, surf_lists = my_lists
        self.balls['verts'], self.balls['edges'], self.balls['surfs'] = ball_lists['verts'], ball_lists['edges'], \
            ball_lists['surfs']
        self.verts['edges'], self.verts['surfs'] = vert_lists['edges'], vert_lists['surfs']
        self.edges = pd.DataFrame(edge_lists)
        self.surfs = pd.DataFrame(surf_lists)
        self.metrics['con'] = now() - self.metrics['start'] - self.metrics['vert']

    def logs_connect(self):
        """
        Connects a network loaded from logs without resolving vertices.
        """
        return net_logs_connect(self)

    def get_real_verts(self):
        my_name = os.getcwd() + '/Data/user_data/' + self.group.sys.name + '_Correct/sys/' + self.group.sys.name + \
                  '_logs.csv'
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
            percentage = min(i / len(self.edges) * 100, 100)
            my_time = time.perf_counter() - self.metrics['start']
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{:02d}:{:2.2f} - Process: building edges: edge {} - {:.2f} %"
                  .format(int(h), int(m), round(s, 2), i, percentage), end="")
            vlocs = [array(self.verts['loc'][_]) for _ in edge['verts']]
            vdist = calc_dist(vlocs[0], vlocs[1])

            if vdist < 1e-8:
                print("\n=== ZERO LENGTH EDGE INPUT ===")
                print(f"edge index = {i}")
                print(f"edge balls = {edge['balls']}")
                print(f"edge verts = {edge['verts']}")
                print(f"vlocs[0] = {vlocs[0]}")
                print(f"vlocs[1] = {vlocs[1]}")
                print(self.verts.iloc[edge['verts']])
                print("==============================\n")
            # Build the edge depending on if it is straight or not
            try:
                edge_points, edge_vals = build_edge(
                    locs=[array(self.balls['loc'][_]) for _ in edge['balls']],
                    rads=[self.balls['rad'][_] for _ in edge['balls']],
                    vlocs=vlocs,
                    blocs=self.balls['loc'],
                    brads=self.balls['rad'],
                    eballs=edge['balls'],
                    res=self.settings['surf_res'],
                    straight=self.settings['net_type'] in {'prm', 'pow'},
                    edub=any([self.verts['dub'][_] in {1, 2} for _ in edge['verts']]),
                    edge_verts=self.verts.iloc[edge['verts']]
                )
            except ValueError:
                print(vlocs, edge['balls'])

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
        analyze(self)

    def build(self, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None,
              calc_verts=None, my_group=None, print_actions=None, print_vert_metrics=False, curr_time=None, verts=None):
        """
        Builds and constructs the complete network structure including vertices, edges, and surfaces.

        This method orchestrates the entire network construction process by:
        1. Sorting and organizing the balls in the network
        2. Finding and verifying vertices
        3. Connecting vertices to form edges
        4. Building surfaces between edges
        5. Analyzing the final network structure

        Parameters
        ----------
        print_actions : bool, optional
            If True, prints detailed progress of network construction steps
        net_type : str, optional
            Specifies the network construction algorithm:
            - 'curv': Curved network
            - 'del': Delaunay network
            - 'pow': Power network
        my_group : list, optional
            List of ball indices to include in network construction
        surf_res : float, optional
            Resolution parameter for surface construction
        max_vert : float, optional
            Maximum allowed vertex size in the network
        box_size : float, optional
            Multiplier for the bounding box size
        build_surfs : bool, optional
            If True, constructs surfaces in the network
        calc_verts : bool, optional
            If False, skips vertex calculations when loading existing network
        """
        # Check for the group name
        if self.group_name is None:
            self.group_name = 'Group'
        # Check to see if the only output for the exports is logs
        limit_mem = False
        if self.settings['build_type'] == 'logs':
            limit_mem = True
        # Sort the balls in the network
        if self.box is None:
            self.sort_balls()
        if verts is not None:
            self.verts = verts
        # Check to see if there are vertices loaded
        if self.verts is None:
            # Find the vertices
            self.find_verts()
            # Check to see if there are vertices
            if self.verts is None or len(self.verts) == 0:
                return
        elif 'vdub' not in self.verts:
            self.metrics['vert'] = 0
            self.verts['dub'] = mark_doublets(self.verts)
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

        # Stop the timer and measure the time
        self.metrics['tot'] = now() - self.metrics['start']
        h, m, s = get_time(self.metrics['tot'])
        num_complete = len([_ for _ in self.balls['complete'] if _])
        print("\r\"{}\" network built - {} complete cell{}, {} verts, {} surfs - {}:{}:{:.2f} s - finished at {}\n"
              .format(self.group_name, num_complete, '' if num_complete == 1 else 's', len(self.verts),
                      len(self.surfs), int(h), int(m), s, datetime.now()), end="")
