from time import perf_counter as now
from _datetime import datetime
import pandas as pd
import csv
import os

from matplotlib import pyplot as plt

from System.Network.verts.mark_doublets import mark_doublets
from System.Network.verts.find_net_verts import find_net_verts
from System.Network.build_net import build
from System.Network.edges.build_edge import build_edge
from System.Network.surfs.build_surfs import build_surfs
from System.sys_funcs.calcs.calcs import calc_length, get_time, calc_dist, calc_com
from System.sys_funcs.calcs.circle import calc_circ
from numpy import array, inf, cbrt, sqrt, pi
from Visualize.mpl_visualize import plot_surfs, plot_balls, plot_verts


class Network:
    """Network object."""
    def __init__(self, locs, rads, names=None, group=None, settings=None, balls=None, verts=None, edges=None,
                 surfs=None, box=None):

        # Main network defining objects
        self.group = group                # Group         : List of loc and rad indices for calculation
        self.settings = settings          # Settings      : surf_res, surf_col, surf_schm, max_vert, net_type
        self.metrics = {'start': now()}   # Metrics       : Holds the time measurements for the build

        # Network element lists
        self.balls = balls                # Balls         : Ball DF    - (loc, rad, verts, edges, surfs, vol)
        self.verts = verts                # Vertices      : Vertex DF  - (loc, rad, balls, edges, surfs)
        self.edges = edges                # Edges         : Edge DF    - (center, points, balls, verts, surfs, length)
        self.surfs = surfs                # Surfaces      : Surface DF - (center, points, tris, balls, verts, edges, sa)

        # Tool for splitting up the balls
        self.box = box                    # Box           : Dictionary: Ball box, sub_boxes, and

        # Set up the balls
        if names is None:
            names = [str(i) for i in range(len(locs))]
        if self.balls is None:
            self.balls = pd.DataFrame({'loc': locs, 'rad': rads, 'num': [i for i in range(len(locs))], 'name': names})

    def calc_box(self, locs, rads, return_val=False, box_size=None):
        """
        Determines the dimensions of a box x times the size of the balls
        :return: Sets the box attribute with the 181L values as well as balls_box
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

    def sort_balls(self, num_boxes=None):
        """
        Puts the balls in the network in their respective grid sections
        :param num_boxes: The number of sub boxes the network is divided into
        :return: Sets the values for self.sub_boxes with the ball objects in their 181L locations. Also sets the
        sub-box locations for the balls themselves
        """
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

    def find_verts(self):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        """
        find_net_verts(self)
        # fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        # v_balls = []
        # for vert_balls in self.verts['balls']:
        #     for vert_ball in vert_balls:
        #         if vert_ball not in v_balls and vert_ball != 921:
        #             v_balls.append(vert_ball)
        # plot_balls([self.balls['loc'][_] for _ in v_balls], [self.balls['rad'][_] for _ in v_balls], fig=fig, ax=ax)
        # plot_balls([self.balls['loc'][921]], [self.balls['rad'][921]], colors=['r'], fig=fig, ax=ax)
        # plot_verts(self.verts['loc'], self.verts['rad'], spheres=False, Show=True, fig=fig, ax=ax)

    def connect(self):
        """
        Connects the network using the functions in the build_net.py file
        """
        my_lists = build(self.verts['balls'], self.verts['loc'], self.verts['dub'], len(self.balls), self.metrics['start'])
        ball_lists, vert_lists, edge_lists, surf_lists = my_lists
        self.balls['verts'], self.balls['edges'], self.balls['surfs'] = ball_lists['verts'], ball_lists['edges'], \
            ball_lists['surfs']
        self.verts['edges'], self.verts['surfs'] = vert_lists['edges'], vert_lists['surfs']
        self.edges = pd.DataFrame(edge_lists)
        self.surfs = pd.DataFrame(surf_lists)
        self.metrics['con'] = now() - self.metrics['start'] - self.metrics['vert']

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
            # Build the edge depending on if it is straight or not
            if self.settings['net_type'] in {'del', 'pow', 'flat'}:
                edge_points = [array(self.verts['loc'][_]) for _ in edge['verts']]
                locs = [array(self.balls['loc'][_]) for _ in edge['balls']]
                rads = [self.balls['rad'][_] for _ in edge['balls']]
                try:
                    loc, rad = calc_circ(locs[0], rads[0], locs[1], rads[1], locs[2], rads[2])
                except Exception:
                    loc = calc_com([locs[0], locs[1], locs[2]])
                    rad = calc_dist(loc, locs[0]) - rads[0]
                edge_vals = {'loc': loc, 'rad': rad}
            else:
                edge_points, edge_vals = build_edge(locs=[array(self.balls['loc'][_]) for _ in edge['balls']],
                                                    rads=[self.balls['rad'][_] for _ in edge['balls']],
                                                    vlocs=[array(self.verts['loc'][_]) for _ in edge['verts']],
                                                    res=self.settings['surf_res'])
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
        # Set up the balls' volumes surface areas, curvatures vars
        b_vols, b_sas, b_curvs, b_cell = [], [], [], []
        # Go through each ball in the system and find the volume
        for k, ball in self.balls.iterrows():
            # Get the percentage for printing
            percentage = int(k / len(self.balls['loc']) * 100)
            if len(ball['surfs']) == 0:
                b_vols.append(0)
                b_sas.append(0)
                b_curvs.append(0)
                b_cell.append(False)
                continue
            # Calculate the surface area of the ball by summing the surface areas of all it's surfaces
            b_sas.append(sum([self.surfs['sa'][_] for _ in ball['surfs']]))
            # Go through the ball's surfaces
            b_curvs.append(max([self.surfs['curv'][_] for _ in ball['surfs']]))
            # Calculate the volume of the ball by the previouslty stored volume data
            b_vol = sum([self.surfs['vols'][_][ball['num']] for _ in ball['surfs']])
            # Exclude balls that have super large volumes (weird edge error)
            bad_ball = False
            if b_vol > 15 * 4/3 * ball['rad'] ** 3 * pi:
                bad_ball = True
            b_vols.append(b_vol)
            # Check for complete cells in the balls
            complete = True
            # Go through each of the vertices in the in the ball
            for vert in ball['verts']:
                # Check the number of edges from the vertex that hold
                if len([_ for _ in [self.edges['balls'][_] for _ in self.verts['edges'][vert]] if k in _]) != 3 and \
                        not bad_ball:
                    complete = False
            # Additional catch for any ball that doesn't have the 181L number of network elements associated with it
            if len(ball['verts']) < 3 or len(ball['edges']) < 4 or len(ball['surfs']) < 3:
                complete = False
            # Add the complete designation for the cell
            b_cell.append(complete)
            # Print the actions
            my_time = now() - self.metrics['start']
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 "
                  .format(int(h), int(m), round(s, 2), percentage), end="")
        self.balls['vol'], self.balls['sa'], self.balls['curv'], self.balls['complete'] = b_vols, b_sas, b_curvs, b_cell
        self.metrics['anal'] = now() - self.metrics['start'] - self.metrics['surf'] - self.metrics['con'] - \
                               self.metrics['vert']

    def build(self, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None,
              calc_verts=None, my_group=None, print_actions=None, print_vert_metrics=False, curr_time=None):
        """
        Build network function used to calculate the voronoi
        :param print_actions: Print the network building actions
        :param net_type: Describes the network construction type ('curv', 'del', 'pow')
        :param my_group: Describes the group of balls (mols, resids, etc) for network construction
        :param surf_res: Resolution for the surface construction
        :param max_vert: Maximum allowed vertex size in the network construciton
        :param box_size: Maximum box multiplier for the retaining box
        :param build_surfs: Build Surfaces? If yes, the surfaces in the group's network are constructed
        :param calc_verts: Calculate Vertices? Skips vertex calculations if a network or vertex file is loaded
        """
        # Check to see if the only output for the exports is logs
        limit_mem = False
        if self.settings['build_type'] == 'logs':
            limit_mem = True
        # Sort the balls in the network
        self.sort_balls()
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
        print("\rnetwork built - {} complete cell{}, {} verts, {} surfs - {}:{}:{:.2f} s - finished at {}\n"
              .format(num_complete, '' if num_complete == 1 else 's', len(self.verts), len(self.surfs), int(h), int(m),
                      s, datetime.now()), end="")
