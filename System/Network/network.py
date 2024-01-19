import time
from _datetime import datetime
import pandas as pd
import csv
import os
from System.Network.verts.find_verts import find_verts
from System.Network.build_net import build, get_time, calc_length
from System.Network.edges.build_edge import build_edge
from System.Network.surfs.build_surf import build_surf
from System.sys_funcs.calcs.calcs import calc_vol, ndx_search, global_vars
from System.sys_funcs.calcs.surf import calc_surf_func, calc_surf_sa, calc_surf_tri_curvs
from numpy import array, inf, cbrt, sqrt, pi


class Network:
    """Network object. Graph that holds the elements of the Voronoi S-Network."""
    def __init__(self, sys, atoms=None, verts=None, edges=None, surfs=None, surf_res=0.2, box_size=1.1, max_vert=40,
                 calc_verts=True, connect_net=True, build_surfs=True, net_type='vor', surf_col='plasma',
                 surf_scheme='curv', sub_net=False, box=None, sub_boxes=None):

        # Main network defining objects
        self.num_splits = None
        self.sys = sys                    # System            :   Route back to outer system
        self.type = net_type              # Network Type      :   String indicating network build type
        self.id = 0                       # Network id #      :
        self.sub_net = sub_net

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
        self.box = box                     # Box              :    Holds a max and min vertex for the retaining box
        self.sub_boxes = None              # Sub boxes        :    3D array holding atoms relative locations
        self.vert_sub_boxes = None         # Vert Sub Boxes   :    Holds the vertices of the network by their location
        self.sub_box_size = None           # Sub box size     :    Holds the size of each sub box
        self.box_max = None                # Box maxes        :    number of x, y, z boxes or rows, columns, aisles
        self.atoms_box = []                # Atoms box        :    min/max vals for the box containing the atoms

        # Diagnostic variables
        self.start_time = time.perf_counter()
        self.my_time = None                # My time          :    Time taken to calculate the network
        self.metrics = {}                  # Build Metrics    :    Holds the time measurements for the build
        self.max_vert_rad = 0              # Max Vertex Rad   :    Maximum real vertex recorded
        self.max_curv = 0

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

    def calc_box(self, locs, rads, return_val=False, box_size=None):
        """
        Determines the dimensions of a box x times the size of the atoms
        :return: Sets the box attribute with the correct values as well as atoms_box
        """
        # Set up the minimum and maximum x, y, z coordinates
        min_vert = array([inf, inf, inf])
        max_vert = array([-inf, -inf, -inf])
        if box_size is None:
            box_size = self.box_size
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
        :return: Sets the values for self.sub_boxes with the atom objects in their correct locations. Also sets the
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

    def sort_verts(self, my_group=None):
        """
        Puts the vertices in the network in their respective grid sections
        :return: Places the vertices into their correct sub_boxes
        """
        # Check to see if a group is provided
        if my_group is not None:
            atom_ndxs = [_['num'] for _ in my_group.atoms]
        else:
            atom_ndxs = [_['num'] for _ in self.atoms]
        atom_ndxs.sort()
        # Instantiate the grid structure of lists is locations representing a grid
        self.vert_sub_boxes = [[[[] for _ in range(self.box_max[2] + 1)]
                               for _ in range(self.box_max[1] + 1)] for _ in range(self.box_max[0] + 1)]

        # Sort the atoms
        drop_verts = []
        for i, vert in enumerate(self.verts):
            # Check that the vertex has at least one atom of interest
            for ndx in vert['ndx']:
                search_ndx = ndx_search(atom_ndxs, ndx)
                if atom_ndxs[search_ndx] == ndx:
                    break
            else:
                drop_verts.append(i)
            # Adjust the maximum radius
            if vert['rad'] > self.max_vert_rad:
                self.max_vert_rad = vert['rad']
            # Find the box they belong to
            box_ndxs = [int((vert['loc'][i] - self.box[0][i]) / self.sub_box_size[i]) for i in range(3)]
            # Add the atom to the box
            try:
                self.vert_sub_boxes[box_ndxs[0]][box_ndxs[1]][box_ndxs[2]].append(vert)
            except IndexError:
                drop_verts.append(i)
            # Add the box to the atom
            vert['box'] = box_ndxs
        # Drop the vertices outside the box
        for vert_num in reversed(drop_verts):
            self.verts.pop(vert_num)
            self.vert_ndxs.pop(vert_num)

    def get_verts(self, cells, reach=0):
        """
        Takes in the cells and the number of additional cells to search and returns an atom list
        :param cells: The initial boxes in the network to stem from
        :param reach: The number of cells out from the initial set of cells to search
        """
        # If a single cell is entered
        if type(cells[0]) is int:
            cells = [cells]
        # Get the min and max of the cells
        ndx_min = [inf, inf, inf]
        ndx_max = [-inf, -inf, -inf]
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

        return verts

    def connect(self):
        """
        Connects the network using the functions in the build_net.py file
        """
        my_lists = build(self.verts['vatoms'], self.verts['vloc'], self.verts['vdub'], len(self.atoms), self.my_time)
        atom_lists, vert_lists, edge_lists, surf_lists = my_lists
        self.atoms['averts'], self.atoms['aedges'], self.atoms['asurfs'] = atom_lists['averts'], atom_lists['aedges'], \
            atom_lists['asurfs']
        self.verts['vedges'], self.verts['vsurfs'] = vert_lists['vedges'], vert_lists['vsurfs']
        self.edges = pd.DataFrame(edge_lists)
        self.surfs = pd.DataFrame(surf_lists)
        self.metrics['con'] = time.perf_counter() - self.my_time - self.metrics['vert']

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

    def find_verts(self, my_group=None, print_metrics=False):
        """
        Using the functions in find_vertices.py finds the vertices in the network
        """
        global_vars(self.sub_boxes, self.box, self.num_splits, self.sys.max_atom_rad, self.sub_box_size)
        # Check to see if a group has been provided
        if my_group is not None:
            atom_nums = my_group.atom_ndxs[:]
        else:
            atom_nums = [i for i in range(len(self.atoms))]
        vert_list_real = self.get_real_verts()
        # Get the indices of the atoms in the network to keep track of the atoms that haven't been visited
        self.atom_ndxs = [_ for _ in atom_nums]
        my_group_atom_ndxs = None
        if my_group is not None:
            my_group_atom_ndxs = my_group.atom_ndxs
        my_guuy = find_verts(alocs=self.atoms['loc'].to_numpy(), arads=self.atoms['rad'].to_numpy(),
                             max_vert=self.max_vert, net_type=self.type, check_atoms=atom_nums,
                             my_group=my_group_atom_ndxs, start_time=self.my_time, print_metrics=print_metrics,
                             vert_box=self.sys.foam_box)
        if my_guuy is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = my_guuy
        # Check for disconnects in the network
        while len(atom_nums) > 0:
            a0 = atom_nums.pop()
            my_guuy = find_verts(a0=a0, alocs=self.atoms['loc'].to_numpy(), arads=self.atoms['rad'].to_numpy(),
                                 max_vert=self.max_vert, net_type=self.type, check_atoms=atom_nums,
                                 my_group=my_group.atom_ndxs, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                                 vloc2s=vloc2s, vrad2s=vrad2s, start_time=self.my_time, print_metrics=print_metrics,
                                 vert_box=self.sys.foam_box, averts=averts)
            if my_guuy is not None:
                vert_ndxs, vlocs, vrads, vloc2s, vrad2s, atom_nums, averts = my_guuy
            if self.sys.foam_box is not None and len(atom_nums) <= 0.25*len(self.atoms['loc']):
                break
        # Create the doublets list
        if vert_list_real is not None and self.type == 'vor':
            missing_verts = [_ for _ in vert_list_real if _ not in vert_ndxs]
            print(missing_verts)
            extra_verts = [_ for _ in vert_ndxs if _ not in vert_list_real]
            print(extra_verts)
        doublets = [0 for _ in range(len(vert_ndxs))]
        # Incorporate the doublets into the vlocs, vatoms, vrads lists and lose the vloc2s and vrad2s
        i = 0
        while i < len(vlocs):
            # Check for doubletness
            if vrad2s[i] is not None:
                # Insert the relevant information into their respective lists
                vert_ndxs.insert(i + 1, vert_ndxs[i])
                vlocs.insert(i + 1, vloc2s[i])
                vrads.insert(i + 1, vrad2s[i])
                doublets.insert(i + 1, 1)
                # Preserve the relational aspects of vrad2s and vloc2s
                vrad2s.insert(i + 1, None)
                vloc2s.insert(i + 1, [None, None, None])
            i += 1

        # Make the dataframe
        self.verts = pd.DataFrame({"vatoms": vert_ndxs, 'vloc': vlocs, 'vrad': vrads, 'vdub': doublets})
        # Clear the print statement
        if self.sys.print_actions:
            print("\r                                                                  ", end="")
        self.metrics['vert'] = time.perf_counter() - self.my_time

    def build_edges(self):
        """
        Builds the edges in the network for use in the surfaces
        """
        # Set the edge points and vals lists
        edges_points, edges_vals = [], []
        # Go through the edges in the network
        for i, edge in self.edges.iterrows():
            # Build the edge depending on if it is straight or not
            straight = True if self.type in ['pow', 'flat', 'del'] else False
            edge_points, edge_vals = build_edge(alocs=[array(self.atoms['loc'][_]) for _ in edge['eatoms']],
                                                arads=[self.atoms['rad'][_] for _ in edge['eatoms']],
                                                vlocs=[array(self.verts['vloc'][_]) for _ in edge['everts']],
                                                res=self.surf_res, straight=straight)
            # Add them to the lists
            edges_points.append(edge_points)
            edges_vals.append(edge_vals)
        # Set the dataframe values
        self.edges['points'], self.edges['vals'] = edges_points, edges_vals

    def build_surfaces(self):
        """
        Takes in a system and returns a fully connected network
        """
        # Make each surface
        points, tris, tri_curvs, curvs, funcs, coms, flats = [], [], [], [], [], [], []
        for i, surf in self.surfs.iterrows():
            # Build the surfaces and print the progress
            my_time = time.perf_counter() - self.my_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {:2}:{:2}:{:.2f} - Process: building surfaces {:.2f} %                                 "
                  .format(int(h), int(m), round(s, 2), min(100.0, 100 * round(i/len(self.surfs), 4))), end="")
            arads = [self.atoms['rad'][_] for _ in surf['satoms']]
            alocs = [self.atoms['loc'][_] for _ in surf['satoms']]
            if arads[0] > arads[1]:
                arads, alocs = [arads[1], arads[0]], [alocs[1], alocs[0]]
            my_surf = build_surf(alocs=alocs, arads=arads, epnts=[self.edges['points'][_] for _ in surf['sedges']],
                                 res=self.surf_res, net_type=self.type)
            surf_points, surf_tris, surf_tri_curvs, surf_curv, surf_func, surf_com, surf_flat = my_surf
            points.append(surf_points)
            tris.append(surf_tris)
            tri_curvs.append(surf_tri_curvs)
            curvs.append(surf_curv)
            funcs.append(surf_func)
            coms.append(surf_com)
            flats.append(surf_flat)
        # Set the dataframe elements
        self.surfs['points'], self.surfs['tris'], self.surfs['tri_curvs'], self.surfs['curv'], self.surfs['func'], \
            self.surfs['com'], self.surfs['flat'] = points, tris, tri_curvs, curvs, funcs, coms, flats
        if self.sys.print_actions:
            print("\r                                                                                             ", end='')
        self.metrics['surf'] = time.perf_counter() - self.my_time - self.metrics['vert'] - self.metrics['con']

    def analyze(self):
        """
        Analyzes the output surfaces, cells and solute vertices for the network for later reference
        """
        # Check to see if my_time has started
        if self.my_time is None:
            self.my_time = time.perf_counter()
        # Get the percentage total number
        tot_num = len(self.edges) + len(self.surfs) + len(self.atoms)
        # Go through the edges in the network
        i = 0
        lengths = []
        for i, edge in self.edges.iterrows():
            percentage = int(i / tot_num * 100)
            # Calculate the length of each edge
            length = calc_length(array(edge['points']))
            lengths.append(length)
            if self.sys.print_actions:
                my_time = time.perf_counter() - self.my_time
                h, m, s = get_time(my_time)
                print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                  "
                      .format(int(h), int(m), round(s, 2), percentage), end="")
        self.edges['length'] = lengths
        # Go through each surface in the system and find the simplices and the surface area
        sas = []
        surfs_tri_curvs, surfs_curvs = [], []
        j = 0
        for j, surf in self.surfs.iterrows():
            percentage = int((i + j + 1) / tot_num * 100)
            sas.append(calc_surf_sa(edges=[self.edges['points'][_] for _ in surf['sedges']], com=array(surf['com']),
                                    tris=surf['tris'], points=surf['points'], flat=surf['flat']))
            # Get the curvature of the surface patch
            if surf['curv'] is None or (surf['curv'] == 0 and not surf['flat']):
                surf_tri_curvs, scurvs = calc_surf_tri_curvs(surf['func'], surf['points'], surf['tris'], surf['curv'])
                surfs_tri_curvs.append(surf_tri_curvs)
                surfs_curvs.append(scurvs)
            else:
                surfs_tri_curvs.append(surf['tri_curvs'])
                surfs_curvs.append(surf['curv'])
            if self.sys.print_actions:
                my_time = time.perf_counter() - self.my_time
                h, m, s = get_time(my_time)
                print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                  "
                      .format(int(h), int(m), round(s, 2), percentage), end="")
        # Get the curvature in the 95th percentile
        my_surf_curvs = surfs_curvs.copy()
        my_surf_curvs.sort()
        try:
            self.max_curv = my_surf_curvs[min(int(0.99 * len(my_surf_curvs)), len(my_surf_curvs) - 1)]
        except IndexError:
            self.max_curv = 0
        # Assign the values
        self.surfs['sa'], self.surfs['curv'], self.surfs['tri_curvs'] = sas, surfs_curvs, surfs_tri_curvs
        # Set up the atoms' volumes surface areas, curvatures vars
        avols, asas, acurvs, acell = [], [], [], []
        asurfs_vols = [[] for _ in range(len(self.surfs))]
        # Go through each atom in the system and find the volume
        for k, atom in self.atoms.iterrows():
            # Get the percentage for printing
            percentage = int((i + j + k + 2) / tot_num * 100)
            avol, asurf_vols = calc_vol(atom['loc'], [self.surfs['points'][_] for _ in atom['asurfs']],
                                        [self.surfs['tris'][_] for _ in atom['asurfs']])
            bad_atom = False
            if avol > 15 * 4/3 * atom['rad'] ** 3 * pi:
                bad_atom = True
            avols.append(avol)
            for i, surf_vol in enumerate(asurf_vols):
                asurfs_vols[atom['asurfs'][i]].append([k, surf_vol])
            # Get/calculate the surface area
            if 'sa' not in atom or atom['sa'] is None or atom['sa'] == 0:
                asas.append(sum([self.surfs['sa'][_] for _ in atom['asurfs']]))
            else:
                asas.append(atom['sa'])
            # Check that the curvature is None
            if 'curv' not in atom or atom['curv'] is None or atom['curv'] == 0:
                atom_curv = 0
                # Go through the atom's surfaces
                for m in atom['asurfs']:
                    surf = self.surfs.iloc[m]
                    if surf['curv'] > atom_curv:
                        atom_curv = surf['curv']
                acurvs.append(atom_curv)
            # Check for complete cells in the atoms
            complete = True
            # Go through each of the vertices in the in the atom
            for vert in atom['averts']:
                # Check the number of edges from the vertex that hold
                if len([_ for _ in [self.edges['eatoms'][_] for _ in self.verts['vedges'][vert]] if k in _]) != 3 and not bad_atom:
                    complete = False
            # Additional catch for any atom that doesn't have the correct number of network elements associated with it
            if len(atom['averts']) < 3 or len(atom['aedges']) < 4 or len(atom['asurfs']) < 3:
                complete = False
            # Add the complete designation for the cell
            acell.append(complete)
            # Print the actions
            if self.sys.print_actions:
                my_time = time.perf_counter() - self.my_time
                h, m, s = get_time(my_time)
                print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 "
                      .format(int(h), int(m), round(s, 2), percentage), end="")
        self.surfs['vols'] = [[asurfs_vols[i][0][1], asurfs_vols[i][1][1]] if asurfs_vols[i][0][0]<asurfs_vols[i][1][0]
                              else [asurfs_vols[i][1][1], asurfs_vols[i][0][1]] for i in range(len(self.surfs))]
        self.atoms['vol'], self.atoms['sa'], self.atoms['curv'], self.atoms['complete'] = avols, asas, acurvs, acell
        self.metrics['anal'] = time.perf_counter() - self.my_time - self.metrics['surf'] - self.metrics['con'] - self.metrics['vert']

    def build(self, surf_res=None, max_vert=None, box_size=None, build_surfs=None, net_type=None,
              calc_verts=None, my_group=None, print_actions=None, print_vert_metrics=False):
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
            self.find_verts(my_group=my_group, print_metrics=print_vert_metrics)
            # Check to see if there are vertices
            if self.verts is None or len(self.verts) == 0:
                return
        else:
            self.metrics['vert'] = 0
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
