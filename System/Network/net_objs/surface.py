import os.path
import matplotlib as mpl
from System.sys_funcs.calcs import calc_tri, calc_dist
from System.Network.net_funcs.build_surf import *
import csv


class Surface:
    """Surface object. Holds the mesh data. Used to analyze interfaces between atoms."""
    def __init__(self, atoms=None, net=None, edges=None, verts=None, doublet=False, points=None, tris=None,
                 perimeter=None, normal=None, sa=0, curvature=None, function=None, file=None, resolution=None,
                 center=None, tri_colors=None, color_scheme=None, color_map=None, ndx=None):

        # Network objects
        self.net = net                  # Network         : Network of the System
        self.atoms = atoms              # Atoms           : List Atom objects for the surface
        self.verts = verts              # Vertices        : List of Vertex objects for surface
        self.edges = edges              # Edges           : List of Edge objects for the surface

        # Main descriptor attributes
        self.ndx = ndx                  # Index           : Indices of the atoms of the surface
        self.file = file                # File            : File address for the reference file holding points and tris
        self.func = function            # Function        : Holds the coefficients of the function describing the surf
        self.res = resolution           # Resolution      : The resolution with which to build the surface

        # Points and tris
        self.points = points            # Points          : The points that make up the surface
        self.flat_points = []           # Flat points     : Points projected into 2d based off of the surface normal
        self.perimeter = perimeter      # Perimeter       : The points around the edges of the surface (IN ORDER)
        self.pflat_points = []          # Flat perimeter  : Flattened points around the perimeter
        self.tris = tris                # Triangles       : A list of connections between the points

        # Coloring values
        self.tri_dists = None           # Distances       : List corresponding to the distance from the center
        self.tri_ins_out = None         # Inside Outside  : List corresponding to inside or outside the atoms
        self.tri_curvs = None           # Curvature       : List corresponding to curvature value for each point

        # Current coloring
        self.tri_dists = None
        self.tri_curvs = None
        self.tri_ins_out = None
        self.tri_colors = tri_colors    # Tri colors      : Holds the color mapped color for each triangle
        self.scheme = color_scheme      # Color Scheme    : Holds the method by which the color map is mapped
        self.color_map = color_map      # Color Map       : Holds the map applied to the triangles on the surface

        # Calculation attributes
        self.sa = sa                    # Surface Area    : The surface area of the
        self.curv = curvature           # Curvature       : The curvature of the surface between the
        self.norm = normal              # Surface Normal  : Normal to the center of the surface
        self.loc = center               # Location        : Center point of the hyperboloid the surface is made from
        self.com = None                 # Center of mass  : The point toward which all building paths travel
        self.doublet = doublet          # Doublet         : Indicates whether a surface is a part of a doublet or not
        self.flat = False               # Is Flat?        : Whether the surface is flat or not

        # Make sure that a0 is the atom with the smaller radius
        if self.atoms is not None:
            if self.atoms[0].rad > self.atoms[1].rad:
                self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]

        # Create and sort the atom's indices
        if atoms is not None:
            self.ndx = [atom.num for atom in atoms]
            self.ndx.sort()

    def calc_func(self):
        """
        Calculates the coefficients for the surface between the two atoms
        :return: The surface has the correct self.func attribute
        """
        # Make sure that a0 is the atom with the smaller radius
        if self.atoms[0].rad > self.atoms[1].rad:
            self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
        # Create a0, a1 variables
        a0, a1 = self.atoms
        # Set the rn vector for the surface since the atoms are sorted
        l0, l1 = np.array(a0.loc), np.array(a1.loc)
        r = l1 - l0
        self.norm = r / np.linalg.norm(r)
        # Grab the centers of the spheres
        x1, y1, z1 = l0
        x2, y2, z2 = l1
        # Calculate the major coefficients (pg. 574 Z. Hu)
        R = a0.rad - a1.rad
        K = (x2 ** 2 - x1 ** 2) + (y2 ** 2 - y1 ** 2) + (z2 ** 2 - z1 ** 2) - R ** 2
        d = x1 - x2, y1 - y2, z1 - z2
        J = 4 * R ** 2 * (x1 ** 2 + y1 ** 2 + z1 ** 2) - K ** 2
        # Instantiate/reset the hyperboloid coefficient vector lists
        ABC, DEF, GHI = [], [], []
        # Calculate hyperboloid coefficients
        for i in range(3):
            ABC.append(4 * R ** 2 - 4 * d[i] ** 2)
            DEF.append(-8 * d[i] * d[(i + 1) % 3])  # The equation asks for D_y, D_z, D_x in that order, hence modulus
            GHI.append(-8 * R ** 2 * l0[i] - 4 * K * d[i])
        # Set the function attribute
        self.func = ABC + DEF + GHI + [J] + [K] + list(d)

    def read_file(self, file=None):
        """
        Reads the file holding the points and the triangles for the surface
        :param file: Specifies the address for the build file
        :return: The surfaces points and triangles are set
        """
        # Check to see if the file exists
        if file is None and self.file is not None:
            file = self.file
        # Check that the provided file works as an address on its own
        if os.path.exists(file):
            file_address = file
        # Check that the file name is a relative location to the system directory
        elif os.path.exists(self.net.sys.dir + file):
            file_address = self.net.sys.dir + file
        # Last brute force a location if the file name is incorrect
        else:
            return
        # Read an off file
        if file_address[-3:].lower() == 'off':
            # Open the file
            with open(file_address, 'r') as my_file:
                # Read the lines
                file_array = my_file.readlines()
                # Get the number of points and triangles
                num_points, num_tris = [int(_) for _ in file_array[1].split()[:2]]
                # Add the points
                self.points = []
                for i in range(4, num_points + 4):
                    line = file_array[i].split()
                    self.points.append([float(_) for _ in line])
                # Add the tris
                self.tris = []
                for i in range(4 + num_points, 4 + num_points + num_tris):
                    line = file_array[i].split()
                    self.tris.append([int(_) for _ in line[1:4]])
        # Read a comma separated file surface file
        elif file_address[-3:].lower() == 'csv':
            # Open the file
            with open(file_address, 'r') as my_file:
                # Get the file element array to read
                read_file = list(csv.reader(my_file, delimiter=","))
                # Get the number of points and triangles
                num_points, num_tris = [int(_) for _ in read_file[1][1:]]
                # Go through the points lines of the file
                self.points = []
                for i in range(3, num_points + 3):
                    self.points.append([float(_) for _ in read_file[i][1:]])
                # Go through the triangles lines of the file
                self.tris = []
                for i in range(4 + num_points, 4 + num_points + num_tris):
                    self.tris.append([int(_) for _ in read_file[i][1:]])

    # Calculate curvature method
    def calc_curv(self):
        """
        Calculates the curvature of the surface
        :return: The curvature attribute is filled
        """
        # Check to see that the function has been calculated or not
        if self.func is None:
            self.calc_func()
        # Made up function to calculate the general curvature of the hyperboloid
        self.curv = np.sqrt(self.func[0]**2 + self.func[1]**2 + self.func[2]**2)

    def calc_sa(self):
        """
        Calculates the surface area of the input surface
        :return: Surface area of the surface
        """
        # Create the surface area variable
        sa = 0
        if self.flat:
            for edge in self.edges:
                if edge.straight:
                    sa += calc_tri([edge.pv0, edge.pv1, self.com])
                else:
                    for i in range(len(edge.points) - 1):
                        p0, p1 = edge.points[i:i + 2]
                        sa += calc_tri([p0, p1, self.com])
        # Go through the triangles in the surface
        for tri in self.tris:
            p0, p1, p2 = self.points[tri[0]], self.points[tri[1]], self.points[tri[2]]
            sa += calc_tri([p0, p1, p2])
        self.sa = sa

    # Build method. Makes the mesh for the surface and calculates the simplices between them
    def build(self, res=None, color=False):
        """
        Main build method for constructing surfaces
        :param res: Specifies the resolution the surface is to be constructed with
        :param color: Bool to color the surf or not
        :return: The surfaces points and triangles are filled
        """
        if self.net.type in {'pow', 'del'} or self.atoms[0].element == self.atoms[1].element:
            self.flat = True
        # Check to see if the file exists
        if self.file is not None:
            self.read_file()
            if self.points is not None and len(self.points) > 1:
                return
        # Set the resolution value that the surface is built with
        if res is None:
            res = self.net.surf_res
        self.res = res
        # Check to see if the function or curvature have been calculated and calculate them if not
        if self.curv is None:
            self.calc_curv()
        # Reset the surface's list of points to empty list and reset the vertex indices list
        self.points = []
        # Build the perimeter of the surface
        build_perimeter(self)
        # Fill the mesh
        fill_mesh(self)
        # Find the simplices of the surface
        find_simps(self)
        # Filter out the bad triangles
        filter_tris(self)
        if color:
            self.color_tris()

    def color_tris(self, color_scheme='dist', color_map='inferno', inverse=False):
        """
        Colors the triangles in the surface based on the specified coloring scheme and map
        :param inverse: Inverts the color of the color map
        :param color_scheme: Determines how the colors will be mapped
        :param color_map: Determines the actual colors of triangles
        :return: The triangles in the surface are colored
        """
        # Set up the color map
        my_cmap = mpl.colormaps.get_cmap(color_map)
        self.color_map = color_map
        # Set the color scheme for the surface
        self.scheme = color_scheme
        # Default is distance based color map
        if color_scheme == 'dist':

            # Check if the tri_dists have been calculated before
            if self.tri_dists is None or len(self.tri_dists) == 0 or len(self.tri_dists) != len(self.tris):
                # Set up the distances
                dists = []
                self.tri_dists = []
                max_dist, min_dist = 0, np.inf
                # Provide value for the points
                for point in self.points:
                    # Calculate the distance
                    my_dist = calc_dist(point, self.loc)
                    dists.append(my_dist)
                    # Record the minimum and maximum distances
                    if my_dist < min_dist:
                        min_dist = my_dist
                    elif my_dist > max_dist:
                        max_dist = my_dist
                # Go through the triangles in the surface
                for i in range(len(self.tris)):
                    # Find the maximum distance point of the triangles
                    self.tri_dists.append(max([dists[_] for _ in self.tris[i]]))
                # Normalize the tri_dists
                self.tri_dists = [(_ - min_dist) / (max_dist - min_dist) for _ in self.tri_dists]
            self.tri_colors = [my_cmap(_) for _ in self.tri_dists]

        elif color_scheme == 'ins_out':

            # Check if the tri_dists have been calculated before
            if self.tri_ins_out is None or len(self.tri_ins_out) == 0 or len(self.tri_ins_out) != len(self.tris):
                # Set up a list of tracking
                inside_array = []
                # Go through the points in the surface
                for point in self.points:
                    # Calculate the distance between the point and the atom
                    my_dist = calc_dist(point, self.atoms[0].loc)
                    if my_dist < self.atoms[0].rad:
                        inside_array.append(True)
                    else:
                        inside_array.append(False)
                # Now add the triangles
                self.tri_ins_out = []
                for tri in self.tris:
                    if inside_array[tri[0]] and inside_array[tri[1]] and inside_array[tri[2]]:
                        self.tri_ins_out.append(my_cmap(0.25))
                    else:
                        self.tri_ins_out.append(my_cmap(0.75))
            self.tri_colors = self.tri_ins_out
        elif color_scheme == 'curv':
            # Check if the surface is flat
            if self.flat:
                self.tri_curvs = [0] * len(self.tris)
            elif self.tri_curvs is None or len(self.tri_curvs) == 0 or len(self.tri_curvs) != len(self.tris):
                # Get the function coeficients
                A, B, C, D, E, F, G, H, I, J = self.func[:10]
                curvs = []
                min_curv, max_curv = np.inf, 0
                # Get the curvature for each point
                for point in self.points:
                    # Label the points
                    x, y, z = point
                    # Get the gradient of the surface at the point
                    delf = [2 * A * x + D * y + F * z + G, 2 * B * y + D * x + E * z + H, 2 * C * z + E * y + F * x + I]
                    # Calculate the norm of the gradient
                    denominator = np.linalg.norm(delf) ** 4
                    # Calculate the determinant of the hessian matrix and the gradient matrix
                    numerator = np.linalg.det([[2 * A, D, F, delf[0]], [D, 2 * B, E, delf[1]], [F, E, 2 * C, delf[2]], delf + [0]])
                    # Get the curvature
                    curv = - numerator / denominator
                    if curv > max_curv:
                        max_curv = curv
                    if curv < min_curv:
                        min_curv = curv
                    curvs.append(curv)
                # Set up the tri_curvs list
                self.tri_curvs = []
                # Go through the curvature values for each point
                for i in range(len(self.tris)):
                    # Get the triangle
                    tri = self.tris[i]
                    # Get the curvatures
                    my_curvs = [curvs[_] for _ in tri]
                    # Find the maximum curvature point
                    curv_val = max(my_curvs)
                    # Add the curve value to the surface's list of curvatures
                    self.tri_curvs.append(curv_val)
                # Normalize the tri_curvs
                self.tri_curvs = [(_ - min_curv) / (max_curv - min_curv) for _ in self.tri_curvs]
                # Set the curvature for the surface
                self.curv = max_curv
            # Set the colors
            self.tri_colors = [my_cmap(_) for _ in self.tri_curvs]
