import numpy as np
import warnings

warnings.filterwarnings("error")


def calc_surf_func(surf):
    """
    Calculates the coefficients for the surface between the two atoms
    :return: The surface has the correct self.func attribute
    """
    # Make sure that a0 is the atom with the smaller radius
    if surf.atoms[0].rad > surf.atoms[1].rad:
        surf.atoms[0], surf.atoms[1] = surf.atoms[1], surf.atoms[0]
    # Create a0, a1 variables
    a0, a1 = surf.atoms
    # Set the rn vector for the surface since the atoms are sorted
    l0, l1 = np.array(a0.loc), np.array(a1.loc)
    r = l1 - l0
    surf.norm = r / np.linalg.norm(r)
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
    surf.func = ABC + DEF + GHI + [J] + [K] + list(d)


def calc_dist(l0, l1):
    """
    Calculate distance function used to simplify code
    :param l0: Point 0 list, array, n-dimensional must match point 1
    :param l1: Point 1 list, array, n-dimensional must match point 0
    :return: float distance between the two points
    """
    # Pythagorean theorem
    return np.sqrt(sum(np.square(np.array(l0) - np.array(l1))))


def calc_angle(p0, p1, p2=None):
    """
    Finds the angle (in rads) between three points. The first being the common point
    :param p0: Point 0 list, array, n-dimensional must match points 1 and 2
    :param p1: Point 1 list, array, n-dimensional must match points 0 and 2
    :param p2: (optional) Point 2 list, array, n-dimensional must match points 0 and 1
    :return: Angle between (p0, O) and (p1, O) or (p0, p1) and (p0, p2)
    """
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = np.array(p0), np.array(p1)
    else:
        v0, v1 = np.array(p1) - np.array(p0), np.array(p2) - np.array(p0)
    # Get the unit vectors
    try:
        n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
        # Calculate the angle between the two vectors with catches for 180 and 0
        angle = np.arccos(np.clip(np.dot(n0, n1), -1.0, 1.0))
    except RuntimeWarning:
        angle = 0.000001
    return angle


# Calculate tetrahedron volume function.
def calc_tetra_vol(p0, p1, p2, p3):
    """
    Calculates the volume of a tetrahedron defined by its vertices
    :param p0:
    :param p1:
    :param p2:
    :param p3:
    :return: Volume of the tetrahedron
    """
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]
    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    vol = (1/6)*abs(np.dot(r03, np.cross(r01, r02)))
    return vol


def calc_tri(points):
    """
    Takes in 3 points in 3 space and returns the area of the triangle created by them
    :param points:
    :return: Area of the triangle made by the three points
    """
    # Get the two triangles vectors
    AB = np.array(points[0]) - np.array(points[1])
    AC = np.array(points[0]) - np.array(points[2])
    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm((np.cross(AB, AC)))


def calc_com(atoms=None, points=None):
    """
    Takes in a set of points and returns the coordinates of the center of mass
    :param atoms: Atom objects
    :param points: lists or arrays
    :return: Center of mass of the inputs
    """
    # If the function was given atoms, get their points
    if atoms:
        points = [atoms[i].loc for i in range(len(atoms))]
    # Set the running sum for the x, y, z values to 0
    tots = [0 for _ in range(len(points[0]))]
    for point in points:
        for i in range(len(points[0])):
            tots[i] += point[i]
    # Return the center of mass of inputs
    return [tots[i]/len(points) for i in range(len(points[0]))]


def calc_circ(atoms):
    """
    Takes in 3 atoms, calculates the center and radius of inscribed circle
    :param atoms: vorPy Atom objects
    :return: Center and radius of the inscribed circle
    """
    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad
    # Get the relevant variables
    R2, R3 = atoms[1].rad, atoms[2].rad
    x2, y2, z2 = atoms[1].loc[0] - l1[0], atoms[1].loc[1] - l1[1], atoms[1].loc[2] - l1[2]
    x3, y3, z3 = atoms[2].loc[0] - l1[0], atoms[2].loc[1] - l1[1], atoms[2].loc[2] - l1[2]
    # Calculate coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (R1 - R2), R1 ** 2 - R2 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (R1 - R3), R1 ** 2 - R3 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3 = y2*z3 - z2*y3, z2*x3 - x2*z3, x2*y3 - y2*x3
    # More coefficients
    F = a3*b2*c1 - a2*b3*c1 - a3*b1*c2 + a1*b3*c2 + a2*b1*c3 - a1*b2*c3
    Fx0 = b3*c2*f1 - b2*c3*f1 - b3*c1*f2 + b1*c3*f2
    Fx1 = b3*c2*d1 - b2*c3*d1 - b3*c1*d2 + b1*c3*d2
    Fy0 = - a3*c2*f1 + a2*c3*f1 + a3*c1*f2 - a1*c3*f2
    Fy1 = - a3*c2*d1 + a2*c3*d1 + a3*c1*d2 - a1*c3*d2
    Fz0 = a3*b2*f1 - a2*b3*f1 - a3*b1*f2 + a1*b3*f2
    Fz1 = a3*b2*d1 - a2*b3*d1 - a3*b1*d2 + a1*b3*d2
    # Catch for F=0 (i.e. no circle exists)
    if F == 0:
        return
    # Find the radius of the tangential circle using the quadratic formula
    a = (Fx1 ** 2 + Fy1 ** 2 + Fz1 ** 2) / F ** 2 - 1
    b = 2 * (Fx0 * Fx1 + Fy0 * Fy1 + Fz0 * Fz1) / F ** 2 - 2 * R1
    c = (Fx0 ** 2 + Fy0 ** 2 + Fz0 ** 2) / F ** 2 - R1 ** 2
    # Calculate the discriminant.
    disc = b ** 2 - 4 * a * c
    # If the discriminant is negative then the tangential circle does not exist.
    if round(disc, 10) > 0:
        # Grab the two roots
        Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
        # If there is only one root return it
        if len(Rs) == 1:
            R = Rs[0]
        # If there are 2 roots choose between them
        else:
            # If the smaller of the two roots is negative return the other root
            if min(Rs) < 0:
                R = max(Rs)
            # If they're both positive, return the smaller of the two
            elif Rs[0] > 0 and Rs[1] > 0:
                R = min(Rs)
            # If they're both negative return
            else:
                return
        # Calculate the vertex based off of our coefficient values and the sphere's radius
        x = Fx0 / F + R * Fx1 / F + l1[0]
        y = Fy0 / F + R * Fy1 / F + l1[1]
        z = Fz0 / F + R * Fz1 / F + l1[2]
        return [[x, y, z], R]


def rotate_points(vec, points, reverse=False):
    """
    Takes in a set of points and a vector and rotates the points and the vector so the v = [0,0,1]
    :param vec: The vector about which the surface is rotated
    :param points: the points of the surface
    :param reverse: Bool for rotating the surface back
    :return: List of rotated points
    """
    # Get the vx, vy, vz vector components
    vx, vy, vz = vec
    # If vy or vz are zero we need a catch for divide by zero error.
    if round(vy, 2) == 0:
        phi = np.pi / 2
    else:
        phi = np.arctan(vx / vy)
    if round(vz, 2) == 0:
        theta = np.pi / 2
    else:
        theta = np.arctan(vy / vz)
    # If the points are to be sent back, provide the negative values for the angles
    if reverse:
        theta, phi = -theta, -phi
    # Get variables for sin(theta), cos(theta), sin(phi), cos(phi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    nps = []
    for p in points:
        px, py, pz = np.round(p[0], 7), np.round(p[1], 7), np.round(p[2], 7)
        # Multiplying the x, y rotation matrices gives the following:
        npx = px * cp - py * sp
        npy = px * ct * sp + py * ct * cp - pz * st
        npz = px * st * sp + py * st * cp + pz * ct
        # Add the new points to the list
        nps.append([npx, npy, npz])
    return nps


def ndx_search(ndxs_list, ndxs):
    """
     Searches a list of indices of atoms sorted by smallest atom and where the vertex would be
    :param ndxs_list: The index for checking
    :param ndxs: The indices to check against
    :return: The vertex index of the vertex or where the vertex should be inserted
    """
    # If the length of the test list is equal to 0 return the next index
    if len(ndxs_list) <= 1:
        # If there exists one vertex already and the new vertex is less than the old vertex return 1
        if len(ndxs_list) > 0 and ndxs > ndxs_list[0]:
            return 1
        # Otherwise, return 0
        return 0
    # Get the middle of the list of vertices
    mid_list_ndx = len(ndxs_list) // 2
    # If the search element (my_list) is greater than the test element (test_lol) search the lower half of test_lol
    if ndxs > ndxs_list[mid_list_ndx]:
        ndxs_ndx = ndx_search(ndxs_list[mid_list_ndx:], ndxs)
        return ndxs_ndx + mid_list_ndx
    # If the search element (my_list) is less than the test element (test_lol) search the upper half of test_lol
    elif ndxs < ndxs_list[mid_list_ndx]:
        ndxs_ndx = ndx_search(ndxs_list[:mid_list_ndx], ndxs)
        return ndxs_ndx
    # If the search element (my_list) is greater than the test element (test_lol) search the lower half of test_lol
    elif ndxs == ndxs_list[mid_list_ndx]:
        return mid_list_ndx


def get_time(seconds):
    """
    Turns seconds into hours, minutes and seconds
    :param seconds:
    :return: hours, minutes, seconds
    """
    # Divide up the values
    hours = seconds // 3600
    minutes = (seconds - (hours * 3600)) // 60
    seconds = seconds - hours * 3600 - minutes * 60
    # Return the values
    return hours, minutes, seconds


def calc_vol(self):
    # Create the volume variable
    vol = 0
    # Go through each surface on the atom
    for surf in self.surfs:
        # Set the surface area
        self.sa += surf.sa
        # Check to see if the surface's volume has been calculated already
        if surf.vols[surf.ndx.index(self.num)] != 0:
            vol += surf.vols[surf.ndx.index(self.num)]
        else:
            # Calculate the volume of the
            for tri in surf.tris:
                p0, p1, p2, p3 = self.loc, surf.points[tri[0]], surf.points[tri[1]], surf.points[tri[2]]
                my_vol = calc_tetra_vol(p0, p1, p2, p3)
                surf.vols[surf.ndx.index(self.num)] = my_vol
                vol += my_vol
    # Return the volume
    self.vol = vol
    return vol


def get_radius(self):
    """
        Finds the radius of the atom from the symbol or vice versa

    :return: The radius of the atom from the symbol or vice versa
    """
    radii = self.sys.radii
    # Get the radius and the element from the name of the atom
    if self.name is not None and self.name in self.sys.special_radii:
        self.element, self.rad = self.sys.special_radii[self.name]
    # If indicated we return the symbol of atom that the radius indicates
    elif self.element is None:
        # Check to see if the radius is in the system
        if self.rad in {radii[_] for _ in radii[1]}:
            self.element = radii[self.rad]
        else:
            # Get the closest atom to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in radii:
                if radii[radius] - self.rad < min_diff:
                    self.element = radii[radius]
    # If we have the type and just want the radius, keep scanning until we find the radius
    elif self.rad is None:
        self.rad = radii[self.element.lower()]


# Calculate curvature method
def calc_surf_curv(self):
    """
    Calculates the curvature of the surface
    :return: The curvature attribute is filled
    """
    # Check to see that the function has been calculated or not
    if self.func is None:
        calc_surf_func(self)
    # Made up function to calculate the general curvature of the hyperboloid
    self.curv = np.sqrt(self.func[0]**2 + self.func[1]**2 + self.func[2]**2)


def calc_surf_sa(self):
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