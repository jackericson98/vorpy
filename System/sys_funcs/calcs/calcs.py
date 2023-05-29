import numpy as np
import warnings
warnings.filterwarnings("error")


def round_func(round_to):
    """
    Nested round function for defining round schemes and rounding multiple values
    :param round_to: int - number of decimal places
    :return: Round function set to round to value
    """
    # Define the inner round function
    def round_(val, new_num=None):
        """
        Inner round function operating on outer defined round to value
        :param val: float/iterable - val(s) to be rounded
        :param new_num: New round to value
        :return: float/list - rounded values
        """
        # Set the new round to number if specified
        if new_num is None:
            new_num = round_to
        # Return the values
        try:
            return round(val, new_num)
        except TypeError:
            return [round(_, new_num) for _ in val]
    # Return the function for the outer function
    return round_


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
    Finds the angle (in rads) between three points
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
    :param p0: Point 0
    :param p1: Point 1
    :param p2: Point 2
    :param p3: Point 3
    :return: Volume of the tetrahedron made by the points
    """
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]

    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    return (1/6)*abs(np.dot(r03, np.cross(r01, r02)))


def calc_tri(points):
    """
    Takes in 3 points and returns the area of the triangle created by them
    :param points: 3D points
    :return: Area of the triangle made by the three points
    """
    # Get the two triangles vectors
    ab = np.array(points[0]) - np.array(points[1])
    ac = np.array(points[0]) - np.array(points[2])

    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm((np.cross(ab, ac)))


def calc_com(points):
    """
    Takes in a set of points and returns the coordinates of the center of mass
    :param points: lists of locations in n-dimensions
    :return: Center of mass of the inputs
    """

    # Set the running sum for the x, y, z values to 0
    tots = [0 for _ in range(len(points[0]))]
    for point in points:
        for i in range(len(points[0])):
            tots[i] += point[i]

    # Return the center of mass of inputs
    return [tots[i]/len(points) for i in range(len(points[0]))]


def calc_length(points):
    """
    Calculates the total length of the points assuming they are in order
    :param points: Points for length calculations
    :return: float total length between consecutive points
    """
    # Reset the length
    length = 0
    # Go through the points in the list
    for m, point in enumerate(points):
        # Make sure not to index error
        if m + 1 < len(points):
            # Add the length to the total
            length += calc_dist(point, points[m + 1])
    return length


def calc_circ(l0, l1, l2, r0, r1, r2):
    """
    Takes in 3 atoms, calculates the center and radius of inscribed circle
    :param : Locations and radii for the circle
    :return: Center and radius of the inscribed circle
    """
    # Move the other atoms to the location of the first
    x2, y2, z2 = l1[0] - l0[0], l1[1] - l0[1], l1[2] - l0[2]
    x3, y3, z3 = l2[0] - l0[0], l2[1] - l0[1], l2[2] - l0[2]
    # Calculate coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (r0 - r1), r0 ** 2 - r1 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (r0 - r2), r0 ** 2 - r2 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
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
    b = 2 * (Fx0 * Fx1 + Fy0 * Fy1 + Fz0 * Fz1) / F ** 2 - 2 * r0
    c = (Fx0 ** 2 + Fy0 ** 2 + Fz0 ** 2) / F ** 2 - r0 ** 2
    # Calculate the discriminant.
    disc = b ** 2 - 4 * a * c
    # If the discriminant is negative then the tangential circle does not exist.
    if round(disc, 10) > 0:
        # Grab the two roots
        rs = [_ for _ in np.roots([a, b, c]) if np.isreal(_)]
        # If there is only one root return it
        if len(rs) == 1:
            r = rs[0]
        # If there are 2 roots choose between them
        else:
            # If the smaller of the two roots is negative return the other root
            if min(rs) < 0:
                r = max(rs)
            # If they're both positive, return the smaller of the two
            elif rs[0] > 0 and rs[1] > 0:
                r = min(rs)
            # If they're both negative return
            else:
                return
        # Calculate the vertex based off of our coefficient values and the sphere's radius
        x = Fx0 / F + r * Fx1 / F + l0[0]
        y = Fy0 / F + r * Fy1 / F + l0[1]
        z = Fz0 / F + r * Fz1 / F + l0[2]
        return [[x, y, z], r]


def calc_surf_norm(l0, l1):
    """
    Calculates the normal between the atoms from small to large
    :param l0: Location of the first atom
    :param l1: Location of the second atom
    :return: The normal to the vector between the atoms
    """
    # Get the vector from l0 to l1
    r = np.array(l1) - np.array(l0)
    r = [_ if _ != 0 else 0.00001 for _ in r]
    return r / np.linalg.norm(r)


def calc_surf_func(l0, r0, l1, r1):
    """
    Calculates the coefficients for the surface between the two atoms
    :return: Returns a function for the hyperboloid between the atoms
    """
    # Check the smaller atom is first
    if r1 < r0:
        l0, r0, l1, r1 = l1, r1, l0, r0
    # Grab the centers of the spheres
    x1, y1, z1 = l0
    x2, y2, z2 = l1
    # Calculate the major coefficients (pg. 574 Z. Hu)
    R = r0 - r1
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
    # Return the function coefficients
    return ABC + DEF + GHI + [J] + [K] + list(d)


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


def calc_surf_sa(edges, com, tris, points, flat):
    """
    Calculates the surface area of the input surface
    :param edges: Edges for the surface if the surface is flat
    :param com: Center of mass point for the surface used for calculations
    :param tris: Triangles for summing
    :param points: Points for the surface
    :param flat: Whether the surface is flat or not
    :return: Surface area of the surface
    """
    # Create the surface area variable
    sa = 0
    if flat:
        for edge in edges:
            for i in range(len(edge.points) - 1):
                sa += calc_tri([*edge.points[i:i + 2], com])
    # Go through the triangles in the surface
    else:
        for tri in tris:
            sa += calc_tri([points[tri[_]] for _ in range(3)])
    # Return the surface area
    return sa


def calc_surf_tri_dists(points, tris, loc, max_val=None):
    """
    Calculate the distances between each triangle and the provided location
    :param points: points from the surface
    :param tris: triangles from the
    :param loc: location for the distance calculations
    :return: List of distance corresponding to the triangles on a surface
    """
    # Set up the distances
    dists = []
    tri_dists = []
    if max_val is not None:
        max_dist = max_val
    max_dist, min_dist = 0, np.inf
    # Provide value for the points
    for point in points:
        # Calculate the distance
        my_dist = calc_dist(point, loc)
        dists.append(my_dist)
        # Record the minimum and maximum distances
        if my_dist < min_dist:
            min_dist = my_dist
        elif my_dist > max_dist:
            max_dist = my_dist
    # Go through the triangles in the surface
    for i in range(len(tris)):
        # Find the maximum distance point of the triangles
        tri_dists.append(max([dists[_] for _ in tris[i]]))
    # Normalize the tri_dists
    return [(_ - min_dist) / (max_dist - min_dist) for _ in tri_dists]


def calc_surf_point_curv(func, point):
    # Get the function coefficients
    A, B, C, D, E, F, G, H, I, J = func[:10]
    # Label the points
    x, y, z = point
    # Get the gradient of the surface at the point
    delf = [2 * A * x + D * y + F * z + G, 2 * B * y + D * x + E * z + H, 2 * C * z + E * y + F * x + I]
    # Calculate the norm of the gradient
    denominator = np.linalg.norm(delf) ** 4
    # Calculate the determinant of the hessian matrix and the gradient matrix
    numerator = np.linalg.det([[2 * A, D, F, delf[0]], [D, 2 * B, E, delf[1]], [F, E, 2 * C, delf[2]],
                               delf + [0]])
    # Get the curvature
    return - numerator / denominator


def calc_surf_tri_curvs(func, points, tris, max_curv):
    """
    Calculates the curvature of the triangles
    :param calc_max_curv:
    :param func: Surface function
    :param points: points for the surface
    :param tris: triangles in the surface
    :return: A list of curvature values for the triangles
    """
    curvs = []
    min_curv = np.inf
    # If the surface normal is within the surface,
    # Get the curvature for each point
    for point in points:
        curv = calc_surf_point_curv(func, point)
        if curv < min_curv:
            min_curv = curv
        elif curv > max_curv:
            max_curv = curv

        curvs.append(curv)
    # Set up the tri_curvs list
    tri_curvs = []
    # Go through the curvature values for each point
    for i in range(len(tris)):
        # Get the triangle
        tri = tris[i]
        # Get the curvatures
        curv_val = sum([curvs[_] for _ in tri])/3
        # Add the curve value to the surface's list of curvatures
        tri_curvs.append(curv_val)
    # Return the values
    return tri_curvs, max_curv


def calc_surf_tri_ins_out(surf):
    """
    Calculates whether the triangles in the surface are inside the overlapping atoms or not
    :param surf: Surface for calculations
    :return: List of bools for if the triangles in the surface are inside or out
    """
    # Set up a list of tracking
    inside_array = []
    # Go through the points in the surface
    for point in surf.points:
        # Calculate the distance between the point and the atom
        my_dist = calc_dist(point, surf.atoms[0].loc)
        # Check if the triangle is inside or not
        if my_dist < surf.atoms[0].rad:
            inside_array.append(True)
        else:
            inside_array.append(False)
    # Now add the triangles
    surf.tri_ins_out = []
    # Color the tris
    for tri in surf.tris:
        if inside_array[tri[0]] and inside_array[tri[1]] and inside_array[tri[2]]:
            surf.tri_ins_out.append(0.25)
        else:
            surf.tri_ins_out.append(0.75)


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
    :param seconds: Number of seconds in the counter
    :return: hours, minutes, seconds
    """
    # Divide up the values
    hours = seconds // 3600
    minutes = (seconds - (hours * 3600)) // 60
    seconds = seconds - hours * 3600 - minutes * 60
    # Return the values
    return hours, minutes, seconds


def calc_vol(atom):
    """
    Calculates the volume of an atom using its surfaces
    :param atom: Atom object for volume calculation
    :return: returns the volume for the atom object
    """
    # Create the volume variable
    vol = 0
    # Go through each surface on the atom
    for surf in atom.surfs:
        # Set the surface area
        atom.sa += surf.sa
        # Check to see if the surface's volume has been calculated already
        if surf.vols[surf.ndx.index(atom.num)] != 0:
            vol += surf.vols[surf.ndx.index(atom.num)]
        else:
            # Calculate the volume of the
            surf_vol = 0
            for tri in surf.tris:
                p0, p1, p2, p3 = atom.loc, surf.points[tri[0]], surf.points[tri[1]], surf.points[tri[2]]
                surf_vol += calc_tetra_vol(p0, p1, p2, p3)
            vol += surf_vol
            surf.vols[surf.ndx.index(atom.num)] = surf_vol
    # Set the volume and return it
    atom.vol = vol
    return vol


def get_radius(atom):
    """
    Finds the radius of the atom from the symbol or vice versa
    :return: The radius of the atom from the symbol or vice versa
    """
    radii, special_radii = atom.sys.radii, atom.sys.special_radii
    # Get the radius and the element from the name of the atom
    if atom.res is not None and atom.res.name in atom.sys.special_radii:
        # Check if no atom name exists or its empty
        if atom.name is not None and atom.name != '':
            for i in range(len(atom.name)):
                name = atom.name[:-i]
                # Check the residue name
                if name in special_radii[atom.res.name]:
                    atom.rad = special_radii[atom.res.name][name]
    # If we have the type and just want the radius, keep scanning until we find the radius
    if atom.rad is None and atom.element.lower() in radii:
        atom.rad = radii[atom.element.lower()]
    # If indicated we return the symbol of atom that the radius indicates
    if atom.rad is None or atom.rad == 0:
        # Check to see if the radius is in the system
        if atom.rad in {radii[_] for _ in radii[1]}:
            atom.element = radii[atom.rad]
        else:
            # Get the closest atom to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in radii:
                if radii[radius] - atom.rad < min_diff:
                    atom.element = radii[radius]
