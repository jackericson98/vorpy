import numpy as np
import warnings
from numba import jit
from numba.core.errors import TypingError
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


def project_to_plane(points, plane_point, plane_normal):
    # Normalize the normal vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Create an orthogonal basis for the plane
    if (plane_normal == np.array([1.0, 0.0, 0.0])).all() or (plane_normal == np.array([-1.0, 0.0, 0.0])).all():
        # Handle the case where the normal is along the x-axis
        u = np.array([0, 1, 0])
    else:
        u = np.cross(plane_normal, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(plane_normal, u)
    v = v / np.linalg.norm(v)

    # Project points onto the plane
    projected_points = []
    for point in points:
        # Vector from point on plane to the point in space
        point_vector = point - plane_point
        # Distance from point to plane
        distance = np.dot(point_vector, plane_normal)
        # Projection of point onto plane
        projection = point - distance * plane_normal
        # Convert projection to 2D coordinates
        u_coord = np.dot(projection - plane_point, u)
        v_coord = np.dot(projection - plane_point, v)
        projected_points.append((u_coord, v_coord))

    return projected_points


def unproject_to_3d(projected_points, plane_point, plane_normal):
    # Normalize the normal vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Create an orthogonal basis for the plane
    if (plane_normal == np.array([1.0, 0.0, 0.0])).all() or (plane_normal == np.array([-1.0, 0.0, 0.0])).all():
        # Handle the case where the normal is along the x-axis
        u = np.array([0, 1, 0])
    else:
        u = np.cross(plane_normal, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(plane_normal, u)
    v = v / np.linalg.norm(v)

    # Map 2D coordinates back to 3D plane
    reconstructed_points = []
    for u_coord, v_coord in projected_points:
        # Reconstruct 3D point on the plane using the basis vectors and plane point
        point_3d = plane_point + u_coord * u + v_coord * v
        reconstructed_points.append(point_3d)

    return reconstructed_points


def map_to_plane(points_2d, plane_point, plane_normal):
    # Normalize the normal vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Create an orthogonal basis for the plane
    if (plane_normal == np.array([1.0, 0.0, 0.0])).all() or (plane_normal == np.array([-1.0, 0.0, 0.0])).all():
        # Handle the case where the normal is along the x-axis
        u = np.array([0, 1, 0])
    else:
        u = np.cross(plane_normal, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(plane_normal, u)
    v = v / np.linalg.norm(v)

    # Map 2D points to the 3D plane
    mapped_points = []
    for point_2d in points_2d:
        u_coord, v_coord = point_2d
        # Calculate the corresponding 3D point
        point_3d = plane_point + u_coord * u + v_coord * v
        mapped_points.append(point_3d)

    return mapped_points


def calc_dist(l0, l1):

    return np.sqrt(sum(np.square(np.array(l0) - np.array(l1))))


@jit(nopython=True)
def calc_dist_numba(l0, l1):
    """
    Calculate distance function used to simplify code
    :param l0: Point 0 list, array, n-dimensional must match point 1
    :param l1: Point 1 list, array, n-dimensional must match point 0
    :return: float distance between the two points
    """
    # Pythagorean theorem
    return np.sqrt(sum(np.square(l0 - l1)))


@jit(nopython=True)
def calc_angle_jit(p0, p1, p2=None):
    """
    Finds the angle (in rads) between three points
    :param p0: Point 0 list, array, n-dimensional must match points 1 and 2
    :param p1: Point 1 list, array, n-dimensional must match points 0 and 2
    :param p2: (optional) Point 2 list, array, n-dimensional must match points 0 and 1
    :return: Angle between (p0, O) and (p1, O) or (p0, p1) and (p0, p2)
    """
    # If no p2 is given, use the origin
    if p2 is None:
        v0, v1 = p0, p1
    else:
        v0, v1 = p1 - p0, p2 - p0
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
    # Calculate the angle between the two vectors with catches for 180 and 0
    my_dot = np.dot(n0, n1)
    if my_dot <= -1.0:
        my_dot = -1.0
    elif my_dot >= 1.0:
        my_dot = 1.0
    angle = np.arccos(my_dot)
    return angle


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
        v0, v1 = p0, p1
    else:
        v0, v1 = p1 - p0, p2 - p0
    n0, n1 = v0/np.linalg.norm(v0), v1/np.linalg.norm(v1)
    # Calculate the angle between the two vectors with catches for 180 and 0
    my_dot = np.dot(n0, n1)
    if my_dot <= -1.0:
        my_dot = -1.0
    elif my_dot >= 1.0:
        my_dot = 1.0
    angle = np.arccos(my_dot)
    return angle


@jit(nopython=True)
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
    r03 = np.array([p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]])

    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    return (1/6)*abs(np.dot(r03, np.cross(r01, r02)))


def calc_tetra_inertia(ps, mass):
    """
    Calculate the moment of inertia tensor of a tetrahedron about its centroid.
    This formula assumes uniform density and calculates the inertia tensor.
    :param ps: List of four vertices of the tetrahedron.
    :param mass: Mass of the tetrahedron.
    :return: 3x3 inertia tensor of the tetrahedron.
    """
    # Placeholder for inertia tensor calculation.
    # For simplicity, this uses an approximate inertia formula for a solid tetrahedron.
    # More accurate calculations can be done by integrating over the volume.
    inertia_tensor = np.zeros((3, 3))

    # Sum contributions from the vertices
    for i in range(4):
        x, y, z = ps[i]
        inertia_tensor[0, 0] += mass * (y ** 2 + z ** 2) / 10.0
        inertia_tensor[1, 1] += mass * (x ** 2 + z ** 2) / 10.0
        inertia_tensor[2, 2] += mass * (x ** 2 + y ** 2) / 10.0
        inertia_tensor[0, 1] -= mass * x * y / 10.0
        inertia_tensor[0, 2] -= mass * x * z / 10.0
        inertia_tensor[1, 2] -= mass * y * z / 10.0

    # Symmetric tensor: fill in the other values
    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    return inertia_tensor


@jit(nopython=True)
def calc_tri(points):
    """
    Takes in 3 points and returns the area of the triangle created by them
    :param points: 3D points
    :return: Area of the triangle made by the three points
    """
    # Get the two triangles vectors
    ab = [points[0][0] - points[1][0], points[0][1] - points[1][1], points[0][2] - points[1][2]]
    ac = [points[0][0] - points[2][0], points[0][1] - points[2][1], points[0][2] - points[2][2]]

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
    return np.array([tots[i]/len(points) for i in range(len(points[0]))])


@jit(nopython=True)
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
            length += calc_dist_numba(point, points[m + 1])
    return length


def calc_sphericity(volume, surface_area):
    """
    Calculate the sphericity of a geometric object based on its volume and surface area.

    Parameters:
    - volume (float): The volume of the object.
    - surface_area (float): The surface area of the object.

    Returns:
    - float: The sphericity of the object.
    """
    if volume <= 0 or surface_area <= 0:
        raise ValueError("Volume and surface area must be positive numbers.")

    # Calculate sphericity using the geometric formula
    sphericity = (np.pi ** (1/3) * (6 * volume) ** (2/3)) / surface_area
    return sphericity


def calc_isoperimetric_quotient(volume, surface_area):
    """
    Calculate the isoperimetric quotient

    Parameters:
        - volume (float): The volume of the object.
        - surface_area (float): The surface area of the object.

    Returns:
        - isoperimetric quotient (float): The isoperimetric quotient of the object
    """
    if volume <= 0 or surface_area <= 0:
        raise ValueError("Volume and surface area must be positive numbers.")

    return (36 * np.pi * volume ** 2) / (surface_area ** 3)


def calc_spikes(ball_loc, surfs):
    # Create the spike list
    spikes = []
    # Loop through the surfaces
    for surf in surfs:
        for point in surf['points']:
            spikes.append(calc_dist(ball_loc, point))

    # Return the minimum, average and maximum spike dist
    return min(spikes), max(spikes)


def calc_cell_box(surfs):
    # Create the mins and maxs varaibles
    mins, maxs = [np.inf, np.inf, np.inf], [-np.inf, -np.inf, -np.inf]
    # Loop through the surfaces
    for surf in surfs:
        for point in surf['points']:
            for i in range(3):
                if point[i] < mins[i]:
                    mins[i] = point[i]
                if point[i] > maxs[i]:
                    maxs[i] = point[i]
    # Return the bounding box for the cell
    return [mins, maxs]


def calc_cell_com(ball_loc, surfs, volume):
    # Create the mass_locs list
    mass_locs = []
    for surf in surfs:
        for tri in surf['tris']:
            # Get the points of the tetrahedron
            ps = [ball_loc, *[surf['points'][_] for _ in tri]]
            # Calculate the centroid of the tetrahedron
            tet_com = [sum([ps[j][i] for j in range(4)]) / 4 for i in range(3)]
            # Calculate the volume of the tetrahedron
            tet_vol = calc_tetra_vol(*ps)
            # Append the volume-weighted centroid
            mass_locs.append([tet_vol * coord for coord in tet_com])

    # Calculate the total center of mass by normalizing with the cell volume
    return np.array([sum(coords) / volume for coords in zip(*mass_locs)])


def calc_cell_moi(ball_loc, surfs, volume, density=1.0):
    """
    Calculate the moment of inertia of a cell with respect to `ball_loc` using tetrahedrons.
    :param ball_loc: Center location of the cell.
    :param surfs: List of surfaces, each containing points and tris (triangles).
    :param volume: Total volume of the cell.
    :param density: Density of the material (default is 1.0).
    :return: 3x3 Moment of inertia tensor of the cell.
    """
    # Create an inertia tensor initialized to zero
    inertia_tensor = np.zeros((3, 3))

    # Iterate through each surface and triangle to calculate the tetrahedron MOI contributions
    for surf in surfs:
        for tri in surf['tris']:
            # Get the points of the tetrahedron
            ps = [ball_loc, *[surf['points'][_] for _ in tri]]

            # Calculate the centroid of the tetrahedron
            tet_com = [sum([ps[j][i] for j in range(4)]) / 4 for i in range(3)]

            # Calculate the volume of the tetrahedron
            tet_vol = calc_tetra_vol(*ps)

            # Calculate the mass of the tetrahedron
            tet_mass = density * tet_vol

            # Calculate the inertia tensor of the tetrahedron about its centroid
            tet_inertia_tensor = calc_tetra_inertia(ps, tet_mass)

            # Calculate the distance vector from the tetrahedron centroid to the cell's center (`ball_loc`)
            r = np.array(tet_com) - np.array(ball_loc)
            r_squared = np.dot(r, r)

            # Use the parallel axis theorem to adjust the inertia tensor to the cell's center
            shift_tensor = tet_mass * (r_squared * np.identity(3) - np.outer(r, r))

            # Add the adjusted tensor to the total inertia tensor
            inertia_tensor += tet_inertia_tensor + shift_tensor

    return inertia_tensor


def combine_inertia_tensors(inertia_tensors, centroids, common_centroid, masses):
    """
    Combine a list of inertia tensors into a single inertia tensor.

    :param inertia_tensors: List of 3x3 inertia tensors (numpy arrays) for each element.
    :param centroids: List of centroid coordinates (numpy arrays) for each element.
    :param common_centroid: The centroid to which all tensors should be shifted (numpy array).
    :param masses: List of masses (or volumes, assuming uniform density) for each element.
    :return: Combined 3x3 inertia tensor.
    """
    # Initialize the total inertia tensor as a zero matrix
    I_total = np.zeros((3, 3))

    # Loop over each element
    for I_i, C_i, m_i in zip(inertia_tensors, centroids, masses):
        # Calculate the displacement vector from the element's centroid to the common centroid
        d = C_i - common_centroid
        d_squared = np.dot(d, d)  # Squared magnitude of the displacement vector

        # Compute the parallel axis theorem adjustment term
        shift_tensor = m_i * (d_squared * np.eye(3) - np.outer(d, d))

        # Shift the inertia tensor of the element to the common centroid and add to total
        I_shifted = I_i + shift_tensor
        I_total += I_shifted

    return I_total


def calc_total_inertia_tensor(spheres, common_point):
    """
    Calculate the combined moment of inertia tensor of a set of spheres.

    :param spheres: List of dictionaries with 'mass', 'radius' (as 'rad'), and 'location' (as 'loc', numpy array).
    :param common_point: The point to which all moments of inertia are shifted (numpy array).
    :return: Total moment of inertia tensor (3x3 numpy array).
    """
    # Initialize the total inertia tensor as a 3x3 zero matrix
    I_total = np.zeros((3, 3))

    # Iterate through each sphere
    for sphere in spheres:
        m = sphere['mass']
        r = sphere['rad']
        loc = sphere['loc']

        # Moment of inertia tensor of the sphere about its own center (3x3 identity scaled by (2/5) * m * r^2)
        I_center = (2 / 5) * m * r**2 * np.eye(3)

        # Calculate the displacement vector from the sphere's center to the common point
        d = loc - common_point
        d_squared = np.dot(d, d)  # Squared magnitude of the displacement vector

        # Calculate the parallel axis shift tensor: m * (d^2 * I3 - d * d^T)
        shift_tensor = m * (d_squared * np.eye(3) - np.outer(d, d))

        # Shift the inertia tensor to the common reference point
        I_shifted = I_center + shift_tensor

        # Add the shifted inertia tensor to the total inertia tensor
        I_total += I_shifted

    return I_total


def calc_contacts(loc, rad, surfs, surf_ndxs):
    """
    Calculate the contact areas and contribution volume for a given sphere.
    :param loc: Center location of the sphere.
    :param rad: Radius of the sphere.
    :param surfs: List of surfaces with points and triangles.
    :return: A list of contact areas for each surface and the total contribution volume.
    """
    # Create the area and volume vals
    contact_areas, contribution_vol = {}, 0

    # Loop through the surfaces
    for i, surf in enumerate(surfs):
        # Initialize contact area for this surface
        contact_area = 0
        new_points = []
        point_inside = []

        # Loop through the points to determine if inside or outside
        for point in surf['points']:
            distance = calc_dist(point, loc)
            if distance <= rad:
                point_inside.append(True)
                new_points.append(point)
            else:
                point_inside.append(False)
                # Get the direction and normalize it
                direction = point - loc
                norm = np.linalg.norm(direction)
                if norm > 0:
                    # Project the point onto the sphere's surface
                    new_points.append(rad * (direction / norm) + loc)
                else:
                    new_points.append(point)  # If the point coincides with the center (rare edge case)

        # Loop through the triangles
        for tri in surf['tris']:
            triangle_points = [surf['points'][index] for index in tri]
            projected_points = [new_points[index] for index in tri]
            inside_flags = [point_inside[index] for index in tri]

            # Determine if the triangle is fully inside, fully outside, or mixed
            all_inside = all(inside_flags)
            all_outside = not any(inside_flags)
            mixed = not all_inside and not all_outside

            if all_inside:
                # Triangle is fully inside the sphere
                contact_area += calc_tri(np.array(triangle_points))
                contribution_vol += calc_tetra_vol(loc, *triangle_points)
            elif all_outside:
                # Triangle is fully outside the sphere
                contribution_vol += calc_tetra_vol(loc, *projected_points)
            elif mixed:
                # Triangle is partially inside and outside
                # We add the volume using a mix of inside and projected points
                mixed_points = [triangle_points[i] if inside_flags[i] else projected_points[i] for i in range(3)]
                contribution_vol += calc_tetra_vol(loc, *mixed_points)
                # Count the triangle as outside for contact area if any point is outside
                if inside_flags.count(True) < 3:
                    contact_area += calc_tri(np.array(triangle_points))

        # Append the contact area for this surface
        contact_areas[surf_ndxs[i]] = contact_area

    return contact_areas, contribution_vol


@jit(nopython=True)
def rotate_points1(vec, points, reverse=False):
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


def rotate_points(vec, points, reverse=False):
    if reverse:
        vec = - vec
    vx, vy, vz = vec
    mag = np.sqrt(vx**2 + vy**2 + vz**2)
    phi = np.arctan2(vy, vx)
    theta = np.arccos(vz / mag)
    if reverse:
        theta, phi = -theta, -phi

    # Forward rotations to align with z-axis
    Rz = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])

    # Combine rotations to align vector with +z direction
    if reverse:
        # Correct sequence for inverse rotation
        rotation_matrix = np.dot(Rz, Ry)
    else:
        # Correct sequence for forward rotation
        rotation_matrix = np.dot(Ry, Rz)

    # Apply rotation to all points
    rotated_points = [np.dot(rotation_matrix, p) for p in points]
    return rotated_points


@jit(nopython=True)
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


def calc_vol(a_loc, surfs_points, surfs_tris):
    """
    Calculates the volume of an ball using its surfaces
    :param ball: ball object for volume calculation
    :return: returns the volume for the ball object
    """
    # Create the volume variable
    surf_vols = []
    # Go through each surface on the ball
    for i in range(len(surfs_points)):
        # Calculate the volume of the
        surf_vol = 0
        for tri in surfs_tris[i]:
            # Calculate the tetrahedron volume between the balls' location and the surface triangle's points
            surf_vol += calc_tetra_vol(np.array(a_loc), surfs_points[i][tri[0]], surfs_points[i][tri[1]], surfs_points[i][tri[2]])
        # Add the surface's volume to the list
        surf_vols.append(surf_vol)
    # Get the total volume by summing the surfaces volumes
    vol = sum(surf_vols)
    # Set the volume and return it
    return vol, surf_vols
