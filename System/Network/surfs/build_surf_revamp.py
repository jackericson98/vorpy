from numpy import array, dot, linalg, float64, allclose, cross, cos, sin, arctan2, c_, zeros, arccos
import matplotlib.pyplot as plt
from System.sys_funcs.calcs.surf import calc_surf_func
from System.sys_funcs.calcs.calcs import calc_dist_numba
from System.Network.surfs.perimeter import build_perimeter, build_perimeter1
from System.Network.surfs.fill import calc_surf_point
from Visualize.mpl_visualize import plot_surfs, plot_balls
import triangle as tr
from shapely.geometry import Point, Polygon


def normalize(v):
    return v / linalg.norm(v)


def rectangle_on_plane(center, normal, width, height):
    normal = normalize(array(normal))
    center = array(center)

    # Find a vector not parallel to the normal
    if allclose(normal, [1, 0, 0]) or allclose(normal, [-1, 0, 0]):
        v = array([0, 1, 0])
    else:
        v = array([1, 0, 0])

    # First vector in the plane
    u1 = normalize(cross(normal, v))
    # Second vector in the plane, perpendicular to the first
    u2 = normalize(cross(normal, u1))

    # Half dimensions offsets
    half_width = (u1 * width / 2)
    half_height = (u2 * height / 2)

    # Rectangle vertices
    p1 = center - half_width - half_height
    p2 = center + half_width - half_height
    p3 = center + half_width + half_height
    p4 = center - half_width + half_height

    return array([p1, p2, p3, p4])


def rotation_matrix(axis, theta):
    axis = normalize(axis)
    a = cos(theta / 2.0)
    b, c, d = -axis * sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def transform_points(points, plane_normal, center_point):
    # Normalize plane normal
    plane_normal = normalize(array(plane_normal))
    center_point = array(center_point)

    # Calculate the rotation
    z_axis = array([0, 0, 1])
    rotation_axis = cross(plane_normal, z_axis)
    if linalg.norm(rotation_axis) != 0:
        rotation_axis = normalize(rotation_axis)
        angle = arccos(dot(plane_normal, z_axis))
        rot_matrix = rotation_matrix(rotation_axis, angle)
    else:
        rot_matrix = np.eye(3) if dot(plane_normal, z_axis) > 0 else rotation_matrix(array([1, 0, 0]), np.pi)

    # Translate points so center_point goes to origin
    translated_points = points - center_point

    # Rotate points to align plane with XY-plane
    aligned_points = dot(translated_points, rot_matrix.T)

    # Drop the z value from the aligned points
    aligned_points = array([_[:2] for _ in aligned_points])

    return aligned_points, rot_matrix, center_point


def triangulate_with_constraints(points):
    # Create segment information for the polygon boundary
    segments = [(i, (i + 1) % len(points)) for i in range(len(points))]

    A = {'vertices': points, 'segments': segments}

    # Triangulate with options: 'p' for Planar Straight Line Graph, 'q' for quality mesh
    B = tr.triangulate(A, 'pq')

    tr.compare(plt, A, B)
    plt.show()

    return B['vertices'], B['triangles']


def reverse_transform_points(points_2d, rot_matrix, center_point):
    # Convert 2D points back to 3D by adding zero as the Z-component
    # Ensure only 3 dimensions (x, y, z where z = 0)
    points_3d = c_[points_2d, zeros(len(points_2d))]  # Now points_3d should be (n, 3)

    # Apply the inverse of the rotation
    inv_rot_matrix = linalg.inv(rot_matrix)
    # Ensure that matrix multiplication can proceed: shapes must align
    rotated_points = dot(points_3d, inv_rot_matrix.T)  # Now should work, as both are (n, 3) @ (3, 3)

    # Translate points back to the original position
    original_points = rotated_points + center_point

    return original_points


def project_points_to_plane(p0, p1, plane_point, plane_normal):

    # Direction vector of the ray
    direction = p1 - p0

    # Calculate t
    numerator = dot(plane_normal, (plane_point - p0))
    denominator = dot(plane_normal, direction)
    if denominator == 0:
        raise ValueError("The ray is parallel to the plane.")

    t = numerator / denominator

    # Calculate the intersection point
    intersection_point = p0 + t * direction
    return intersection_point


def build_surf(locs, rads, edge_points, plotting=False):
    """
    1. Calc_surf_func
    2. Build Perimeter
    3. Calculate plane
    4. Project perimeter onto plane using main point (smaller radius ball_center)
    5. Rotate and transform Points, so they are in the xy plane
    6. Make 2D Constrained Delaunay Triangulation using 2D perimeter Points
    7. Rotate and translate back to the
    7. Project these points through the plane onto the surface
    8. Return the solved surface
    """
    # Step 0: set up and check all inputs
    locs, rads = array(locs), array(rads)
    # Make sure the smaller of the two is ball 0
    if rads[0] > rads[1]:
        locs, rads = [locs[1], locs[0]], [rads[1], [rads[0]]]

    # Step 1: Calculate the surface function
    sfunc = calc_surf_func(locs[0], rads[0], locs[1], rads[1])

    # Step 2: Build the perimeter
    perimeter = build_perimeter1(edge_points)

    # Step 3: Calculate the plane between the two balls
    p_norm = (locs[1] - locs[0]) / linalg.norm(locs[1] - locs[0])
    p_point = locs[0] + 0.5 * (calc_dist_numba(locs[0], locs[1]) + rads[0] - rads[1]) * p_norm

    # Step 4: Project the perimeter points toward the ball 1 loc onto the plane
    plane_points = array([project_points_to_plane(_, locs[0], p_point, p_norm) for _ in perimeter])

    # Step 5: translate and rotate onto the x_y plane
    twoD_points, rot_mtx, ctr_pnt = transform_points(plane_points, p_norm, p_point)

    # Step 6: Fill in with Delaunay triangulation
    fill_2d_points = triangulate_with_constraints(twoD_points)

    # Step 7: Rotate and translate back to the plane
    new_plane_points = reverse_transform_points(twoD_points, rot_mtx, ctr_pnt)

    # Step 8: Project onto the surface
    surface_points = array([calc_surf_point(locs, _, sfunc) for _ in new_plane_points])


    # PLot if needed:
    if plotting:
        # Create the figure and the axis
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        # Plot the initial perimeter
        ax.scatter([_[0] for _ in perimeter], [_[1] for _ in perimeter], [_[2] for _ in perimeter])
        # Plot the plane points
        ax.scatter([_[0] for _ in plane_points], [_[1] for _ in plane_points], [_[2] for _ in plane_points])
        # Plot the rotated points
        ax.scatter([_[0] for _ in new_plane_points], [_[1] for _ in new_plane_points], [_[2] for _ in new_plane_points], marker='x')
        # Plot the rotated points
        ax.scatter([_[0] for _ in surface_points], [_[1] for _ in surface_points], [_[2] for _ in surface_points], marker='x')
        # Plot the plane between the two balls
        plot_surfs([rectangle_on_plane(p_point, p_norm, 10, 5)], [array([[0, 1, 2], [2, 3, 0]])],
                   fig=fig, ax=ax, alpha=0.5, simps=True, colors=['k'])

        # Plot the two balls
        plot_balls(locs, rads, fig=fig, ax=ax, colors=['k', 'k'], alpha=0.2)
        # Show everything
        plt.show()



######### Test Test Test Test Test Test Test Test Test Test Test Test Test Test Test Test Test Test Test Test ##########

if __name__ == '__main__':
    locs = [array([67.06, 55.12, 9.29]), array([65.29, 53.08, 9.49])]
    rads = [float64(1.5), float64(1.8)]
    epnts = [[array([67.20873939, 53.43244362, 8.8953191]), array([67.20822024, 53.43284059, 8.89546888]),
              array([67.01350033, 53.58296371, 8.94946707]), array([66.83848062, 53.7204175, 8.99353309]),
              array([66.67826622, 53.84892583, 9.02911298]), array([66.52897971, 53.9715037, 9.0572365]),
              array([66.38773625, 54.09042886, 9.07860851]), array([66.25227132, 54.20749363, 9.09377461]),
              array([66.12066289, 54.32420655, 9.1032204]), array([65.9911148, 54.44197393, 9.10740577]),
              array([65.86178581, 54.56226913, 9.10674586]), array([65.73064625, 54.68679343, 9.10156137]),
              array([65.59533483, 54.81764213, 9.09201957]), array([65.4529778, 54.95750599, 9.07807429]),
              array([65.29992124, 55.10995776, 9.05939806]), array([65.13130223, 55.27990382, 9.03528632]),
              array([65.03913856, 55.37351182, 9.02083172])],
             [array([65.10329148, 55.27870307, 8.56303136]), array([65.10370876, 55.27828959, 8.56325535]),
              array([65.26328532, 55.12083727, 8.64579738]), array([65.4127934, 54.97474465, 8.71655514]),
              array([65.55594055, 54.83646892, 8.77691374]), array([65.69592215, 54.70309072, 8.82744886]),
              array([65.83524145, 54.57248527, 8.86786263]), array([65.97578927, 54.44321884, 8.89714016]),
              array([66.11884144, 54.31449124, 8.91382496]), array([66.26508679, 54.18601296, 8.91647324]),
              array([66.41481565, 54.05772967, 8.90416386]), array([66.56830603, 53.92942756, 8.87676172]),
              array([66.72628086, 53.80037483, 8.83472587]), array([66.89027128, 53.66911784, 8.77858485]),
              array([67.06285852, 53.53339303, 8.70836657]), array([67.10776927, 53.4984138, 8.68852939])],
             [array([67.20771228, 53.43276231, 8.8873286]), array([67.20771629, 53.43276106, 8.88735985]),
              array([67.20873939, 53.43244362, 8.8953191]), array([67.20873939, 53.43244362, 8.8953191])],
             [array([67.20771228, 53.43276231, 8.8873286]), array([67.20732115, 53.43301179, 8.8865559]),
              array([67.10776927, 53.4984138, 8.68852939]), array([67.10776927, 53.4984138, 8.68852939])],
             [array([65.03913856, 55.37351182, 9.02083172]), array([65.03905101, 55.37352768, 9.01990397]),
              array([65.01688773, 55.37923023, 8.78203989]), array([65.01688773, 55.37923023, 8.78203989])],
             [array([65.10329148, 55.27870307, 8.56303136]), array([65.10301514, 55.27901676, 8.56372144]),
              array([65.03789273, 55.35433388, 8.72819172]), array([65.01688773, 55.37923023, 8.78203989])]]
    build_surf(locs, rads, epnts)