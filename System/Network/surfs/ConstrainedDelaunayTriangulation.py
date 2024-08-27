from shapely import Polygon, Point, LinearRing
from shapely.plotting import plot_polygon
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist


def plot_points_and_tris(pnts=None, trs=None, pcol=None, tcol=None, plot_points=True, Show=False):

    if trs is not None:
        for tri in trs:
            p0, p1, p2 = [pnts[_] for _ in tri]
            plt.plot([p0[0], p1[0], p2[0], p0[0]], [p0[1], p1[1], p2[1], p0[1]], c=tcol)
    if pnts is not None and plot_points:
        plt.scatter([_[0] for _ in pnts], [_[1] for _ in pnts], c=pcol)
    if Show:
        plt.show()


def generate_spiderweb(box, res, center=None, ring_scaler=None):
    # Set up the ring scaler variable
    if ring_scaler is None:
        ring_scaler = 1
    # Set up th minimum and maximum values
    min_x, max_x, min_y, max_y = box[0][0] - 2 * res, box[1][0] + 2 * res, box[0][1] - 2 * res, box[1][1] + 2 * res
    # Check if center is None
    if center is None:
        center = [min_x + 0.5 * (max_x - min_x), min_y + 0.5 * (max_y - min_y)]
    # Get the center points variable
    cx, cy = center
    # Get the corners
    corners = [box[0], [min_x, max_y], [max_x, min_y], box[1]]
    # Find the maximum possible radius based on the distance from the center to the corners
    max_radius = max([calc_dist(center, _) for _ in corners])
    # Get the number of rings based on the
    num_rings = int(max_radius / res) + 1
    # Create concentric circles of points
    points = [center]  # Start with the center point
    for i in range(1, num_rings + 1):
        # Create the new radius for the next ring
        radius = max_radius * (i / num_rings)
        # Increase the number of ring points as we go out
        num_points_per_ring = int(2 * np.pi * radius / res) + 1
        # Loop through the ring points adding if needed
        for j in range(num_points_per_ring):
            # Find the angle to place the next point
            angle = 2 * np.pi * j / num_points_per_ring
            # Get the x and y values
            x, y = cx + radius * np.cos(angle), cy + radius * np.sin(angle)
            # Check the location of the point and if it is outside the given box
            if min_x > x or x > max_x or min_y > y or y > max_y:
                continue
            points.append((x, y))

    # Convert list to numpy array for Delaunay triangulation
    points = np.array(points)

    return points


def sort_tris(perimeter, tris, polygon, points):
    """
    Sorts the triangles into different groups of inside and outside
    """
    # Set up the different triangles lists
    in_, out, mid, mid_designations = [], [], [], []
    # Loop through the triangles
    for tri in tris:
        # Create the list of true and false
        tri_points = [True if (polygon.contains(points[_]) or _ < len(perimeter)) else False for _ in tri]
        # If the entire triangle is contained add it to the inside triangles list
        if all(tri_points):
            in_.append(tri)
        # If half of the points are in add it to the mid_list and add the designations
        elif any(tri_points):
            good_tri_points = [i for i in range(3) if tri_points[i]]
            if len(good_tri_points) == 1 and tri[good_tri_points[0]] < len(perimeter):
                out.append(tri)
                # plot_polygon(polygon)
                # plot_points_and_tris([[_.x, _.y] for _ in points], [tri], plot_points=False, Show=True, tcol='k')
                continue

            mid.append(tri)
            mid_designations.append(tri_points)
        # Otherwise it is outside and we don't need it
        else:
            out.append(tri)
    # Return the lists
    return in_, out, mid, mid_designations


def find_shared_edge(triangle1, triangle2):
    """Helper function to find a shared edge between two triangles."""
    edges1 = {tuple(sorted((triangle1[i], triangle1[(i + 1) % 3]))) for i in range(3)}
    edges2 = {tuple(sorted((triangle2[i], triangle2[(i + 1) % 3]))) for i in range(3)}
    shared_edge = edges1.intersection(edges2)
    return shared_edge


def order_triangles(triangles, tri_designations):
    # Make sure the triangle list is actually included
    if not triangles:
        return []

    # Start with the first triangle
    ordered = [triangles.pop(0)]
    ordered_designations = [tri_designations[0]]

    # Try to find a matching triangle that shares an edge with the last in the ordered list
    while triangles:
        last_triangle = ordered[-1]
        found_next = False

        for i, triangle in enumerate(triangles):
            if find_shared_edge(last_triangle, triangle):
                ordered.append(triangle)
                ordered_designations.append(tri_designations[i])
                triangles.pop(i)
                found_next = True
                break

        if not found_next:
            # If no next triangle is found, the list may be broken or not all triangles are connected
            print("Warning: Not all triangles are connected in a loop.")
            break

    # Optionally, verify if the first and last triangles share an edge to close the loop
    if not find_shared_edge(ordered[0], ordered[-1]):
        print("Warning: The first and last triangles do not share an edge. Loop is not closed.")

    return ordered, ordered_designations


def normalize(v):
    """ Normalize a vector. """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v  # Handle the zero vector case gracefully
    return v / norm


def direction_toward_side(A, B, C):
    """
    Finds a normalized direction vector from vertex A towards side BC in 3D.

    :param A: Iterable (list or tuple), coordinates of vertex A.
    :param B: Iterable (list or tuple), coordinates of vertex B.
    :param C: Iterable (list or tuple), coordinates of vertex C.
    :return: Numpy array, normalized direction vector.
    """
    # Ensure inputs are numpy arrays and convert to 3D by adding a zero z-component if necessary
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    C = np.array(C, dtype=float)
    if A.shape[0] == 2:
        A = np.append(A, 0)
    if B.shape[0] == 2:
        B = np.append(B, 0)
    if C.shape[0] == 2:
        C = np.append(C, 0)

    # Vectors in the plane of the triangle
    vec_BC = B - C
    vec_BA = B - A

    # Normal to the triangle
    normal = np.cross(vec_BC, vec_BA)

    # Ensure normal is not a zero vector
    if np.linalg.norm(normal) == 0:
        raise ValueError("The points are collinear or invalid input for a normal vector.")

    # Direction vector from A towards BC, perpendicular to BC in the plane of the triangle
    direction = np.cross(normal, vec_BC)


    # plt.scatter([A[0], B[0], C[0]], )
    # Normalize the direction vector
    return normalize(direction[:2])


def find_intersection(p0, d, l1, l2):
    """
    Find intersection of a vector with a line segment.

    :param p0: np.array, starting point of the vector.
    :param d: np.array, direction vector.
    :param l1: np.array, one endpoint of the line segment.
    :param l2: np.array, other endpoint of the line segment.
    :return: np.array or None, the intersection point if it exists.
    """
    # Convert to numpy arrays for vector operations
    p0 = np.array(p0)
    d = np.array(d)
    l1 = np.array(l1)
    l2 = np.array(l2)

    # Vector from l1 to l2
    segment_vector = l2 - l1

    # Solve the system of linear equations to find t and u
    matrix = np.vstack([d, -segment_vector]).T
    try:
        result = np.linalg.solve(matrix, l1 - p0)
        t, u = result
    except np.linalg.LinAlgError:
        # The matrix is singular and the system cannot be solved
        return None

    # Check if the solution is within the bounds of the line segment
    if 0 <= u <= 1:
        # Intersection point is within the segment
        intersection = p0 + t * d
        return intersection
    else:
        # Intersection occurs outside the segment
        return None


def fit_edge_triangles(polygon, perimeter, triangles, tri_designations, points):
    """
    Move all triangles from the outside of the polygon to the inside
    Step 1: sort all triangles so they are next to their neigboring triangles
    """
    # Step 1: Put the triangles in order
    # ordered_triangles, ordered_desigs = order_triangles(triangles, tri_designations)

    # Make the points copy list
    points_copy = points.copy()

    # Move every other triangle's points
    for i, tri in enumerate(triangles):

        new_ordered_desigs, outside_indices, inside_indices = [], [], []
        for t in tri:
            my_point = Point(points[t])
            if polygon.contains(my_point) or t < len(perimeter):
                new_ordered_desigs.append(True)
                inside_indices.append(t)
            else:
                new_ordered_desigs.append(False)
                outside_indices.append(t)
        # # Check if there is one or two points outside of the polygon
        # # if len([0 for _ in new_ordered_desigs if _]) == 2:
        # #     continue
        # # Get the index of the point on the outside of the polygon
        # inside_indices = [tri[j] for j in range(3) if new_ordered_desigs[j]]
        # # Get the other two indices
        # outside_indices = [tri[j] for j in range(3) if not new_ordered_desigs[j]]
        # # Get the direction the outside point needs to move
        # move_dir = - direction_toward_side(points[inside_index], points[outside_indices[0]], points[outside_indices[1]])
        # Check the current perimeter segment to see if there is an overlap
        for outside_point_index in outside_indices:
            move_dir = points[inside_indices[0]] - points[outside_point_index]

            overlaps = []
            for p_seg in range(len(perimeter)):
                overlap = find_intersection(points[outside_point_index], move_dir,
                                            perimeter[p_seg % len(perimeter)], perimeter[(p_seg + 1) % len(perimeter)])
                if overlap is not None:
                    overlaps.append(overlap)

            # Find the closest overlap point to the previous point
            curr_dist, main_point = np.inf, None
            for overlap_point in overlaps:
                dist = calc_dist(overlap_point, points[outside_point_index])
                if dist < curr_dist:
                    curr_dist, main_point = dist, overlap_point

            # Assign the new point value
            if main_point is not None:
                points_copy[outside_point_index] = main_point

        # plot_polygon(polygon)
        # plot_points_and_tris(points, [tri], tcol='r', plot_points=False)
        #
        # plot_points_and_tris(points_copy, [tri], tcol='b', plot_points=False)
        # plt.plot([points[outside_point_index][0], points[outside_point_index][0] + move_dir[0]], [points[outside_point_index][1], points[outside_point_index][1] + move_dir[1]])
        #


    # Return the new points
    return points_copy


def filter_points_and_tris(points, triangles):
    # Find all unique indices used in triangles
    used_indices = set(idx for triangle in triangles for idx in triangle)

    # Map old indices to new indices
    index_map = {old_index: new_index for new_index, old_index in enumerate(sorted(used_indices))}

    # Create a new list of points that are actually used
    new_points = [points[idx] for idx in sorted(used_indices)]

    # Update triangles with new indices
    new_triangles = [(index_map[idx1], index_map[idx2], index_map[idx3]) for idx1, idx2, idx3 in triangles]

    return new_points, new_triangles


def triangulate_2D_Surface(perimeter, all_points=None, res=0.2, center=None):
    """
    takes in 2d perimeter points and returns an evenly filled and triangulated 2d surface
    1. Get the maximum and minimum possible x and y values for the perimeter
    2. Make a pre_triangulated grid based of a set resolution and these parameters
    3. Make the Polygon, LinearRing, and Point objects for the grid
    4. Record the points within the polygon and map their original indices to the new point_indices
    """
    if all_points is None:
        # Step 1: Get the maximum and minimum values for the perimeter with an additional cushion
        px, py = [_[0] for _ in perimeter], [_[1] for _ in perimeter]
        box = [[min(px), min(py)], [max(px), max(py)]]

        # Step 2: Create the grid for mapping to the surface with the given triangles
        grid_points = generate_spiderweb(box, res, center)

        # Step 3: Add the perimeter points to the grid points
        all_points = np.concatenate((perimeter, grid_points), axis=0)

    # Step 4: Triangulate the points
    triangles = Delaunay(all_points).simplices

    # Step 5: Make the points and polygon objects
    poly, my_points = Polygon(perimeter), [Point(_) for _ in all_points]

    # Step 6: Check the triangles for the ones that are in and the ones that are out
    in_tris, out_tris, mid_tris, mid_tri_designations = sort_tris(perimeter, triangles, poly, my_points)
    # plot_polygon(poly)
    # plot_points_and_tris(all_points, mid_tris, plot_points=False, Show=True)

    # Step 7: Move the points for the mid triangles that are on the outside of the polygon to the border
    new_points = fit_edge_triangles(poly, perimeter, mid_tris.copy(), mid_tri_designations, all_points)

    # plt.show()
    # Step 8: Filter out the points we don't want
    spoints, stris = filter_points_and_tris(new_points, mid_tris + in_tris)

    plot_polygon(poly)
    plot_points_and_tris(new_points, in_tris, tcol='b', plot_points=False, Show=False)
    plot_points_and_tris(new_points, mid_tris, tcol='r', plot_points=False, Show=True)
    # plot_points_and_tris(spoints, stris, tcol='r', plot_points=True, Show=True)

    # plot_points_and_tris(new_points, in_tris, tcol='k', plot_points=False, Show=True)
    return spoints, stris




if __name__ == '__main__':
    # Real perimeter points we are testing
    # perimeter = np.array([[ 0.63663669, -1.22153212], [ 0.63619087, -1.22091145], [ 0.47012759, -0.99400362], [ 0.32167774, -0.7997529 ], [ 0.18484762, -0.62957204], [ 0.05502084, -0.47715853], [-0.07129593, -0.33795481], [-0.19689486, -0.2084431 ], [-0.32416444, -0.08569257], [-0.45534415, 0.03294308], [-0.59277716, 0.15004052], [-0.73918503, 0.26828829], [-0.89799568, 0.3906871 ], [-1.07378565,  0.52082571], [-1.27295209,  0.66331013], [-1.50483439,  0.82448172], [-1.63723295,  0.91482777], [-1.63723295,  0.91482777], [-1.63799014,  0.91414042], [-1.83616539,  0.73860756], [-1.83616539,  0.73860756], [-1.83616539,  0.73860756], [-1.84656425,  0.67033029], [-1.87993191,  0.46297227], [-1.88007721,  0.46210507], [-1.88007721,  0.46210507], [-1.87934367,  0.461845  ], [-1.60765386,  0.36232849], [-1.36896132,  0.26817206], [-1.15470852,  0.17596575], [-0.95871494,  0.08255767], [-0.77693182, -0.0149893 ], [-0.60685917, -0.11951755], [-0.44708807, -0.23372636], [-0.29678835, -0.36000575], [-0.1551099 , -0.50031552], [-0.02065755, -0.65633335], [ 0.10876597, -0.82995915], [ 0.23616031, -1.02405134], [ 0.36528648, -1.24326412], [ 0.39836064, -1.30278741], [ 0.39836064, -1.30278741], [ 0.39836064, -1.30278741], [ 0.62952189, -1.22766827], [ 0.63042709, -1.2273798 ], [ 0.63042709, -1.2273798 ], [ 0.63045137, -1.22735692], [ 0.63663669, -1.22153212], [ 0.63663669, -1.22153212]])
    perimeter = np.array([[-0.32031466, -2.6643524 ], [-0.3205994 , -2.66304939], [-0.42536602, -2.21684867], [-0.52055895, -1.87382853], [-0.61456839, -1.59510605], [-0.7142733 , -1.35898797], [-0.82636876, -1.15315994], [-0.95813227, -0.97030988], [-1.11784694, -0.80535512], [-1.3155826 , -0.65327181], [-1.5655073 , -0.50748314], [-1.71682049, -0.43416263], [-1.71682049, -0.43416263], [-1.71767751, -0.43511787], [-1.95447763, -0.69418377], [-2.22736201, -0.98407975], [-2.22736201, -0.98407975], [-2.22736201, -0.98407975], [-2.30481144, -1.08777617], [-2.57414254, -1.42852264], [-2.82556298, -1.7293156 ], [-3.12569491, -2.07584162], [-3.61492918, -2.62377672], [-3.61668081, -2.62571273], [-3.61668081, -2.62571273], [-3.38520617, -2.82632159], [-3.12411656, -3.09746262], [-3.12313616, -3.09859058], [-3.12313616, -3.09859058], [-2.9499568 , -3.26150044], [-2.71828124, -3.51081639], [-2.52042276, -3.75737855], [-2.26048103, -4.1366458 ], [-2.25961652, -4.13802131], [-2.25961652, -4.13802131], [-2.25961652, -4.13802131], [-2.1989681 , -4.16364231], [-2.19873135, -4.16374444], [-2.19873135, -4.16374444], [-2.19873135, -4.16374444], [-1.7982394 , -4.27739202], [-1.79667726, -4.27788227], [-1.79667726, -4.27788227], [-1.50951021, -3.92630864], [-1.16060074, -3.51540185], [-0.91051971, -3.23654014], [-0.66390494, -2.97968507], [-0.32148554, -2.66532207], [-0.32031466, -2.6643524 ]])
    minx, miny, maxx, maxy = min([_[0] for _ in perimeter]), min([_[1] for _ in perimeter]), max([_[0] for _ in perimeter]), max([_[1] for _ in perimeter])
    points = np.concatenate((perimeter, generate_spiderweb([[minx, miny], [maxx, maxy]], res=0.2, center=[0, 0])), axis=0)
    tris = Delaunay(points).simplices
    plot_points_and_tris(pnts=points, trs=tris, plot_points=False, tcol='r')
    poly = Polygon(perimeter)
    plot_polygon(Polygon(perimeter), add_points=False)
    plt.show()

    in_tris, out_tris, mid_tris, mid_tri_designations = sort_tris(perimeter, tris, Polygon(perimeter), [Point(_) for _ in points])
    plot_points_and_tris(pnts=points, trs=in_tris, plot_points=False, tcol='g')
    plot_points_and_tris(pnts=points, trs=out_tris, plot_points=False, tcol='r')
    plot_points_and_tris(pnts=points, trs=mid_tris, plot_points=False, tcol='y')
    plot_polygon(Polygon(perimeter), add_points=False)
    plt.show()


    # Step 7: Move the points for the mid triangles that are on the outside of the polygon to the border
    new_points = fit_edge_triangles(poly, perimeter, mid_tris.copy(), mid_tri_designations, points)



    plot_points_and_tris(pnts=new_points, trs=in_tris, plot_points=False, tcol='g')
    plot_points_and_tris(pnts=new_points, trs=out_tris, plot_points=False, tcol='r')
    plot_points_and_tris(pnts=new_points, trs=mid_tris, plot_points=False, tcol='y')
    plot_polygon(Polygon(perimeter), add_points=False)
    plt.show()


    plot_polygon(Polygon(perimeter), add_points=False)
    plot_points_and_tris(new_points, in_tris + mid_tris, plot_points=False, Show=True, tcol='b')



    # triangulate_2D_Surface(perimeter, 0.1)
    #
    # # Test the spider web
    # points = generate_spiderweb([5, 5], [[-10, -10], [10, 10]], 2)
    # tris = Delaunay(points).simplices
    # plot_points_and_tris(points, tris, pcol='k', tcol='b', Show=True)




############################################### Junk Junk Junk Junk ####################################################


# def create_equilateral_triangle_grid(min_x, max_x, min_y, max_y, res, plotting=True):
#
#     # Calculate the number of columns and rows based on the bounding box and resolution
#     num_cols = int((max_x - min_x) / res) + 1
#     num_rows = int((max_y - min_y) / res) + 1
#
#     # Generate points in a staggered grid
#     points = []
#     for i in range(num_rows):
#         for j in range(num_cols):
#             # Staggering even and odd rows
#             offset_x = res / 2 if i % 2 == 1 else 0
#             x = min_x + j * res + offset_x
#             y = min_y + i * res
#             points.append((x, y))
#
#     # Triangulate the points into equilateral triangles
#     # Creating two triangles for each rectangular cell in the grid, except for possibly the last row/column if incomplete
#     triangles = []
#     for i in range(num_rows - 1):
#         for j in range(num_cols - 1):
#             base_index = i * num_cols + j
#             if i % 2 == 0:
#                 # For even rows: downward pointing triangle
#                 triangles.append((base_index, base_index + num_cols, base_index + 1))
#                 if j < num_cols - 2:
#                     # Upward pointing triangle, skip this on the last column of even rows
#                     triangles.append((base_index + 1, base_index + num_cols, base_index + num_cols + 1))
#             else:
#                 # For odd rows: upward pointing triangle
#                 triangles.append((base_index, base_index + num_cols + 1, base_index + num_cols))
#                 if j < num_cols - 2:
#                     # Downward pointing triangle, skip this on the last column of odd rows
#                     triangles.append((base_index + 1, base_index + num_cols + 1, base_index + num_cols + 2))
#     if plotting:
#         plt.scatter([_[0] for _ in points], [_[1] for _ in points])
#         for tri in triangles:
#             p1, p2, p3 = [points[_] for _ in tri]
#             plt.plot([p1[0], p2[0], p3[0], p1[0]], [p1[1], p2[1], p3[1], p1[1]])
#         plt.show()
#     return np.array(points), np.array(triangles)
#
#
# def fill_points(perimeter, resolution):
#     min_x, max_x, min_y, max_y = min([_[0] for _ in perimeter]), max([_[0] for _ in perimeter]), min([_[1] for _ in perimeter]), max([_[1] for _ in perimeter])
#     my_range = max_x - min_x, max_y - min_y
#     my_poly = Polygon(perimeter)
#     my_lin_rin = LinearRing(perimeter)
#     grid_points = []
#     for i in range(int(my_range[0] / resolution) + 1):
#         for j in range(int(my_range[1] / resolution) + 1):
#
#             my_p = Point((min_x + i * resolution, min_y + j * resolution))
#             grid_points.append(my_p)
#     good_indices = []
#     good_points = []
#     for i, point in enumerate(grid_points):
#         if my_poly.contains(point) and point.distance(my_lin_rin) > 0.5 * resolution:
#             good_indices.append(i)
#             good_points.append([point.x, point.y])
#     plt.scatter([_[0] for _ in good_points], [_[1] for _ in good_points])
#     plot_polygon(my_poly)
#     plt.show()