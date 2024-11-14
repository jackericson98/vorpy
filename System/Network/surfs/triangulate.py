from shapely import Polygon, Point
from shapely.plotting import plot_polygon
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist, calc_tri, calc_com, project_to_plane
from System.sys_funcs.calcs.surf import calc_2d_surf_sa
from scipy.spatial._qhull import QhullError


def plot_points_and_tris(pnts=None, trs=None, pcol=None, tcol=None, plot_points=True, Show=False):

    if trs is not None:
        for tri in trs:
            p0, p1, p2 = [pnts[_] for _ in tri]
            plt.plot([p0[0], p1[0], p2[0], p0[0]], [p0[1], p1[1], p2[1], p0[1]], c=tcol)
    if pnts is not None and plot_points:
        plt.scatter([_[0] for _ in pnts], [_[1] for _ in pnts], c=pcol)
    plt.grid(False)
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


def is_within(perimeter, point, surf_loc, surf_norm):
    # First see if the perimeter is a list
    if type(perimeter) is list:
        # Check if we need to project to the plane
        if len(perimeter[0]) == 3:
            # We need to map to plane
            perimeter = project_to_plane(perimeter, surf_loc, surf_norm)
            # Same with the point
            point = project_to_plane([point], surf_loc, surf_norm)[0]
        # Set the shapely objects
        try:
            perimeter, point = Polygon(perimeter), Point(point)
        except TypeError:
            return False
    # Return the result
    return perimeter.contains(point)


def sort_tris(perimeter, tris, polygon, points):
    """
    Sorts the triangles into different groups of inside and outside
    """
    # Set up the different triangles lists
    in_, out, mid = [], [], []
    # Point Designation dictionary
    point_desigs = {}
    # Loop through the points
    for i, point in enumerate(points):
        # Check if the point is on the perimeter
        if i < len(perimeter):
            # Assign as an edge point
            point_desigs[i] = 'e'
        # Check if the point is inside
        elif polygon.contains(point):
            # Assign as inside
            point_desigs[i] = 'i'
        # Otherwise it is outside
        else:
            # Assign the outside
            point_desigs[i] = 'o'

    # Loop through the triangles
    for tri in tris:
        # Create the list of designations
        tri_point_desigs = [point_desigs[_] for _ in tri]

        # First check that all 3 points are within the perimeter
        if tri_point_desigs == ['i', 'i', 'i']:
            # Add to the in list
            in_.append(tri)
        # If all three points are edges this will need to be checked
        elif tri_point_desigs == ['e', 'e', 'e']:
            # Check the center of mass
            com = calc_com([[points[_].x, points[_].y] for _ in tri])
            # Check if the polygon contains this center of mass
            if polygon.contains(Point(com)):
                # Add to the ins
                in_.append(tri)
            # Otherwise add to the out list
            else:
                # Add to the outs
                out.append(tri)
        # Contain at least 1 point inside the perimeter means middle
        elif 'i' in tri_point_desigs:
            # Check if the other two are sides
            if 'o' in tri_point_desigs:
                mid.append(tri)
            else:
                # Add to the mids
                in_.append(tri)
        # Last if there is a mixture of edges and outs it is an out triangle
        else:
            # Check the center of mass
            com = calc_com([[points[_].x, points[_].y] for _ in tri])
            # Check if the polygon contains this center of mass
            if polygon.contains(Point(com)):
                # Add to the ins
                mid.append(tri)
            # Otherwise add to the out list
            else:
                # Add to the outs
                out.append(tri)
    # Return the lists
    return in_, out, mid, point_desigs


def reassign_tri_points(perimeter, mid_tris, polygon, points):
    """
    Three goals: 1. Filter out mid triangles that go to zero, move the points of the mid triangles that lie outside the
    perimeter to the closest perimeter point, 3.
    """
    point_mapping = {}
    new_tris = []
    # Loop through each triangle and assign the
    for tri in mid_tris:
        # Create the new triangle variable that will store the new triangle indices
        new_tri = []
        # Go through the points in the triangle and reassign the triangle indices to the closest perimeter point
        for tri_point in tri:
            # Check the index first:
            if tri_point in point_mapping:
                # Assign the triangle point to the new point mapping
                new_tri.append(point_mapping[tri_point])
            # Next check if the point is inside the polygon (we have to check for if it is in the perimeter as well
            elif polygon.contains(Point(points[tri_point])) or tri_point < len(perimeter):
                # Add the same value to the tri point so that it does not move
                new_tri.append(tri_point)
                # Add the tri_point mapping to the dictionary so we can easily loop next time it comes up
                point_mapping[tri_point] = tri_point
            # Now if the point has not been found and it is outside the perimeter, we need to assign a perimeter point
            else:
                # Create the distance and perimeter closest point variables
                close_point, dist, my_point = None, np.inf, points[tri_point]
                # First loop through the perimeter points so we can test for closeness
                for i, perim_point in enumerate(perimeter):
                    # Calculate the distance to the perimeter point
                    new_dist = calc_dist(perim_point, my_point)
                    # Check if it is closer than our current point
                    if new_dist < dist:
                        # Assign the new value
                        dist, close_point = new_dist, i
                # Check that the perimeter point is not None
                if close_point is not None:
                    # Assign the new triangle point
                    new_tri.append(close_point)
                    # Assign the point mapping
                    point_mapping[tri_point] = close_point
        # Check that triangle is something we actually want and/or is complete
        if len(new_tri) < 3:
            continue
        # Calculate the area of the triangle
        if round(calc_tri(np.array([points[_] for _ in new_tri])), 10) > 0:
            # Check to make sure it is actually inside the polygon
            if polygon.contains(Point(calc_com([points[_] for _ in new_tri]))):
                # Add it to the triangles for return
                new_tris.append(new_tri)
    # Return the new set opf triangles
    return new_tris


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


def triangulate_2D_Surface(perimeter, all_points=None, res=0.2, center=None, timer=False, plotting=False,
                           filter_hard=False):
    """
    takes in 2d perimeter points and returns an evenly filled and triangulated 2d surface
    1. Get the maximum and minimum possible x and y values for the perimeter
    2. Make a pre_triangulated grid based of a set resolution and these parameters
    3. Make the Polygon, LinearRing, and Point objects for the grid
    4. Record the points within the polygon and map their original indices to the new point_indices
    """
    # start = time.perf_counter()
    # If no points were provided we need to map the points to the surface
    if all_points is None:

        # Step 1: Get the maximum and minimum values for the perimeter with an additional cushion
        px, py = [_[0] for _ in perimeter], [_[1] for _ in perimeter]
        box = [[min(px), min(py)], [max(px), max(py)]]

        # Step 2: Create the grid for mapping to the surface with the given triangles
        grid_points = generate_spiderweb(box, res, center)

        # Step 3: Add the perimeter points to the grid points
        all_points = np.concatenate((perimeter, grid_points), axis=0)

    # if timer:
    #     spider_time = time.perf_counter() - start
    #     start = time.perf_counter()
    # Step 4: Triangulate the points
    # First Check to see if we need to filter really specifically for points outside the perimeter

    poly, my_points = Polygon(perimeter), [Point(_) for _ in all_points]
    if filter_hard:
        re_ass_dict = {}
        good_points = []
        good_point_points = []
        skip_num = 0
        for i, point in enumerate(my_points):
            if not poly.contains(point) and i >= len(perimeter):
                skip_num += 1
            else:
                re_ass_dict[i - skip_num] = i
                good_points.append(all_points[i])
                good_point_points.append(point)

    else:
        good_points, good_point_points = all_points, my_points
    try:
        triangles = Delaunay(good_points).simplices
    except QhullError as e:
        try:
            triangles = Delaunay(good_points, qhull_options='QJ').simplices
        except QhullError as e2:
            return all_points, []

    if filter_hard:
        new_triangles = []
        for tri in triangles:
            new_triangles.append([re_ass_dict[_] for _ in tri])
        triangles = new_triangles
    # if timer:
    #     delaunay_time = time.perf_counter() - start
    #     start = time.perf_counter()
    # Step 5: Make the points and polygon objects
    # if plotting:
    #     plot_points_and_tris(all_points, triangles, tcol='r', plot_points=True, Show=False)
    #     plot_polygon(poly)
    #     plt.show()
    # if timer:
    #     make_polygon_time = time.perf_counter() - start
    #     start = time.perf_counter()
    # Step 6: Check the triangles for the ones that are in and the ones that are out
    in_tris, out_tris, mid_tris, mid_tri_designations = sort_tris(perimeter, triangles, poly, my_points)

    # if plotting:
    #     plot_points_and_tris(all_points, in_tris, tcol='g')
    #     plot_points_and_tris(all_points, mid_tris, tcol='y')
    #     plot_polygon(poly)
    #     plt.show()

    # if timer:
    #     tri_desig_time = time.perf_counter() - start
    #     start = time.perf_counter()
    # If the points were provided
    mid_tris = reassign_tri_points(perimeter, mid_tris, poly, all_points)

    # if timer:
    #     tri_reassignment = time.perf_counter() - start
    #
    # spoints, stris = filter_points_and_tris(all_points, in_tris + mid_tris)
    # if plotting:
    #     plot_polygon(poly)
    #     plot_points_and_tris(all_points, in_tris, tcol='b', plot_points=False, Show=False)
    #     plot_points_and_tris(all_points, mid_tris, tcol='r', plot_points=False, Show=True)
    # plot_points_and_tris(spoints, stris, tcol='r', plot_points=True, Show=True)
    # if timer:
    #     time_vals = {'spider': spider_time, 'Delaunay': delaunay_time, 'designations': tri_desig_time, 'reassign': tri_reassignment}
    # plot_points_and_tris(new_points, in_tris, tcol='k', plot_points=False, Show=True)
    # Final check to see if we have triangulated correctly
    my_area, poly_area = calc_2d_surf_sa(in_tris + mid_tris, all_points), poly.area
    if round(my_area, 4) != round(poly_area, 4) and not filter_hard:
        return triangulate_2D_Surface(perimeter, all_points, res, center, timer, False, filter_hard=True)
    return all_points, in_tris + mid_tris


def triangulate_2D_Surface1(perimeter, res=0.2, center=None, filter_hard=False):
    """
    takes in 2d perimeter points and returns an evenly filled and triangulated 2d surface
    1. Get the maximum and minimum possible x and y values for the perimeter
    2. Make a pre_triangulated grid based of a set resolution and these parameters
    3. Make the Polygon, LinearRing, and Point objects for the grid
    4. Record the points within the polygon and map their original indices to the new point_indices
    """

    # Step 1: Get the maximum and minimum values for the perimeter with an additional cushion
    px, py = [_[0] for _ in perimeter], [_[1] for _ in perimeter]
    box = [[min(px), min(py)], [max(px), max(py)]]

    # Step 2: Create the grid for mapping to the surface with the given triangles
    grid_points = generate_spiderweb(box, res, center)

    # Step 3: Set up the shapely objects and test for insideness
    poly, all_ppoints = Polygon(perimeter), [Point(_) for _ in perimeter]
    # Create a list of all points
    all_points = perimeter.copy()
    # Loop through the grid points
    for point in grid_points:
        # Create the shapely point object
        test_point = Point(point)
        # Check for insideness of the point and add the objects if it is
        if poly.contains(test_point):
            all_points.append(point)
            all_ppoints.append(test_point)


    # Step 5: Create the triangulation of the points
    try:
        triangles = Delaunay(all_points).simplices
    except QhullError as e:
        try:
            triangles = Delaunay(all_points, qhull_options='QJ').simplices
        except QhullError as e2:
            return all_points, []

    # Step 6: Sort the triangles and reassign the points
    in_tris, out_tris, mid_tris, mid_tri_designations = sort_tris(perimeter, triangles, poly, all_ppoints)
    # plot_polygon(poly)
    # plot_points_and_tris(all_points, mid_tris, tcol='y', plot_points=False)
    # plot_points_and_tris(all_points, out_tris, tcol='r', plot_points=False)
    # plot_points_and_tris(all_points, in_tris, tcol='g', plot_points=False, Show=True)

    mid_tris = reassign_tri_points(perimeter, mid_tris, poly, all_points)

    # Step 7: Return the values
    return all_points, in_tris + mid_tris



if __name__ == '__main__':
    # Real perimeter points we are testing
    perimeter = np.array([[ 0.63663669, -1.22153212], [ 0.63619087, -1.22091145], [ 0.47012759, -0.99400362], [ 0.32167774, -0.7997529 ], [ 0.18484762, -0.62957204], [ 0.05502084, -0.47715853], [-0.07129593, -0.33795481], [-0.19689486, -0.2084431 ], [-0.32416444, -0.08569257], [-0.45534415, 0.03294308], [-0.59277716, 0.15004052], [-0.73918503, 0.26828829], [-0.89799568, 0.3906871 ], [-1.07378565,  0.52082571], [-1.27295209,  0.66331013], [-1.50483439,  0.82448172], [-1.63723295,  0.91482777], [-1.63723295,  0.91482777], [-1.63799014,  0.91414042], [-1.83616539,  0.73860756], [-1.83616539,  0.73860756], [-1.83616539,  0.73860756], [-1.84656425,  0.67033029], [-1.87993191,  0.46297227], [-1.88007721,  0.46210507], [-1.88007721,  0.46210507], [-1.87934367,  0.461845  ], [-1.60765386,  0.36232849], [-1.36896132,  0.26817206], [-1.15470852,  0.17596575], [-0.95871494,  0.08255767], [-0.77693182, -0.0149893 ], [-0.60685917, -0.11951755], [-0.44708807, -0.23372636], [-0.29678835, -0.36000575], [-0.1551099 , -0.50031552], [-0.02065755, -0.65633335], [ 0.10876597, -0.82995915], [ 0.23616031, -1.02405134], [ 0.36528648, -1.24326412], [ 0.39836064, -1.30278741], [ 0.39836064, -1.30278741], [ 0.39836064, -1.30278741], [ 0.62952189, -1.22766827], [ 0.63042709, -1.2273798 ], [ 0.63042709, -1.2273798 ], [ 0.63045137, -1.22735692], [ 0.63663669, -1.22153212], [ 0.63663669, -1.22153212]])
    # perimeter = np.array([[-0.32031466, -2.6643524 ], [-0.3205994 , -2.66304939], [-0.42536602, -2.21684867], [-0.52055895, -1.87382853], [-0.61456839, -1.59510605], [-0.7142733 , -1.35898797], [-0.82636876, -1.15315994], [-0.95813227, -0.97030988], [-1.11784694, -0.80535512], [-1.3155826 , -0.65327181], [-1.5655073 , -0.50748314], [-1.71682049, -0.43416263], [-1.71682049, -0.43416263], [-1.71767751, -0.43511787], [-1.95447763, -0.69418377], [-2.22736201, -0.98407975], [-2.22736201, -0.98407975], [-2.22736201, -0.98407975], [-2.30481144, -1.08777617], [-2.57414254, -1.42852264], [-2.82556298, -1.7293156 ], [-3.12569491, -2.07584162], [-3.61492918, -2.62377672], [-3.61668081, -2.62571273], [-3.61668081, -2.62571273], [-3.38520617, -2.82632159], [-3.12411656, -3.09746262], [-3.12313616, -3.09859058], [-3.12313616, -3.09859058], [-2.9499568 , -3.26150044], [-2.71828124, -3.51081639], [-2.52042276, -3.75737855], [-2.26048103, -4.1366458 ], [-2.25961652, -4.13802131], [-2.25961652, -4.13802131], [-2.25961652, -4.13802131], [-2.1989681 , -4.16364231], [-2.19873135, -4.16374444], [-2.19873135, -4.16374444], [-2.19873135, -4.16374444], [-1.7982394 , -4.27739202], [-1.79667726, -4.27788227], [-1.79667726, -4.27788227], [-1.50951021, -3.92630864], [-1.16060074, -3.51540185], [-0.91051971, -3.23654014], [-0.66390494, -2.97968507], [-0.32148554, -2.66532207], [-0.32031466, -2.6643524 ]])
    minx, miny, maxx, maxy = min([_[0] for _ in perimeter]), min([_[1] for _ in perimeter]), max([_[0] for _ in perimeter]), max([_[1] for _ in perimeter])
    points = np.concatenate((perimeter, generate_spiderweb([[minx, miny], [maxx, maxy]], res=0.1, center=[0, 0])), axis=0)
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
    mid_tris = reassign_tri_points(perimeter, mid_tris, poly, points)



    plot_points_and_tris(pnts=points, trs=in_tris, plot_points=False, tcol='g')
    plot_points_and_tris(pnts=points, trs=out_tris, plot_points=False, tcol='r')
    plot_points_and_tris(pnts=points, trs=mid_tris, plot_points=False, tcol='y')
    plot_polygon(Polygon(perimeter), add_points=False)
    plt.show()


    plot_polygon(Polygon(perimeter), add_points=False)
    plot_points_and_tris(points, in_tris + mid_tris, plot_points=False, Show=True, tcol='b')



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