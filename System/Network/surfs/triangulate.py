from shapely import Polygon, Point, LineString
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
        else:
            # Assign as inside
            point_desigs[i] = 'i'

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
        else:
            # Add to the mids
            in_.append(tri)
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


def triangulate_2D_Surface(perimeter, res=0.2, center=None):
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
    poly, linestring, all_ppoints = Polygon(perimeter), LineString(perimeter), [Point(_) for _ in perimeter]
    # Check for points close to the edge
    buffer = linestring.buffer(res / 2)
    # Create a list of all points
    all_points = perimeter.copy()
    # Loop through the grid points
    for point in grid_points:
        # Create the shapely point object
        test_point = Point(point)
        # Check for insideness of the point and add the objects if it is
        if poly.contains(test_point):
            # Remove any of the points that are too close to the perimeter to prevent bad triangles
            if not buffer.contains(test_point):
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

    # Step 7: Check if the mid tris exist and hard fix them
    if len(mid_tris) > 0:
        mid_tris = reassign_tri_points(perimeter, mid_tris, poly, all_points)

    # Step 8: Verify the surface
    # my_sa = calc_2d_surf_sa(in_tris + mid_tris, all_points)
    # poly_sa = poly.area
    # if round(my_sa, 5) != round(poly_sa, 5):
    #     print(my_sa, poly_sa)
    #     plot_polygon(poly, add_points=False)
    #     plot_points_and_tris(all_points, in_tris + mid_tris, tcol='grey', plot_points=False, Show=True)

    # Step 7: Return the values
    return all_points, in_tris + mid_tris
