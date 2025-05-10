from shapely import Polygon, Point, LineString
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import calc_tri
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import project_to_plane
from scipy.spatial._qhull import QhullError


def plot_points_and_tris(pnts=None, trs=None, pcol=None, tcol=None, plot_points=True, Show=False):
    """
    Plots points and triangles in a 2D space with customizable colors and display options.

    This function visualizes a set of points and their triangulation by:
    1. Drawing triangles between points if provided
    2. Scattering points if enabled
    3. Supporting custom colors for both points and triangles
    4. Offering display control through the Show parameter

    Parameters
    ----------
    pnts : list of numpy.ndarray, optional
        List of 2D point coordinates to plot
    trs : list of tuples, optional
        List of triangle indices referencing points in pnts
    pcol : str or list, optional
        Color specification for points
    tcol : str or list, optional
        Color specification for triangles
    plot_points : bool, optional
        Whether to display the points
    Show : bool, optional
        Whether to immediately display the plot

    Notes
    -----
    - Points are plotted as scatter points
    - Triangles are drawn as closed polygons
    - Grid lines are disabled by default
    - Supports both single color and color list specifications
    """

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
    """
    Generates a spiderweb-like pattern of points for surface triangulation.

    This function creates a set of points arranged in concentric circles (rings) around a center point,
    with the number of points per ring increasing with radius to maintain approximately uniform spacing.
    The points are constrained to lie within a specified bounding box.

    Parameters
    ----------
    box : list of list
        Bounding box defined as [[min_x, min_y], [max_x, max_y]]
    res : float
        Approximate resolution (spacing) between points
    center : list, optional
        Center point [x, y] for the spiderweb pattern
    ring_scaler : float, optional
        Scaling factor for ring spacing

    Returns
    -------
    numpy.ndarray
        Array of 2D points arranged in a spiderweb pattern

    Notes
    -----
    - Points are arranged in concentric circles around the center
    - Number of points per ring increases with radius to maintain uniform spacing
    - Points outside the bounding box are excluded
    - The pattern starts with a single point at the center
    - Ring spacing is determined by the resolution parameter
    """
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
    """
    Determines if a point lies within a given perimeter using geometric containment checks.

    This function checks whether a point is contained within a perimeter by:
    1. Handling both 2D and 3D perimeters by projecting 3D points to a plane
    2. Converting the perimeter and point to Shapely geometric objects
    3. Using Shapely's contains() method to determine containment

    Parameters
    ----------
    perimeter : list or numpy.ndarray
        List of points defining the perimeter boundary
    point : numpy.ndarray
        The point to check for containment
    surf_loc : numpy.ndarray
        Surface location for 3D projection
    surf_norm : numpy.ndarray
        Surface normal vector for 3D projection

    Returns
    -------
    bool
        True if the point is within the perimeter, False otherwise

    Notes
    -----
    - Handles both 2D and 3D perimeters by projecting to a plane when needed
    - Uses Shapely's geometric operations for robust containment checking
    - Returns False if the perimeter or point cannot be converted to valid geometric objects
    """
    # First see if the perimeter is a list
    if type(perimeter) is list or isinstance(perimeter, np.ndarray):
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
    Sorts triangles into groups based on their position relative to a perimeter.

    This function categorizes triangles into three groups:
    1. Inside triangles (completely within the perimeter)
    2. Outside triangles (completely outside the perimeter)
    3. Middle triangles (spanning the perimeter boundary)

    Parameters
    ----------
    perimeter : list of numpy.ndarray
        List of points defining the perimeter boundary
    tris : list of list of int
        List of triangles, where each triangle is a list of three point indices
    polygon : shapely.geometry.Polygon
        Shapely polygon object representing the perimeter
    points : list of shapely.geometry.Point
        List of Shapely point objects representing all points in the triangulation

    Returns
    -------
    tuple
        A tuple containing:
        - in_ : list of list of int
            List of triangles completely inside the perimeter
        - out : list of list of int
            List of triangles completely outside the perimeter
        - mid : list of list of int
            List of triangles spanning the perimeter boundary
        - point_desigs : dict
            Dictionary mapping point indices to their designation ('e' for edge, 'i' for inside)

    Notes
    -----
    - Triangles are categorized based on their vertices' positions relative to the perimeter
    - Edge triangles (all vertices on perimeter) are checked using their center of mass
    - Middle triangles (spanning boundary) are included in the inside group
    - Point designations are used to track which points are on the perimeter vs inside
    """
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
    Reassigns points in middle triangles to create valid triangles within the perimeter.

    This function processes middle triangles (those with points both inside and outside the perimeter) by:
    1. Filtering out invalid triangles with zero area
    2. Moving points outside the perimeter to their closest perimeter points
    3. Ensuring all resulting triangles are valid and contained within the perimeter

    Parameters
    ----------
    perimeter : list of numpy.ndarray
        List of points defining the perimeter boundary
    mid_tris : list of list of int
        List of middle triangles, where each triangle is a list of 3 point indices
    polygon : shapely.geometry.Polygon
        Polygon object representing the perimeter boundary
    points : list of numpy.ndarray
        List of all points in the triangulation

    Returns
    -------
    list of list of int
        List of valid triangles after point reassignment, where each triangle is a list of 3 point indices

    Notes
    -----
    - Maintains a mapping of original points to their new perimeter point assignments
    - Only includes triangles that have positive area and are fully contained within the perimeter
    - Preserves points that are already on the perimeter or inside the polygon
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
    """
    Filters points and triangles to remove unused points and reindex the remaining points.

    This function:
    1. Identifies all points that are actually used in the triangles
    2. Creates a mapping from old indices to new indices
    3. Creates a new list of points containing only the used points
    4. Updates triangle indices to reference the new point list

    Parameters
    ----------
    points : list of numpy.ndarray
        List of point coordinates in 2D space
    triangles : list of tuple
        List of triangles, where each triangle is a tuple of three point indices

    Returns
    -------
    tuple
        A tuple containing:
        - new_points : list of numpy.ndarray
            List of filtered point coordinates
        - new_triangles : list of tuple
            List of triangles with updated indices

    Notes
    -----
    - Removes any points that are not referenced by any triangle
    - Maintains the relative ordering of points in the new list
    - Preserves triangle connectivity while updating indices
    - Useful for cleaning up triangulation results
    """
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
    Triangulates a 2D surface defined by a perimeter of points.

    This function creates a triangulated mesh of a 2D surface by:
    1. Determining the bounding box of the perimeter points
    2. Generating a uniform grid of points within the bounding box
    3. Creating geometric objects (Polygon, LineString, Point) for spatial analysis
    4. Filtering grid points to only those inside the perimeter
    5. Performing Delaunay triangulation on the filtered points
    6. Sorting and validating triangles to ensure proper surface coverage

    Parameters
    ----------
    perimeter : list of numpy.ndarray
        List of 2D points defining the perimeter of the surface
    res : float, optional
        Resolution for grid point generation, defaults to 0.2
    center : numpy.ndarray, optional
        Center point for grid generation, defaults to None

    Returns
    -------
    tuple
        A tuple containing:
        - all_points : list of numpy.ndarray
            List of all points used in the triangulation
        - triangles : list of tuple
            List of triangles as tuples of point indices

    Notes
    -----
    - Uses Delaunay triangulation for mesh generation
    - Includes points near the perimeter to ensure proper edge coverage
    - Handles Qhull errors by attempting alternative triangulation options
    - Validates surface area against the original polygon area
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
    # plot_polygon(poly, add_points=False)
    # plot_points_and_tris(all_points, in_tris + mid_tris, tcol='grey', plot_points=False, Show=True)

    # Step 7: Return the values
    return all_points, in_tris + mid_tris
