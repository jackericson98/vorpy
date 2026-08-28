from shapely import Polygon, Point, LineString
try:
    from shapely import contains_xy
except ImportError:
    contains_xy = None
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import calc_tri
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import project_to_plane
from scipy.spatial._qhull import QhullError
from time import perf_counter


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
    Generate concentric ring points used for surface triangulation.

    The original spiderweb geometry and point ordering are preserved while
    coordinate generation and bounds filtering are vectorized within each ring.
    """
    if ring_scaler is None:
        ring_scaler = 1

    min_x, max_x, min_y, max_y = (
        box[0][0] - 2 * res,
        box[1][0] + 2 * res,
        box[0][1] - 2 * res,
        box[1][1] + 2 * res
    )

    if center is None:
        center = [
            min_x + 0.5 * (max_x - min_x),
            min_y + 0.5 * (max_y - min_y)
        ]

    cx, cy = center
    corners = [box[0], [min_x, max_y], [max_x, min_y], box[1]]
    max_radius = max(calc_dist(center, corner) for corner in corners)
    num_rings = int(max_radius / res) + 1

    point_chunks = [np.asarray(center, dtype=float).reshape(1, 2)]

    for i in range(1, num_rings + 1):
        radius = max_radius * (i / num_rings)
        num_points_per_ring = int(2 * np.pi * radius / res) + 1

        j = np.arange(num_points_per_ring, dtype=np.float64)
        angles = 2 * np.pi * j / num_points_per_ring

        x = cx + radius * np.cos(angles)
        y = cy + radius * np.sin(angles)

        mask = (
            (x >= min_x) &
            (x <= max_x) &
            (y >= min_y) &
            (y <= max_y)
        )

        if np.any(mask):
            point_chunks.append(np.column_stack((x[mask], y[mask])))

    if len(point_chunks) == 1:
        return point_chunks[0]

    return np.concatenate(point_chunks, axis=0)

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


def sort_tris(perimeter, tris, polygon, numeric_points, timing=None):
    """
    Production triangle sorting.

    Keeps the validated optimizations:
      - indexed perimeter/interior classification
      - numeric-coordinate centroid calculation
      - vectorized all-perimeter containment via contains_xy
      - exact reconstruction in original Delaunay triangle order

    Detailed per-triangle profiling has been removed.
    """
    n_perimeter = len(perimeter)

    point_desigs = {
        i: ('e' if i < n_perimeter else 'i')
        for i in range(len(numeric_points))
    }

    in_, out, mid = [], [], []

    tri_modes = np.zeros(len(tris), dtype=np.uint8)
    perimeter_centroids = []

    all_interior = 0
    all_perimeter = 0
    mixed = 0

    for tri_idx, tri in enumerate(tris):
        i0, i1, i2 = tri

        e0 = i0 < n_perimeter
        e1 = i1 < n_perimeter
        e2 = i2 < n_perimeter

        if not e0 and not e1 and not e2:
            all_interior += 1

        elif e0 and e1 and e2:
            all_perimeter += 1
            tri_modes[tri_idx] = 1

            p0 = numeric_points[i0]
            p1 = numeric_points[i1]
            p2 = numeric_points[i2]

            cx = (p0[0] + p1[0] + p2[0]) / 3.0
            cy = (p0[1] + p1[1] + p2[1]) / 3.0

            perimeter_centroids.append((cx, cy))

        else:
            mixed += 1

    perimeter_inside = None

    if perimeter_centroids:
        if contains_xy is not None:
            centroid_arr = np.asarray(perimeter_centroids, dtype=float)
            perimeter_inside = np.asarray(
                contains_xy(
                    polygon,
                    centroid_arr[:, 0],
                    centroid_arr[:, 1]
                ),
                dtype=bool
            )
        else:
            perimeter_inside = np.asarray(
                [
                    polygon.contains(Point(com))
                    for com in perimeter_centroids
                ],
                dtype=bool
            )

    perimeter_result_pos = 0

    for tri_idx, tri in enumerate(tris):
        if tri_modes[tri_idx] == 0:
            in_.append(tri)
        else:
            if perimeter_inside[perimeter_result_pos]:
                in_.append(tri)
            else:
                out.append(tri)
            perimeter_result_pos += 1

    if timing is not None:
        timing['tri_all_interior'] = timing.get('tri_all_interior', 0) + all_interior
        timing['tri_all_perimeter'] = timing.get('tri_all_perimeter', 0) + all_perimeter
        timing['tri_mixed'] = timing.get('tri_mixed', 0) + mixed
        timing['tri_perimeter_containment_tests'] = (
            timing.get('tri_perimeter_containment_tests', 0) + all_perimeter
        )

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


def triangulate_2D_Surface(perimeter, res=0.2, center=None, timing=None):
    """
    Triangulates a 2D surface defined by a perimeter of points.

    This function creates a triangulated mesh of a 2D surface by:
    1. Determining the bounding box of the perimeter points
    2. Generating a uniform grid of points within the bounding box
    3. Creating Polygon and LineString objects for spatial analysis
    4. Filtering grid points using vectorized containment predicates
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

    tri_start = perf_counter()

    def _add(key, elapsed):
        if timing is not None:
            timing[key] = timing.get(key, 0.0) + elapsed

    # Step 1: Bounding box. Calculation is unchanged.
    t = perf_counter()
    px, py = [_[0] for _ in perimeter], [_[1] for _ in perimeter]
    box = [[min(px), min(py)], [max(px), max(py)]]
    _add('tri_bbox', perf_counter() - t)

    # Step 2: Historical spiderweb generation.
    t = perf_counter()
    grid_points = generate_spiderweb(box, res, center)
    _add('tri_spiderweb', perf_counter() - t)

    # Step 3: Historical Shapely setup.
    t = perf_counter()
    poly = Polygon(perimeter)
    linestring = LineString(perimeter)
    buffer = linestring.buffer(res / 2)
    all_points = perimeter.copy()
    _add('tri_shapely_setup', perf_counter() - t)

    # Step 4: Filter the exact same spiderweb points, preserving their order.
    #
    # Shapely 2.x exposes contains_xy(), which evaluates the same GEOS
    # containment predicate without constructing a Python Point object for
    # every rejected candidate. This is intentionally limited to the
    # filtering implementation; candidate coordinates, predicates,
    # accepted-point ordering, and downstream triangulation are unchanged.
    t = perf_counter()
    accepted = 0

    if contains_xy is not None and len(grid_points) > 0:
        grid_arr = np.asarray(grid_points)

        inside_mask = np.asarray(
            contains_xy(poly, grid_arr[:, 0], grid_arr[:, 1]),
            dtype=bool
        )
        inside_points = grid_arr[inside_mask]

        if len(inside_points) > 0:
            buffer_mask = np.asarray(
                contains_xy(buffer, inside_points[:, 0], inside_points[:, 1]),
                dtype=bool
            )
            accepted_points = inside_points[~buffer_mask]
        else:
            accepted_points = inside_points

        accepted = len(accepted_points)

        # Preserve historical accepted-point ordering exactly without
        # constructing Shapely Point objects. Delaunay and downstream
        # classification use numeric coordinates only.
        all_points.extend(accepted_points)

    else:
        # Compatibility fallback: original scalar implementation.
        for point in grid_points:
            test_point = Point(point)
            if poly.contains(test_point):
                if not buffer.contains(test_point):
                    all_points.append(point)
                    accepted += 1

    _add('tri_point_filter', perf_counter() - t)

    # Step 5: Historical Delaunay + QJ fallback.
    t = perf_counter()
    qhull_fallback = 0
    try:
        triangles = Delaunay(all_points).simplices
    except QhullError:
        qhull_fallback = 1
        try:
            triangles = Delaunay(all_points, qhull_options='QJ').simplices
        except QhullError:
            _add('tri_delaunay', perf_counter() - t)
            if timing is not None:
                timing['tri_calls'] = timing.get('tri_calls', 0) + 1
                timing['tri_qhull_fallbacks'] = timing.get('tri_qhull_fallbacks', 0) + qhull_fallback
                timing['tri_grid_points'] = timing.get('tri_grid_points', 0) + len(grid_points)
                timing['tri_accepted_points'] = timing.get('tri_accepted_points', 0) + accepted
                timing['tri_input_points'] = timing.get('tri_input_points', 0) + len(all_points)
                timing['tri_total'] = timing.get('tri_total', 0.0) + (perf_counter() - tri_start)
            return all_points, []
    _add('tri_delaunay', perf_counter() - t)

    # Step 6: Historical triangle classification.
    t = perf_counter()
    in_tris, out_tris, mid_tris, mid_tri_designations = sort_tris(
        perimeter, triangles, poly, all_points, timing=timing
    )
    _add('tri_sort', perf_counter() - t)

    # Step 7: Historical mid-triangle reassignment.
    t = perf_counter()
    mid_before = len(mid_tris)
    if len(mid_tris) > 0:
        mid_tris = reassign_tri_points(perimeter, mid_tris, poly, all_points)
    _add('tri_reassign', perf_counter() - t)

    if timing is not None:
        timing['tri_calls'] = timing.get('tri_calls', 0) + 1
        timing['tri_qhull_fallbacks'] = timing.get('tri_qhull_fallbacks', 0) + qhull_fallback
        timing['tri_grid_points'] = timing.get('tri_grid_points', 0) + len(grid_points)
        timing['tri_accepted_points'] = timing.get('tri_accepted_points', 0) + accepted
        timing['tri_input_points'] = timing.get('tri_input_points', 0) + len(all_points)
        timing['tri_raw_triangles'] = timing.get('tri_raw_triangles', 0) + len(triangles)
        timing['tri_inside_triangles'] = timing.get('tri_inside_triangles', 0) + len(in_tris)
        timing['tri_outside_triangles'] = timing.get('tri_outside_triangles', 0) + len(out_tris)
        timing['tri_mid_before'] = timing.get('tri_mid_before', 0) + mid_before
        timing['tri_mid_after'] = timing.get('tri_mid_after', 0) + len(mid_tris)
        timing['tri_total'] = timing.get('tri_total', 0.0) + (perf_counter() - tri_start)

    return all_points, in_tris + mid_tris