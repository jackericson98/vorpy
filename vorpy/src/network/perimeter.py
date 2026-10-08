import numpy as np


def _sort_key(value):
    return repr(value)


def stable_vertex_identity(row):
    """Return a table-index-independent identity for one AW vertex."""
    memberships = row.get('balls')
    balls = tuple(sorted(int(value) for value in memberships)) if memberships is not None else ()
    location = tuple(float(value) for value in row.get('loc'))
    return balls, location


def _reverse_metadata(metadata):
    if metadata is None or not hasattr(metadata, 'items'):
        return metadata
    result = dict(metadata)
    for key in ('vnorm', 'dnorm', 'vnorm0', 'dnorm0', 'vnorm1', 'dnorm1'):
        if key in result and result[key] is not None:
            result[key] = -np.asarray(result[key], dtype=float)
    for left, right in (('loc', 'loc2'), ('rad', 'rad2')):
        if left in metadata and right in metadata:
            result[left], result[right] = metadata[right], metadata[left]
    for suffix, other in (('0', '1'),):
        for prefix in ('vmid', 'vdist', 'pa'):
            left, right = f'{prefix}{suffix}', f'{prefix}{other}'
            if left in metadata and right in metadata:
                result[left], result[right] = metadata[right], metadata[left]
        for prefix in ('vnorm', 'dnorm'):
            left, right = f'{prefix}{suffix}', f'{prefix}{other}'
            if left in metadata and right in metadata:
                result[left] = -np.asarray(metadata[right], dtype=float)
                result[right] = -np.asarray(metadata[left], dtype=float)
    return result


def canonicalize_edge_orientation(points, endpoint_keys, metadata=None):
    """Orient an edge by stable endpoint identities and reverse its metadata."""
    if len(endpoint_keys) != 2:
        raise ValueError("An edge must have exactly two endpoint identities.")
    if _sort_key(endpoint_keys[0]) <= _sort_key(endpoint_keys[1]):
        return list(points), metadata, False
    return list(reversed(points)), _reverse_metadata(metadata), True


def _surface_normal(locs):
    direction = np.asarray(locs[1], dtype=float) - np.asarray(locs[0], dtype=float)
    direction = np.asarray([value if value != 0 else 0.0001 for value in direction])
    return direction / np.linalg.norm(direction)


def _surface_location(locs, rads, net_type):
    distance_vector = np.asarray(locs[1], dtype=float) - np.asarray(locs[0], dtype=float)
    distance = np.sqrt(sum(np.square(distance_vector)))
    normal = _surface_normal(locs)
    if net_type == 'aw':
        return np.asarray(locs[0]) + (
            rads[0] + 0.5 * (distance - (rads[0] + rads[1]))
        ) * normal
    if net_type == 'prm':
        return np.asarray(locs[0]) + 0.5 * distance * normal
    if net_type == 'pow':
        offset = 0.5 * (distance ** 2 + rads[0] ** 2 - rads[1] ** 2) / distance
        return np.asarray(locs[0]) + offset * normal
    return None


def _polygon_normal(points):
    values = np.asarray(points, dtype=float)
    if len(values) < 3:
        return np.zeros(3, dtype=float)
    return np.sum(np.cross(values, np.roll(values, -1, axis=0)), axis=0)


def _canonical_edge_record(points, endpoints):
    if len(endpoints) != 2:
        raise ValueError("Each boundary edge must have exactly two endpoints.")
    start, end = tuple(endpoints)
    values = [np.asarray(point, dtype=float) for point in points]
    if len(values) < 2 or any(point.shape != (3,) for point in values):
        raise ValueError("Boundary edge samples must contain 3D points.")
    values, _, _ = canonicalize_edge_orientation(values, (start, end))
    if _sort_key(start) <= _sort_key(end):
        return start, end, values
    return end, start, values


def _canonical_boundary_perimeter(locs, rads, epnts, edge_endpoints, net_type):
    if len(epnts) != len(edge_endpoints):
        raise ValueError("Boundary edge endpoints must align with edge samples.")

    records = [
        _canonical_edge_record(points, endpoints)
        for points, endpoints in zip(epnts, edge_endpoints)
    ]
    incident = {}
    for edge_index, (start, end, _) in enumerate(records):
        incident.setdefault(start, []).append(edge_index)
        incident.setdefault(end, []).append(edge_index)
    if any(len(edges) < 2 for edges in incident.values()):
        raise ValueError("Boundary connectivity is not a collection of closed cycles.")
    if any(len(edges) > 2 for edges in incident.values()):
        raise ValueError("Boundary connectivity has a higher-order junction.")

    def edge_key(edge_index):
        start, end, points = records[edge_index]
        return (
            _sort_key(start),
            _sort_key(end),
            tuple(tuple(float(value) for value in point) for point in points),
        )

    def trace(start_vertex, first_edge):
        current = start_vertex
        edge_index = first_edge
        used = set()
        perimeter = []
        while True:
            if edge_index in used:
                raise ValueError("Boundary cycle repeats an edge before closing.")
            start, end, points = records[edge_index]
            if current == start:
                next_vertex, ordered_points = end, points
            elif current == end:
                next_vertex, ordered_points = start, list(reversed(points))
            else:
                raise ValueError("Boundary edge is not incident to the current vertex.")
            perimeter.extend(ordered_points)
            used.add(edge_index)
            current = next_vertex
            if current == start_vertex:
                return used, perimeter
            candidates = [
                candidate
                for candidate in incident[current]
                if candidate not in used
            ]
            if len(candidates) != 1:
                raise ValueError("Boundary traversal is ambiguous or does not close.")
            edge_index = candidates[0]

    remaining = set(range(len(records)))
    cycles = []
    while remaining:
        first_edge = min(remaining, key=edge_key)
        start_vertex = min(records[first_edge][:2], key=_sort_key)
        first_edges = sorted(incident[start_vertex], key=edge_key)
        forward_used, forward = trace(start_vertex, first_edges[0])
        reverse_used, reverse = trace(start_vertex, first_edges[-1])
        if forward_used != reverse_used:
            raise ValueError("Boundary traversal produced inconsistent cycles.")
        remaining.difference_update(forward_used)
        forward_score = float(np.dot(_polygon_normal(forward), _surface_normal(locs)))
        reverse_score = float(np.dot(_polygon_normal(reverse), _surface_normal(locs)))
        if reverse_score > forward_score:
            perimeter = reverse
        elif forward_score > reverse_score:
            perimeter = forward
        else:
            perimeter = min(forward, reverse, key=lambda values: repr(values))
        cycles.append(perimeter)

    if len(cycles) != 1:
        raise ValueError(
            f"Multiple boundary cycles ({len(cycles)}) are not supported by "
            "the single-polygon surface triangulator."
        )
    return cycles[0], _surface_location(locs, rads, net_type), _surface_normal(locs)


def build_perimeter(
    locs, rads, epnts, net_type='aw', edge_endpoints=None,
    allow_legacy_fallback=False,
):
    """
    Builds a perimeter by sorting and connecting edge points in a continuous loop around a surface.

    This function takes a set of edge points and organizes them into a continuous perimeter by:
    1. Starting with the first edge's points
    2. Finding the closest remaining edge points to continue the perimeter
    3. Connecting edges in either forward or reverse order to maintain continuity

    Parameters
    ----------
    locs : list of numpy.ndarray
        List of ball locations in 3D space
    rads : list of float
        List of ball radii
    epnts : list of list of numpy.ndarray
        List of edge points, where each edge is a list of points
    net_type : str, optional
        Type of network ('aw', 'pow', or 'prm'), defaults to 'aw'

    Returns
    -------
    tuple
        A tuple containing:
        - perimeter : list of numpy.ndarray
            Ordered list of points forming the continuous perimeter
        - surf_loc : numpy.ndarray
            Location of the surface center
        - surf_norm : numpy.ndarray
            Normal vector of the surface

    Notes
    -----
    - The function ensures the perimeter forms a continuous loop by connecting edges at their closest points
    - Edge points can be added in reverse order if it provides a better connection
    - The surface location and normal are calculated based on the network type
    - For 'aw' networks, the surface is positioned at the average of the ball radii
    - For 'prm' networks, the surface is positioned at the midpoint between balls
    - For 'pow' networks, the surface position is calculated using the power diagram formula
    """
    if not epnts:
        raise ValueError("Cannot build a perimeter without edge points.")

    # A failed edge projection can produce an array of NaNs.  Let the caller
    # discard that surface rather than allowing the traversal below to keep a
    # ``None`` index and fail with an opaque list.pop error.
    for edge_index, edge_points in enumerate(epnts):
        points = np.asarray(edge_points, dtype=float)
        if points.ndim != 2 or points.shape[0] < 2 or points.shape[1] != 3:
            raise ValueError(
                f"Invalid edge {edge_index}: expected at least two 3D points."
            )
        if not np.all(np.isfinite(points)):
            raise ValueError(f"Invalid edge {edge_index}: non-finite point detected.")

    if net_type == 'aw' and edge_endpoints is not None:
        try:
            return _canonical_boundary_perimeter(
                locs, rads, epnts, edge_endpoints, net_type
            )
        except ValueError as error:
            message = str(error)
            unsupported_boundary = (
                message == "Boundary connectivity has a higher-order junction."
                or message.startswith("Multiple boundary cycles (")
            )
            if not allow_legacy_fallback or not unsupported_boundary:
                raise
            # shortcut: unsupported AW boundary topology retains the historical
            # walk until the single-polygon triangulator supports such cycles.
            records = [
                _canonical_edge_record(points, endpoints)
                for points, endpoints in zip(epnts, edge_endpoints)
            ]
            records.sort(
                key=lambda record: (
                    _sort_key(record[0]),
                    _sort_key(record[1]),
                    tuple(tuple(float(value) for value in point) for point in record[2]),
                )
            )
            epnts = [record[2] for record in records]

    # Add the first edge's vertex location and set of points to the perimeter points list
    perimeter = epnts[0][:]
    # Make a copy of the edges to organize excluding the first edge
    edges_points = epnts[1:]

    # Keep looping while we haven't gone through the edges
    while edges_points:

        # Set the max distance to infinity, the index for the intended edge to None and the reverse bool to False
        d, ndx, reverse = np.inf, None, False

        # Go through each of the remaining edges in the list
        for i in range(len(edges_points)):
            # Calculate the distance between the most recently recorded point and the first/last points in the edge
            d0 = np.sqrt(sum(np.square(np.array(perimeter[-1]) - np.array(edges_points[i][0]))))
            d1 = np.sqrt(sum(np.square(np.array(perimeter[-1]) - np.array(edges_points[i][-1]))))
            # If the first edge point is closer to the last perimeter point and the last isn't closer add that edge
            if d0 < d and d0 < d1:
                d, ndx, reverse = d0, i, False
            # Otherwise, if the last edge point is the closest add the edge in reverse
            elif d1 < d:
                d, ndx, reverse = d1, i, True
        if ndx is None:
            raise ValueError(
                "Unable to connect edge points into a perimeter; "
                "the surface geometry is disconnected or non-finite."
            )

        # Pull the edge from the list of edges
        my_edge_points = edges_points.pop(ndx)
        # Add the edge's point in the right order and then add the 181L vertex
        if reverse:  # In reverse order
            my_edge_points = my_edge_points[::-1]
        perimeter += my_edge_points
    d = np.sqrt(sum(np.square(np.array(locs[1]) - np.array(locs[0]))))
    # Get the center of the surface
    r = np.array(locs[1]) - np.array(locs[0])
    r = np.array([_ if _ != 0 else 0.0001 for _ in r])
    surf_norm = r / np.linalg.norm(r)
    surf_loc = None
    if net_type == 'aw':
        surf_loc = np.array(locs[0]) + (rads[0] + 0.5 * (d - (rads[0] + rads[1]))) * surf_norm
    elif net_type == 'prm':
        surf_loc = np.array(locs[0]) + 0.5 * d * surf_norm
    elif net_type == 'pow':
        d0 = 0.5 * (d ** 2 + rads[0] ** 2 - rads[1] ** 2) / d
        surf_loc = np.array(locs[0]) + d0 * surf_norm
    return perimeter, surf_loc, surf_norm


def build_perimeter1(epnts):
    # Add the first edge's vertex location and set of points to the perimeter points list
    perimeter = epnts[0][:]
    # Make a copy of the edges to organize excluding the first edge
    edges_points = epnts[1:]

    # Keep looping while we haven't gone through the edges
    while edges_points:

        # Set the max distance to infinity, the index for the intended edge to None and the reverse bool to False
        d, ndx, reverse = np.inf, None, False

        # Go through each of the remaining edges in the list
        for i in range(len(edges_points)):
            # Calculate the distance between the most recently recorded point and the first/last points in the edge
            d0 = np.sqrt(sum(np.square(np.array(perimeter[-1]) - np.array(edges_points[i][0]))))
            d1 = np.sqrt(sum(np.square(np.array(perimeter[-1]) - np.array(edges_points[i][-1]))))
            # If the first edge point is closer to the last perimeter point and the last isn't closer add that edge
            if d0 < d and d0 < d1:
                d, ndx, reverse = d0, i, False
            # Otherwise, if the last edge point is the closest add the edge in reverse
            elif d1 < d:
                d, ndx, reverse = d1, i, True
        # Pull the edge from the list of edges
        my_edge_points = edges_points.pop(ndx)
        # Add the edge's point in the right order and then add the 181L vertex
        if not reverse:  # In order
            perimeter += my_edge_points
        else:  # Reverse order
            perimeter += my_edge_points[::-1]
    # Return the perimeter
    return perimeter
