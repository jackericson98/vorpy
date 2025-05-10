import numpy as np


def build_perimeter(locs, rads, epnts, net_type='aw'):
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
