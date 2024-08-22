import numpy as np


def build_perimeter(locs, rads, epnts, net_type='aw'):
    """
    Sorts the edges of the surface to create a list of points in order around the perimeter
    :param net_type:
    :param surf: Surface object for perimeter building
    :return: None
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
        if not reverse:  # In order
            perimeter += my_edge_points
        else:  # Reverse order
            perimeter += my_edge_points[::-1]
    d = np.sqrt(sum(np.square(np.array(locs[1]) - np.array(locs[0]))))
    # Get the center of the surface
    r = np.array(locs[1]) - np.array(locs[0])
    r = [_ if _ != 0 else 0.0001 for _ in r]
    surf_norm = r / np.linalg.norm(r)
    surf_loc = None
    if net_type == 'aw':
        surf_loc = np.array(locs[0]) + (rads[0] + 0.5 * (d - (rads[0] + rads[1]))) * surf_norm
    elif net_type == 'del':
        surf_loc = np.array(locs[0]) + 0.5 * d * surf_norm
    elif net_type == 'pow':
        surf_loc = np.array(locs[0]) + 0.5 * (surf_norm ** 2 + rads[0] ** 2 - rads[1] ** 2) / surf_norm
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
