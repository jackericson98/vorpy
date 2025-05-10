import numpy as np


def draw_line(points, radius=0.02, color=None, edge_org=None):
    """
    Creates a 3D cylindrical line segment between points using triangular mesh representation.

    This function generates a triangular mesh representation of a cylindrical line segment
    connecting a series of points in 3D space. The cylinder is constructed using triangular
    faces and can be customized with different radii and colors.

    Args:
        points (list of numpy.ndarray): List of 3D points to connect with the line segment
        radius (float, optional): Radius of the cylindrical line segment. Defaults to 0.02
        color (tuple, optional): RGB color tuple for the line segment. Defaults to None
        edge_org (numpy.ndarray, optional): Vector defining the orientation of the edge.
            Defaults to [0, 0, 1] if not specified

    Returns:
        tuple: A tuple containing:
            - draw_points (list): List of 3D points defining the vertices of the triangular mesh
            - draw_tris (list): List of triangular faces connecting the vertices
    """
    if edge_org is None:
        edge_org = [0, 0, 1]
    # Initiate the draw attributes
    draw_points, draw_tris = [], []
    r = None
    # Go through the points
    for i in range(len(points)):
        # If we are at the end of the points list, use the previous point for calibration
        p0 = np.array(points[i])
        if i < len(points) - 1:
            p1 = np.array(points[i + 1])
            r = p1 - p0
        # Find the vector and its normal between the two points
        rn = r / np.linalg.norm(r)
        # In the case that the vector between the points is in the z direction only, move it
        if rn[0] == 0 and rn[1] == 0:
            r = r + np.array([0.001, 0.001, 0])
            rn = r / np.linalg.norm(r)
        # Take the cross product with the +z direction and normalize it
        v0_0x = np.cross(rn, np.array(p0 - edge_org))
        v0_0n = v0_0x / np.linalg.norm(v0_0x)
        # Calculate the location of the first point
        p0_0 = v0_0n * radius + p0
        # Take the cross product of the edge vector and the vector to the first point and normalize it
        v0_1x = np.cross(rn, v0_0n)
        v0_1nx = v0_1x / np.linalg.norm(v0_1x)
        # Find the vectors for the other two points (30/60/90 triangle)
        v0_1 = - 0.5 * radius * v0_0n + 0.5 * np.sqrt(3) * radius * v0_1nx
        v0_2 = - 0.5 * radius * v0_0n - 0.5 * np.sqrt(3) * radius * v0_1nx
        # Get the points and add them to the list of draw points
        p0_1, p0_2 = v0_1 + p0, v0_2 + p0
        draw_points += [p0_0, p0_1, p0_2]
    # Go through the points
    for i in range(len(points) - 1):
        # List the points
        p0_0, p0_1, p0_2, p1_0, p1_1, p1_2 = range(3 * i, 3 * (i + 2))
        # Create the triangles
        draw_tris += [[p0_0, p0_1, p1_0], [p1_0, p1_1, p0_1],
                           [p0_1, p0_2, p1_1], [p1_1, p1_2, p0_2],
                           [p0_2, p0_0, p1_2], [p1_2, p1_0, p0_0]]
    # Return the points and triangles
    return draw_points, draw_tris


# Draw Edge Function. Takes in an edge and updates its attributes draw_points, draw_tris
def draw_edge(edge, radius=0.02, color=None):
    """
    Draws an edge in triangles and points
    :param edge: Edge object for exporting
    :param radius: Radius of the edge to be drawn
    :param color: Color for the edge drawing
    :return: None
    """
    # # Get the edge direction to point away from
    # rads = [_.rad for _ in edge.balls]
    # min_ball = edge['balls'][rads.index(min(rads))]
    # if edge.points is None or len(edge.points) <= 1:
    #     edge.points, edge.vals = build_edge(edge['balls'], edge['verts'], edge.net.surf_res)
    # Calculate the lines
    return draw_line(edge.points, radius, color=color)
