import numpy as np


def draw_line(points, radius=0.05, color="#000000", edge_org=None):
    """
    Creates a 3D cylindrical line segment between points using triangular mesh representation.

    Args:
        points (list of array-like): List of 3D points.
        radius (float): Radius of the tube.
        color (tuple, optional): Unused here, kept for compatibility.
        edge_org (array-like, optional): Optional reference vector for initializing the
            first frame only.

    Returns:
        tuple:
            - draw_points: list of 3D vertex positions
            - draw_tris: list of triangle indices
    """
    if len(points) < 2:
        return [], []

    pts = [np.asarray(p, dtype=float) for p in points]
    draw_points = []
    draw_tris = []
    eps = 1e-12

    tangents = []
    for i in range(len(pts)):
        if i == 0:
            t = pts[1] - pts[0]
        elif i == len(pts) - 1:
            t = pts[-1] - pts[-2]
        else:
            t = pts[i + 1] - pts[i - 1]

        norm_t = np.linalg.norm(t)
        if norm_t < eps:
            raise ValueError(f"Degenerate segment near point index {i}.")

        tangents.append(t / norm_t)

    t0 = tangents[0]

    if edge_org is not None:
        ref = np.asarray(edge_org, dtype=float)
    else:
        ref = np.array([0.0, 0.0, 1.0])

    if abs(np.dot(ref, t0)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])

    n0 = np.cross(t0, ref)
    norm_n0 = np.linalg.norm(n0)
    if norm_n0 < eps:
        ref = np.array([0.0, 1.0, 0.0])
        n0 = np.cross(t0, ref)
        norm_n0 = np.linalg.norm(n0)
        if norm_n0 < eps:
            raise ValueError("Could not construct an initial normal vector.")

    n0 = n0 / norm_n0
    b0 = np.cross(t0, n0)
    b0 = b0 / np.linalg.norm(b0)

    normals = [n0]
    binormals = [b0]

    for i in range(1, len(pts)):
        t_prev = tangents[i - 1]
        t_curr = tangents[i]
        n_prev = normals[-1]

        axis = np.cross(t_prev, t_curr)
        axis_norm = np.linalg.norm(axis)

        if axis_norm < eps:
            n_curr = n_prev.copy()
        else:
            axis = axis / axis_norm
            angle = np.arccos(np.clip(np.dot(t_prev, t_curr), -1.0, 1.0))
            n_curr = (
                n_prev * np.cos(angle)
                + np.cross(axis, n_prev) * np.sin(angle)
                + axis * np.dot(axis, n_prev) * (1.0 - np.cos(angle))
            )

        n_curr = n_curr - np.dot(n_curr, t_curr) * t_curr
        norm_n_curr = np.linalg.norm(n_curr)
        if norm_n_curr < eps:
            n_curr = normals[-1].copy()
            n_curr = n_curr - np.dot(n_curr, t_curr) * t_curr
            norm_n_curr = np.linalg.norm(n_curr)
            if norm_n_curr < eps:
                raise ValueError(f"Could not stabilize frame at point index {i}.")

        n_curr = n_curr / norm_n_curr
        b_curr = np.cross(t_curr, n_curr)
        b_curr = b_curr / np.linalg.norm(b_curr)

        normals.append(n_curr)
        binormals.append(b_curr)

    for i, p in enumerate(pts):
        n = normals[i]
        b = binormals[i]

        p0 = p + radius * n
        p1 = p + radius * (-0.5 * n + 0.5 * np.sqrt(3.0) * b)
        p2 = p + radius * (-0.5 * n - 0.5 * np.sqrt(3.0) * b)

        draw_points.extend([p0, p1, p2])

    for i in range(len(pts) - 1):
        p0_0, p0_1, p0_2 = 3 * i, 3 * i + 1, 3 * i + 2
        p1_0, p1_1, p1_2 = 3 * (i + 1), 3 * (i + 1) + 1, 3 * (i + 1) + 2

        draw_tris.extend([
            [p0_0, p0_1, p1_0], [p1_0, p1_1, p0_1],
            [p0_1, p0_2, p1_1], [p1_1, p1_2, p0_2],
            [p0_2, p0_0, p1_2], [p1_2, p1_0, p0_0],
        ])

    return draw_points, draw_tris


# Draw Edge Function. Takes in an edge and updates its attributes draw_points, draw_tris
def draw_edge(edge, radius=0.05, color=None):
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
