import os
import time
from vorpy.src.output.draw import draw_edge
from vorpy.src.output.colors import color_dict


def write_edges(net, edges, file_name, color=None, directory=None, profile=True):
    """Write selected network edges to an OFF file."""

    if directory is not None:
        os.chdir(directory)

    if edges is None or len(edges) == 0:
        return

    edges = list(edges)

    if color is None:
        color = 'gray'

    if color in color_dict:
        color = color_dict[color]

    # ------------------------------------------------------------------
    # Ensure cache columns exist
    # ------------------------------------------------------------------

    if 'draw_tris' not in net.edges:
        net.edges['draw_tris'] = [[] for _ in range(len(net.edges))]

    if 'draw_points' not in net.edges:
        net.edges['draw_points'] = [[] for _ in range(len(net.edges))]

    # ------------------------------------------------------------------
    # Gather rows
    # ------------------------------------------------------------------

    edge_rows = [net.edges.iloc[ndx] for ndx in edges]

    # ------------------------------------------------------------------
    # Generate/retrieve drawing geometry
    # ------------------------------------------------------------------

    edges_draw_points = []
    edges_draw_tris = []
    newly_drawn = 0
    cached = 0

    for ndx, edge in zip(edges, edge_rows):
        draw_points = edge['draw_points']
        draw_tris = edge['draw_tris']

        if draw_points is None or draw_tris is None or len(draw_points) == 0 or len(draw_tris) == 0:
            draw_points, draw_tris = draw_edge(edge)

            # Store directly in the corresponding DataFrame cell.
            net.edges.at[net.edges.index[ndx], 'draw_points'] = draw_points
            net.edges.at[net.edges.index[ndx], 'draw_tris'] = draw_tris

            newly_drawn += 1
        else:
            cached += 1

        edges_draw_points.append(draw_points)
        edges_draw_tris.append(draw_tris)

    # ------------------------------------------------------------------
    # Count output geometry
    # ------------------------------------------------------------------

    count_start = time.perf_counter()

    num_points = sum(len(points) for points in edges_draw_points)
    num_tris = sum(len(tris) for tris in edges_draw_tris)

    # ------------------------------------------------------------------
    # Write OFF
    # ------------------------------------------------------------------

    with open(file_name + ".off", 'w', buffering=1024 * 1024) as file:
        file.write(f"OFF\n{num_points} {num_tris} 0\n\n\n")

        # --------------------------------------------------------------
        # Points
        # --------------------------------------------------------------

        buffer = []
        chunk_size = 10000

        for draw_points in edges_draw_points:
            for point in draw_points:
                buffer.append(
                    f"{round(float(point[0]), 4)} "
                    f"{round(float(point[1]), 4)} "
                    f"{round(float(point[2]), 4)}\n"
                )

                if len(buffer) >= chunk_size:
                    file.write(''.join(buffer))
                    buffer.clear()

        if buffer:
            file.write(''.join(buffer))
            buffer.clear()

        # --------------------------------------------------------------
        # Faces
        # --------------------------------------------------------------

        vertex_offset = 0
        buffer = []

        for draw_points, draw_tris in zip(edges_draw_points, edges_draw_tris):
            for tri in draw_tris:
                buffer.append(
                    f"3 {tri[0] + vertex_offset} "
                    f"{tri[1] + vertex_offset} "
                    f"{tri[2] + vertex_offset} "
                    f"{color[0]} {color[1]} {color[2]}\n"
                )

                if len(buffer) >= chunk_size:
                    file.write(''.join(buffer))
                    buffer.clear()

            vertex_offset += len(draw_points)

        if buffer:
            file.write(''.join(buffer))



def write_edges1(edges, file_name, color=None, directory=None):
    """
    Writes an off file for the edges specified
    :param edges: Edges to be output
    :param file_name: Name for the output file
    :param color: Color for the edges
    :param directory: Output directory
    :return: None
    """
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no surfaces are provided return
    if edges is None or len(edges) == 0:
        return
    # If no color is given, make the color random
    if color is None:
        color = [0.5, 0.5, 0.5]

    edges_draw_points, edges_draw_tris = [], []
    for i, edge in edges.iterrows():
        draw_points, draw_tris = draw_edge(edge)
        edges_draw_points.append(draw_points)
        edges_draw_tris.append(draw_tris)

    num_verts, num_tris = 0, 0
    # Go through and create each edge
    for i, edge in edges.iterrows():
        num_verts += len(edge['points']) * 3
        num_tris += (len(edge['points']) - 1) * 6
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for edge_draw_points in edges_draw_points:
            # Go through the points on the surface
            for point in edge_draw_points:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for edge_draw_tris in edges_draw_tris:
            # Go through the triangles in the surface
            for tri in edge_draw_tris:
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(edge_draw_tris)
