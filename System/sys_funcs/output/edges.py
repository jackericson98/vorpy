import os
from System.sys_funcs.draw.draw import draw_edge


def write_edges(edges, file_name, color=None, directory=None):
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no surfaces are provided return
    if edges is None or len(edges) == 0:
        return
    # If no color is given, make the color random
    if color is None:
        color = [0.5, 0.5, 0.5]
    # Check that the edge has been drawn
    for edge in edges:
        if edge.draw_points is None or edge.draw_tris is None:
            draw_edge(edge)
    num_verts, num_tris = 0, 0
    # Go through and create each edge
    for i in range(len(edges)):
        num_verts += len(edges[i].points) * 3
        num_tris += (len(edges[i].points) - 1) * 6
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i in range(len(edges)):
            # Go through the points on the surface
            for point in edges[i].draw_points:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i in range(len(edges)):
            edge = edges[i]
            # Go through the triangles in the surface
            for j in range(len(edge.draw_tris)):
                # Get the triangle and colors
                tri = edge.draw_tris[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(edge.draw_points)
