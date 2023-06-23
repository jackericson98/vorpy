import os
from System.sys_funcs.draw.draw import draw_edge


def write_edges(net, edges, file_name, color=None, directory=None):
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
    # Check that the edge has been drawn
    if 'draw_tris' not in net.edges:
        net.edges['draw_tris'] = [[] for _ in range(len(net.edges))]
    if 'draw_points' not in net.edges:
        net.edges['draw_points'] = [[] for _ in range(len(net.edges))]
    edges_draw_points, edges_draw_tris = [], []
    for ndx in edges:
        edge = net.edges.iloc[ndx]
        if 'draw_points' not in edge or edge['draw_points'] is None or edge['draw_tris'] is None or edge['draw_tris'] == []:
            draw_points, draw_tris = draw_edge(edge)
            edges_draw_points.append(draw_points)
            edges_draw_tris.append(draw_tris)
        else:
            edges_draw_points.append(edge['draw_points'])
            edges_draw_tris.append(edge['draw_tris'])
    j = 0
    net_edges_draw_points, net_edges_draw_tris = [], []
    for i in range(len(net.edges)):
        if i in edges:
            net_edges_draw_tris.append(edges_draw_tris[j])
            net_edges_draw_points.append(edges_draw_points[j])
            j += 1
        else:
            net_edges_draw_tris.append(net.edges['draw_tris'][i])
            net_edges_draw_points.append(net.edges['draw_points'][i])
    net.edges['draw_tris'], net.edges['draw_points'] = net_edges_draw_tris, net_edges_draw_points
    my_edges = [net.edges.iloc[_] for _ in edges]
    num_verts, num_tris = 0, 0
    # Go through and create each edge
    for edge in my_edges:
        num_verts += len(edge['points']) * 3
        num_tris += (len(edge['points']) - 1) * 6
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for edge in my_edges:
            # Go through the points on the surface
            for point in edge['draw_points']:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for edge in my_edges:
            # Go through the triangles in the surface
            for tri in edge['draw_tris']:
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(edge['draw_points'])
