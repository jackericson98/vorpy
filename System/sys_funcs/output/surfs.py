import os
import numpy as np
from System.sys_funcs.output.draw import color_tris


def write_surfs(net, surfs, file_name, color=False, directory=None):
    """
    Writes files given a list of surfaces into the current directory or the given one
    :param surfs: Surface object
    :param file_name: Name of the output file for the surfaces
    :param color: Color of the output surface
    :param directory:
    :return:
    """
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no surfaces are provided return
    if surfs is None or len(surfs) == 0:
        return
    # If no color is given, make the color random
    if color is False:
        color = (1, 0, 0)
    if 'tri_colors' not in net.surfs:
        net.surfs['tri_colors'] = [[] for _ in range(len(net.surfs))]
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        num_verts, num_tris = 0, 0
        for ndx in surfs:
            surf = net.surfs.iloc[ndx]
            num_verts += len(surf['points'])
            num_tris += len(surf['tris'])
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for ndx in surfs:
            surf = net.surfs.iloc[ndx]
            # Go through the points on the surface
            for point in surf['points']:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for ndx in surfs:
            surf = net.surfs.iloc[ndx]
            tri_colors = [color for _ in range(len(surf['tris']))]
            if net.settings['net_type'] == 'aw':
                if net.settings['surf_scheme'] == 'gauss':
                    max_val = max(net.surfs['gauss_curv'])
                else:
                    max_val = max(net.surfs['mean_curv'])
                tri_colors = color_tris(surf=surf, color_map=net.settings['surf_col'],
                                        color_scheme=net.settings['surf_scheme'],
                                        color_factor=net.settings['scheme_factor'], max_val=max_val)

            # Go through the triangles in the surface
            for j, tri in enumerate(surf['tris']):
                # Get the triangle and colors
                color = tri_colors[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(surf['points'])


def write_surfs1(surfs, file_name, settings, color=False, directory=None):
    """
    Writes files given a list of surfaces into the current directory or the given one
    :param surfs: Surface object
    :param file_name: Name of the output file for the surfaces
    :param color: Color of the output surface
    :param directory:
    :return:
    """
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no surfaces are provided return
    if surfs is None or len(surfs) == 0:
        return
    # If no color is given, make the color random
    if color is False:
        color = np.random.rand(3)
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        num_verts, num_tris = 0, 0
        for i, surf in surfs.iterrows():
            num_verts += len(surf['points'])
            num_tris += len(surf['tris'])
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i, surf in surfs.iterrows():
            # Go through the points on the surface
            for point in surf['points']:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i, surf in surfs.iterrows():
            tri_colors = [color for _ in range(len(surf['tris']))]
            if settings['net_type'] == 'aw':
                max_val = max(surfs['curv'])
                tri_colors = color_tris(surf=surf, color_map=settings['surf_col'],
                                        color_scheme=settings['surf_scheme'],
                                        color_factor=settings['scheme_factor'], max_val=max_val)

            # Go through the triangles in the surface
            for j, tri in enumerate(surf['tris']):
                # Get the triangle and colors
                color = tri_colors[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(surf['points'])
