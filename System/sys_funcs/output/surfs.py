import os
import numpy as np
from System.sys_funcs.draw.draw import color_tris


def write_surfs(surfs, file_name, color=False, directory=None):
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
        # Get the maximum curvature of the
        # Count the number of triangles and vertices there are
        num_verts, num_tris = 0, 0
        for i, surf in enumerate(surfs):
            if surf.points is None:
                continue
            if surf.tri_colors is None:
                color_tris(surf=surf, color_map=surf.color_map, color_scheme=surf.scheme, max_val=surfs[0].net.max_curv)
            num_verts += len(surf.points)
            num_tris += len(surf.tris)
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i in range(len(surfs)):
            # Go through the points on the surface
            for point in surfs[i].points:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i in range(len(surfs)):
            surf = surfs[i]
            # Go through the triangles in the surface
            for j in range(len(surfs[i].tris)):
                # Get the triangle and colors
                tri = surf.tris[j]
                colors = surf.tri_colors
                if colors is not None:
                    # If the surface is flat, average out the colors
                    color = colors[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(surfs[i].points)

