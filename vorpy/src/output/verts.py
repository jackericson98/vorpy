import os
import numpy as np
from vorpy.src.output.colors import color_dict


def write_off_verts(net, verts, file_name, atom_type=None, directory=None, color=None, vert_rad=0.05):
    """
    Exports vertex data to an OFF (Object File Format) file for visualization.

    This function creates an OFF file containing geometric data for specified vertices from a network.
    Each vertex is represented as a small octahedron centered at its location, with the size controlled
    by the vertex radius parameter. The vertices can be colored uniformly using the provided color values.

    Args:
        net: Network object containing vertex data
        verts: List of vertex indices to export
        file_name: Base name for the output OFF file
        atom_type: Optional atom type for vertex representation (default: 'He')
        directory: Optional output directory path. If None, uses current directory
        color: Optional RGB color tuple for vertex coloring. If None, uses default red (1,0,0)
        vert_rad: Radius of the octahedron representing each vertex (default: 0.05)

    Returns:
        None: Creates an OFF file with the specified vertex data
    """
    # If no color is given, make the color random
    if color is None:
        color = [1, 0, 0]
    if color in color_dict:
        color = color_dict[color]
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    if atom_type is None:
        atom_type = 'He'
    # If no surfaces are provided return
    if verts is None or len(verts) == 0:
        return
    loc_points, loc_tris = [], []
    for vert in verts:
        loc = net.verts['loc'][vert]
        # Draw the point
        xp, xn = loc + np.array([vert_rad, 0, 0]), loc - np.array([vert_rad, 0, 0])
        yp, yn = loc + np.array([0, vert_rad, 0]), loc - np.array([0, vert_rad, 0])
        zp, zn = loc + np.array([0, 0, vert_rad]), loc - np.array([0, 0, vert_rad])
        # Connect the points
        loc_points.append([xp, xn, yp, yn, zp, zn])
        loc_tris.append([[0, 2, 4], [0, 2, 5], [0, 3, 4], [0, 3, 5], [1, 2, 4], [1, 2, 5], [1, 3, 4], [1, 3, 5]])
    num_verts, num_tris = 8 * len(loc_points), 6 * len(loc_tris)
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i in range(len(verts)):
            # Go through the points on the surface
            for point in loc_points[i]:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i in range(len(verts)):
            tri_list = loc_tris[i]
            # Go through the triangles in the surface
            for j in range(len(tri_list)):
                # Get the triangle and colors
                tri = tri_list[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(loc_points[i])


def write_off_verts1(verts, file_name, atom_type=None, directory=None, color=None, vert_rad=0.05):
    """
    Creates a pdb file for vertex representation
    :param vert_rad:
    :param color:
    :param pdb:
    :param verts:
    :param file_name:
    :param atom_type:
    :param directory:
    :return:
    """
    # If no color is given, make the color random
    if color is None:
        color = [1, 0, 0]
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    if atom_type is None:
        atom_type = 'He'
    # If no surfaces are provided return
    if verts is None or len(verts) == 0:
        return
    loc_points, loc_tris = [], []
    for i, vert in verts.iterrows():
        # Get the location of the vert
        loc = vert['loc']
        # Draw the point
        xp, xn = loc + np.array([vert_rad, 0, 0]), loc - np.array([vert_rad, 0, 0])
        yp, yn = loc + np.array([0, vert_rad, 0]), loc - np.array([0, vert_rad, 0])
        zp, zn = loc + np.array([0, 0, vert_rad]), loc - np.array([0, 0, vert_rad])
        # Connect the points
        loc_points.append([xp, xn, yp, yn, zp, zn])
        loc_tris.append([[0, 2, 4], [0, 2, 5], [0, 3, 4], [0, 3, 5], [1, 2, 4], [1, 2, 5], [1, 3, 4], [1, 3, 5]])
    # Set the number of triangles and points for the off file record
    num_verts, num_tris = 8 * len(loc_points), 6 * len(loc_tris)
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i in range(len(verts)):
            # Go through the points on the surface
            for point in loc_points[i]:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i in range(len(verts)):
            tri_list = loc_tris[i]
            # Go through the triangles in the surface
            for j in range(len(tri_list)):
                # Get the triangle and colors
                tri = tri_list[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(loc_points[i])
