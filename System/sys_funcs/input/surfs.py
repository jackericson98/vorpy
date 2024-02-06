import numpy as np
from os import path
import csv


def read_surf_file(surf, file=None):
    """
    Reads the file holding the points and the triangles for the surface
    :param surf:
    :param file: Specifies the address for the build file
    :return: The surfaces points and triangles are set
    """
    # Check to see if the file exists
    if file is None and surf.file is not None:
        file = surf.file
    # Check that the provided file works as an address on its own
    if path.exists(file):
        file_address = file
    # Check that the file name is a relative location to the system directory
    elif path.exists(surf.net.sys.dir + file):
        file_address = surf.net.sys.dir + file
    # Last brute force a location if the file name is incorrect
    else:
        return
    # Read an off file
    if file_address[-3:].lower() == 'off':
        # Open the file
        with open(file_address, 'r') as my_file:
            # Read the lines
            file_array = my_file.readlines()
            # Get the number of points and triangles
            num_points, num_tris = [int(_) for _ in file_array[1].split()[:2]]
            # Add the points
            surf['points'] = np.array([])
            for i in range(4, num_points + 4):
                line = file_array[i].split()
                np.append(surf['points'], [float(_) for _ in line])
            # Add the tris
            surf['tris'] = np.array([])
            for i in range(4 + num_points, 4 + num_points + num_tris):
                line = file_array[i].split()
                np.append(surf['tris'], [int(_) for _ in line[1:4]])
    # Read a comma separated file surface file
    elif file_address[-3:].lower() == 'csv':
        # Open the file
        with open(file_address, 'r') as my_file:
            # Get the file element array to read
            read_file = list(csv.reader(my_file, delimiter=","))
            # Get the number of points and triangles
            num_points, num_tris = [int(_) for _ in read_file[1][1:]]
            # Go through the points lines of the file
            surf['points'] = np.array([])
            for i in range(3, num_points + 3):
                np.append(surf['points'], [float(_) for _ in read_file[i][1:]])
            # Go through the triangles lines of the file
            surf['tris'] = np.array([])
            for i in range(4 + num_points, 4 + num_points + num_tris):
                np.append(surf['tris'], [int(_) for _ in read_file[i][1:]])
