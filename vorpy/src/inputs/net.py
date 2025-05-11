import csv
import os.path
import numpy as np
from vorpy.src.inputs.integrate_net import integrate_net


def read_net(net, file_name):
    """
    Read and process a network file into a network object.

    This function parses a network file containing vertex, edge, and surface data,
    and integrates it into a network object. The file format includes:
    - Network settings (type, resolution, vertex limits, etc.)
    - Vertex data with coordinates and atom indices
    - Edge data connecting vertices
    - Surface data with triangulation information
    - Connection data linking vertices to edges and surfaces

    Parameters:
    -----------
    net : Network
        The network object to populate with the parsed data
    file_name : str
        Path to the network file to be read

    Returns:
    --------
    Network
        The modified network object containing the integrated data

    Notes:
    ------
    - Supports CSV format with specific section headers (n, v, e, s, c)
    - Integrates vertices, edges, surfaces, and their connections
    - Updates network settings and metrics
    """
    # Check the file_name
    if os.path.exists(file_name):
        file = file_name
    else:
        return
    # Set up the data lists
    verts, edges, surfs, cons = np.array([]), np.array([]), np.array([]), np.array([])
    # Open the file
    with open(file, 'r') as net_file:
        # Get the file element array to read
        nt_fl = csv.reader(net_file, delimiter=",")
        # Set the read type to header
        reading = ""
        # Go through the lines in the read_file
        for i, line in enumerate(nt_fl):
            # Read the first line
            if line[0] in {'n', 'v', 'e', 's', 'c'}:
                reading = line[0]
                continue
            if reading == 'net':
                net.settins['net_type'], net.settings['surf_res'], net.settings['max_vert'], net.settings['ox_sizeb'] = \
                    [int(line[0]), line[1], float(line[2]), float(line[3]), float(line[4])]
                # Read the verts
            elif reading == 'v':
                # Add the data
                np.append(verts, line)
            # Read the edges
            elif reading == 'e':
                # Add the edge data
                np.append(edges, line)
            # Read the surfaces
            elif reading == 's':
                # Surface points
                if line[0] == 'pts':
                    np.append(surfs, ({"atoms": {*line[1:3]}}))
                    surfs[-1]['points'] = line[3:]
                # Triangles
                elif line[0] == 'tris':
                    surfs[-1]['tris'] = line[3:]
            # Read the connections
            elif reading == 'c':
                # Add the connections
                np.append(cons, line)
    # Integrate the data
    integrate_net(net, verts, edges, surfs, cons)
    # Return the network
    return net


# Input index function. Takes in an index file and loads it into the list of indices
def read_ndx(sys, file=None):
    """
    Read and process an index file into a system object.

    This function parses index files, which contain atom group definitions and metadata.
    The function converts the index data into a standardized format for further processing.

    The index format includes:
    - Group names in square brackets
    - Atom indices for each group
    - Multiple groups can be defined

    Parameters:
    -----------
    sys : System
        The system object to populate with index data
    file : str, optional
        Path to the index file. If None, uses sys.ndx_file

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom groups from index data
        - Storing group names in sys.ndx_names
        - Storing atom indices in sys.ndxs
    """
    # If no file is provided, check the system
    if file is None:
        file = sys.ndx_file
    # Get the file information and make sure to close the file when done
    try:
        with open(file, 'r') as f:
            my_file = f.readlines()
    except FileNotFoundError:
        return
    # Set up the indices lists and the current index
    curr_ndx = -1
    indices = []
    names = []
    # Go through the lines in the file
    for line in my_file:
        # Split the line into
        line = line.split()
        # Add the
        if line[0] == "[":
            curr_ndx += 1
            names.append([line[1]])
        else:
            for i in range(len(line)):
                indices[curr_ndx].append(line[i])
    # Set the systems indices
    sys.ndx_names = names
    # Set the systems indices
    sys.ndxs = [[sys.atoms[ndx] for ndx in indices[i]] for i in range(len(indices))]
