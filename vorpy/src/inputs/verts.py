import pandas as pd


def read_verts(group, file):
    """
    Reads vertex data from a file and creates a DataFrame containing vertex information.

    This function processes a file containing vertex data for a network, where each vertex
    represents a connection between atoms (balls). The file format expects:
    - First line: Network type verification
    - Subsequent lines: Vertex data with:
        * 4 ball indices
        * 3D coordinates
        * Radius
        * Additional coordinates and radius for certain vertices
    - 'END' marker to indicate completion

    Parameters:
    -----------
    group : Group
        The group object containing network settings and type information
    file : str
        Path to the vertex data file

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing vertex information with columns:
        - balls: List of 4 ball indices
        - loc: 3D coordinates
        - rad: Radius
        - loc2: Optional second set of coordinates
        - rad2: Optional second radius

    Notes:
    ------
    - Verifies network type matches the first line of the file
    - Combines related vertex data when indicated by flag (line[8] == 1)
    - Returns None if file doesn't end with 'END' marker
    """

    # Initialize an empty list to store vertex data
    verts = []
    # Initialize a flag to indicate if the file has been fully processed
    finished = False
    # Open the file for reading
    with open(file, 'r') as my_file:
        # Read each line of the file
        for i, line in enumerate(my_file.readlines()):
            # Split the line into a list of strings
            line = line.split(' ')
            # If the line is the first line, check if the network type matches
            if i == 0:
                if len(line) == 15 and group.settings['net_type'].lower() != line[-1].lower()[:-1]:
                    print("\nWarning - Loaded Vertices do not match the set network type\n\n")
                continue
            # If the line is the end marker, set the finished flag to True
            if line[0] == 'END':
                finished = True
                continue
            # If the line indicates a combined vertex, update the last vertex
            if int(line[8]) == 1 and [int(_) for _ in line[:4]] == verts[-1]['balls']:
                verts[-1]['loc2'] = [float(_) for _ in line[4:7]]
                verts[-1]['rad2'] = float(line[7])
            else:
                # Otherwise, add a new vertex to the list
                verts.append(
                    {'balls': [int(_) for _ in line[:4]], 'loc': [float(_) for _ in line[4:7]], 'rad': float(line[7]),
                     'loc2': None, 'rad2': None})
    # If the file has been fully processed, return the DataFrame
    if finished:
        return pd.DataFrame(verts)
    # Otherwise, return None
    return None
