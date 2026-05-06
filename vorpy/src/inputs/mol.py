import numpy as np
from vorpy.src.objects.atom import make_atom
from pandas import DataFrame


def read_mol(sys, file=None):
    """
    Read and process a MOL format file into a system object.

    This function parses MOL format files, which contain molecular structure data including:
    - Atom coordinates
    - Element information
    - Bond connectivity

    Parameters:
    -----------
    sys : System
        The system object to populate with MOL data
    file : str, optional
        Path to the MOL file. If None, uses sys.files['base_file']

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom objects from MOL data
        - Storing atoms in a pandas DataFrame
        - Storing bond information
        - Initializing empty lists for chains and residues
    """

    # Check the file variable and if it is none get the systems base file
    if file is None:
        file = sys.files['base_file']

    # Create the dictionary that holds the information from the
    file_dict = {'balls': [], 'Additional Lines': [], 'bonds': []}
    # Open the file
    with open(file, 'r') as rf:
        # Create the index for counting the atoms
        index = 0

        # Loop through the lines
        for line in rf.readlines():

            # Split the line
            line_info = line.split()

            # Check for if it is an atom dood
            if len(line_info) >= 10:

                # Pull the location
                location = np.array([float(_) for _ in line_info[:3]])

                # Create the ball
                ball = make_atom(sys, location=location, element=line[3], index=index)

                # Add the ball
                file_dict['balls'].append(ball)

                # Increment the index
                index += 1

            # If the length of the line is 4 it is the bonds
            elif len(line_info) == 4:

                # Add the bond to the
                file_dict['bonds'].append([int(_) for _ in line_info])

            # Otherwise add the
            else:
                # Add the line to the extra lines list
                file_dict['Additional Lines'].append(line)
    # Return the dataframe
    sys.balls = DataFrame(file_dict['balls'])
    sys.data = file_dict['Additional Lines']
    sys.chains, sys.residues = [], []
