import pandas as pd
import numpy as np
from vorpy.src.objects import make_atom


def read_txt(sys, file=None):
    """
    Reads atom data from a text file and populates the system with atom objects.

    This function processes a text file containing atom coordinates and radii, creating
    atom objects that are added to the system. The file can be either space-delimited
    or comma-delimited, with each line containing:
    - X coordinate
    - Y coordinate
    - Z coordinate
    - Radius

    Parameters:
    -----------
    sys : System
        The system object to populate with atoms
    file : str, optional
        Path to the text file containing atom data. If None, uses sys.files['base_file']

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom objects from the file data
        - Adding atoms to sys.balls
        - Initializing empty lists for sys.residues and sys.chains
    """
    # If no file is specified add the base file from the system
    if file is None:
        # If the base file is not set
        if sys.files['base_file'] is None:
            # Print an error message
            print("No base file set")
            # Return
            return
        # Otherwise set the file to the base file
        file = sys.files['base_file']
    # Open the txt file and read it
    with open(file, 'r') as read_file:
        # Create the balls list
        balls = []
        # Loop through the lines in the file
        for i, line in enumerate(read_file.readlines()):
            # Split the file by commas if the file is comma delimited
            if ',' in line:
                line = line.split(',')
            else:
                line = line.split()
            # Remove any of the blank entries for the line
            line = [_ for _ in line if _ != ""]
            # Get the location
            loc = np.array([float(_) for _ in line[:3]])
            # Get the radius
            rad = float(line[3])
            # Add the ball
            balls.append(make_atom(sys, loc, rad, i))
    # Create the balls, residues and the chains
    sys.balls, sys.residues, sys.chains = pd.DataFrame(balls), [], []
