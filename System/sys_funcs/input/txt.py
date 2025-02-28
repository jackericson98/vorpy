import pandas as pd
import numpy as np
from System.sys_objs.atom import make_atom


def read_txt(sys, file=None):
    """
    Read text for atoms
    """
    # If no file is specified add the base file from the system
    if file is None:
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

