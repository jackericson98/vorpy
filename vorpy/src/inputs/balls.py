import numpy as np
import pandas as pd
from vorpy.src.objects import make_atom


def read_balls(sys):
    """
    Function to read and process ball data into a system object.
    
    Parameters:
    -----------
    sys : System
        The system object to populate with ball data
        
    Returns:
    --------
    None
    
    Notes:
    ------
    - Converts input atom data into ball objects with locations and radii
    - Creates a pandas DataFrame to store the ball objects
    - Initializes empty lists for residues and chains
    - Sets the system name to 'balls'
    """
    # Initialize empty list to store ball objects
    balls = []

    # Iterate over each atom in the system
    for i, ball in enumerate(sys.atoms):
        # Extract location and radius from atom data
        loc = np.array([float(_) for _ in ball[0]])
        rad = float(ball[1])

        # Create ball object using make_atom function
        balls.append(make_atom(sys, loc, rad, i))

    # Create pandas DataFrame to store ball objects
    sys.balls = pd.DataFrame(balls)

    # Initialize empty lists for residues and chains
    sys.residues = []
    sys.chains = []

    # Set system name to 'balls'
    sys.name = 'balls'
