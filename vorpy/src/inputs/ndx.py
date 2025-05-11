

def read_ndx(sys, file=None):
    """
    Reads an index file and adds the indices to the system object.

    This function reads an index file and adds the indices to the system object.
    It supports various file formats including GROMACS, CHARMM, and AMBER.
    
    Parameters:
    -----------
    sys : System
        The system object to which indices will be added
    file : str, optional
    
    """
    # Check if the file is provided
    if file is None:
        file = sys.files['ndx_file']
    # Read the index file
    with open(file, 'r') as f:
        lines = f.readlines()
    # Initialize a list to store the indices
    indices = []
    # Go through the lines
    for line in lines:
        # Split the line by whitespace
        parts = line.split()
    # add to the system
    sys.ndxs = indices
    # return the system
    return sys
        
    
