import numpy as np
from pandas import DataFrame
from vorpy.src.objects.atom import make_atom

try:
    import tkinter as tk
    from tkinter import filedialog
    HAS_TKINTER = True
except ImportError:
    HAS_TKINTER = False


def read_gro(sys, file=None):
    """
    Read and process a GROMACS (.gro) format file into a system object.

    This function parses GROMACS coordinate files, which contain atom positions and metadata
    in a fixed-width format. The function converts the GROMACS data into a standardized
    vorpy ball dataframe format for further processing.

    The GROMACS format includes:
    - Residue sequence numbers
    - Residue names
    - Atom names
    - Atom indices
    - Cartesian coordinates (x, y, z)

    Parameters:
    -----------
    sys : System
        The system object to populate with GROMACS data
    file : str, optional
        Path to the GROMACS file. If None, uses sys.files['base_file']

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom objects from GROMACS data
        - Storing atoms in a pandas DataFrame
        - Initializing empty lists for chains and residues
    """

    # Get the file if the file is not specified
    if file is None:
        file = sys.files['base_file']
    # Create the dictionary that holds the balls and the additional information
    file_dict = {'balls': [], 'Additional Lines': []}
    # Line splits
    line_splits = [0, 5, 8, 15, 20, 28, 36, 44]
    # Value types
    val_types = [int, str, str, int, float, float, float]
    # Line values
    line_vals = ['res_seq', 'res_name', 'atom_name', 'index', 'x', 'y', 'z']
    # Open the file
    with open(file, 'r') as read_file:
        # Loop through the lines
        for line in read_file.readlines():
            try:
                # Split the line into its constituent parts
                ball = {line_vals[j]: val_types[j](line[line_splits[j]: line_splits[j + 1]].strip()) for j in range(7)}
                # Make an atom
                ball = make_atom(sys, location=np.array([ball['x'], ball['y'], ball['z']]), index=ball['index'],
                                 name=ball['atom_name'], res_name=ball['res_name'])
                # Add the atom to the list
                file_dict['balls'].append(ball)
            # If the line is not a valid GROMACS line, add it to the additional lines list
            except ValueError:
                file_dict['Additional Lines'].append(line)
    # Add the information to the system
    sys.balls, sys.data = DataFrame(file_dict['balls']), file_dict['Additional Lines']
    # Initialize empty lists for chains and residues
    sys.chains, sys.residues = [], []


if __name__ == '__main__':
    # Create a root window
    root = tk.Tk()
    # Hide the root window
    root.withdraw()
    # Make the root window topmost
    root.wm_attributes('-topmost', 1)
    # Get the file
    my_file = filedialog.askopenfilename()
    # Read the GROMACS file
    read_gro(sys=None, file=my_file)

