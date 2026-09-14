import numpy as np
from pandas import DataFrame
from vorpy.src.input_progress import iter_input_lines
from vorpy.src.objects.atom import make_atom

try:
    import tkinter as tk
    from tkinter import filedialog
    HAS_TKINTER = True
except ImportError:
    HAS_TKINTER = False


def read_gro(sys, file=None, progress=None):
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
    file_dict = {'balls': [], 'Additional Lines': []}
    lines = iter_input_lines(file, progress)
    title = next(lines, "")
    count_line = next(lines, "")
    try:
        atom_count = int(count_line.strip())
    except ValueError as error:
        raise ValueError("Invalid GRO atom count") from error
    if atom_count < 0:
        raise ValueError("Invalid GRO atom count")
    file_dict['Additional Lines'].extend((title, count_line))
    for index in range(atom_count):
        line = next(lines, "")
        try:
            # GRO uses five-character residue/atom names and nanometres.
            # Infer coordinate precision from the decimal-point spacing.
            decimals = [i for i, char in enumerate(line[20:]) if char == '.']
            width = decimals[1] - decimals[0] if len(decimals) >= 3 else 8
            location = np.array([float(line[20 + axis * width:20 + (axis + 1) * width])
                                 for axis in range(3)]) * 10.0
            ball = make_atom(sys, location=location, index=index,
                             name=line[10:15].strip(), res_name=line[5:10].strip(),
                             res_seq=int(line[:5]))
            ball['gro_id'] = int(line[15:20])
        except (ValueError, IndexError) as error:
            raise ValueError(f"Invalid GRO atom record {index + 1}: {line.rstrip()}") from error
        file_dict['balls'].append(ball)
    file_dict['Additional Lines'].extend(lines)
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

