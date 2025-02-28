import numpy as np
from System.sys_objs.atom import make_atom
from pandas import DataFrame


# Read mol method. Interprets the data from a .mol file type
def read_mol1(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'mol':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    atoms = []
    # Go through the lines in the file
    for i in range(len(my_file)):
        # Get the line
        line = my_file[i]
        # If the line is an atom line add the data
        if len(line) > 6:
            # Add the data
            atoms.append(make_atom([line[0], line[1], line[2]], element=line[3], index=i))
    sys.atoms = DataFrame(atoms)


def read_mol(sys, file):
    """
    Read .mol file


    """

    # Check the file variable and if it is none get the systems base file
    if file is None:
        file = sys.files['base_file']

    # Create the dictionary that holds the information from the
    file_dict = {'balls': [], 'Additional Lines': [], 'bonds': []}

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
