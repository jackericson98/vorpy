from objects import System, Atom, Vertex
import numpy as np


def read_pdb(file):
    """Function to read the info from a pdb file and create a molecule object"""
    # Open and read the file
    file = open(file).readlines()
    # Van der Waals radii for later translation
    vdw_rads = {'H': 1.1, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'P': 1.8, 'S': 1.8, 'Na': 2.27}
    # Create molecule object
    mySys = System()
    # Get non-pertinent info from the pdb file. This will likely be excluded for space and time once in production
    mySys.info = {
        "header": get_pdb_data(file, 'HEADER'),
        "title": get_pdb_data(file, 'TITLE'),
        "compound": get_pdb_data(file, 'COMPOUND'),
        "source": get_pdb_data(file, 'SOURCE'),
        "key_words": get_pdb_data(file, 'KEYWDS'),
        "exp_data": get_pdb_data(file, 'EXPDTA'),
        "author": get_pdb_data(file, 'AUTHOR'),
        "revisions": get_pdb_data(file, 'REVDAT'),
        "journal": get_pdb_data(file, 'JRNL'),
        "remarks": get_pdb_data(file, 'REMARK'),
        "debreif": get_pdb_data(file, 'DBREF'),
        "seqadv": get_pdb_data(file, 'SEQADV'),
        "formul": get_pdb_data(file, 'FORMUL'),
        "residues": get_pdb_data(file, 'SEQRES'),
        "helix": get_pdb_data(file, 'HELIX'),
        "sheet": get_pdb_data(file, 'SHEET'),
        "crystal": get_pdb_data(file, 'CRYST'),
        "origin": get_pdb_data(file, 'ORIG'),
        "scale": get_pdb_data(file, 'SCALE'),
        "terminals": get_pdb_data(file, 'TER'),
        "het_atom": get_pdb_data(file, 'HETATM'),
        "master": get_pdb_data(file, 'MASTER'),
    }
    # Grab all the lines that start with ATOM
    atoms = get_pdb_data(file, 'ATOM')
    # Make a sphere object in our Molecule's atoms list for each line in atoms
    for atom in atoms:
        # Create sphere object with radius grabbed from van der waals dictionary and coordinates from the data
        mySys.atoms.append(Atom([float(atom[-6]), float(atom[-5]), float(atom[-4])], vdw_rads[atom[-1]]))
    # Return the molecule we created
    return mySys


def get_pdb_data(file, word):
    """Function that goes through each line of the input file and returns the lines corresponding to the word. The
    word must be an exact match."""
    # Grab the test_data from the specific word
    data = []
    # Go through each line in the file and check if the first word is the word we are looking for
    for i in range(len(file)):
        if word in file[i][:len(word)]:  # check the first len(word) letters
            # Add the split test_data to our list and remove the word at the beginning of the list
            data.append(file[i].split()[1:])
    # Return all matching data
    return data


# Random system function. Creates a system with atoms placed in random locations with random radii
def random_system(anums=30, dmax=15, rmax=1):
    # Instantiate the system
    mySys = System()
    # Create the atoms
    for i in range(anums):
        # Choose a random set of 3 numbers between dmax and -dmax. Choose a random radius between 0 and rmax
        mySys.atoms.append(Atom(np.random.rand(3)*2*dmax - dmax, np.random.rand()*rmax))
    # Return the system
    return mySys
