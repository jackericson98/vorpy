from System.sys_objs.atom import make_atom
from pandas import DataFrame


# Read gro method. Interprets the data from a .cif file type
def read_gro(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'gro':
        file = sys.base_file
    atoms = []
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
        for i, line in enumerate(my_file):
            line.split()
            if 2 <= i < len(my_file):
                atoms.append(make_atom(location=[line[3], line[4], line[5]], system=sys, index=i, name=line[1]))
    sys.atoms = DataFrame(atoms)
