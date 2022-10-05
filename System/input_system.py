from System.Network.network import *
from System.atom import Atom, get_radius
import os


# Get name method. Extracts the name from the file name
def get_name(file):
    if not file:
        return
    filename = ""
    i = -1
    # Go through each char in the path from the back and stop at the first slash
    while file[i] != "/":
        filename = filename + file[i]
        i -= 1
    # Reverse to normal and trim the extension and the dot
    return filename[::-1][:-4]


# Read PDB atom data function. Gets atom information from pdb data line
def read_pdb_atom(line):
    # Make sure that the line starts with ATOM
    if line[:4].lower() != 'atom':
        return
    # Create the atom
    atom = Atom([float(line[30:38]), float(line[38:46]), float(line[46:54])], get_radius(line[76:78]),
                symbol=line[76:78], res=line[17:20], chain=line[21], res_seq=line[22:26])
    # If no chain is specified, set the chain to 'None'
    if atom.chain == ' ':
        atom.chain = 'Mol'
    # Return the atom
    return atom


# Get pdb data method. Finds the lines of the file with prefixes and returns them as a list
def get_pdb_data(sys, word):
    # Get the file information
    sys.file = open(sys.file_address).readlines()
    sys.sys_file_name = get_name(sys.file_address)
    # Special case for Atom lines
    if word.lower() == 'atom':
        atoms = []
    # Go through each line in the file and check if the first word is the word we are looking for
        for i in range(len(sys.file)):
            line = sys.file[i]
            if line and line[:4].lower() == 'atom':  # Check if the line starts with atom
                atom = read_pdb_atom(line)
                atoms.append(atom)
        return atoms
    # Standard case
    else:
        data = []
        for i in range(len(sys.file)):
            if word == sys.file[i][:len(word)]:  # check the first len(word) letters
                # Add the split test_data to our list and remove the word at the beginning of the list
                data.append(sys.file[i].split()[1:])
    return data


# Get pdb method. Finds the atoms and
def get_pdb(sys):
    # .PDB file type standards.
    keys = ['HEADER', 'TITLE', 'COMPOUND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK',
            'DBREF', 'SEQADV', 'FORMUL', 'SEQRES', 'HELIX', 'SHEET', 'CRYST', 'ORIG', 'SCALE', 'TER', 'HETATM',
            'MASTER']
    # Define the keys
    pdb_stds = ['header', 'title', 'compound', 'source', 'key_words', 'exp_data', 'author', 'revisions', 'journal',
                'remarks', 'debrief', 'seq_adv', 'formula', 'residues', 'helix', 'sheet', 'crystal', 'origin',
                'scale', 'terminals', 'het_atom', 'master']
    # Set the keys
    for i in range(len(pdb_stds)):
        sys.info[keys[i]] = get_pdb_data(sys, pdb_stds[i])
    # Grab the lines that start with ATOM and create Atom objects
    sys.atoms = get_pdb_data(sys, 'ATOM')


# Get cif function. Finds the data in a cif file
def get_cif(sys):
    # Get the system file
    sys.file = open(sys.file_address).readlines()
    num = int(sys.file[0][4:])
    # Go through each line of the file
    for i in range(len(sys.file)):
        # Split the line
        sys.file[i] = sys.file[i].split()
        # Add the atoms
        if sys.file[i] == int(num) and len(sys.file[i]) >= 7:
            sys.atoms.append(Atom([sys.file[i][9], sys.file[i][10], sys.file[i][11]], get_radius(sys.file[i][3]),
                                  symbol=sys.file[i][3]))


# Get gro method. Finds data in a gro file
def get_gro(sys):
    sys.file = open(sys.file_address).readlines()
    sys.info['header'] = sys.file[0]
    # Go through each line in the file and create an atom object
    for line in sys.file[2:-2]:
        sys.atoms.append(Atom([line[3], line[4], line[5]], get_radius(line[1][0]), symbol=line[1][0]))


# Get mol method. Finds data in a mol file
def get_mol(sys):
    sys.file = open(sys.file_address).readlines()
    for line in sys.file:
        if len(line) > 6:
            sys.atoms.append(Atom([line[0], line[1], line[2]], get_radius(line[3]), symbol=line[3]))


# Add vertices function. Takes in a system and a file with vertices in it and adds the verts to the system
def add_verts(sys, file_address):
    # Reset the network and open the network
    sys.net.verts, sys.net.surfs, sys.net.edges = [], [], []
    vert_file = open(file_address).readlines()
    # Go through each of the vertices file
    for i in range(len(vert_file)):
        # Set up the line variable and split it
        line = vert_file[i]
        line = line.split()
        # Set uo the line2 variable
        line2 = None
        # If there is another line after this one, check it for the same atoms
        if i + 1 < len(vert_file):
            line2 = vert_file[i + 1]
            line2 = line2.split()
        atoms = [sys.atoms[int(line[_])] for _ in range(4, 8)]
        # Check if the next line has the same atom indices as the current line
        if line2 is not None and atoms == [sys.atoms[int(line2[_])] for _ in range(4, 8)]:
            print("Doublet")
            # Doublet vertex
            my_vert = Vertex(atoms, location=[float(line[0]), float(line[1]), float(line[2])], radius=float(line[3]),
                             net=sys.net, doublet=True, loc2=[float(line2[0]), float(line2[1]),
                             float(line2[2])], rad2=float(line2[3]))
            # Skip the next line
            i += 1
        else:
            # Regular vertex
            my_vert = Vertex(atoms, location=[float(line[0]), float(line[1]), float(line[2])], radius=float(line[3]),
                             net=sys.net)
        # Add the vertex to the system
        sys.net.verts.append(my_vert)


# Add vertices function. Takes in a system and a file with vertices in it and adds the verts to the system
def add_grant_verts(sys, file_address):
    # Reset the network and open the network
    sys.net.verts, sys.net.surfs, sys.net.edges = [], [], []
    vert_file = open(file_address).readlines()
    # Go through each of the vertices file
    for line in vert_file[1:]:
        line = line.split(',')
        atoms = [sys.atoms[int(line[i]) - 1] for i in range(4)]
        loc = [float(line[-3][1:]), float(line[-2]), float(line[-1][:-2])]
        rad = calc_dist(atoms[0].loc, loc) - atoms[0].rad
        vert = Vertex(atoms, location=loc, radius=rad, net=sys.net)
        sys.net.verts.append(vert)


# Add Voronota data method. Takes in voronota data and adds it to the System
def add_vta_data(sys, ball_file, vert_file):
    # Set the voronota system indicator to True
    sys.net.flat_faces = True
    # Create the System and load the files
    vert_file = open(vert_file).readlines()
    ball_file = open(ball_file).readlines()
    # Interpret the balls
    balls = []
    for i in range(len(ball_file)):
        # Split the data
        data = ball_file[i].split(" ")
        # Grab the data reference for the atoms
        balls.append(sys.atoms[int(data[5])])
    # Interpret the vertices
    for i in range(len(vert_file)):
        # Split the data
        data = vert_file[i].split(" ")
        # Add the vertex data
        loc, rad = [float(data[4]), float(data[5]), float(data[6])], float(data[7])
        atoms = [balls[int(data[0])], balls[int(data[1])], balls[int(data[2])], balls[int(data[3])]]
        myVert = Vertex(atoms=atoms, net=sys.net, location=loc, radius=rad)
        sys.net.verts.append(myVert)

