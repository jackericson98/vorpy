import os.path as path
import os
import csv
from System.sys_objs.atom import Atom
from System.Network.network import Network
from System.Network.net_objs.vertex import Vertex
from System.sys_objs.residue import Residue
from System.sys_objs.chain import Chain, Sol


def read_pdb(sys, file=None):
    """
    Interprets pdb data into a system of atom objects
    :param sys: System to add the pdb information to
    :param file: .pdb file to be added to the system
    :return: list of tuples of locations and radii
    """
    # Check to see if the file is provided and use the base file if not
    if file is None and sys.base_file[-3:] == 'pdb':
        file = sys.base_file
    if path.exists(file) and file[0] == '.' and sys.vpy_dir is not None:
        file_address = sys.vpy_dir + file[1:]
    elif path.exists(file):
        file_address = file
    elif sys.vpy_dir is not None and path.exists(sys.vpy_dir + file):
        file_address = sys.vpy_dir + file
    elif sys.dir is not None and path.exists(sys.dir + file):
        file_address = sys.dir + file
    elif sys.dir is not None and path.exists(sys.dir + file[1:]):
        file_address = sys.dir + file[1:]
    else:
        return
    # Get the file information and make sure to close the file when done
    with open(file_address, 'r') as f:
        my_file = f.readlines()
    # Add the system name and reset the atoms and data lists
    sys.name = path.basename(sys.base_file)[:-4]
    # Set up the atom and the data lists
    atoms, data, atom_count = [], [], 0
    sys.chains, sys.residues = [], []
    chains, resids = {}, {}
    # Go through each line in the file and check if the first word is the word we are looking for
    for i in range(len(my_file)):
        # Check to make sure the line isn't empty
        if len(my_file[i]) == 0:
            continue
        # Pull the file line and first word
        line = my_file[i]
        word = line[:4].lower()
        # Check to see if the line is an atom line
        if line and word == 'atom':  # Check if the line starts with atom
            # Check for the "m" situation
            if line[76:78] == ' M':
                continue
            # Create the atom
            atom = Atom(location=[float(line[30:38]), float(line[38:46]), float(line[46:54])], system=sys,
                        element=line[76:78].strip(), res_seq=int(line[22:26]), name=line[12:16], seg_id=line[72:76],
                        index=atom_count)
            # Add the atom to the atoms list
            atoms.append(atom)
            atom_count += 1
            # If no chain is specified, set the chain to 'None'
            res_str, chain_str = line[17:20], line[21]
            if chain_str == ' ':
                if res_str.lower() in {'sol', 'hoh'}:
                    chain_str = 'SOL'
                elif res_str.lower() in {'cl', 'mg', 'na', 'k'} and 'SOL' in chains:
                    chain_str = 'SOL'
                else:
                    chain_str = 'A'
            # Create the chain and residue dictionaries
            res_name, chn_name = line[17:20] + str(atom.res_seq), chain_str
            # If the chain has been made before
            if chn_name in chains:
                # Get the chain from the dictionary and add the atom
                my_chn = chains[chn_name]
                my_chn.add_atom(atom)
                atom.chn = my_chn
            # Create the chain
            else:
                # If the chain is the sol chain
                if res_str.lower() == 'sol':
                    my_chn = Sol(atoms=[atom], residues=[], name=chn_name)
                    sys.sol = my_chn
                # If the chain is not sol create a regular chain object
                else:
                    my_chn = Chain(atoms=[atom], residues=[], name=chn_name)
                    sys.chains.append(my_chn)
                # Set the chain in the dictionary and give the atom it's chain
                chains[chn_name] = my_chn
                atom.chn = my_chn

            # Assign the atoms and create the residues
            if res_name in resids:
                my_res = resids[res_name]
                my_res.atoms.append(atom)
            else:
                my_res = Residue(atoms=[atom], name=res_str, sequence=atom.res_seq, chain=atom.chn)
                atom.chn.residues.append(my_res)
                resids[res_name] = my_res
                if res_str.lower() != 'sol':
                    sys.residues.append(my_res)
                else:
                    sys.sol.residues.append(my_res)
            # Assign the residue to the atom
            atom.res = my_res
        # If the line is not an atom line store the other data
        else:
            data.append(my_file[i].split())
    # Set the atoms and the data
    sys.atoms, sys.data = atoms, data


# Read cif function. Interprets the data in a cif file
def read_cif(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'cif':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Get the starting number for the line
    num = int(my_file[0].split()[0])
    # Go through each line of the file
    for i in range(len(my_file)):
        # Split the line
        line = my_file[i].split()
        # Add the atoms
        if line == int(num) and len(line) >= 7:
            sys.atoms.append(Atom([line[9], line[10], line[11]], element=line[3], index=i))


# Read gro method. Interprets the data from a .cif file type
def read_gro(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'gro':
        file = sys.base_file
    # Check that the system atoms list has been created
    if sys.atoms is None:
        sys.atoms = []
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
        for i, line in enumerate(my_file):
            line.split()
            if 2 <= i < len(my_file):
                sys.atoms.append(Atom(location=[line[3], line[4], line[5]], system=sys, index=i, name=line[1]))


# Read mol method. Interprets the data from a .mol file type
def read_mol(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'mol':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Go through the lines in the file
    for i in range(len(my_file)):
        # Get the line
        line = my_file[i]
        # If the line is an atom line add the data
        if len(line) > 6:
            # Add the data
            sys.atoms.append(Atom([line[0], line[1], line[2]], element=line[3],
                                  index=i))


# Add Voronota data method. Takes in voronota data and adds it to the System
def read_vta_data(sys, ball_file, vert_file):
    # If no network has been created, make one
    if sys.net is None:
        sys.net = Network(sys, sys.atoms, verts=[], edges=[], surfs=[])
    if sys.net.verts is None:
        sys.net.verts = []
    # Create the System and load the files
    with open(ball_file, 'r') as b:
        b_file = b.readlines()
    with open(vert_file, 'r') as v:
        v_file = v.readlines()
    # Interpret the balls
    balls = []
    for i in range(len(b_file)):
        # Split the data
        data = b_file[i].split(" ")
        # Grab the data reference for the atoms
        balls.append(sys.atoms[int(data[5]) - 1])
    # Interpret the vertices
    for i in range(len(v_file)):
        # Split the data
        data = v_file[i].split(" ")
        # Add the vertex data
        loc, rad = [float(data[4]), float(data[5]), float(data[6])], float(data[7])
        atoms = [balls[int(data[0])], balls[int(data[1])], balls[int(data[2])], balls[int(data[3])]]
        ndx = [sys.atoms.index(atom) for atom in atoms]
        ndx.sort()
        my_vert = Vertex(atoms=atoms, net=sys.net, ndx=ndx, location=loc, radius=rad)
        sys.net.verts.append(my_vert)


# Input index function. Takes in an index file and loads it into the list of indices
def read_ndx(sys, file=None):
    # If no file is provided, check the system
    if file is None:
        file = sys.ndx_file
    # Get the file information and make sure to close the file when done
    try:
        with open(file, 'r') as f:
            my_file = f.readlines()
    except FileNotFoundError:
        return
    # Set up the indices lists and the current index
    curr_ndx = -1
    indices = []
    names = []
    # Go through the lines in the file
    for line in my_file:
        # Split the line into
        line = line.split()
        # Add the
        if line[0] == "[":
            curr_ndx += 1
            names.append([line[1]])
        else:
            for i in range(len(line)):
                indices[curr_ndx].append(line[i])
    # Set the systems indices
    sys.ndx_names = names
    sys.ndxs = [[sys.atoms[ndx] for ndx in indices[i]] for i in range(len(indices))]


# Import vertices function.
def read_verts(net, file=None):
    # If file is None use the system's vertex file
    if file is None:
        file = net.sys.verts_file
    # Open the file
    try:
        with open(file) as f:
            my_file = f.readlines()
    except FileNotFoundError:
        print("\r No such file exists", end="")
        return
    # Set up the vertices list
    verts = []
    last_vert = None
    # Go through the lines in the file
    for line in my_file[1:]:
        line = line.split()
        if line[0].lower() != 'vert':
            continue
        new_vert = Vertex(atoms=[net.atoms[int(_)] for _ in line[1:5]], location=[float(_) for _ in line[5:8]],
                          radius=float(line[8]), ndx=[int(_) for _ in line[1:5]])
        verts.append(new_vert)
        if last_vert is not None and last_vert.ndx == new_vert.ndx:
            # Link the doublets
            last_vert.doublet, last_vert.loc2, last_vert.rad2 = new_vert, new_vert.loc, new_vert.rad
            new_vert.doublet, new_vert.loc2, new_vert.rad2 = last_vert, last_vert.loc, last_vert.rad
        # Assign the vertex
        last_vert = new_vert
    # Set the network's vertices
    net.verts = verts


def read_surf_file(surf, file=None):
    """
    Reads the file holding the points and the triangles for the surface
    :param surf:
    :param file: Specifies the address for the build file
    :return: The surfaces points and triangles are set
    """
    # Check to see if the file exists
    if file is None and surf.file is not None:
        file = surf.file
    # Check that the provided file works as an address on its own
    if os.path.exists(file):
        file_address = file
    # Check that the file name is a relative location to the system directory
    elif os.path.exists(surf.net.sys.dir + file):
        file_address = surf.net.sys.dir + file
    # Last brute force a location if the file name is incorrect
    else:
        return
    # Read an off file
    if file_address[-3:].lower() == 'off':
        # Open the file
        with open(file_address, 'r') as my_file:
            # Read the lines
            file_array = my_file.readlines()
            # Get the number of points and triangles
            num_points, num_tris = [int(_) for _ in file_array[1].split()[:2]]
            # Add the points
            surf.points = []
            for i in range(4, num_points + 4):
                line = file_array[i].split()
                surf.points.append([float(_) for _ in line])
            # Add the tris
            surf.tris = []
            for i in range(4 + num_points, 4 + num_points + num_tris):
                line = file_array[i].split()
                surf.tris.append([int(_) for _ in line[1:4]])
    # Read a comma separated file surface file
    elif file_address[-3:].lower() == 'csv':
        # Open the file
        with open(file_address, 'r') as my_file:
            # Get the file element array to read
            read_file = list(csv.reader(my_file, delimiter=","))
            # Get the number of points and triangles
            num_points, num_tris = [int(_) for _ in read_file[1][1:]]
            # Go through the points lines of the file
            surf.points = []
            for i in range(3, num_points + 3):
                surf.points.append([float(_) for _ in read_file[i][1:]])
            # Go through the triangles lines of the file
            surf.tris = []
            for i in range(4 + num_points, 4 + num_points + num_tris):
                surf.tris.append([int(_) for _ in read_file[i][1:]])
