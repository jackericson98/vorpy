import csv
import os.path

from System.sys_objs.atom import Atom, get_radius
from System.Network.network import Network
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge
from System.Network.net_objs.surface import Surface


# Get name method. Strips the location and extension from the file
def get_name(file):
    # Set up the file name variable
    filename = ""
    i = -1
    # Go through each char in the path from the back and stop at the first slash
    while file[i] != "/":
        filename = filename + file[i]
        i -= 1
    # Reverse to normal and trim the extension and the dot
    return filename[::-1][:-4]


# Read pdb function. Interprets pdb data into a system of atom objects
def read_pdb(sys, file=None):

    # Check to see if the file is provided and use the base file if not
    if file is None and sys.base_file[-3:] == 'pdb':
        file = sys.base_file
    if os.path.exists(file) and file[0] == '.' and sys.vpy_dir is not None:
        file_address = sys.vpy_dir + file[1:]
    elif os.path.exists(file):
        file_address = file
    elif sys.vpy_dir is not None and os.path.exists(sys.vpy_dir + file):
        file_address = sys.vpy_dir + file
    elif sys.dir is not None and os.path.exists(sys.dir + file):
        file_address = sys.dir + file
    elif sys.dir is not None and os.path.exists(sys.dir + file[1:]):
        file_address = sys.dir + file[1:]
    else:
        return
    # Get the file information and make sure to close the file when done
    with open(file_address, 'r') as f:
        my_file = f.readlines()
    # Add the system name and reset the atoms and data lists
    sys.name = get_name(sys.base_file)
    # Set up the atom and the data lists
    atoms, data, atom_count = [], [], 0
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
            atom = Atom(location=[float(line[30:38]), float(line[38:46]), float(line[46:54])], radius=get_radius(line[76:78], sys), system=sys,
                        element=line[76:78], mol_class=line[17:20], chain=line[21], res_seq=line[22:26], name=line[12:16],
                        ocp=line[54:60], t_fact=line[60:66], seg_id=line[72:76], charge=line[78:80], index=atom_count)
            # If no chain is specified, set the chain to 'None'
            if atom.mol == ' ' and atom.mol_class.lower() != 'sol' and atom.res_seq.lower() != 'sol':
                atom.mol = 'MOL'
            # Add the atom to the
            atoms.append(atom)
            atom_count += 1
        # If the line is not an atom line store the other data
        else:
            data.append(my_file[i].split())
    # Return the atoms and the data
    return atoms, data


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
            sys.atoms.append(Atom([line[9], line[10], line[11]], get_radius(line[3], system=sys), element=line[3],
                                  index=i))


# Read gro method. Interprets the data from a .cif file type
def read_gro(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'gro':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Go through each line in the file and create an atom object
    for i in range(2, len(my_file) - 2):
        line = my_file[i]
        sys.atoms.append(Atom([line[3], line[4], line[5]], get_radius(line[1][0], system=sys), element=line[1][0],
                              index=i))


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
            sys.atoms.append(Atom([line[0], line[1], line[2]], get_radius(line[3], system=sys), element=line[3],
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
        myVert = Vertex(atoms=atoms, net=sys.net, ndx=ndx, location=loc, radius=rad)
        sys.net.verts.append(myVert)


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


def read_net(sys, file=None):
    # Open the file
    if file is None:
        file = sys.net_file
        if sys.net_file is None:
            return
    # Get the directory for the surfaces
    net_dir = os.path.dirname(file)
    # Keep using the same directory, this will cut down on clutter
    sys.dir = net_dir
    # Open the file
    with open(file, 'r') as my_file:
        # Get the file element array to read
        read_file = list(csv.reader(my_file, delimiter=","))
        # Get the network information
        net_verts, net_edges, net_surfs = [int(_) for _ in read_file[1][5:8]]
        # Create the network if needed
        if sys.net is None:
            sys.net = Network(sys=sys, atoms=sys.atoms)
        # Create the blank objects
        sys.net.verts = [Vertex(net=sys.net) for _ in range(net_verts)]
        sys.net.edges = [Edge(net=sys.net) for _ in range(net_edges)]
        sys.net.surfs = [Surface(net=sys.net) for _ in range(net_surfs)]
        # Add the settings
        sys.net.surf_res, sys.net.max_vert, sys.net.box_size = [float(_) for _ in read_file[1][1:4]]
        sys.net.build_surfs = bool(read_file[1][4])
        sys.net.calc_box()
        # Add the vertices
        for i in range(3, 3 + net_verts):
            print("\rloading vertices - {}%".format(round(100 * (i - 3) / net_verts, 2)), end="")
            vert = sys.net.verts[i - 3]
            vert.loc = [float(_) for _ in read_file[i][1:4]]
            vert.rad = float(read_file[i][4])
            vert.atoms = [sys.atoms[int(_)] for _ in read_file[i][5:9]]
            vert.ndx = [int(_) for _ in read_file[i][5:9]]
            vert.edges = [sys.net.edges[int(_)] for _ in read_file[i][9:14] if _ != '']
            surf_ndxs = [int(_) for _ in read_file[i][14:] if _ != '']
            vert.surfs = [sys.net.surfs[_] for _ in surf_ndxs]
            if i >= 3 and vert.ndx == sys.net.verts[i - 2].ndx:
                vert.doublet = sys.net.verts[i - 2]
                sys.net.verts[i - 2].doublet = vert
            for atom in vert.atoms:
                atom.verts.append(vert)
            for surf in vert.surfs:
                if surf.verts is None:
                    surf.verts = []
                surf.verts.append(vert)
        # Add the edges
        for i in range(4 + net_verts, 4 + net_verts + net_edges):
            print("\rloading edges - {}%".format(round(100 * (i - 4 - net_verts) / net_edges, 2)), end="")
            edge = sys.net.edges[i - 4 - net_verts]
            edge.point_refs = [int(_) for _ in read_file[i][1:4] if _ != '']
            edge.atoms = [sys.atoms[int(_)] for _ in read_file[i][4:7]]
            edge.ndx = [int(_) for _ in read_file[i][4:7]]
            edge.verts = [sys.net.verts[int(_)] for _ in read_file[i][7:9]]
            edge.surfs = [sys.net.surfs[int(_)] for _ in read_file[i][9:] if _ != '']
            for atom in edge.atoms:
                atom.edges.append(edge)
            for surf in edge.surfs:
                if surf.edges is None:
                    surf.edges = []
                surf.edges.append(edge)
        # Add the surfaces
        # noinspection PyTypeChecker
        for i in range(5 + net_verts + net_edges, 5 + net_verts + net_edges + net_surfs):
            print("\rloading surfaces - {}%".format(round(100 * (i - 5 - net_verts - net_edges) / net_surfs, 2)), end="")
            surf = sys.net.surfs[i - 5 - net_verts - net_edges]
            surf.atoms = [sys.atoms[int(_)] for _ in read_file[i][5:7]]
            if surf.atoms[0].rad > surf.atoms[1].rad:
                surf.atoms[0], surf.atoms[1] = surf.atoms[1], surf.atoms[0]
            surf.ndx = [int(_) for _ in read_file[i][5:7]]
            if read_file[i][1] != '':
                surf.file = read_file[i][1]
            if read_file[i][2] != '':
                surf.res = float(read_file[i][2])
            if read_file[i][3] != '':
                surf.sa = float(read_file[i][3])
            if read_file[i][4].isdigit():
                surf.curv = float(read_file[i][4])
            if isinstance(read_file[i][16], tuple):
                surf.func = [float(_) for _ in read_file[i][7:16]] + [float(_) for _ in read_file[i][16:]]
            else:
                surf.func = [float(_) for _ in read_file[i][7:]]
            for atom in surf.atoms:
                atom.surfs.append(surf)
        if surf.atoms[0].rad > surf.atoms[1].rad:
            surf.atoms[0], surf.atoms[1] = surf.atoms[1], surf.atoms[0]
    # Set the network to connected
    sys.net.connect_net = False
