from System.Network.network import *
from System.atom import Atom, get_radius
from System.Network.connect_network import Edge, Surface


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


# Get pdb data method. Finds the lines of the file with prefixes and returns them as a list
def read_pdb(sys):
    # Get the file information
    sys.file = open(sys.file_address).readlines()
    sys.sys_file_name = get_name(sys.file_address)
    atoms = []
    data = []
    # Go through each line in the file and check if the first word is the word we are looking for
    for i in range(len(sys.file)):
        line = sys.file[i]
        word = line[:4].lower()
        if line and word == 'atom':  # Check if the line starts with atom
            # Create the atom
            atom = Atom([float(line[30:38]), float(line[38:46]), float(line[46:54])], get_radius(line[76:78]),
                        symbol=line[76:78], res=line[17:20], chain=line[21], res_seq=line[22:26], name=line[12:16],
                        ocp=line[54:60], t_fact=line[60:66], seg_id=line[72:76], charge=line[78:80])
            # If no chain is specified, set the chain to 'None'
            if atom.chain == ' ':
                atom.chain = 'Mol'
            # Add the atom to the
            atoms.append(atom)
        else:

            for i in range(len(sys.file)):
                if word == sys.file[i][:len(word)]:  # check the first len(word) letters
                    # Add the split test_data to our list and remove the word at the beginning of the list
                    data.append(sys.file[i].split()[1:])
    return atoms, data


# Get cif function. Finds the data in a cif file
def read_cif(sys):
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
def read_gro(sys):
    sys.file = open(sys.file_address).readlines()
    sys.info['header'] = sys.file[0]
    # Go through each line in the file and create an atom object
    for line in sys.file[2:-2]:
        sys.atoms.append(Atom([line[3], line[4], line[5]], get_radius(line[1][0]), symbol=line[1][0]))


# Get mol method. Finds data in a mol file
def read_mol(sys):
    sys.file = open(sys.file_address).readlines()
    for line in sys.file:
        if len(line) > 6:
            sys.atoms.append(Atom([line[0], line[1], line[2]], get_radius(line[3]), symbol=line[3]))


# Add Voronota data method. Takes in voronota data and adds it to the System
def add_vta_data(sys, ball_file, vert_file):
    # Set the voronota system indicator to True
    sys.myNet.flat_faces = True
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
        myVert = Vertex(atoms=atoms, net=sys.myNet, location=loc, radius=rad)
        sys.myNet.verts.append(myVert)


# Import network function. Imports vorpy-created text document and creates network objects
def import_net(net, filename):

    # Open the file
    file = open(filename).readlines()

    # Go through the file, line by line
    for i in range(len(file)):

        ####################################### Prepare the line #######################################################

        # Get the line
        line = file[i]
        # Check for empty lines
        if len(line) == 0:
            continue
        # Split the information in the file
        line.split()
        # If there is another line after this one, check it for the same atoms
        if len(file) > i + 1 and len(file[i + 1]) > 0:
            line2 = file[i + 1]
            line2.split()

        ######################################### Get Objects ##########################################################

        # Network

        if line[0].lower() == 'netw':

            # Load the network information
            net.min_dist = float(line[1])
            net.beta_val = float(line[2])
            net.box_size = float(line[3])
            net.sol_verts = bool(line[4])
            net.curved_faces = bool(line[5])
            net.flat_faces = bool(line[6])

        # Atoms

        # Check to see if the line is an atom line
        elif line[0].lower() == "atom":
            # Get the atom from the network
            myAtom = net.atoms[int(line[1])]
            # Set attributes for the atom
            myAtom.box = [int(_) for _ in line[2:5]]
            myAtom.cell_vol = float(line[5])
            # Get the indices of the objects by checking going through the next 3 lines
            myAtom.load_ndxs = [[int(_) for _ in line[1:]] for line in file[i + 1: i + 4]]
            # Skip the next 3 lines
            i += 3

        # Vertices

        # Check for if the line is a vertex line
        elif line[0].lower() == "vert":
            # Get the indices of the atoms in the vertex and then the atoms themselves
            ndxs = [int(_) for _ in line[1:5]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Get the location and radius of the vertex
            loc, rad = [float(_) for _ in line[5:8]], float(line[8])
            # Get the doublet information
            dub, loc2, rad2 = bool(line[9]), [float(_) for _ in line[-4:-1]], float(line[-1])
            # Set up the default vertex
            myVert = Vertex(atoms=atoms, net=net, location=loc, radius=rad, doublet=dub, loc2=loc2, rad2=rad2)
            # Get the indices of the objects by checking going through the next 2 lines
            myVert.load_ndxs = [[int(_) for _ in line[1:]] for line in file[i + 1: i + 3]]
            # Skip the next 2 lines
            i += 2
            # Add the vertex to the network
            net.verts.append(myVert)

        # Edges

        # Check for if the line is an edge line or not
        elif line[0].lower() == "edge":
            # Get the indices and atoms for the edges
            ndxs = [int(_) for _ in line[1:4]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Create the edge
            myEdge = Edge(atoms, net)
            # Get the indices of the objects by checking going through the next 2 lines
            myEdge.load_ndxs = [[int(_) for _ in line[1:]] for line in file[i + 1: i + 3]]
            # Skip the next 2 lines
            i += 2
            # Set up the points list
            points = []
            # Get the points for the edge
            while file[i + 1][0].lower() == 'edpt':
                points.append([float(_) for _ in file[i + 1][1:]])
                i += 1
            # Add the points, location and radius to the edge
            myEdge.points, myEdge.loc, myEdge.rad = points, [float(_) for _ in line[4:7]], float(line[7])
            myEdge.pv0, myEdge.pv1 = [float(_) for _ in line[8:11]], [float(_) for _ in line[11:14]]
            myEdge.doublet = bool(line[14])
            # Add the edge to the network
            net.edges.append(myEdge)

        # Surfaces

        # Check for if the line is an edge line or not
        elif line[0].lower() == "surf":

            # Get the indices and atoms for the edges
            ndxs = [int(_) for _ in line[1:3]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Create the surface and add it to the network
            mySurf = Surface(atoms=atoms, net=net)
            # Get the surface normal and the surface area
            mySurf.rn, mySurf.sa = [float(_) for _ in line[3:6]], float(line[6])
            # Get the indices of the objects by checking going through the next 2 lines
            mySurf.load_ndxs = [[int(_) for _ in line[1:]] for line in file[i + 1: i + 3]]
            # Set up the list of points and perimeter points
            points, tris = [], []
            # Get the surface points
            while file[i + 1][0].lower() == 'supt':
                # Add the points to the perimeter and increment the counter
                points.append([float(_) for _ in file[i + 1][1:]])
                i += 1

            # Get the triangles
            while file[i + 1][0].lower() == 'sutr':
                # Add the points to the perimeter and increment the counter
                tris.append([int(_) for _ in file[i + 1][1:]])
                i += 1
            # Add the points, triangles and perimeter points
            mySurf.points, mySurf.tris, mySurf.perimeter = points, tris, points[:int(line[7])]
            # Add the surface to the network
            net.surfs.append(mySurf)
    # Connect the objects
    connect_input_net(net)


# Connect input network function. Connect the network objects using their load index lists
def connect_input_net(net):

    # Connect the atom objects
    for atom in net.atoms:
        # Connect the atom's verts
        for i in range(len(atom.load_ndxs[0])):
            atom.verts.append(net.verts[atom.load_ndxs[0][i]])
        # Connect the atom's edges
        for i in range(len(atom.load_ndxs[1])):
            atom.edges.append(net.edges[atom.load_ndxs[1][i]])
        # Connect the atom's surfaces
        for i in range(len(atom.load_ndxs[2])):
            atom.surfs.append(net.surfs[atom.load_ndxs[2][i]])

    # Connect the vertex objects
    for vert in net.verts:
        # Connect the atom's verts
        for i in range(len(vert.load_ndxs[0])):
            vert.edges.append(net.edges[vert.load_ndxs[0][i]])
        # Connect the atom's edges
        for i in range(len(vert.load_ndxs[1])):
            vert.surfs.append(net.surfs[vert.load_ndxs[1][i]])

    # Connect the edge objects
    for edge in net.edges:
        # Connect the edge's vertices
        for i in range(len(edge.load_ndxs[0])):
            edge.verts.append(net.verts[edge.load_ndxs[0][i]])
        # Connect the edge's surfaces
        for i in range(len(edge.load_ndxs[1])):
            edge.surfs.append(net.surfs[edge.load_ndxs[1][i]])

    # Connect the surface objects
    for surf in net.surfs:
        # Connect the surface's vertices
        for i in range(len(surf.load_ndxs[0])):
            surf.verts.append(net.verts[surf.load_ndxs[0][i]])
        # Connect the surface's edges
        for i in range(len(surf.load_ndxs[1])):
            surf.edges.append(net.edges[surf.load_ndxs[1][i]])
