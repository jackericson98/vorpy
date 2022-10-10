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
    file = open(sys.base_file).readlines()
    sys.sys_file_name = get_name(sys.base_file)
    atoms = []
    data = []
    # Go through each line in the file and check if the first word is the word we are looking for
    for i in range(len(file)):
        # Check to make sure the line isn't empty
        if len(file[i]) == 0:
            continue
        # Pull the file line and first word
        line = file[i]
        word = line[:4].lower()
        # Check to see if the line is an atom line
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
        # If the line is not an atom line store the other data
        else:
            data.append(file[i].split())
    # Return the atoms and the data
    return atoms, data


# Get cif function. Finds the data in a cif file
def read_cif(sys):
    # Get the system file
    sys.base_file = open(sys.base_file).readlines()
    num = int(sys.base_file[0][4:])
    # Go through each line of the file
    for i in range(len(sys.base_file)):
        # Split the line
        sys.base_file[i] = sys.base_file[i].split()
        # Add the atoms
        if sys.base_file[i] == int(num) and len(sys.base_file[i]) >= 7:
            sys.atoms.append(Atom([sys.base_file[i][9], sys.base_file[i][10], sys.base_file[i][11]], get_radius(sys.base_file[i][3]),
                                  symbol=sys.base_file[i][3]))


# Get gro method. Finds data in a gro file
def read_gro(sys):
    sys.base_file = open(sys.base_file).readlines()
    sys.info['header'] = sys.base_file[0]
    # Go through each line in the file and create an atom object
    for line in sys.base_file[2:-2]:
        sys.atoms.append(Atom([line[3], line[4], line[5]], get_radius(line[1][0]), symbol=line[1][0]))


# Get mol method. Finds data in a mol file
def read_mol(sys):
    sys.base_file = open(sys.base_file).readlines()
    for line in sys.base_file:
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
    net.verts, net.edges, net.surfs = [], [], []

    # Go through the file, line by line
    for i in range(len(file)):

        ####################################### Prepare the line #######################################################

        # Get the line
        line = file[i]
        # Split the information in the file
        line = line.split()
        # Check for empty lines
        if len(line) == 0:
            continue

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
            # Get atom index lines and split them
            acon_lines = file[i + 1: i + 4]
            acon_lines = [line.split() for line in acon_lines]
            # Get the indices of the objects by checking going through the next 3 lines
            myAtom.load_ndxs = [[int(_) for _ in line[1:]] for line in acon_lines]
            # Skip the next 3 lines
            i += 3

        # Vertices

        # Check for if the line is a vertex line
        elif line[0].lower() == "vert":
            # Get the indices of the atoms in the vertex and then the atoms themselves
            ndxs = [int(_) for _ in line[2:6]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Get the location and radius of the vertex
            loc, rad, loc2, rad2 = [float(_) for _ in line[6:9]], float(line[9]), None, None
            # Get the doublet information
            dub = False
            if line[10].lower() == 'true':
                dub, loc2, rad2 = True, [float(_) for _ in line[11:14]], float(line[14])
            # Set up the default vertex
            myVert = Vertex(atoms=atoms, net=net, location=loc, radius=rad, doublet=dub, loc2=loc2, rad2=rad2)
            # Get the vertex index lines
            vcon_lines = file[i + 1: i + 3]
            vcon_lines = [line.split() for line in vcon_lines]
            # Get the indices of the objects by checking going through the next 2 lines
            myVert.load_ndxs = [[int(_) for _ in line[1:]] for line in vcon_lines]
            # Add the vertex to the network
            net.verts.append(myVert)
            # Skip the next 2 lines
            i += 2

        # Edges

        # Check for if the line is an edge line or not
        elif line[0].lower() == "edge":
            # Get the indices and atoms for the edges
            ndxs = [int(_) for _ in line[2:5]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Create the edge
            myEdge = Edge(atoms, net)
            # Get the edge index lines
            econ_lines = file[i + 1: i + 3]
            econ_lines = [line.split() for line in econ_lines]
            # Get the indices of the objects by checking going through the next 2 lines
            myEdge.load_ndxs = [[int(_) for _ in line[1:]] for line in econ_lines]
            # Skip the next 2 lines
            i += 2
            # Set up the points list
            points = []
            # Get the points for the edge
            while file[i + 1][0].lower() == 'edpt':
                points.append([float(_) for _ in file[i + 1][1:]])
                i += 1
            # Add the points, location and radius to the edge
            myEdge.points, myEdge.loc, myEdge.rad = points, [float(_) for _ in line[5:8]], float(line[8])
            myEdge.pv0, myEdge.pv1 = [float(_) for _ in line[9:12] if _ != 'None'], [float(_) for _ in line[12:15] if _ != 'None']
            # Check for doubletness
            if line[15] == 'True':
                myEdge.doublet = True
            # Add the edge to the network
            net.edges.append(myEdge)

        # Surfaces

        # Check for if the line is an edge line or not
        elif line[0].lower() == "surf":

            # Get the indices and atoms for the edges
            ndxs = [int(_) for _ in line[2:4]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Create the surface and add it to the network
            mySurf = Surface(atoms=atoms, net=net)
            # Get the surface normal and the surface area
            mySurf.rn, mySurf.sa = [float(_) for _ in line[4:7]], float(line[7])
            # Get the surface index lines
            scon_lines = file[i + 1: i + 3]
            scon_lines = [line1.split() for line1 in scon_lines]
            # Get the indices of the objects by checking going through the next 2 lines
            mySurf.load_ndxs = [[int(_) for _ in line1[1:]] for line1 in scon_lines]
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
            mySurf.points, mySurf.tris, mySurf.perimeter = points, tris, points[:int(line[8])]
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
        # Reset the vertices edges and surfaces
        vert.edges, vert.surfs = [], []
        # Connect the atom's verts
        for i in range(len(vert.load_ndxs[0])):
            vert.edges.append(net.edges[vert.load_ndxs[0][i]])
        # Connect the atom's edges
        for i in range(len(vert.load_ndxs[1])):
            vert.surfs.append(net.surfs[vert.load_ndxs[1][i]])

    # Connect the edge objects
    for edge in net.edges:
        # Reset the edge's vertices and surfaces
        edge.verts, edge.surfs = [], []
        # Connect the edge's vertices
        for i in range(len(edge.load_ndxs[0])):
            edge.verts.append(net.verts[edge.load_ndxs[0][i]])
        # Connect the edge's surfaces
        for i in range(len(edge.load_ndxs[1])):
            edge.surfs.append(net.surfs[edge.load_ndxs[1][i]])

    # Connect the surface objects
    for surf in net.surfs:
        # Reset the surfaces and edges
        surf.verts, surf.edges = [], []
        # Connect the surface's vertices
        for i in range(len(surf.load_ndxs[0])):
            surf.verts.append(net.verts[surf.load_ndxs[0][i]])
        # Connect the surface's edges
        for i in range(len(surf.load_ndxs[1])):
            surf.edges.append(net.edges[surf.load_ndxs[1][i]])
