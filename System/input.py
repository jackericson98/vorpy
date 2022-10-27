from System.atom import Atom, get_radius
from System.Network.network import Vertex, Edge, Surface


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
def read_pdb(sys):
    # Get the file information and make sure to close the file when done
    with open(sys.base_file, 'r') as f:
        file = f.readlines()
    # Add the system name and reset the atoms and data lists
    sys.name = get_name(sys.base_file)
    # Set up the atom and the data lists
    atoms, data = [], []
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
                        element=line[76:78], residue=line[17:20], chain=line[21], res_seq=line[22:26], name=line[12:16],
                        ocp=line[54:60], t_fact=line[60:66], seg_id=line[72:76], charge=line[78:80])
            # If no chain is specified, set the chain to 'None'
            if atom.chain == ' ':
                atom.chain = 'ZZ'
            # Add the atom to the
            atoms.append(atom)
        # If the line is not an atom line store the other data
        else:
            data.append(file[i].split())
    # Return the atoms and the data
    return atoms, data


# Read cif function. Interprets the data in a cif file
def read_cif(sys):
    # Get the file information and make sure to close the file when done
    with open(sys.base_file, 'r') as f:
        file = f.readlines()
    # Get the starting number for the line
    num = int(file)
    # Go through each line of the file
    for i in range(len(file)):
        # Split the line
        file[i] = file[i].split()
        # Add the atoms
        if file[i] == int(num) and len(file[i]) >= 7:
            sys.atoms.append(Atom([file[i][9], file[i][10], file[i][11]],
                                  get_radius(file[i][3]), element=file[i][3]))


# Read gro method. Interprets the data from a .cif file type
def read_gro(sys):
    # Get the file information and make sure to close the file when done
    with open(sys.base_file, 'r') as f:
        file = f.readlines()
    # Go through each line in the file and create an atom object
    for line in file[2:-2]:
        sys.atoms.append(Atom([line[3], line[4], line[5]], get_radius(line[1][0]), element=line[1][0]))


# Read mol method. Interprets the data from a .mol file type
def read_mol(sys):
    # Get the file information and make sure to close the file when done
    with open(sys.base_file, 'r') as f:
        file = f.readlines()
    # Go through the lines in the file
    for line in file:
        # If the line is an atom line add the data
        if len(line) > 6:
            # Add the data
            sys.atoms.append(Atom([line[0], line[1], line[2]], get_radius(line[3]), element=line[3]))


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


# Input index function. Takes in an index file and loads it into the list of indices
def input_index(sys):
    # Get the file information and make sure to close the file when done
    with open(sys.index_file, 'r') as f:
        file = f.readlines()
    # Set up the indices lists and the current index
    curr_ndx = -1
    indices = []
    names = []
    # Go through the lines in the file
    for line in file:
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


# Import network function. Imports vorpy-created text document and creates network objects
def import_net(net, filename, verts_only=False):

    # Open the file
    file = open(filename).readlines()
    # Instantiate the lists
    net.verts, net.edges, net.surfs = [], [], []
    # Instantiate the current objects
    curr_atom, curr_vert, curr_edge, curr_surf = Atom(), Vertex(), Edge(), Surface()
    perim_len = 0
    # Go through the file, line by line
    for i in range(len(file)):
        # Get the line
        line = file[i]
        # Check for empty lines
        if len(line) == 0:
            continue
        # Split the line
        line = line.split()

        # Network
        if len(line) == 0:
            continue
        # Check for the network signifier
        elif line[0].lower() == 'netw':
            # Load the network information
            net.min_dist = float(line[1])
            net.beta_val = float(line[2])
            net.box_size = float(line[3])
            net.my_time = float(line[4])
            net.cpu_time = float(line[5])
            net.sol_verts = bool(line[6])
            net.curved_faces = bool(line[7])
            net.flat_faces = bool(line[8])

        # Atoms
        # Check to see if the line is an atom line
        elif line[0].lower() == "atom":
            # Get the atom from the network
            myAtom = net.atoms[int(line[1])]
            # Set attributes for the atom
            myAtom.box = [int(_) for _ in line[2:5]]
            myAtom.cell_vol = float(line[5])
            curr_atom = myAtom
        # Add the connections for the atom
        elif line[0].lower() == "acon":
            # Get the indices of the objects by checking going through the next 3 lines
            curr_atom.load_ndxs.append([int(_) for _ in line[1:]])

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
            myVert = Vertex(atoms=atoms, net=net, location=loc, radius=rad, doublet=dub, loc2=loc2, rad2=rad2, ndx=ndxs)
            curr_vert = myVert
            # Add the vertex to the network
            net.verts.append(myVert)
        # Get the vertex connections
        elif line[0].lower() == 'vcon':
            # Get the indices of the objects by checking going through the next 2 lines
            curr_vert.load_ndxs.append([int(_) for _ in line[1:]])

        # Edges
        # Check for if the line is an edge line or not
        elif line[0].lower() == "edge":
            # Quick check for verts_only
            if verts_only:
                net.build()
                return
            # Get the indices and atoms for the edges
            ndxs = [int(_) for _ in line[2:5]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Create the edge
            myEdge = Edge(atoms, net)
            # Add the points, location and radius to the edge
            myEdge.loc, myEdge.rad = [float(_) for _ in line[5:8]], float(line[8])
            myEdge.pv0, myEdge.pv1 = [float(_) for _ in line[9:12] if _ != 'None'], \
                                     [float(_) for _ in line[12:15] if _ != 'None']
            myEdge.points = []
            # Check for doubletness
            if line[15] == 'True':
                myEdge.doublet = True
            # Add the edge to the network
            net.edges.append(myEdge)
            curr_edge = myEdge
        # Edge connections
        elif line[0].lower() == 'econ':
            # Get the indices of the objects by checking going through the next 2 lines
            curr_edge.load_ndxs.append([int(_) for _ in line[1:]])
        # Edge points
        elif line[0].lower() == 'epnt':
            curr_edge.points.append([float(_) for _ in line[1:]])

        # Surfaces
        # Check for if the line is for a surface
        elif line[0].lower() == "surf":
            # Get the indices and atoms for the edges
            ndxs = [int(_) for _ in line[2:4]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Create the surface and add it to the network
            mySurf = Surface(atoms=atoms, net=net)
            # Get the surface normal and the surface area
            mySurf.rn, mySurf.sa, mySurf.points, mySurf.tris = [float(_) for _ in line[4:7]], float(line[7]), [], []
            # Set the current surface
            curr_surf = mySurf
            perim_len = int(line[8])
            net.surfs.append(mySurf)
        # Add the surface connections
        elif line[0].lower() == 'scon':
            # Get the indices of the objects by checking going through the next 2 lines
            curr_surf.load_ndxs.append([int(_) for _ in line[1:]])
        # Surface point
        elif line[0].lower() == 'spnt':
            curr_surf.points.append([float(_) for _ in line[1:]])
            if len(curr_surf.points) == perim_len:
                curr_surf.perimeter = curr_surf.points
        # Surface triangles
        elif line[0].lower() == 'stri':
            curr_surf.tris.append([int(_) for _ in line[1:]])

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
