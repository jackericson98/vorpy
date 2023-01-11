import csv
from System.sys_objs.atom import Atom, get_radius
from System.Network.network import Network, Vertex, Edge, Surface


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
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()

    # Add the system name and reset the atoms and data lists
    sys.name = get_name(sys.base_file)
    # Set up the atom and the data lists
    atoms, data = [], []
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
                        element=line[76:78], residue=line[17:20], chain=line[21], res_seq=line[22:26], name=line[12:16],
                        ocp=line[54:60], t_fact=line[60:66], seg_id=line[72:76], charge=line[78:80])
            # If no chain is specified, set the chain to 'None'
            if atom.mol == ' ' and atom.res.lower() != 'sol' and atom.res_seq.lower() != 'sol':
                atom.mol = 'MOL'
            # Add the atom to the
            atoms.append(atom)
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
            sys.atoms.append(Atom([line[9], line[10], line[11]], get_radius(line[3], system=sys), element=line[3]))


# Read gro method. Interprets the data from a .cif file type
def read_gro(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'gro':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Go through each line in the file and create an atom object
    for line in my_file[2:-2]:
        sys.atoms.append(Atom([line[3], line[4], line[5]], get_radius(line[1][0], system=sys), element=line[1][0]))


# Read mol method. Interprets the data from a .mol file type
def read_mol(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'mol':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Go through the lines in the file
    for line in my_file:
        # If the line is an atom line add the data
        if len(line) > 6:
            # Add the data
            sys.atoms.append(Atom([line[0], line[1], line[2]], get_radius(line[3], system=sys), element=line[3]))


# Add Voronota data method. Takes in voronota data and adds it to the System
def read_vta_data(sys, ball_file, vert_file):
    # If no network has been created, make one
    if sys.net is None:
        sys.net = Network(sys, sys.atoms, verts=[], edges=[], surfs=[], flat_faces=True)
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
    sys.net.flat_faces = True


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


# Import network function. Imports vorpy-created text document and creates network objects
def read_old_net(net, file, verts_only=False):
    # Open the file
    try:
        with open(file) as f:
            my_file = f.readlines()
    except FileNotFoundError:
        print("\r No such file exists", end="")
        return
    # Instantiate the lists
    net.verts, net.edges, net.surfs = [], [], []
    # Instantiate the current objects
    curr_atom, curr_vert, curr_edge, curr_surf = Atom(), Vertex(), Edge(), Surface()
    perim_len = 0
    num_verts, num_edges, num_surfs = 0, 0, 0
    print("\rLoading Network: ", end="")
    # Go through the file, line by line
    for i in range(len(my_file)):
        # Get the line
        line = my_file[i]
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
            net.surf_res = float(line[1])
            net.max_vert = float(line[2])
            net.box_size = float(line[3])
            net.my_time = float(line[4])
            net.cpu_time = float(line[5])
            net.sol_verts = bool(line[6])
            net.curved_faces = bool(line[7])
            net.flat_faces = bool(line[8])

            num_verts, num_edges, num_surfs = [int(_) for _ in line[9:12]]


        # Atoms
        # Check to see if the line is an atom line
        elif line[0].lower() == "atom":
            # Running print statement tracking the progress of loading the atoms
            print("\rLoading Atoms: {:.2f}%       ".format(min(100, 100 * (int(line[1]) + 1) / len(net.atoms))), end="")
            # Get the atom from the network
            myAtom = net.atoms[int(line[1])]
            # Set attributes for the atom
            myAtom.box = [int(_) for _ in line[2:5]]
            myAtom.vol = float(line[5])
            curr_atom = myAtom
        # Add the connections for the atom
        elif line[0].lower() == "acon":
            # Get the indices of the objects by checking going through the next 3 lines
            curr_atom.load_ndxs.append([int(_) for _ in line[1:]])

        # Vertices
        # Check for if the line is a vertex line
        elif line[0].lower() == "vert":
            # Running print statement tracking the progress of loading the vertices
            print("\rLoading Vertices: {:.2f}%         ".format(min(100, 100 * (int(line[1]) + 1) / num_verts)), end="")
            # Get the indices of the atoms in the vertex and then the atoms themselves
            ndxs = [int(_) for _ in line[2:6]]
            atoms = [net.atoms[ndx] for ndx in ndxs]
            # Get the location and radius of the vertex
            loc, rad, loc2, rad2 = [float(_) for _ in line[6:9]], float(line[9]), None, None
            # Set up the default vertex
            myVert = Vertex(atoms=atoms, net=net, location=loc, radius=rad, loc2=loc2, rad2=rad2, ndx=ndxs)
            curr_vert = myVert
            # Get the doublet information
            if line[10].lower() == 'true':
                loc2, rad2 = [float(_) for _ in line[11:14]], float(line[14])
                if loc2 == net.verts[-1].loc:
                    net.verts[-1].doublet = myVert
                    myVert.doublet = net.verts[-1].doublet
            # Add the vertex to the network
            net.verts.append(myVert)
        # Get the vertex connections
        elif line[0].lower() == 'vcon':
            # Get the indices of the objects by checking going through the next 2 lines
            curr_vert.load_ndxs.append([int(_) for _ in line[1:]])

        # Edges
        # Check for if the line is an edge line or not
        elif line[0].lower() == "edge":
            # Running print statement tracking the progress of loading the edges
            print("\rLoading Edges: {:.2f}%            ".format(min(100, 100 * (int(line[1]) + 1) / num_edges)), end="")
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
            # Running print statement tracking the progress of loading the surfaces
            print("\rLoading Surfaces: {:.2f}%         ".format(min(100, 100 * (int(line[1]) + 1) / num_surfs)), end="")
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
    connect_read_net_old(net)


# Connect input network function. Connect the network objects using their load index lists
def connect_read_net_old(net):
    # Connect the atom objects
    print("\rConnecting Atoms            ", end="")
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
    print("\rConnecting Vertices               ", end="")
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
    print("\rConnecting Edges                  ", end="")
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
    print("\rConnecting Surfaces               ", end="")
    for surf in net.surfs:
        # Reset the surfaces and edges
        surf.verts, surf.edges = [], []
        # Connect the surface's vertices
        for i in range(len(surf.load_ndxs[0])):
            surf.verts.append(net.verts[surf.load_ndxs[0][i]])
        # Connect the surface's edges
        for i in range(len(surf.load_ndxs[1])):
            surf.edges.append(net.edges[surf.load_ndxs[1][i]])
    print("\rNetwork Loaded                    ", end="")


def read_net(sys, file=None):
    # Open the file
    if file is None:
        file = sys.net_file
    if sys.dir is None:
        sys.set_output_directory()
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
        sys.net.sol_verts = bool(read_file[1][4])
        # Add the vertices
        for i in range(3, 3 + net_verts):
            vert = sys.net.verts[i - 3]
            vert.loc = [float(_) for _ in read_file[i][1:4]]
            vert.rad = float(read_file[i][4])
            vert.atoms = [sys.atoms[int(_)] for _ in read_file[i][5:9]]
            vert.ndx = [int(_) for _ in read_file[i][5:9]]
            vert.edges = [sys.net.edges[int(_)] for _ in read_file[i][9:14] if _ != '']
            surf_ndxs = [int(_) for _ in read_file[i][14:] if _ != '']
            vert.surfs = [sys.net.surfs[_] for _ in surf_ndxs]
            for atom in vert.atoms:
                atom.verts.append(vert)
            for surf in vert.surfs:
                if surf.verts is None:
                    surf.verts = []
                surf.verts.append(vert)
        # Add the edges
        for i in range(4 + net_verts, 4 + net_verts + net_edges):
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
            surf = sys.net.surfs[i - 5 - net_verts - net_edges]
            surf.atoms = [sys.atoms[int(_)] for _ in read_file[i][4:6]]
            surf.ndx = [int(_) for _ in read_file[i][4:6]]
            surf.file = read_file[i][1]
            surf.sa = float(read_file[i][2])
            surf.curv = float(read_file[i][3])
            surf.func = [float(_) for _ in read_file[i][5:16]] + [[float(_) for _ in read_file[i][16:]]]
            for atom in surf.atoms:
                atom.surfs.append(surf)
    # Go through and add the surfaces if they have files
    for surf in sys.net.surfs:
        surf.read_file()
