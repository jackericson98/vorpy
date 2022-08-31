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


# Get pdb data method. Finds the lines of the file with prefixes and returns them as a list
def get_pdb_data(sys, word):
    # Get the file information
    sys.file = open(sys.file_address).readlines()
    sys.name = get_name(sys.file_address)
    # Split each line in the file
    for i in range(len(sys.file)):
        sys.file[i] = sys.file[i].split()
    # Special case for Atom lines
    if word.lower() == 'atom':
        atoms = []
    # Go through each line in the file and check if the first word is the word we are looking for
        for i in range(len(sys.file)):
            line = sys.file[i]
            if line and line[0].lower() == 'atom':  # Check if the line starts with atom
                atom = Atom([float(line[5]), float(line[6]), float(line[7])], get_radius(line[-1]),
                            symbol=line[-1], res=line[2], chain=line[3], res_seq=line[4])
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


# Create pdb method. Creates a pdb file type in the current working directory
def create_pdb(sys, directory=None):
    # Create the output file
    file = open(sys.name + "_structure.pdb", 'w')
    # If the file exists, copy it over
    if sys.file is not None:
        for line in sys.file:
            file.write(line)
        return
    # Make sure each atom has a type before we change directories
    for atom in sys.atoms:
        # Give each atom a type if not indicated
        if atom.type == "":
            atom.type = sys.get_radius(atom.rad, return_symbol=True)
    # Move to the indicated directory
    if directory:
        os.chdir(directory)

    # Go through each atom in the system
    for i in range(len(sys.atoms)):
        a = sys.atoms[i]
        loc = [str(round(a.loc[0], 3)), str(round(a.loc[1], 3)), str(round(a.loc[2], 3))]
        # Write the lines for the atom
        file.write("ATOM" + " " * (7 - len(str(i+1))) +
                   str(i + 1) + "  " +
                   a.type + " " * (4 - len(a.type)) +
                   a.res + " " * (4 - len(a.res)) +
                   a.chain + " " * (5 - len(a.chain) - len(a.res_seq)) +
                   " " * 4 + " " * (8 - len(loc[0])) +
                   loc[0] + " " * (8 - len(loc[1])) +
                   loc[1] + " " * (8 - len(loc[2])) +
                   loc[2] + " " * 2 +
                   "1.00  0.00" + " " * (12 - len(a.type)) + a.type + "\n")


# Get gro method. Finds data in a gro file
def get_gro(sys):
    sys.info['header'] = sys.file[0]
    # Go through each line in the file and create an atom object
    for line in sys.file[2:-2]:
        sys.atoms.append(Atom([line[3], line[4], line[5]], sys.get_radius(line[1][0]), symbol=line[1][0]))


# Get mol method. Finds data in a mol file
def get_mol(sys):
    for line in sys.file:
        if len(line) > 6:
            sys.atoms.append(Atom([line[0], line[1], line[2]], sys.get_radius(line[3]), symbol=line[3]))


# Add Voronota data method. Takes in voronota data and adds it to the System
def add_vta_data(sys, ball_file, vert_file):
    # Set the voronota system indicator to True
    sys.net.vta = True
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


def export_sys(sys):
    os.chdir(sys.dir)
    # If the file is none create a pdb for the file
    if sys.file is None:
        create_pdb(sys)
    # Set the name of the file to be created if no name exists
    if sys.name is None:
        sys.name = "mySystem"
    # Set the counters to 0
    tot_verts, tot_tris = 0, 0
    # Get the total number of vertices and tris
    for surf in sys.net.surfs:
        tot_verts += len(surf.points)
        tot_tris += len(surf.tris)
    # System file
    sys_file = open(sys.name + "_System.off", 'w')
    sys_file.write("OFF\n" + str(tot_verts) + " " + str(tot_tris) + " 0\n\n\n")
    # Go through the surfaces and add the points
    for i in range(len(sys.net.surfs)):
        # Go through the points on the surface
        for point in sys.net.surfs[i].points:
            # Add the point to the system file and the surface's file (rounded to 4 decimal points)
            str_point = [str(round(point[_], 4)) for _ in range(3)]
            sys_file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
    num_verts = 0
    # Go through each surface and add the faces
    for i in range(len(sys.net.surfs)):
        for tri in sys.net.surfs[i].tris:
            # Add the triangle to the system file and the surface's file
            str_tri = [str(tri[_] + num_verts) for _ in range(3)]
            sys_file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " 1 0 0\n")
        # Keep counting triangles for the system file
        num_verts += len(sys.net.surfs[i].points)
    os.chdir('..')


def export_surfs(sys):
    # Surfaces Folder
    os.mkdir(sys.dir + "/Surfaces")
    os.chdir(sys.dir + "/Surfaces")
    # Go through each surface and create a file for each adding the vertex points
    surf_ndxs = []
    for i in range(len(sys.net.surfs)):
        # Find the relative surface index and add it to the list
        surf_ndxs.append(str(sys.atoms.index(sys.net.surfs[i].atoms[0]) + 1) + "_" +
                         str(sys.atoms.index(sys.net.surfs[i].atoms[1]) + 1))
        # Name the file with the surface's index
        file = open(str("surf_" + surf_ndxs[i] + ".off"), 'w')
        # Start the file with "OFF", the number of points for the surface and the number of triangles
        file.write("OFF\n" +
                   str(len(sys.net.surfs[i].points)) + " " + str(len(sys.net.surfs[i].tris)) + " 0\n\n\n")
        # Go through the points on the surface
        for point in sys.net.surfs[i].points:
            # Add the point to the system file and the surface's file (rounded to 4 decimal points)
            str_point = [str(round(point[_], 4)) for _ in range(3)]
            file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
    num_verts = 0
    # Go through each surface opening the previously created file and add the faces
    for i in range(len(sys.net.surfs)):
        file = open(str("surf_" + surf_ndxs[i] + ".off"), 'a')
        for tri in sys.net.surfs[i].tris:
            # Add the triangle to the system file and the surface's file
            file.write("3 " + str(tri[0]) + " " + str(tri[1]) + " " + str(tri[2]) + " 1 0 0\n")
        # Keep counting triangles for the system file
        num_verts += len(sys.net.surfs[i].points)
    os.chdir("..")


def export_atoms(sys):
    # Atoms Folder
    os.mkdir(sys.dir + "/Atoms")
    os.chdir(sys.dir + "/Atoms")
    # Add the vertices and triangles for each surface of each atom
    for i in range(len(sys.atoms)):
        # Create a file for each atom
        atom_file = open(str("atom_" + str(i + 1) + "_cell.off"), 'w')
        # Calculate the number of points and triangles
        tot_verts, tot_tris = 0, 0
        for surf in sys.atoms[i].surfs:
            tot_verts += len(surf.points)
            tot_tris += len(surf.tris)
        # Create the off header
        atom_file.write("OFF\n" + str(tot_verts) + " " + str(tot_tris) + " 0\n\n\n")
        # Go through each surface of the atom and add the vertices
        for j in range(len(sys.atoms[i].surfs)):
            for point in sys.atoms[i].surfs[j].points:
                str_point = [str(round(point[_], 4)) for _ in range(3)]
                atom_file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts = 0
        # Go through each surface opening the previously created file and add the faces
        for j in range(len(sys.atoms[i].surfs)):
            atom_file = open(str("atom_" + str(i + 1) + "_cell.off"), 'a')
            for tri in sys.atoms[i].surfs[j].tris:
                atom_file.write("3 " + str(tri[0] + num_verts) + " " + str(tri[1] + num_verts) + " " +
                                str(tri[2] + num_verts) + " 1 0 0\n")
            num_verts += len(sys.atoms[i].surfs[j].points)
    os.chdir('..')


def export_mols(sys):
    # Create the molecules folder
    os.mkdir(sys.dir + '/Molecules')
    os.chdir(sys.dir + '/Molecules')

    os.chdir('..')


def export_analysis(sys):
    os.chdir(sys.dir)
    # Create the Atoms Folder
    atom_info = open("atom_info.txt", 'w')
    # Add the vertices and triangles for each surface of each atom
    for i in range(len(sys.atoms)):
        atom_info.write("Atom " + str(i + 1) + "\n")
        atom_info.write(" Cell Volume = {}\n".format(sys.atoms[i].vol))
        atom_info.write(" Surfaces:\n")
        for j in range(len(sys.atoms[i].surfs)):
            if sys.atoms[i] == sys.atoms[i].surfs[j].atoms[0]:
                a1 = sys.atoms[i].surfs[j].atoms[1]
            else:
                a1 = sys.atoms[i].surfs[j].atoms[0]
            atom_info.write("  Surface " + str(j + 1) + ", Made with Atom " + str(sys.atoms.index(a1) + 1) +
                            ", Surface Area = " + str(sys.atoms[i].surfs[j].sa) + "\n")
        atom_info.write("\n")
    os.chdir("..")
