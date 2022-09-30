import os

from System.calcs import *


# Set output directory function. prevents the system from making duplicate output directories
def set_output_dir(sys, dir_name=None):
    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        dir_name = sys.vorpy_directory + "/Data/User_data/" + sys.name
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = str(i)
            # If no file with the system name exists change the string to empty
            if i == 0:
                i_str = ""
            # Try to create the directory
            os.mkdir(dir_name + i_str)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    sys.output_directory = dir_name + i_str


# Export vertices function. Creates a vertex file for future reference
def export_myVerts(sys):
    # Change to the output directory
    os.chdir(sys.output_directory)
    # Create the vertices file
    file = open(os.getcwd() + "/" + sys.name + "Vertices.txt", 'w')
    # Go through the vertices in the system
    for i in range(len(sys.net.verts)):
        # Get the vertex and the vertexes' atoms' indices
        vert = sys.net.verts[i]
        # Write the vertex into the file (x y z r a0 a1 a2 a3 a4)
        file.write(str(vert.loc[0]) + " " + str(vert.loc[1]) + " " + str(vert.loc[2]) + " " + str(vert.rad) + " "
                   + str(vert.ndx[0]) + " " + str(vert.ndx[1]) + " " + str(vert.ndx[2]) + " " + str(vert.ndx[3]) + '\n')


# Create pdb method. Creates a pdb file type in the current working directory
def create_pdb(sys, directory=None):
    # Create the output file
    file = open(sys.name + "_structure.pdb", 'w')
    # If the file exists, copy it over
    if sys.file is not None:
        for line in open(sys.file_address):
            file.write(str(line))
        return
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
                   a.element + " " * (4 - len(a.element)) +
                   a.res + " " * (4 - len(a.res)) +
                   a.chain + " " * (5 - len(a.chain) - len(a.res_seq)) +
                   " " * 4 + " " * (8 - len(loc[0])) +
                   loc[0] + " " * (8 - len(loc[1])) +
                   loc[1] + " " * (8 - len(loc[2])) +
                   loc[2] + " " * 2 +
                   "1.00  0.00" + " " * (12 - len(a.element)) + a.element + "\n")


# Write surfaces function. Writes files given a list of surfaces
def write_surfs(surfs, file_name, color=None):
    # If no color is given, make the color white
    if color is None:
        color = [1, 1, 1]
    # Create the file
    file = open(file_name + ".off", 'w')
    # Count the number of triangles and vertices there are
    num_verts, num_tris = 0, 0
    for i in range(len(surfs)):
        num_verts += len(surfs[i].points)
        num_tris += len(surfs[i].tris)
    # Write the numbers into the file
    file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
    # Go through the surfaces and add the points
    for i in range(len(surfs)):
        # Go through the points on the surface
        for point in surfs[i].points:
            # Add the point to the system file and the surface's file (rounded to 4 decimal points)
            str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
            file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
    num_verts, tri_count = 0, 0
    # Go through each surface and add the faces
    for i in range(len(surfs)):
        for tri in surfs[i].tris:
            # Add the triangle to the system file and the surface's file
            str_tri = [str(tri[_] + num_verts) for _ in range(3)]
            file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                       str(color[1]) + " " + str(color[2]) + "\n")
        # Keep counting triangles for the system file
        num_verts += len(surfs[i].points)


# Export my system function. Used to create and export the surfaces of a system as one file
def export_mySys(sys, n, max_num):
    # Get the percentage and update the print statement
    percentage = int((n + 1) / max_num * 100)
    print("\rExporting System: ",
          '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
    # If the file is none create a pdb for the file
    create_pdb(sys, sys.file_address)
    # Set the name of the file to be created if no name exists
    if sys.name is None:
        sys.name = "mySystem"
    # Write the surfaces
    write_surfs(sys.net.surfs, sys.name + "_system")


# Export my surfaces function. Used to create and export the surfaces of a system as separate files
def export_mySurfs(sys, n, max_num):
    # Surfaces Folder
    os.mkdir(sys.output_directory + "/Surfaces")
    os.chdir(sys.output_directory + "/Surfaces")
    # Go through each surface and create a file for each adding the vertex points
    surf_ndxs = []
    for i in range(len(sys.net.surfs)):
        percentage = int((n + (i + 1) / 2) / max_num * 100)
        print("\rExporting System: ",
              '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
        # Find the relative surface index and add it to the list
        surf_ndxs.append(str(sys.atoms.index(sys.net.surfs[i].atoms[0]) + 1) + "_" +
                         str(sys.atoms.index(sys.net.surfs[i].atoms[1]) + 1))
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        # Write the surface
        write_surfs([sys.net.surfs[i]], "surf_" + surf_ndxs[i], color)
    export_info(sys, sys.name + "_surface_info", "surfaces", interfaces=[[surf] for surf in sys.net.surfs])
    os.chdir(sys.output_directory)


# Export my atoms function. Used to create and export the surfaces surrounding each atom of a system as separate files
def export_myAtoms(sys, n, max_num):
    # Atoms Folder
    os.mkdir(sys.output_directory + "/Atoms")
    os.chdir(sys.output_directory + "/Atoms")
    # Add the vertices and triangles for each surface of each atom
    for i in range(len(sys.atoms)):
        percentage = int((n + (i + 1) / 2) / max_num * 100)
        print("\rExporting System: ",
              '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        # Write the surfaces
        write_surfs(sys.atoms[i].surfs, "atom_" + str(i + 1) + "_cell", color)
    # Export the information file
    export_info(sys, sys.name + "_Atom_info", "", sys.atoms)
    os.chdir(sys.output_directory)


# Export my mols function. Used to create and export the surfaces the interfaces between molecules of the system  and
# the cells of the atoms of each molecule as separate files
def export_myMols(sys, n, max_num, export_residues=False):
    # Create the molecules folder
    os.mkdir(sys.output_directory + '/Molecules')
    chains = []
    chain_lists = []
    # Create the chains
    for atom in sys.atoms:
        # If the chain hasn't been found create it and add the atom to it
        if atom.chain not in chains:
            os.mkdir(sys.output_directory + '/Molecules/' + atom.chain)
            chains.append(atom.chain)
            chain_lists.append([atom])
        # If the chain has been found add the atom to the chain's list of atoms
        else:
            chain_lists[chains.index(atom.chain)].append(atom)
    # Go through each of the chains
    for i in range(len(chains)):
        # Percentage print statement
        percentage = int((n + (i + 1)) / max_num * 100)
        print("\rExporting System: ",
              '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
        chain_file_names = []
        # Go through the other chains and create a file for their interfaces
        for j in range(len(chains)):
            if chains[j] == chains[i]:
                continue
            chain_file_names.append(chains[i] + '_' + chains[j] + '_interface')
        # Get the surfaces for the interfaces of each molecule
        chain_iface_surfs = [[] for _ in range(len(chain_file_names))]
        chain_surfs = []
        # Find the file that each surf belongs to
        for surf in sys.net.surfs:
            if surf.atoms[0].chain == chains[i] != surf.atoms[1].chain:
                chain_iface_surfs[chain_file_names.index(chains[i] + '_' + surf.atoms[1].chain + '_interface')].append(surf)
                chain_surfs.append(surf)
            elif surf.atoms[1].chain == chains[i] != surf.atoms[0].chain:
                chain_iface_surfs[chain_file_names.index(chains[i] + '_' + surf.atoms[0].chain + '_interface')].append(surf)
                chain_surfs.append(surf)
            else:
                continue
        # Make a file for the whole chain
        os.chdir(sys.output_directory + "/Molecules")
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        write_surfs(chain_surfs, chains[i], color)
        # Make a file for each of the interfaces
        os.chdir(sys.output_directory + "/Molecules/" + chains[i])
        for j in range(len(chain_iface_surfs)):
            # Give the surfaces random colors
            color = [np.random.random_sample() for _ in range(3)]
            write_surfs(chain_iface_surfs[j], chain_file_names[j], color)
        # Export the residues for the molecule
        if export_residues:
            export_myResidues(sys, chain_lists[i], chains[i])
        # Change back to the chain directory
        os.chdir(sys.output_directory + "/Molecules/" + chains[i])
        # Get the information for the molecules in the system
        export_info(sys, sys.name + "_mol_" + chains[i] + "_info", chains[i], chain_lists[i], chain_iface_surfs, chain_file_names)
        # Change back to the main directory
        os.chdir(sys.output_directory)


# Export residues function.
def export_myResidues(sys, mol_atoms, chain):
    # Create a 'residues' folder for each chain
    os.mkdir(sys.output_directory + "/Molecules/" + chain + "/Residues")
    os.chdir(sys.output_directory + "/Molecules/" + chain + "/Residues")
    # Create the lists of residues and their respective atoms
    residues = []
    res_lists = []
    # Create the residues
    for atom in mol_atoms:
        # Check to see that the atom has a residue
        if atom.res == "":
            atom.res = "None"
        # If the chain hasn't been found create it and add the atom to it
        if atom.res not in residues:
            os.mkdir(sys.output_directory + '/Molecules/' + chain + "/" + atom.res)
            residues.append(atom.res)
            res_lists.append([atom])
        # If the chain has been found add the atom to the chain's list of atoms
        else:
            res_lists[residues.index(atom.res)].append(atom)
    # Write the residue files
    for i in range(len(residues)):
        res_surfs = []
        for atom in res_lists[i]:
            res_surfs += atom.surfs
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        write_surfs(res_surfs, residues[i], color)


# Export analysis function. Creates an analysis file for the user.
def export_info(sys, file_name, set_name, atoms=None, interfaces=None, interface_names=None):
    # Create the Atoms Folder
    info_file = open(file_name + ".txt", 'w')
    info_file.write(set_name + "\n")
    # Go through each atom in the list providing contributing information to the set of atoms
    if atoms is not None:
        # Get the total volume for the atoms in the given system
        tot_vol = 0
        info_file.write("\nAtoms\n")
        # Go through each atom adding to the running total volume
        for i in range(len(atoms)):
            tot_vol += atoms[i].cell_vol
        # Write out the total volume for the set of atoms
        info_file.write("\nTotal Volume for " + set_name + " = " + str(tot_vol) + "\n\n")
        # Go through each atom recording volume and surface areas for the respective surfaces
        for i in range(len(atoms)):
            # Write the header for each atom
            info_file.write("Atom " + str(i) + ": Chain - " + atoms[i].chain + "\n")
            info_file.write(" Cell Volume = {}\n".format(atoms[i].cell_vol))
            # Write the surface information
            info_file.write(" Surfaces:\n")
            for j in range(len(atoms[i].surfs)):
                # Grab the other atom in the surface
                if atoms[i] == atoms[i].surfs[j].atoms[0]:
                    a1 = atoms[i].surfs[j].atoms[1]
                else:
                    a1 = atoms[i].surfs[j].atoms[0]
                # Write the information for the surface
                info_file.write("  Surface " + str(j + 1) + ", Made with Atom " + str(sys.atoms.index(a1) + 1) +
                                ", Surface Area = " + str(atoms[i].surfs[j].sa) + "\n")
            info_file.write("\n\n")

    # Go through each interface given
    if interfaces is not None:
        # Set the interface names to nothing if they have not been given
        if interface_names is None:
            interface_names = ["" for i in range(len(interfaces))]
        # Get the total surface area for the atoms in the given system
        total_sa = 0
        info_file.write("\nInterfaces\n")
        surface_areas = []
        # Go through each interface listed
        for i in range(len(interfaces)):
            iface = interfaces[i]
            iface_sa = 0
            # Go through each surface in the interface
            for j in range(len(iface)):
                iface_sa += iface[j].sa
            surface_areas.append(iface_sa)
            total_sa += iface_sa
        # Write out the total volume for the set of atoms
        info_file.write("\nTotal Surface Area for " + set_name + " = " + str(total_sa) + "\n\n")
        # Go through each interface writing the information for it
        for i in range(len(interfaces)):
            info_file.write("\nInterface " + interface_names[i] + ":\n")
            info_file.write("\nTotal Surface Area = " + str(surface_areas[i]))
            info_file.write("\nCurvature = \n\n")
