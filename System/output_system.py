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
        # If the vertex is a doublet, write down the radius and location
        if vert.doublet:
            file.write(str(vert.loc2[0]) + " " + str(vert.loc2[1]) + " " + str(vert.loc2[2]) + " " + str(vert.rad2) +
                       " " + str(vert.ndx[0]) + " " + str(vert.ndx[1]) + " " + str(vert.ndx[2]) + " " + str(vert.ndx[3])
                       + '\n')


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
def export_myMols(sys, export_residues=False):
    # Create the molecules folder
    os.mkdir(sys.output_directory + '/Molecules')
    # Go through the chains in the system
    for mol in sys.mols:
        # Set up the list of molecule surfaces
        mol_surfs = []
        # Get the name of the molecule
        mol_name = mol[0].chain
        # Make the output directory for the molecule
        os.mkdir(sys.output_directory + '/Molecules/' + mol_name)
        # Change to the directory of the molecule
        os.chdir(sys.output_directory + '/Molecules/' + mol_name)
        # Get the surrounding surfaces for the molecule
        for surf in sys.net.surfs:
            # Grab the surface's atoms
            a0, a1 = surf.atoms
            # Check to see if the surface is a part of the current molecule and another molecule or not
            if (a0.chain == mol_name and a1.chain != mol_name) or (a0.chain != mol_name and a1.chain == mol_name):
                mol_surfs.append(surf)
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        # Create the molecule file
        write_surfs(mol_surfs, "Molecule_" + mol_name, color)

        # Create the interfaces folder and change to it
        os.mkdir(sys.output_directory + '/Molecules/' + mol_name + '/Interfaces')
        # Change to the directory of the molecule
        os.chdir(sys.output_directory + '/Molecules/' + mol_name + '/Interfaces')
        # Create a liat to hold the interfaces for the molecule
        interfaces = [[] for _ in range(len(sys.mols))]
        # Create the interfaces of the molecule with the other molecules
        for i in range(len(sys.mols)):
            # Get the second molecule
            mol2 = sys.mols[i]
            # Make sure the 2 molecules aren't the same
            if mol == mol2:
                continue
            # Get the second molecule's name
            mol2_name = mol2[0].chain
            # Find the surfaces in the 2 molecules' interface
            for surf in sys.net.surfs:
                # Grab the surface's atoms
                a0, a1 = surf.atoms
                # Check to see if the surface is a part of the current molecule and another molecule or not
                if (a0.chain == mol2_name and a1.chain == mol_name) or (a0.chain == mol_name and a1.chain == mol2_name):
                    # Add the surface to the list of surfaces in the interface
                    interfaces[i].append(surf)
            # Only create the file if one or more surface exists in the interface between the molecules
            if len(interfaces[i]) >= 1:
                write_surfs(interfaces[i], mol_name + "_" + mol2_name + "_Interface")
        # Change back out of the directory
        os.chdir(sys.output_directory + '/Molecules/' + mol_name)
        # If the residues have been requested to be exported, do so
        if export_residues:
            export_myResidues(sys, mol, mol_name)
        # Export the information for the molecule
        export_info(sys, mol_name + "_Information", mol_name + "_Information", interfaces=interfaces)


# Export residues function.
def export_myResidues(sys, mol, mol_name):
    # Create a 'residues' folder for each chain
    os.mkdir(sys.output_directory + "/Molecules/" + mol_name + "/Residues")
    os.chdir(sys.output_directory + "/Molecules/" + mol_name + "/Residues")
    # Go through each residue adding the surrounding surfaces
    for res in sys.residues:
        # Check the residue to see if it belongs to the molecule
        if res[0].chain != mol[0].chain:
            continue
        # Get the residue sequence to check against
        res_seq = res[0].res_seq
        surfs = []
        for atom in res:
            surfs += atom.surfs
        # Set up a list of surfaces for the residue's outer surfaces
        res_surfs = []
        # Go through the surfaces in the
        for surf in surfs:
            # Get the surface's atoms
            a0, a1 = surf.atoms
            # Check to see if the surface is a part of the current molecule and another molecule or not
            if (a0.res_seq == res_seq and a1.res_seq != res_seq) or (a0.res_seq != res_seq and a1.res_seq == res_seq):
                res_surfs.append(surf)
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        # Create the molecule file
        write_surfs(surfs=res_surfs, file_name=res[0].res, color=color)


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
