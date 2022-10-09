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


# Add vertices function. Takes in a system and a file with vertices in it and adds the verts to the system
def add_verts(sys, file_address):
    # Reset the network and open the network
    sys.myNet.verts, sys.myNet.surfs, sys.myNet.edges = [], [], []
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
                             net=sys.myNet, doublet=True, loc2=[float(line2[0]), float(line2[1]),
                                                                float(line2[2])], rad2=float(line2[3]))
            # Skip the next line
            i += 1
        else:
            # Regular vertex
            my_vert = Vertex(atoms, location=[float(line[0]), float(line[1]), float(line[2])], radius=float(line[3]),
                             net=sys.myNet)
        # Add the vertex to the system
        sys.myNet.verts.append(my_vert)


# Add vertices function. Takes in a system and a file with vertices in it and adds the verts to the system
def add_grant_verts(sys, file_address):
    # Reset the network and open the network
    sys.myNet.verts, sys.myNet.surfs, sys.myNet.edges = [], [], []
    vert_file = open(file_address).readlines()
    # Go through each of the vertices file
    for line in vert_file[1:]:
        line = line.split(',')
        atoms = [sys.atoms[int(line[i]) - 1] for i in range(4)]
        loc = [float(line[-3][1:]), float(line[-2]), float(line[-1][:-2])]
        rad = calc_dist(atoms[0].loc, loc) - atoms[0].rad
        vert = Vertex(atoms, location=loc, radius=rad, net=sys.myNet)
        sys.myNet.verts.append(vert)


# Export my surfaces function. Used to create and export the surfaces of a system as separate files
def export_mySurfs(sys, n, max_num):
    # Surfaces Folder
    os.mkdir(sys.output_directory + "/Surfaces")
    os.chdir(sys.output_directory + "/Surfaces")
    # Go through each surface and create a file for each adding the vertex points
    surf_ndxs = []
    for i in range(len(sys.myNet.surfs)):
        percentage = int((n + (i + 1) / 2) / max_num * 100)
        print("\rExporting System: ",
              '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
        # Find the relative surface index and add it to the list
        surf_ndxs.append(str(sys.atoms.index(sys.myNet.surfs[i].atoms[0]) + 1) + "_" +
                         str(sys.atoms.index(sys.myNet.surfs[i].atoms[1]) + 1))
        # Give the surfaces random colors
        color = [np.random.random_sample() for _ in range(3)]
        # Write the surface
        write_surfs([sys.myNet.surfs[i]], "surf_" + surf_ndxs[i], color)
    export_info(sys, sys.name + "_surface_info", "surfaces", interfaces=[[surf] for surf in sys.myNet.surfs])
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
        for surf in sys.myNet.surfs:
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
            for surf in sys.myNet.surfs:
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


# Export vertices function. Creates a vertex file for future reference
def export_myVerts(sys):
    # Change to the output directory
    os.chdir(sys.output_directory)
    # Create the vertices file
    file = open(os.getcwd() + "/" + sys.name + "Vertices.txt", 'w')
    # Go through the vertices in the system
    for i in range(len(sys.myNet.verts)):
        # Get the vertex and the vertexes' atoms' indices
        vert = sys.myNet.verts[i]
        # Write the vertex into the file (x y z r a0 a1 a2 a3 a4)
        file.write(str(vert.loc[0]) + " " + str(vert.loc[1]) + " " + str(vert.loc[2]) + " " + str(vert.rad) + " "
                   + str(vert.ndx[0]) + " " + str(vert.ndx[1]) + " " + str(vert.ndx[2]) + " " + str(vert.ndx[3]) + '\n')
        # If the vertex is a doublet, write down the radius and location
        if vert.doublet:
            file.write(str(vert.loc2[0]) + " " + str(vert.loc2[1]) + " " + str(vert.loc2[2]) + " " + str(vert.rad2) +
                       " " + str(vert.ndx[0]) + " " + str(vert.ndx[1]) + " " + str(vert.ndx[2]) + " " + str(vert.ndx[3])
                       + '\n')



    # Export method. Takes in an export type: 'Atoms', 'surfs'
    def export(self, export_all=True, export_sys=False, export_atoms=False, export_mols=False,
               export_surfs=False, export_sys_pdb=False, export_reses=False):
        if self.output_directory is None:
            set_output_dir(self)
        os.chdir(self.output_directory)
        # Go through all the possible user inputs and choose the correct export function
        n = 0
        lengths = [1, len(self.myNet.surfs), len(self.atoms), len(self.mols)]
        exports = [export_sys, export_surfs, export_atoms, export_mols]
        # Find the total number of things being exported for the loading bar
        max_num_arr = [lengths[i] for i in range(len(lengths)) if exports[i] or export_all]
        max_num = sum(max_num_arr)
        # Export the system
        if export_sys or export_all:
            export_mySys(self, n, max_num)
            n += 1
        # Export the pdb of the system
        if export_sys_pdb or export_all:
            create_pdb(self, os.getcwd())
        # Export the individual surfaces
        if export_surfs or export_all:
            export_mySurfs(self, n, max_num)
            n += len(self.myNet.surfs)
        # Export the atom cells
        if export_atoms or export_all:
            export_myAtoms(self, n, max_num)
            n += len(self.atoms)
        # Export the molecules
        if export_mols or export_all:
            export_myMols(self, n, max_num, export_residues=export_reses or export_all)
            n += len(self.mols)
        print("\rExporting System:  ########## 100 %")
        print("\rSystem Exported")


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
                                ", Surface Area = " + str(atoms[i].surfs[j].interface_sa) + "\n")
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
                iface_sa += iface[j].interface_sa
            surface_areas.append(iface_sa)
            total_sa += iface_sa
        # Write out the total volume for the set of atoms
        info_file.write("\nTotal Surface Area for " + set_name + " = " + str(total_sa) + "\n\n")
        # Go through each interface writing the information for it
        for i in range(len(interfaces)):
            info_file.write("\nInterface " + interface_names[i] + ":\n")
            info_file.write("\nTotal Surface Area = " + str(surface_areas[i]))
            info_file.write("\nCurvature = \n\n")