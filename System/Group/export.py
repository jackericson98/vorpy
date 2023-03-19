from System.sys_funcs.output import *
from System.Group.sort import get_iface


def export_iface(grp, bff=None, info_file=False, interface_atoms=False):
    """
    Exports the information from the given interface as a txt file
    :param bff: Group for the interface with self
    :param info_file: Whether to export a txt file with info on the interface or not
    :param interface_atoms: Whether to export a pdb file with the atoms around the interface or not
    :return:
    """
    # Check to see if there is a bff or not
    if bff is not None:
        grp.bff = bff
    # Check to see that the interface has been calculated
    if grp.iface_surfs is None or grp.iface_atoms is None:
        get_iface(grp)
    # Set the interface name
    interface_name = grp.name + "_" + grp.bff.name + "_interface"
    # Move to the output directory
    os.chdir(grp.sys.dir)
    # Create and move to the interface directory
    os.mkdir(os.getcwd() + "/" + interface_name)
    os.chdir(os.getcwd() + "/" + interface_name)
    # Write the surfaces for the interface
    write_surfs(grp.iface_surfs, interface_name)
    # Check to see of the user wants to export the interface's atoms
    if interface_atoms:
        # Get the two sets of interface atoms
        write_pdb(grp.iface_atoms, interface_name + "_" + grp.name + "_atoms", grp.sys)
        write_pdb(grp.bff.iface_atoms, interface_name + "_" + grp.bff.name + "_atoms", grp.bff.sys)
    # Check to see if the user wants to export the interface's information
    if info_file:
        info = open(interface_name + "_info.txt", 'w')
        info.write("Interface between " + grp.name + " and " + grp.bff.name + " : \n")
        info.write("Number of Surfaces: " + str(len(grp.iface_surfs)))
        info.write("Surface Area: " + str(grp.iface_sa))
        info.close()


def group_exports(grp, all_=False, atoms=False, shell=False, fill=False, surfaces=False, layers=False, num_layers=50,
            info=False, iface=False, verts=False, surr_atoms=False, ext_atoms=False, shell_edges=False,
            shell_verts=False, edges=False):
    """
    Exports specified export types for the group
    :param grp:
    :param all_: All possible exports for the group will be exported to the group directory
    :param atoms: Exports a new pdb file contasining only the atoms of the group
    :param shell: Exports the outer surfaces of the group
    :param fill: Exports all surfaces in the group as one object
    :param surfaces: Exports all surfaces in the group as seperate files, named by their atoms
    :param layers: Exports all layers surrounding the group, unless num_layers is specified
    :param num_layers: Controls the number of exported layers for the group
    :param info: Exports the information for the group
    :param iface: Exports the interface for the group, bff must be specified first
    :param verts: Exports the vertices of the group as an off file
    :param surr_atoms: Exports the atoms directly surrounding the group (residues intact)
    :param ext_atoms: Exports the outermost atoms in the group's set of atoms (must be a part of shell)
    :param shell_edges: Exports only the outermost edges for the group as an OFF file
    :param shell_verts: Exports the outermost vertices for the group
    :param edges: Exports all edges for the group
    :return: The specified export is placed in the group's directory
    """
    # Set the surface colors and scheme
    if grp.surf_color is None:
        grp.surf_color = grp.sys.net.surf_col
    if grp.surf_scheme is None:
        grp.surf_scheme = grp.sys.net.surf_scm
    # Get the surfaces if they haven't been got
    if grp.surfs is None or len(grp.surfs) == 0:
        grp.build_surfs()
    # Get the surface coloring scheme
    scheme = ""
    if grp.surf_scheme is not None:
        scheme = "_" + grp.surf_scheme
    # Create the output directory inside the system's directory
    if grp.dir is None:
        i = 1
        my_dir = grp.sys.dir + "/" + grp.name
        first = True
        while os.path.exists(my_dir):
            if first:
                my_dir += "__"
                first = False
            my_dir = my_dir[:-(1 + len(str(i)))] + '_' + str(i)
            i += 1
        grp.dir = my_dir
        os.mkdir(grp.dir)
    os.chdir(grp.dir)
    # If the user wants to export the atoms for the group
    if atoms or all_:
        write_pdb(atoms=grp.atoms, name=grp.name + "_atoms", sys=grp.sys)
    # If the user wants to export the shell for the group
    if shell or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # noinspection PyUnresolvedReferences
        if grp.layer_surfs is not None and len(grp.layer_surfs) > 0:
            write_surfs(surfs=grp.layer_surfs[0], file_name=grp.name + "_shell" + scheme, directory=grp.dir, color_map=grp.surf_color, color_scheme=grp.surf_scheme)
    # If the user wants a filled shell for the group
    if fill or all_:
        grp.build_surfs()
        write_surfs(surfs=grp.surfs, file_name=grp.name + "_fill" + scheme, directory=grp.dir, color_map=grp.surf_color, color_scheme=grp.surf_scheme)
    # If the user wants separate surfaces for the group
    if surfaces or all_:
        i = 1
        my_dir = grp.dir + "/surfaces"
        while os.path.exists(my_dir):
            if my_dir[-1] == 's':
                my_dir += '__'
            my_dir  = my_dir[:-2] + str(i)
            i += 1
        os.mkdir(my_dir)
        os.chdir(my_dir)
        for surf in grp.surfs:
            write_surfs([surf], file_name="_".join([str(_) for _ in surf.ndx]), directory=my_dir, color_map=grp.surf_color, color_scheme=grp.surf_scheme)
        os.chdir(grp.dir)
    # If the user wants layers
    if layers or all_:
        # First check to see if the number of layers is greater than 1
        if grp.layer_atoms is None or len(grp.layer_atoms) <= 1:
            grp.get_layers(max_layers=num_layers)
        # Create the layers directory
        i = 1
        my_dir = os.getcwd() + "/layers"
        while os.path.exists(my_dir):
            if my_dir[-1] == 's':
                my_dir += '__'
            my_dir = my_dir[:-2] + str(i)
            i += 1
        os.mkdir(my_dir)
        os.chdir(my_dir)
        # Create the layer and atoms files
        for i in range(len(grp.layer_surfs)):
            write_pdb(grp.layer_atoms[i + 1], name=str(i) + "_atoms", sys=grp.sys)
            write_surfs(grp.layer_surfs[i], file_name=str(i) + "_surfs", color_map=grp.surf_color, color_scheme=grp.surf_scheme)
        # If the user wants info and layers create a layers info file
        if info or all_:
            grp.get_info()
            # Create the information file
            info = open(grp.name + "_layer_info.txt", 'w')
            info.write(grp.name + " body: \n")
            # Go through the layers in the group's layers
            for i in range(len(grp.layer_surfs)):
                info.write("Number of atoms: " + str(len(grp.layer_atoms[i])) + "\n")
                info.write("Volume: " + str(grp.layer_info[i][0]) + "\n")
                info.write("Surface Area: " + str(grp.layer_info[i][1]) + "\n")
            info.close()
        # Change back to the group directory
        os.chdir(grp.dir)
    # If the user wants to export the interface
    if (iface or all_) and grp.bff is not None:
        get_iface(grp)
        export_iface(grp, [grp, grp.bff], info_file=info)
    # If the user wants a full information file on the group
    if info or all_:
        os.chdir(grp.dir)
        grp.get_info()
        info = open("cell_" + grp.name + "_info.txt", 'w')
        info.write(grp.name + " body: \n")
        info.write("Number of atoms: " + str(len(grp.atoms)) + "\n")
        info.write("Volume: " + str(grp.vol) + "\n")
        info.write("Surface Area: " + str(grp.sa) + "\n")
        info.close()
    if verts or all_:
        if grp.verts is None:
            grp.get_verts()
        write_verts(verts=grp.verts, file_name=grp.name + "_verts", directory=grp.dir)
    if surr_atoms or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        write_pdb(atoms=grp.layer_atoms[1], name=grp.name + "_surr_atoms", directory=grp.dir)
    if ext_atoms or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        write_pdb(atoms=grp.layer_atoms[0], name=grp.name + "_ext_atoms", directory=grp.dir)
    if shell_verts or all_:
        if grp.layer_verts is None:
            # Get the first layer
            grp.get_layers(max_layers=1, build_surfs=False)
        write_verts(grp.layer_verts[0], file_name=grp.name + "_shell_verts", directory=grp.dir)
    if edges or all_:
        if grp.edges is None:
            grp.get_edges()
        write_edges(edges=grp.edges, file_name=grp.name + "_edges", directory=grp.dir)
    if shell_edges or all_:
        if grp.edges is None:
            grp.get_edges()
        if grp.layer_edges is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        write_edges(grp.layer_edges[0], file_name=grp.name + "_shell_edges", directory=grp.dir)
    os.chdir("..")
    # Change back to the system directory
    os.chdir(grp.sys.dir)