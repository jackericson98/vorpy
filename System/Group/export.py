from System.sys_funcs.output.output import *
from System.Group.sort import get_iface
from System.sys_funcs.output.verts import write_verts
from System.sys_funcs.output.edges import write_edges


def export_iface_verts(grp, directory=None):
    """
    Exports the interfacial vertices between the group and its bff
    :param grp: Group object for exporting
    :param directory: Directory to export to
    """
    # Move to the directory
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)
    # write the vertices
    write_verts(grp.iface_verts, directory=directory, file_name=grp.sys.net.type + "_verts")


def export_iface_edges(grp, directory=None):
    """
    Exports the edges of the interface
    :param grp: Group to pull the interface from
    :param directory: Output directory for the interface
    """
    # Move to the directory
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)
    # write the vertices
    write_edges(grp.iface_edges, directory=directory, file_name=grp.sys.net.type + "_edges")


def export_iface(grp, info_file=True, interface_atoms=False, directory=None):
    """
    Exports the information from the given interface as a txt file
    :param grp: Group object for interface
    :param info_file: Whether to export a txt file with info on the interface or not
    :param interface_atoms: Whether to export a pdb file with the atoms around the interface or not
    :param directory: Output directory for the interface stuff
    """
    # Move to the directory
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)
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
    write_verts(grp.iface_verts, file_name=interface_name + "_verts")
    write_edges(grp.iface_edges, file_name=interface_name + "_edges")
    # Check to see of the user wants to export the interface's atoms
    if interface_atoms:
        # Get the two sets of interface atoms
        write_pdb(grp.iface_atoms, interface_name + "_" + grp.name + "_atoms", grp.sys)
        write_pdb(grp.bff.iface_atoms, interface_name + "_" + grp.bff.name + "_atoms", grp.bff.sys)
    # Check to see if the user wants to export the interface's information
    if info_file:
        export_iface_info(grp=grp, directory=directory)


def export_iface_info(grp, directory=None):
    """
    Exports the information for an interface
    :param grp: The group that holds the interface information
    :param directory: Output directory for the group interface info
    """
    # Move to the directory
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)
    # Create the file
    with open("info.txt", 'w', encoding='utf-8') as info:
        # Write the main header
        info.write(grp.name + " - " + grp.bff.name + " interface \n\n")
        # Information sub header
        info.write("Interface:\n\n")
        # Write the information
        info.write("  {} Surfaces, {} {} atoms, {} {} atoms\n".format(len(grp.iface_surfs), len(grp.atoms), grp.name,
                                                                      len(grp.bff.atoms), grp.bff.name))
        # Network counts
        info.write("  {} Vertices, {} Edges\n\n".format(len(grp.iface_verts), len(grp.iface_edges),
                                                        len(grp.iface_surfs)))
        # Write the analysis header
        info.write("\nAnalysis:\n\n")
        # Write the analysis
        info.write(u"  Surface Area = {:.5f} \u212B\u00B2, Average Curvature = {:.5}\n\n"
                   .format(grp.iface_sa, grp.iface_curv))
        # Surfaces header
        info.write("\nSurfaces:\n\n")
        # Go through each of the surfaces in the group
        for surf in grp.iface_surfs:
            info.write("  Surface {} - \n".format(surf.ndx))
            info.write("    Surface Area = {:.5f} \u212B\u00B2\n".format(surf.sa))
            info.write("    Volume contributions = {:.5f}, {:.5f} \u212B\u00B3\n".format(surf.vols[0],
                                                                                                        surf.vols[1]))
            info.write("    Gaussian Curvature = {:.5f}\n".format(surf.curv))


def export_info(grp, directory=None):
    """
    Exports the information for a group
    :param grp: The group to have information exported
    :param directory: Output directory for the group
    """
    # Move to the directory
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)
    # Change to the directory of the group
    os.chdir(grp.dir)
    # Get the information for the group
    grp.get_info()
    # Open the export information file
    with open("info.txt", 'w', encoding="utf-8") as info:
        # Write the main header
        info.write("{} - {}\n\n".format(grp.name, grp.sys.name))
        # System counts header
        info.write("Group system information:\n")
        # System counts
        info.write("  {} Atoms, {} Residues, {} Chains\n\n".format(len(grp.atoms), len(grp.residues),
                                                                   len(grp.chains)))
        # Network counts header
        info.write("Group Network information:\n")
        # Network counts
        info.write("  {} Vertices, {} Edges, {} Surfaces\n\n".format(len(grp.verts), len(grp.edges),
                                                                     len(grp.surfs)))
        # Analysis header
        info.write("Analysis:\n")
        # Analysis information
        info.write(u"  Surface Area: {:.5f} \u212B\u00B2, Volume: {:.5f} \u212B\u00B3\n\n".format(grp.sa,
                                                                                                  grp.vol))


def group_exports(grp, all_=False, iface=False, atoms=False, surfs=False, sep_surfs=False, edges=False,
                  sep_edges=False, verts=False, sep_verts=False, layers=-1, info=False, surr_atoms=False,
                  ext_atoms=False, shell=False):
    """
    Exports specified export types for the group using bools
    :param grp: Group object for export
    :param all_: All possible exports for the group will be exported to the group directory
    :param iface: Exports the interface for the group, bff must be specified first
    :param atoms: Exports a new pdb file containing only the atoms of the group
    :param surfs: Exports all surfaces in the group as one object
    :param sep_surfs: Exports all surfaces in the group as separate files, named by their atoms
    :param edges: Exports all edges for the group
    :param sep_edges: Exports separate edges for the group
    :param verts: Exports the vertices of the group as an off file
    :param sep_verts: Exports the separate vertices for the group
    :param layers: Exports all layers surrounding the group, unless num_layers is specified
    :param info: Exports the information for the group
    :param surr_atoms: Exports the atoms directly surrounding the group (residues intact)
    :param ext_atoms: Exports the outermost atoms in the group's set of atoms (must be a part of shell)
    :param shell: Whether to do the shell objects or the full object
    """
    # Set the surface colors and scheme
    if grp.surf_color is None:
        grp.surf_color = grp.sys.net.surf_col
    if grp.surf_scheme is None:
        grp.surf_scheme = grp.sys.net.surf_scm
    # Get the surfaces if they haven't been got
    if grp.surfs is None or len(grp.surfs) == 0:
        grp.build_surfs()
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
        write_pdb(atoms=grp.atoms, file_name="atoms", sys=grp.sys)
    # If the user wants to export the shell for the group
    if shell or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # noinspection PyUnresolvedReferences
        if grp.layer_surfs is not None and len(grp.layer_surfs) > 0:
            write_surfs(surfs=grp.layer_surfs[0], file_name="shell", directory=grp.dir)
        if edges:
            if grp.edges is None:
                grp.get_edges()
            if grp.layer_edges is None:
                grp.get_layers(max_layers=1, build_surfs=False)
            write_edges(grp.layer_edges[0], file_name="shell_edges", directory=grp.dir)
        if verts:
            if grp.layer_verts is None:
                # Get the first layer
                grp.get_layers(max_layers=1, build_surfs=False)
            write_verts(grp.layer_verts[0], file_name="shell_verts", directory=grp.dir)
    if edges or all_:
        if grp.edges is None:
            grp.get_edges()
        write_edges(edges=grp.edges, file_name="edges", directory=grp.dir)
    if sep_verts:
        os.mkdir(grp.dir + "/verts")
        if grp.layer_verts is None:
            # Get the first layer
            grp.get_layers(max_layers=1, build_surfs=False)
        for i, vert in enumerate(grp.layer_verts[0]):
            write_verts([vert], str(i), directory=grp.dir + "/verts")
    if sep_edges:
        os.mkdir(grp.dir + "/edges")
        if grp.layer_edges is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        for i, edge in enumerate(grp.layer_edges[0]):
            write_edges([edge], str(i), directory=grp.dir + "/edges")
    # If the user wants a filled shell for the group
    if surfs or all_:
        grp.build_surfs()
        write_surfs(surfs=grp.surfs, file_name="fill", directory=grp.dir)
    # If the user wants separate surfaces for the group
    if sep_surfs or all_:
        i = 1
        my_dir = grp.dir + "/surfaces"
        while os.path.exists(my_dir):
            if my_dir[-1] == 's':
                my_dir += '__'
            my_dir = my_dir[:-2] + str(i)
            i += 1
        os.mkdir(my_dir)
        for surf in grp.surfs:
            write_surfs([surf], file_name="_".join([str(_) for _ in surf.ndx]), directory=my_dir)
    # If the user wants layers
    if layers > 0 or all_:
        # First check to see if the number of layers is greater than 1
        if grp.layer_atoms is None or len(grp.layer_atoms) <= 1:
            grp.get_layers(max_layers=layers)
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
            write_pdb(grp.layer_atoms[i + 1], file_name=str(i) + "_atoms", sys=grp.sys)
            write_surfs(grp.layer_surfs[i], file_name=str(i) + "_surfs")
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
    if (iface or all_) and grp.bff is not None and grp.bff != []:
        get_iface(grp)
        if len(grp.iface_surfs) > 0:
            export_iface(grp, info_file=info)
    # If the user wants a full information file on the group
    if info or all_:
        export_info(grp)
    if verts or all_:
        if grp.verts is None:
            grp.get_verts()
        write_verts(verts=grp.verts, file_name="verts", directory=grp.dir)
    if surr_atoms or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        write_pdb(atoms=grp.layer_atoms[1], file_name="surr_atoms", directory=grp.dir, sys=grp.sys)
    if (ext_atoms or all_) and len(grp.atoms) > 15:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        write_pdb(atoms=grp.layer_atoms[0], file_name="ext_atoms", directory=grp.dir)

    os.chdir("..")
    # Change back to the system directory
    os.chdir(grp.sys.dir)
