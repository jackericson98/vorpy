from System.sys_funcs.output.output import *
from System.sys_funcs.output.verts import write_off_verts
from System.sys_funcs.output.edges import write_edges
from System.sys_funcs.output.net import write_logs


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
        info.write("  {} Atoms, {} Residues, {} Chains\n\n".format(len(grp.ball_ndxs), len(grp.rsds), len(grp.chns)))
        # Network counts header
        info.write("Group Network information:\n")
        # Network counts
        info.write("  {} Vertices, {} Edges, {} Surfaces\n\n".format(len(grp.net.verts), len(grp.net.edges), len(grp.net.surfs)))
        # Analysis header
        info.write("Analysis:\n")
        # Analysis information
        info.write(u"  Surface Area: {:.5f} \u212B\u00B2, Volume: {:.5f} \u212B\u00B3, Density: {:.5f}\n\n"
                   .format(grp.sa, grp.vol, grp.density))


def group_exports(grp, all_=False, atoms=False, atom_surfs=False, atom_edges=False, atom_verts=False, surfs=False,
                  sep_surfs=False, shell_surfs=False, edges=False, sep_edges=False, shell_edges=False,
                  verts=False, sep_verts=False, shell_verts=False, layers=-1, info=False, surr_atoms=False, logs=False,
                  ext_atoms=False):
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
    if grp.settings['surf_col'] is None:
        grp.settings['surf_col'] = grp.net.settings['surf_col']
    if grp.settings['surf_scheme'] is None:
        grp.settings['surf_scheme'] = grp.net.settings['surf_scheme']
    # Get the surfaces if they haven't been got
    if grp.net.surfs is None or len(grp.net.surfs) == 0:
        return
    # Create the output directory inside the system's directory
    if grp.dir is None:
        i = 1
        my_dir = grp.sys.files['dir'] + "/" + grp.name
        first = True
        while os.path.exists(my_dir):
            if first:
                my_dir += "__"
                first = False
            my_dir = my_dir[:-(1 + len(str(i)))] + '_' + str(i)
            i += 1
        grp.dir = my_dir
        os.mkdir(grp.dir)
    # Go back to the group directory
    os.chdir(grp.dir)
    # If the user wants to export the atoms for the group
    if atoms or all_:
        if grp.sys.files['base_file'][-3:] == 'txt':
            pass
        else:
            write_pdb(atoms=grp.ball_ndxs, file_name="group_atoms", sys=grp.sys)
    # If the atoms surfaces are selected go for it
    if atom_verts or atom_edges or atom_surfs or all_:
        if not path.exists(grp.dir + '/atoms'):
            os.mkdir(grp.dir + '/atoms')
        write_atom_cells(grp.net, atoms=grp.ball_ndxs, directory=grp.dir + '/atoms', surfs=atom_surfs or all_,
                         edges=atom_edges or all_, verts=atom_verts or all_)
        os.chdir(grp.dir)
    # Export the log file
    if logs or all_:
        write_logs(grp)
    # If the user wants to export the shell for the group
    if shell_surfs or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # noinspection PyUnresolvedReferences
        if grp.layer_surfs is not None and len(grp.layer_surfs) > 0:
            write_surfs(net=grp.net, surfs=grp.layer_surfs[0], file_name="shell_surfs", directory=grp.dir)
    # If the user wants all of the surfaces in one file
    if surfs or all_:
        write_surfs(grp.net, [i for i in range(len(grp.net.surfs))], 'surfs')
    # Separate surfaces
    if sep_surfs or all_:
        # Make the surfaces directory
        if not os.path.exists(grp.dir + '/surfs'):
            os.mkdir(grp.dir + '/surfs')
        # Create the surfaces' files
        for j, my_surf in grp.net.surfs.iterrows():
            write_surfs(grp.net, [j], file_name='b{}_b{}'.format(*my_surf['balls']), directory=grp.dir + '/surfs')
    # Shell edges
    if shell_edges or all_:
        if grp.layer_edges is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        write_edges(grp.net, grp.layer_edges[0], file_name="shell_edges", directory=grp.dir, color=[1, 0, 0])
    # All one big edge file
    if edges or all_:
        write_edges(grp.net, edges=[i for i in range(len(grp.net.edges))], file_name="edges", directory=grp.dir, color=[0, 1, 0])
    # If the separate edges are called
    if sep_edges or all_:
        # Make the edges directory
        if not os.path.exists(grp.dir + '/edges'):
            os.mkdir(grp.dir + '/edges')
        for j, my_edge in grp.net.edges.iterrows():
            write_edges(grp.net, [j], 'b{}_b{}_b{}'.format(*my_edge['balls']), directory=grp.dir + '/edges')
    # Run the separate vertices
    if sep_verts:
        # Make the vertices directory
        if not path.exists(grp.dir + '/verts'):
            os.mkdir(grp.dir + "/verts")
        for j, vert in grp.net.verts.iterrows():
            write_off_verts(grp.net, [j], 'b{}_b{}_b{}_b{}'.format(*vert['balls']), directory=grp.dir + "/verts")
    # Export all the vertices in one file
    if verts or all_:
        write_off_verts(grp.net, [i for i in range(len(grp.net.verts))], directory=grp.dir, file_name='verts', color=[0, 0, 1])
    # Export the shell vertices
    if shell_verts or all_:
        if grp.layer_verts is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        write_off_verts(grp.net, grp.layer_verts[0], file_name="shell_verts", directory=grp.dir)
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
            write_surfs(grp.net, grp.layer_surfs[i], file_name=str(i) + "_surfs")
        # If the user wants info and layers create a layers info file
        if info or all_:
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
    # If the user wants a full information file on the group
    if info or all_:
        export_info(grp)
    # Surrounding atoms
    if surr_atoms or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        try:
            write_pdb(atoms=grp.layer_atoms[1], file_name="surr_atoms", directory=grp.dir, sys=grp.sys)
        except IndexError:
            pass
    if (ext_atoms or all_) and len(grp.atoms) > 15:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        write_pdb(sys=grp.sys, atoms=grp.layer_atoms[0], file_name="ext_atoms", directory=grp.dir)

    os.chdir("..")
    # Change back to the system directory
    os.chdir(grp.sys.files['dir'])
