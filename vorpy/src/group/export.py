import os
from os import path
from vorpy.src.output import write_pdb
from vorpy.src.output import write_atom_cells
from vorpy.src.output import write_logs
from vorpy.src.output import write_surfs
from vorpy.src.output import write_edges
from vorpy.src.output import write_off_verts


def export_info(grp, directory=None):
    """
    Exports comprehensive information about a group to a text file.
    
    This function generates a detailed report containing:
    - Basic group and system identification
    - System statistics (atom, residue, and chain counts)
    - Network topology information (vertices, edges, surfaces)
    - Geometric analysis (surface area, volume, density)
    
    The information is written to 'info.txt' in the specified directory.
    
    Parameters
    ----------
    grp : Group
        The Group object containing the data to be exported
    directory : str, optional
        Output directory path. If not specified, uses the group's default directory
    
    Returns
    -------
    None
    
    Notes
    -----
    - The function automatically changes to the specified directory before writing
    - If the directory doesn't exist, it will use the group's default directory
    - The output file is encoded in UTF-8 to support special characters
    - The function calls grp.get_info() to ensure all metrics are up to date
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
                  ext_atoms=False, concave_colors=False):
    """
    Exports various components of a Group object to files based on specified parameters.
    This function provides flexible export options for different aspects of a molecular group,
    including atoms, surfaces, edges, vertices, and surrounding elements.

    Parameters
    ----------
    grp : Group
        The Group object containing the data to be exported
    all_ : bool, optional
        If True, exports all possible components of the group. Default is False
    atoms : bool, optional
        If True, exports a PDB file containing only the atoms of the group. Default is False
    atom_surfs : bool, optional
        If True, exports the surfaces associated with each atom. Default is False
    atom_edges : bool, optional
        If True, exports the edges associated with each atom. Default is False
    atom_verts : bool, optional
        If True, exports the vertices associated with each atom. Default is False
    surfs : bool, optional
        If True, exports all surfaces in the group as a single object. Default is False
    sep_surfs : bool, optional
        If True, exports each surface as a separate file, named by their constituent atoms. Default is False
    shell_surfs : bool, optional
        If True, exports all surfaces for the group's shell. Default is False
    edges : bool, optional
        If True, exports all edges in the group as a single object. Default is False
    sep_edges : bool, optional
        If True, exports each edge as a separate file. Default is False
    shell_edges : bool, optional
        If True, exports all edges for the group's shell. Default is False
    verts : bool, optional
        If True, exports all vertices as a single OFF file. Default is False
    sep_verts : bool, optional
        If True, exports each vertex as a separate file. Default is False
    shell_verts : bool, optional
        If True, exports all vertices for the group's shell. Default is False
    layers : int, optional
        Number of layers to export around the group. If -1, exports all layers. Default is -1
    info : bool, optional
        If True, exports group information to info.txt. Default is False
    surr_atoms : bool, optional
        If True, exports atoms directly surrounding the group with intact residues. Default is False
    logs : bool, optional
        If True, exports log files. Default is False
    ext_atoms : bool, optional
        If True, exports the outermost atoms in the group's shell. Default is False
    concave_colors : bool, optional
        If True, exports the concave colors for the surfaces. Default is False
    Returns
    -------
    None

    Notes
    -----
    - All exports are written to the group's directory (grp.dir)
    - If the directory doesn't exist, it will be created
    - Surface colors and schemes are inherited from network settings if not specified
    - For atom-related exports, a subdirectory 'atoms' is created if needed
    - Layer exports require the group to have calculated layers first

    Examples
    --------
    # Export all components of a group
    >>> group_exports(my_group, all_=True)
    # This will create a comprehensive export of the group, including:
    # - A PDB file of all atoms
    # - All surfaces, edges, and vertices
    # - Individual atom components (surfaces, edges, vertices)
    # - Layer information
    # - Group information
    # - Log files
    
    # Export only surfaces and edges
    >>> group_exports(my_group, surfs=True, edges=True)
    # This will create:
    # - A single file containing all surfaces
    # - A single file containing all edges
    # Useful for visualization of the group's structure without atom details
    
    # Export atom-related components
    >>> group_exports(my_group, atoms=True, atom_surfs=True, atom_edges=True)
    # This will create:
    # - A PDB file of the group's atoms
    # - Individual surface files for each atom in the 'atoms' subdirectory
    # - Individual edge files for each atom in the 'atoms' subdirectory
    # Useful for detailed analysis of individual atoms
    
    # Export surrounding atoms and information
    >>> group_exports(my_group, surr_atoms=True, info=True)
    # This will create:
    # - A file containing atoms surrounding the group
    # - An info.txt file with group statistics
    # Useful for analyzing the group's environment and properties
    """
    # Set the surface colors and scheme
    if grp.settings['surf_col'] is None:
        grp.settings['surf_col'] = grp.net.settings['surf_col']
    # Set the surface scheme
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
                         edges=atom_edges or all_, verts=atom_verts or all_, concave_colors=concave_colors)
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
            write_surfs(net=grp.net, surfs=grp.layer_surfs[0], file_name="shell_surfs", directory=grp.dir, concave_colors=concave_colors, ref_surfs=grp.ball_ndxs)
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
        write_edges(grp.net, grp.layer_edges[0], file_name="shell_edges", directory=grp.dir, color=grp.settings['edge_col'])
    # All one big edge file
    if edges or all_:
        write_edges(grp.net, edges=[i for i in range(len(grp.net.edges))], file_name="edges", directory=grp.dir, color=grp.settings['edge_col'])
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
        write_off_verts(grp.net, [i for i in range(len(grp.net.verts))], directory=grp.dir, file_name='verts', color=grp.settings['vert_col'])
    # Export the shell vertices
    if shell_verts or all_:
        if grp.layer_verts is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        write_off_verts(grp.net, grp.layer_verts[0], file_name="shell_verts", directory=grp.dir, color=grp.settings['vert_col'])
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
    # Check to see if there is verts file in the system directory
    for file in os.listdir(grp.sys.files['dir']):
        # Move the verts file to the group directory
        if file.endswith('_verts.txt'):
            os.rename(grp.sys.files['dir'] + "/" + file, grp.dir + "/" + file)
    os.chdir("..")
    # Change back to the system directory
    os.chdir(grp.sys.files['dir'])
