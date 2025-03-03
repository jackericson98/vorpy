from System.sys_funcs.output.atoms import write_pdb, write_atom_cells
from System.sys_funcs.output.surfs import write_surfs
from System.sys_objs.atom import special_radii
import os
from os import path
import shutil

###################################################### Export Functions ################################################


def export_micro(sys):
    """
    Smallest output function. Outputs the information for the system, the groups, and the system's interfaces.
    """
    # Export the information for the system.
    sys.exports(info=True)
    # Loop through the groups in the system
    for group in sys.group:
        # Set up the group directory
        if group.dir is None:
            group.dir = sys.files['dir'] + '/' + group.name
            os.mkdir(group.dir)
        # Export the information for the group
        group.exports(info=True)
    # Loop through the interfaces for the groups.
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            # Export the interface information
            iface.export(info=True)


def export_tiny(sys):
    """
    Second smallest of the exports. Outputs are:

    System:
        1. General Information
        2. Set balls script for pymol
        4. The PDB file
        5. The balls file
    Groups:
        1. General Information
        2. Shell for the group
        3. Logs for the group
    Interfaces:
        1.
    """
    sys.exports(info=True, set_atoms=True, pbd=True, balls=True)
    for group in sys.groups:
        group.dir = sys.files['dir'] + '/' + group.name
        os.mkdir(group.dir)
        group.export(info=True, shell=True, logs=True)
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(info=True)


def export_med(sys):
    """
    Medium export. Exports the pdb, the set atoms script and the general information for the system. The group gets the
    logs, the shell for the group, the surfaces for the group, the full set of edges, the shell edges, and the vertices
    """
    # Export the system exports
    sys.exports(pdb=True, set_atoms=True, info=True)

    # Loop through the groups and give their exports
    for group in sys.groups:
        # Set and make the group directory
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name + '_' + group.settings['net_type']):
            group.dir = sys.files['dir'] + '/' + group.name + '_' + group.settings['net_type']
            os.mkdir(group.dir)
        # Do the group exports
        group.exports(shell_surfs=True, surfs=True, shell_edges=True, edges=True, shell_verts=True, verts=True,
                      logs=True, atoms=True, surr_atoms=True)
        # Check to see if the verts are in the system directory and if so move them to the group folder
        if os.path.exists(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt'):
            shutil.move(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt',
                        group.dir + '/' + group.settings['net_type'] + '_verts.txt')
    # Export the interfaces
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(surfs=True, atoms=True, info=True)


def export_large(sys):
    """
    Large group exports. Exports the basic system files and the shell vertices, the shell surfaces, the information,
    the edges, the vertices, the atosm the surrounding atoms, the logs, the atom surfaces, the atom edges, and the
    atom vertices for each group
    """
    # Export the system exports
    sys.exports(pdb=True, set_atoms=True, info=True)
    # Loop through the groups and export the listed items
    for group in sys.groups:
        # Set and make the group directory
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name + '_' + group.settings['net_type']):
            group.dir = sys.files['dir'] + '/' + group.name + '_' + group.settings['net_type']
            os.mkdir(group.dir)
        # Export the group exports
        group.exports(shell_verts=True, shell_edges=True, shell_surfs=True, info=True, edges=True, verts=True,
                      atoms=True, surr_atoms=True, logs=True, atom_surfs=True, atom_edges=True, atom_verts=True)
        # Check to see if the verts are in the system directory and if so move them to the group folder
        if os.path.exists(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt'):
            shutil.move(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt',
                        group.dir + '/' + group.settings['net_type'] + '_verts.txt')
    # Export the interfaces
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(balls=True, surfs=True, edges=True, verts=True, info=True)


def export_all(sys):
    """
    Export all. Exports everything there is to export and makes a massive comprehensive set of files that will take a
    lot of space
    """
    # Export the system stuff
    sys.exports(pdb=True, info=True, set_atoms=True)
    # For each group in the system export the
    for group in sys.groups:
        # Set and make the group directory
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name + '_' + group.settings['net_type']):
            group.dir = sys.files['dir'] + '/' + group.name + '_' + group.settings['net_type']
            os.mkdir(group.dir)
        group.dir = sys.files['dir'] + '/' + group.name
        os.mkdir(group.dir)
        group.exports(atoms=True, shell=True, surfs=True, info=True, ext_atoms=True, sep_surfs=True, sep_edges=True,
                      sep_verts=True, verts=True, edges=True, surr_atoms=True, logs=True)

        # Check to see if the verts are in the system directory and if so move them to the group folder
        if os.path.exists(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt'):
            shutil.move(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt',
                        group.dir + '/' + group.settings['net_type'] + '_verts.txt')
    # Make the
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(all=True)


################################################ Other Exports #########################################################


def other_exports(sys, usr_npt):
    """

    :param sys:
    :param usr_npt:
    :return:
    """
    # If the first word is atom
    if usr_npt.lower() in {"a", "atoms"}:
        write_atom_cells(sys.net.atoms['num'], sys.files['dir'])
    # If the first word is logs
    elif usr_npt.lower() in {'logs', 'lgs'}:
        for group in sys.groups:
            group.exports(logs=True)
        sys.exports(pdb=True, set_atoms=True)
    # If the first word is shell
    elif usr_npt.lower() in {'shell', 'shl'}:
        for grp in sys.groups:
            grp.exports(shell_surfs=True)
    # If the first word is network
    elif usr_npt.lower() in {'net', 'network'}:
        sys.exports(network=True)


####################################################### Main Funcs #####################################################


def set_sys_dir(sys, dir_name=None):
    """
    Sets the directory for the output data. If the directory exists add 1 to the end number
    :param sys: System to assign the output directory to
    :param dir_name: Name for the directory
    :return:
    """

    # Make sure a user_data path exists
    if sys.files['root_dir'] is not None and not os.path.exists(sys.files['root_dir'] + "/Data/user_data"):
        os.mkdir(sys.files['root_dir'] + "/Data/user_data")
    elif sys.files['root_dir'] is None and not os.path.exists("./Data/user_data"):
        if not os.path.exists('./Data'):
            os.mkdir(os.path.abspath('.') + '/Data/user_data')
        else:
            os.mkdir(os.path.abspath('./Data') + '/user_data')

    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        if sys.files['dir'] is not None:
            dir_name = sys.files['dir'] + '/' + sys.name
        elif sys.files['root_dir'] is not None:

            dir_name = sys.files['root_dir'] + "/Data/user_data/" + sys.name
        else:
            dir_name = os.getcwd() + "/Data/user_data/" + sys.name
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = '_' + str(i)
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
    sys.files['dir'] = dir_name + i_str


def export_sys(sys, all_=False, pdb=False, full_network_object=False, alter_atoms_script=False, info=False):
    """
        Prepares the output directory and system for output. Keeps things consistent
        :return:
        """
    # Check to see if the pdb directory is suitable
    if sys.files['dir'] is None:
        sys.set_output_directory()
    # If the information is requested, export it
    if info or all_:
        os.chdir(sys.files['dir'])
        export_sys_info(sys)
    if pdb or all_:
        os.chdir(sys.files['dir'])
        # Export a pdb file for the system
        write_pdb([_ for i, _ in sys.balls.iterrows()], sys.name, sys)
        os.chdir(sys.files['dir'])
    if (full_network_object or all_) and sys.net.build_surfs:
        # Export a full system
        write_surfs(sys.net.surfs, "full_sys", directory=sys.files['dir'])
    # Write the alter atoms script
    if alter_atoms_script or all_:
        pass
        os.chdir(sys.files['dir'])
        set_pymol_atoms(sys)


def set_pymol_atoms(sys):

    """
    Creates a script to set the radii of the spheres in pymol
    :param sys:
    :return:
    """
    # If we have special circumstances for the atoms in our base file, output the already created set pymol atoms
    if sys.type == 'foam' or sys.type == 'coarse':
        # Get the directory for the base_file and copy the set atoms file
        try:
            shutil.copyfile(path.dirname(sys.files['base_file']) + '/set_atoms.pml', sys.files['dir'] + '/sys/set_atoms.pml')
        except FileNotFoundError:
            # Create the file
            with open('set_atoms.pml', 'w') as file:
                for i, ball in sys.balls.iterrows():
                    file.write(
                        "alter r. {} and n. {}, vdw={}\n".format(ball['res_name'], ball['name'], ball['rad']))
                file.write("\nrebuild")
        return
    # Check to see if the atoms in the system are all accounted for
    for i, res in enumerate(sys.residues):
        if res.name not in special_radii:
            special_radii[res.name] = {sys.balls['name'][j]: round(sys.balls['rad'][j], 2) for j in res.atoms}
    # Create the file
    with open('set_atoms.pml', 'w') as file:
        # Write the change radii script for the system's set atomic radii
        for radius in sys.element_radii:
            if radius != '':
                file.write("alter {} and e. {}, vdw={}\n".format(sys.name, radius, sys.element_radii[radius]))
        # Change the radii for special atoms
        for res in special_radii:
            for atom in special_radii[res]:
                res_str = "r. {} ".format(res) if res != "" else ""
                file.write("alter {} and {}and n. {}, vdw={}\n".format(sys.name, res_str, atom, special_radii[res][atom]))
        # Rebuild the system
        file.write("\nrebuild")


def export_sys_info(sys):
    # Open the file
    with open(sys.name + "_info.txt", 'w') as info:
        # Write the header
        info.write(sys.name + " Network")
        # Write the chain header
        info.write("\n\n++++++++++++++++++++++++  Chains  +++++++++++++++++++++++++++++++\n\n")
        # Go through the chains in the system
        if sys.chains is not None:
            for chain in sys.chains:
                # Write the chain header
                info.write("Chain {} - {} atoms, {} residues\n\n".format(chain.name, len(chain.atoms), len(chain.residues)))
                # Quick check to see if the chain has been calculated
                if chain.vol is not None and chain.vol < 0:
                    # Write the chain information
                    info.write("  Volume = {}, Surface Area = {}\n\n\n".format(chain.vol, chain.sa))
        # Draw a separating line
        info.write("\n\n++++++++++++++++++++++++  Groups  +++++++++++++++++++++++++++++++\n\n")
        for group in sys.groups:
            # Write the group header
            info.write("Group {} - {} atoms, {} residues, {} chains\n\n".format(group.name, len(group.atms), len(group.rsds), len(group.chns)))
            # Write the group info
            info.write("  Volume = {}, Surface Area = {}\n\n\n".format(group.vol, group.sa))

