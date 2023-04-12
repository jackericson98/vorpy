import os
from System.sys_funcs.output.atoms import write_pdb, write_atom_cells
from System.sys_funcs.output.surfs import write_surfs
from System.sys_funcs.output.net import export_net_logs

###################################################### Export Functions ################################################


def export_min1(sys):
    sys.exports(info=True)
    for group in sys.group:
        group.exports(info=True)


def export_min2(sys):
    sys.exports(info=True, set_atoms=True)
    for group in sys.groups:
        group.export(info=True, shell=True)


def export_med(sys):
    sys.exports(pdb=True, set_atoms=True, info=True, network=True, logs=True)
    for group in sys.groups:
        group.exports(shell=True, info=True, edges=True, atoms=True)


def export_large(sys):
    sys.exports(pdb=True, set_atoms=True, info=True, logs=True, network=True)
    for group in sys.groups:
        group.exports(shell=True, info=True, edges=True, verts=True, atoms=True, surr_atoms=True)
        os.mkdir(group.dir + "/atoms")
        write_atom_cells(group.atoms, directory=group.dir + "/atoms")


def export_all(sys):
    sys.exports(pdb=True, info=True, network=True, logs=True, set_atoms=True)
    for group in sys.groups:
        group.exports(atoms=True, shell=True, surfs=True, info=True, ext_atoms=True, sep_surfs=True, sep_edges=True, sep_verts=True, verts=True, edges=True)
    os.mkdir(sys.dir + "/atoms")
    write_atom_cells(sys.atoms, directory=sys.dir + "/atoms", verts=True, edges=True)

################################################ Other Exports #########################################################


def other_exports(sys, usr_npt):
    """

    :param sys:
    :param usr_npt:
    :return:
    """
    # If the first word is atom
    if usr_npt.lower() in {"a", "atoms"}:
        write_atom_cells(sys.atoms, sys.dir)
    # If the first word is logs
    elif usr_npt.lower() in {'logs', 'lgs'}:
        sys.exports(logs=True)


####################################################### Main Funcs #####################################################


def set_sys_dir(sys, dir_name=None):
    """
    Sets the directory for the output data. If the directory exists add 1 to the end number
    :param sys: System to assign the output directory to
    :param dir_name: Name for the directory
    :return:
    """
    if not os.path.exists("./Data/user_data"):
        os.mkdir("./Data/user_data")
    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        if sys.vpy_dir is not None:
            dir_name = sys.vpy_dir + "/Data/user_data/" + sys.name
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
    sys.dir = dir_name + i_str


def export_sys(sys, all_=False, network=False, pdb=False, surfaces=False, full_network_object=False,
               alter_atoms_script=False, info=False, logs=False):
    """
        Prepares the output directory and system for output. Keeps things consistent
        :return:
        """
    # Check to see if the pdb directory is suitable
    if sys.dir is None:
        if os.path.dirname(sys.base_file)[-9:] != 'test_data':
            sys.dir = os.path.dirname(sys.base_file)
        else:
            sys.set_output_directory()
    if network or all_:
        os.chdir(sys.dir)
        # Export the network
        sys.export_net()
    if pdb or all_:
        if not os.path.exists(sys.dir + '/sys'):
            os.mkdir(sys.dir + "/sys")
        os.chdir(sys.dir + "/sys")
        # Export a pdb file for the system
        write_pdb(sys.atoms, sys.name, sys)
        os.chdir(sys.dir)
    if surfaces or all_:
        if not os.path.exists(sys.dir + '/surfs'):
            os.mkdir(sys.dir + "/surfs")
        # Export a pdb file for the system
        for surf in sys.net.surfs:
            write_surfs(surfs=[surf], file_name="_".join([str(_) for _ in surf.ndx]), directory=sys.dir + "/surfs")
        os.chdir(sys.dir)
    if (full_network_object or all_) and sys.net.build_surfs:
        if not os.path.exists(sys.dir + '/sys'):
            os.mkdir(sys.dir + "/sys")
        # Export a full system
        write_surfs(sys.net.surfs, "full_sys", directory=sys.dir + "/sys")
    # Write the alter atoms script
    if alter_atoms_script or all_:
        if not os.path.exists(sys.dir + '/sys'):
            os.mkdir(sys.dir + "/sys")
        os.chdir(sys.dir + "/sys")
        set_pymol_atoms(sys)
    # If the information is requested, export it
    if info or all_:
        if not os.path.exists(sys.dir + "/sys"):
            os.mkdir(sys.dir + "/sys")
        os.chdir(sys.dir + "/sys")
        export_sys_info(sys)
    # Export the log file
    if logs or all_:
        if not os.path.exists(sys.dir + "/sys"):
            os.mkdir(sys.dir + "/sys")
        os.chdir((sys.dir + "/sys"))
        export_net_logs(sys.net)
    os.chdir(sys.dir)


def set_pymol_atoms(sys):

    """
    Creates a script to set the radii of the spheres in pymol
    :param sys:
    :return:
    """
    # Create the file
    with open('set_atoms.pml', 'w') as file:
        # Write the change radii script for the system's set atomic radii
        for radius in sys.radii:
            if radius != '':
                file.write("alter (elem {}), vdw={}\n".format(radius, sys.radii[radius]))
        # Change the radii for special atoms
        for res in sys.special_radii:
            for atom in sys.special_radii[res]:
                res_str = "residue {} ".format(res) if res != "" else ""
                file.write("alter ({}name {}), vdw={}\n".format(res_str, atom, sys.special_radii[res][atom]))
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
            info.write("Group {} - {} atoms, {} residues, {} chains\n\n".format(group.name, len(group.atoms), len(group.residues), len(group.chains)))
            # Write the group info
            info.write("  Volume = {}, Surface Area = {}\n\n\n".format(group.vol, group.sa))

