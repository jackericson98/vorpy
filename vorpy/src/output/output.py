import os
import shutil
from vorpy.src.output import write_atom_cells


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
            os.makedirs(group.dir, exist_ok=True)
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
        os.makedirs(group.dir, exist_ok=True)
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
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name):
            group.dir = sys.files['dir'] + '/' + group.name
            # Catch for if the group name is too long
            try:
                os.makedirs(group.dir, exist_ok=True)
            except FileNotFoundError:
                group.dir = sys.files['dir'] + '/group'
        # Do the group exports
        group.exports(shell_surfs=True, surfs=True, shell_edges=True, edges=True, shell_verts=True, verts=True,
                      logs=True, atoms=True, surr_atoms=True)
        # # Check to see if the verts are in the system directory and if so move them to the group folder
        # if os.path.exists(sys.files['dir'] + '/verts.txt'):
        #     shutil.move(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt',
        #                 group.dir + '/' + group.settings['net_type'] + '_verts.txt')
    # Export the interfaces
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(surfs=True, atoms=True, edges=True, logs=True, verts=True, info=True)


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
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name):
            group.dir = sys.files['dir'] + '/' + group.name
            os.makedirs(group.dir, exist_ok=True)
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
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name):
            group.dir = sys.files['dir'] + '/' + group.name
            os.makedirs(group.dir, exist_ok=True)
        group.dir = sys.files['dir'] + '/' + group.name
        os.makedirs(group.dir, exist_ok=True)
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

