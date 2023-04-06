

def argv_export(my_sys, usr_npt, interfaces=False):
    """
    Exports the specified elements from the system.

    Specifications:

        'default' - (optional) exports network, information, no Sol shell, pdb, set atoms from the system and atoms,
                    surrounding atoms, shell, info from each group
        'all' - exports everything from the full system and each group, needs a warning first
        'info' - exports only the information files for the system and each group
        'surfs' - exports the built surfaces individually from the system and all group surfaces + verts and edges


    :param interfaces:
    :param my_sys:
    :param usr_npt:
    :return:
    """
    # Check for interfaces
    if interfaces:
        for group in my_sys.groups:
            group.exports(iface=True)
        if len(usr_npt) == 0:
            return
    # Go through each of the inputs in the exports
    if len(usr_npt) == 0:
        usr_npt.append(['default'])
    # Export the specified exports
    for npt in usr_npt:
        export_npt(my_sys, npt[0])


def export_npt(my_sys, usr_npt):
    """
    Exports the selection from the usr_npts
    :param my_sys:
    :param usr_npt:
    :return:
    """

    """
    _____________________________________Default____________________________________________
    """
    # If nothing is specified export the defaults
    if usr_npt.lower() == 'default':
        # Export the system's network, information file, the
        my_sys.exports(network=True, pdb=True, set_atoms=True, info=True)
        # Go through each of the system's groups
        for grouping in my_sys.groups:
            # Export the atoms, shell and info file
            grouping.exports(shell=True, atoms=True, info=True, surr_atoms=True, shell_edges=True, shell_verts=True)

    """
    ____________________________________All__________________________________________________________
    """

    # If nothing is specified export the defaults
    if usr_npt.lower() == 'all':
        # Export everything from the system
        my_sys.exports(all_=True)
        # Go through each of the system's groups
        for grouping in my_sys.groups:
            # Export everything from the group
            grouping.exports(all_=True)

    """
    ____________________________________Medium_____________________________________________________________
    """
    # if medium export is specified
    if usr_npt.lower() == 'med':
        # Export everything from the system
        my_sys.exports(network=True, pdb=True, set_atoms=True, info=True)
        # Go through each of the system's groups
        for grouping in my_sys.groups:
            # Export a medium amount of things from the group
            grouping.exports(shell=True, edges=True, verts=True, info=True)


    """
    ____________________________________info_________________________________________________________
    """

    # If nothing is specified export the defaults
    if usr_npt.lower() == 'info':
        # Export the system's network, information file, the
        my_sys.exports(network=True, info=True)
        # Go through each of the system's groups
        for grouping in my_sys.groups:
            # Export the atoms, shell and info file
            grouping.exports(info=True)

    """
    ___________________________________surfs___________________________________________________________
    """
    # If nothing is specified export the defaults
    if usr_npt.lower() == 'surfs':
        # Export the system's network, information file, the
        my_sys.exports(surfaces=True, set_atoms=True, pdb=True, network=True)
        # Go through each of the system's groups
        for grouping in my_sys.groups:
            # Export the atoms, shell and info file
            grouping.exports(shell=True, atoms=True, fill=True, surfaces=True, verts=True, edges=True)

    """
    _______________________________________edges_________________________________________________________
    """
    if usr_npt.lower() == 'e':
        for grouping in my_sys.groups:
            grouping.exports(edges=True)
    """
    _______________________________________verts_________________________________________________________
    """
    if usr_npt.lower() == 'v':
        for grouping in my_sys.groups:
            grouping.exports(verts=True)

    """
    ________________________________simple export______________________________________________________
    """
    if usr_npt.lower() == 'simple':
        for grouping in my_sys.groups:
            grouping.exports(atoms=True, info=True, shell=True)
    """
    _____________________________________________heavy export__________________________________________
    """
