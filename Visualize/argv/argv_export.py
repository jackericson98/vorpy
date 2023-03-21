

def argv_export(my_sys, usr_npt):
    """
    Exports the specified elements from the system.

    Specifications:

        'default' - (optional) exports network, information, no Sol shell, pdb, set atoms from the system and atoms,
                    surrounding atoms, shell, info from each group
        'all' - exports everything from the full system and each group, needs a warning first
        'info' - exports only the information files for the system and each group
        'surfs' - exports the built surfaces individually from the system and all group surfaces + verts and edges


    :param my_sys:
    :param usr_npt:
    :return:
    """
    # Go through each of the inputs in the exports
    if len(usr_npt) == 0:
        usr_npt.append(['default'])
    # Export the specified exports
    for npt in usr_npt:
        export_npt(my_sys, npt[0])


def export_npt(my_sys, usr_npt):

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
            grouping.exports(shell=True, atoms=True, info=True, surr_atoms=True, shell_edges=True)

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

    """
    _______________________________________verts_________________________________________________________
    """

    """
    ________________________________simple export______________________________________________________
    """
    """
    _____________________________________________heavy export__________________________________________
    """


