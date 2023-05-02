import os.path
from System.sys_funcs.output.output import export_min1, export_min2, export_med, export_large, export_all, other_exports, set_sys_dir


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
        if npt[0].lower() in {'dir', 'directory'} and len(npt) == 2 and (os.path.isdir(npt[1]) or npt[1] == 'gsu_logs'):
            if npt[1] == 'gsu_logs':
                set_sys_dir(my_sys, "C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/{}/logs".format(my_sys.name))
            elif npt[1] == 'gsu_logs1':
                set_sys_dir(my_sys, "/Users/jackericson/Library/CloudStorage/OneDrive-GeorgiaStateUniversity/GSU NSC/Jack/Vorpy/test_data/{}/logs".format(my_sys.name))
            else:
                my_sys.dir = npt[1]
    for npt in usr_npt:
        if npt[0].lower() in {'dir', 'directory'}:
            continue
        export_npt(my_sys, npt[0])


def export_npt(my_sys, usr_npt=None):
    """
    Exports the selection from the usr_npts
    :param my_sys:
    :param usr_npt:
    :return:
    """

    # If nothing is specified export the defaults
    if usr_npt is None or usr_npt.lower() in {'default', '2', 'medium', '', 'med'}:
        export_med(sys=my_sys)

    # Small export
    elif usr_npt.lower() in {"tiny", "i", "info", "0", "smallest"}:
        export_min1(my_sys)

    # Medium small export
    elif usr_npt.lower() in {"small", "s", "1"}:
        export_min2(my_sys)

    # Large Export
    elif usr_npt.lower() in {"large", "l", "3"}:
        export_large(my_sys)

    # Export all
    elif usr_npt.lower() in {'all', 'a', 'everything'}:
        export_all(my_sys)

    else:
        other_exports(my_sys, usr_npt)
