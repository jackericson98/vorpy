import os
from vorpy.src.output import export_micro
from vorpy.src.output import export_tiny
from vorpy.src.output import export_med
from vorpy.src.output import export_large
from vorpy.src.output import export_all
from vorpy.src.output import other_exports


FORMAT_COMMANDS = {'ft', 'file_type', 'file_format', 'format', 'mesh_format'}
FORMAT_VALUES = {'off': 'off', 'ply': 'ply', 'vtp': 'vtp', 'vtk': 'vtp'}
ONLY_COMMANDS = {'only', 'just'}


def _set_mesh_format(my_sys, value):
    """Apply a geometry-format modifier to the system and all current groups."""
    value = str(value).strip().lower().lstrip('.')
    file_type = FORMAT_VALUES.get(value)
    if file_type is None:
        raise ValueError(
            f"Unsupported geometry file format {value!r}. "
            "Choose off, ply, or vtp."
        )

    my_sys.file_type = file_type
    for group in getattr(my_sys, 'groups', []) or []:
        if getattr(group, 'settings', None) is None:
            group.settings = {}
        group.settings['file_type'] = file_type

    print(f"Geometry file format: {file_type.upper()}")



def argv_export(my_sys, usr_npt, add_on=None):
    """
    Command-line export handler.

    Default behavior:
    - If no export type is specified, export the large preset.
    - Format and directory commands modify the export but do not replace its preset.
    - If only modifiers are specified, export the large preset.
    - ``only`` explicitly restricts output to the export names that follow it.
    - If explicit exports are specified, only export those.
    """

    if usr_npt is None:
        usr_npt = []

    # First pass: handle export modifiers and collect actual export commands.
    export_commands = []
    skip_archive = False

    for npt in usr_npt:
        if len(npt) == 0:
            continue

        command = npt[0].lower()

        if command in {'dir', 'directory'}:
            if len(npt) < 2:
                continue

            out_dir = " ".join(npt[1:])

            if add_on is not None:
                out_dir = out_dir + add_on

            os.makedirs(out_dir, exist_ok=True)

            my_sys.dir = out_dir
            my_sys.files['dir'] = out_dir

            print(f"Export directory set to: {out_dir}")

        elif command in FORMAT_COMMANDS:
            if len(npt) < 2:
                raise ValueError(
                    "A geometry format is required after the format command. "
                    "Example: -e ft ply"
                )
            _set_mesh_format(my_sys, npt[1])

        elif command in {'no_archive', 'no-archive', 'no_network', 'no-network'}:
            skip_archive = True

        elif command in ONLY_COMMANDS:
            skip_archive = True
            if len(npt) < 2:
                raise ValueError(
                    "At least one export type is required after 'only'. "
                    "Example: -e only logs"
                )
            export_commands.extend([[export_name] for export_name in npt[1:]])

        else:
            export_commands.append(npt)

    # The interface CLI has a deliberately small default bundle.  Legacy
    # non-interface commands retain the historical Large preset.
    if len(export_commands) == 0:
        if getattr(my_sys, '_compact_interface_workflow', False):
            from vorpy.src.output.visualization import export_compact_visualization_bundle
            export_compact_visualization_bundle(my_sys)
            return
        export_commands.append(['large'])

    skip_archive = skip_archive or any(npt[0].lower() in {'none', 'no', 'skip'} for npt in export_commands)
    # Keep archive policy scoped to this invocation, including repeated presets.
    previous = {name: getattr(my_sys, name) for name in
                ('_export_skip_archive', '_export_archive_written') if hasattr(my_sys, name)}
    my_sys._export_skip_archive = skip_archive
    my_sys._export_archive_written = False
    try:
        for npt in export_commands:
            export_npt(my_sys, npt[0])
    finally:
        for name in ('_export_skip_archive', '_export_archive_written'):
            if name in previous:
                setattr(my_sys, name, previous[name])
            else:
                delattr(my_sys, name)


def export_npt(my_sys, usr_npt=None):
    """
    Handles the export of system data based on user-specified export type.

    Parameters:
    -----------
    my_sys : System
        The system object containing the data to be exported
    usr_npt : str, optional
        The export type specification. If None or 'default', performs a large export.
        Valid options include:
        - 'default': Large export (default)
        - '2'/'medium'/'med': Medium export
        - 'tiny'/'i'/'info'/'0'/'smallest': Small export
        - 'small'/'s'/'1': Medium-small export
        - 'large'/'l'/'3': Large export
        - 'all'/'a'/'everything': Full export
        - Other custom export types

    Returns:
    --------
    None
        The function performs exports but does not return any values.
    """

    # print("\n=== EXPORT DEBUG ===")
    # print(f"export type      = {usr_npt}")
    # print(f"system name      = {my_sys.name}")
    # print(f"system dir       = {getattr(my_sys, 'dir', None)}")
    # print(f"files['dir']     = {my_sys.files.get('dir', None)}")
    # print(f"round_to         = {getattr(my_sys, 'round_to', None)}")
    # print(f"groups           = {len(getattr(my_sys, 'groups', []))}")
    #
    # for i, grp in enumerate(getattr(my_sys, 'groups', [])):
    #     print(f"\nGROUP {i}")
    #     print(f"  name           = {grp.name}")
    #     print(f"  dir            = {getattr(grp, 'dir', None)}")
    #
    #     if hasattr(grp, 'net'):
    #         net = getattr(grp, "net", None)
    #
    #         print(f"  net            = {net}")
    #
    #         if net is None:
    #             print("  verts          = None")
    #             print("  edges          = None")
    #             print("  surfs          = None")
    #         else:
    #             print(
    #                 f"  verts          = "
    #                 f"{None if net.verts is None else len(net.verts)}"
    #             )
    #             print(
    #                 f"  edges          = "
    #                 f"{None if net.edges is None else len(net.edges)}"
    #             )
    #             print(
    #                 f"  surfs          = "
    #                 f"{None if net.surfs is None else len(net.surfs)}"
    #             )
    #
    # print("====================\n")

    # If nothing is specified export the large default.
    if usr_npt is None or usr_npt.lower() in {'default', ''}:

        export_large(my_sys)

    elif usr_npt.lower() in {'2', 'medium', 'med'}:

        # print("RUNNING export_med()\n")

        export_med(sys=my_sys)

    elif usr_npt.lower() in {"tiny", "i", "info", "0", "smallest"}:

        # print("RUNNING export_micro()\n")

        export_micro(my_sys)

    elif usr_npt.lower() in {"small", "s", "1"}:

        # print("RUNNING export_tiny()\n")

        export_tiny(my_sys)

    elif usr_npt.lower() in {"large", "l", "3"}:

        # print("RUNNING export_large()\n")

        export_large(my_sys)

    elif usr_npt.lower() in {'all', 'a', 'everything'}:

        # print("RUNNING export_all()\n")

        export_all(my_sys)

    else:

        # print(f"RUNNING other_exports({usr_npt})\n")

        other_exports(my_sys, usr_npt)
