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

        elif command in ONLY_COMMANDS:
            if len(npt) < 2:
                raise ValueError(
                    "At least one export type is required after 'only'. "
                    "Example: -e only logs"
                )
            export_commands.extend([[export_name] for export_name in npt[1:]])

        else:
            export_commands.append(npt)

    # No preset, or only modifiers such as ``dir``/``ft``: export Large.
    if len(export_commands) == 0:
        export_commands.append(['large'])

    # Second pass: run exports
    for npt in export_commands:
        if len(npt) == 0:
            continue

        export_npt(my_sys, npt[0])


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
