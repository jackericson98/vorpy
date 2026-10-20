import os
from vorpy.src.command.calculation_plan import build_export_plan
from vorpy.src.output import export_micro
from vorpy.src.output import export_tiny
from vorpy.src.output import export_med
from vorpy.src.output import export_large
from vorpy.src.output import export_all
from vorpy.src.output import other_exports


FORMAT_COMMANDS = {'ft', 'file_type', 'file_format', 'format', 'mesh_format'}
FORMAT_VALUES = {'off': 'off', 'ply': 'ply', 'vtp': 'vtp', 'vtk': 'vtp'}
ONLY_COMMANDS = {'only', 'just'}


def _export_scope(system):
    """Infer a scope only for callers that bypass Command.parse_commands."""
    if getattr(system, '_compact_interface_workflow', False) or getattr(system, 'ifaces', None):
        return 'interface'
    if getattr(system, 'groups', None):
        return 'group'
    return 'system'


def _run_component_plan(system, plan):
    """Dispatch the normalized plan through the existing writer entry points."""
    if not plan.included_components:
        return
    from vorpy.src.output import output

    # Interface runs reserve their canonical directories in the CLI after
    # collision validation.  Direct callers still need the normal root setup.
    if (
        not getattr(system, '_cli_interface_output_prepared', False)
        and getattr(system, 'balls', None) is not None
    ):
        output.prepare_canonical_root(system)

    group_options = {
        'edges': {'edges': True},
        'surfaces': {'surfs': True},
        'vertices': {'verts': True},
        'shell': {'shell_surfs': True},
        'shell_edges': {'shell_edges': True},
        'shell_vertices': {'shell_verts': True},
        'group_atoms': {'atoms': True},
        'surrounding_atoms': {'surr_atoms': True},
        'surrounding_residues': {'surr_resids': True},
        'external_atoms': {'ext_atoms': True},
        'separate_surfaces': {'sep_surfs': True},
        'separate_edges': {'sep_edges': True},
        'separate_vertices': {'sep_verts': True},
        'atom_surfaces': {'atom_surfs': True},
        'atom_edges': {'atom_edges': True},
        'atom_vertices': {'atom_verts': True},
        'logs': {'logs': True},
        'info': {'info': True},
    }
    interface_options = {
        'edges': {'edges': True},
        'surfaces': {'surfs': True},
        'vertices': {'verts': True},
        'interface_atoms': {'interface_atoms': True},
        'surrounding_atoms': {'atoms': True},
        'surrounding_residues': {'atoms': True},
        'contact_tables': {'contact_tables': True},
        'logs': {'logs': True},
        'info': {'info': True},
        'results': {},
        'provenance': {},
        'legacy_results': {},
    }
    system_options = {
        'info': {'info': True},
        'source_structure': {'pdb': True},
        'pymol_atoms': {'set_atoms': True},
    }
    legacy_components = {
        'dual': 'dual',
        'visualization': 'visualize',
        'water': 'water',
        'mol': 'mol',
        'cif': 'cif',
        'xyz': 'xyz',
        'txt': 'txt',
    }

    groups = [
        group for group in (getattr(system, 'groups', None) or [])
        if getattr(group, 'net', None) is not None
        or getattr(group, 'source_view', None) is not None
    ]
    ifaces = [
        iface for iface in (getattr(system, 'ifaces', None) or [])
        if getattr(iface, 'net', None) is not None
    ]
    components = plan.included_components
    unsupported_interface = {
        component: 'NOT_SUPPORTED'
        for component in components
        if component in {'curvature', 'topology', 'shell', 'shell_edges', 'shell_vertices'}
    }
    for iface in ifaces:
        iface.export_component_status = dict(unsupported_interface)
    system_calls = [
        (component, system_options[component])
        for component in components
        if component in system_options and callable(getattr(system, 'exports', None))
    ]
    group_call_count = sum(
        sum(component in group_options for component in components)
        for group in groups if callable(getattr(group, 'exports', None))
    )
    interface_call_count = sum(
        sum(component in interface_options for component in components)
        for iface in ifaces if callable(getattr(iface, 'export', None))
    )
    checkpoint = (
        plan.checkpoint_policy in {'shared', 'required'}
        or (
            plan.checkpoint_policy == 'explicit-only'
            and 'checkpoint' in components
        )
    ) and not getattr(system, '_export_skip_archive', False)
    total = len(system_calls) + group_call_count + interface_call_count + int(checkpoint)
    progress = output.ExportProgress(total, system)

    for component in components:
        legacy = legacy_components.get(component)
        if legacy is not None:
            other_exports(system, legacy)

    for component, kwargs in system_calls:
        output._run_export(progress, f'system {component}', system.exports, **kwargs)

    for group in groups:
        calls = [
            (component, group_options[component])
            for component in components if component in group_options
        ]
        if not calls or not callable(getattr(group, 'exports', None)):
            continue
        output._set_group_directory(system, group)
        with output._group_export_cache(group):
            for component, kwargs in calls:
                output._run_export(
                    progress, f'{group.name}: {component}', group.exports, **kwargs
                )

    for iface in ifaces:
        calls = [
            (component, interface_options[component])
            for component in components if component in interface_options
        ]
        if not calls or not callable(getattr(iface, 'export', None)):
            continue
        if getattr(iface, 'dir', None) is None and hasattr(iface, 'set_dir'):
            iface.set_dir()
        with output.export_color_cache(iface.net):
            for component, kwargs in calls:
                if component == 'results':
                    func, kwargs = output._export_canonical_interface_log, {'iface': iface}
                elif component == 'provenance':
                    func, kwargs = output._export_interface_provenance, {'iface': iface}
                elif component == 'legacy_results':
                    func, kwargs = output._export_legacy_interface_results, {'iface': iface}
                else:
                    func = iface.export
                output._run_export(
                    progress, f'{iface.name}: {component}', func, **kwargs
                )

    if checkpoint:
        output._run_export(
            progress, 'network archive', output._export_network_archive, sys=system
        )
    progress.finish()


def _plan_dispatches_components(plan):
    """Whether the normalized component set, rather than legacy names, wins."""
    return any(
        operation and operation[0] in {'preset', 'add', 'no', 'remove', 'only', 'none'}
        for operation in plan.operations
    )


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
    - If no export type is specified, interface workflows use the first-class
      interface preset; other workflows use the large preset.
    - Format and directory commands modify the export but do not replace its preset.
    - If only modifiers are specified, export the large preset.
    - ``only`` explicitly restricts output to the export names that follow it.
    - If explicit exports are specified, only export those.
    """

    if usr_npt is None:
        usr_npt = []

    plan = getattr(my_sys, '_export_plan', None)
    if plan is None:
        plan = build_export_plan(usr_npt, scope=_export_scope(my_sys))
        my_sys._export_plan = plan

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

            if getattr(my_sys, 'verbose', False):
                print(f"Export directory set to: {out_dir}")

        elif command in FORMAT_COMMANDS:
            if len(npt) < 2:
                raise ValueError(
                    "A geometry format is required after the format command. "
                    "Example: -e ft ply"
                )
            _set_mesh_format(my_sys, npt[1])

        elif command in {
            'checkpoint', 'no_checkpoint', 'no_archive', 'no-archive',
            'no_network', 'no-network',
        }:
            if command in {'no_checkpoint', 'no_archive', 'no-archive', 'no_network', 'no-network'}:
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

    # Interface results have a first-class default bundle. Legacy non-interface
    # commands retain the historical Large preset.
    if len(export_commands) == 0:
        if getattr(my_sys, '_compact_interface_workflow', False):
            from vorpy.src.output.output import export_interface_default
            export_interface_default(my_sys)
            return
        export_commands.append(['large'])

    skip_archive = (
        skip_archive
        or plan.checkpoint_policy == 'none'
        or any(npt[0].lower() in {'none', 'no', 'skip'} for npt in export_commands)
    )
    # Keep archive policy scoped to this invocation, including repeated presets.
    previous = {name: getattr(my_sys, name) for name in
                ('_export_skip_archive', '_export_archive_written') if hasattr(my_sys, name)}
    my_sys._export_skip_archive = skip_archive
    my_sys._export_archive_written = False
    normalized_dispatch = _plan_dispatches_components(plan)
    legacy_commands = [
        operation[1]
        for operation in plan.operations
        if operation and operation[0] == 'legacy'
    ]
    try:
        if normalized_dispatch:
            # Keep explicit legacy/API exports additive, but route every
            # preset and component selection through the canonical plan.
            for command in legacy_commands:
                export_npt(my_sys, command)
            _run_component_plan(my_sys, plan)
        else:
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
        - 'interface': First-class physical interface export bundle
        - 'nonpolar': Opt-in nonpolar interface geometry tables
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

    elif usr_npt.lower() in {'interface', 'interfaces'}:

        # Keep the explicit preset equivalent to ordinary ``-i`` defaults.
        from vorpy.src.output.output import export_interface_default
        export_interface_default(my_sys)

    else:

        # print(f"RUNNING other_exports({usr_npt})\n")

        other_exports(my_sys, usr_npt)
