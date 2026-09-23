"""Export a live solved group using the legacy GUI's presets."""
import os
from pathlib import Path

PRESETS = {
    'Small': ('info', 'shell_surfs', 'logs'),
    'Medium': ('info', 'shell_surfs', 'surfs', 'shell_edges', 'edges', 'shell_verts',
               'verts', 'logs', 'atoms', 'surr_atoms', 'surr_resids'),
    'Large': ('shell_verts', 'shell_edges', 'shell_surfs', 'info', 'edges',
              'verts', 'atoms', 'surr_atoms', 'surr_resids', 'logs',
              'atom_surfs', 'atom_edges', 'atom_verts'),
}


def export_result(result, directory, size, custom, round_to, file_type,
                  atom_format=None, color_overrides=None):
    group = result.export_group
    if group is None or result.defaults_stale:
        raise ValueError('Solve the current structure before exporting.')
    if atom_format is not None and atom_format not in {'pdb', 'xyz'}:
        raise ValueError('Atom structure format must be PDB or XYZ')
    parent = Path(directory).expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)
    # A new directory keeps repeated exports from overwriting existing files.
    from tempfile import mkdtemp
    destination = Path(mkdtemp(prefix=Path(group.name).name + '_', dir=parent))
    cwd, old_dir, old_parent = Path.cwd(), group.dir, group.sys.files['dir']
    options = dict(custom) if size == 'Custom' else dict.fromkeys(PRESETS[size], True)
    atomic_options = {key: options.pop(key, False) for key in
                      ('atoms', 'surr_atoms', 'surr_resids', 'full_molecule')} if atom_format else {}
    # None preserves the legacy calling convention; the Workbench passes its
    # atom format and viewer colors explicitly.
    full_molecule = options.pop('full_molecule', False)
    network = getattr(group, 'net', None)
    missing = object()
    previous_colors = getattr(network, '_export_color_provider', missing)
    try:
        if color_overrides is not None and network is not None:
            from vorpy.workbench.services.export_colors import make_color_provider
            network._export_color_provider = make_color_provider(result, color_overrides)
        group.dir = str(destination)
        group.sys.files['dir'] = str(parent)
        group.exports(**options, round_to=round_to, file_type=file_type)
        if atom_format or full_molecule:
            from vorpy.workbench.services.atomic_export import export_atomic_selections
            if full_molecule:
                atomic_options['full_molecule'] = True
            export_atomic_selections(result, atomic_options, destination, atom_format or 'pdb')
    finally:
        if network is not None:
            if previous_colors is missing:
                if hasattr(network, '_export_color_provider'):
                    del network._export_color_provider
            else:
                network._export_color_provider = previous_colors
        os.chdir(cwd)
        group.dir = old_dir
        group.sys.files['dir'] = old_parent
    return destination
