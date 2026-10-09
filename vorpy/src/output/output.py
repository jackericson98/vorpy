"""Data-driven output presets and standalone export commands for VorPy."""

import json
import os
import shutil
import time
from contextlib import contextmanager
from pathlib import Path

from vorpy.src.output.atoms import write_atom_cells
from vorpy.src.output.curvature_colors import export_color_cache
from vorpy.src.analyze.nonpolar_interface import export_geometry_tables, build_geometry_tables
from vorpy.src.output.layout import unique_output_directory


SYSTEM_PRESETS = {
    'micro': [('info', {'info': True})],
    'tiny': [('info', {'info': True})],
    'small': [('info', {'info': True})],
    'medium': [
        ('PDB', {'pdb': True}),
        ('PyMOL atoms', {'set_atoms': True}),
        ('info', {'info': True}),
    ],
    'large': [
        ('PDB', {'pdb': True}),
        ('PyMOL atoms', {'set_atoms': True}),
        ('info', {'info': True}),
    ],
    'all': [
        ('PDB', {'pdb': True}),
        ('PyMOL atoms', {'set_atoms': True}),
        ('info', {'info': True}),
    ],
}

GROUP_PRESETS = {
    'micro': [('info', {'info': True})],
    'tiny': [('info', {'info': True})],
    'small': [
        ('info', {'info': True}),
        ('surfaces', {'surfs': True}),
        ('logs', {'logs': True}),
    ],
    'medium': [
        ('info', {'info': True}),
        ('shell surfaces', {'shell_surfs': True}),
        ('surfaces', {'surfs': True}),
        ('shell edges', {'shell_edges': True}),
        ('edges', {'edges': True}),
        ('shell vertices', {'shell_verts': True}),
        ('vertices', {'verts': True}),
        ('logs', {'logs': True}),
        ('atoms', {'atoms': True}),
        ('surrounding atoms', {'surr_atoms': True}),
        ('surrounding residues', {'surr_resids': True}),
    ],
    'large': [
        ('shell vertices', {'shell_verts': True}),
        ('shell edges', {'shell_edges': True}),
        ('shell surfaces', {'shell_surfs': True}),
        ('surfaces', {'surfs': True}),
        ('info', {'info': True}),
        ('edges', {'edges': True}),
        ('vertices', {'verts': True}),
        ('atoms', {'atoms': True}),
        ('surrounding atoms', {'surr_atoms': True}),
        ('surrounding residues', {'surr_resids': True}),
        ('logs', {'logs': True}),
        ('atom surfaces', {'atom_surfs': True}),
        ('atom edges', {'atom_edges': True}),
        ('atom vertices', {'atom_verts': True}),
    ],
    'all': [
        ('atoms', {'atoms': True}),
        ('shell surfaces', {'shell_surfs': True}),
        ('surfaces', {'surfs': True}),
        ('separate surfaces', {'sep_surfs': True}),
        ('shell edges', {'shell_edges': True}),
        ('edges', {'edges': True}),
        ('separate edges', {'sep_edges': True}),
        ('shell vertices', {'shell_verts': True}),
        ('vertices', {'verts': True}),
        ('separate vertices', {'sep_verts': True}),
        ('surrounding atoms', {'surr_atoms': True}),
        ('surrounding residues', {'surr_resids': True}),
        ('external atoms', {'ext_atoms': True}),
        ('logs', {'logs': True}),
        ('info', {'info': True}),
        ('atom surfaces', {'atom_surfs': True}),
        ('atom edges', {'atom_edges': True}),
        ('atom vertices', {'atom_verts': True}),
    ],
}

INTERFACE_PRESETS = {
    'micro': [('info', {'info': True, 'buried_water': False})],
    'tiny': [('info', {'info': True, 'buried_water': False})],
    'small': [
        ('surfaces', {'surfs': True}),
        ('info', {'info': True, 'buried_water': False}),
    ],
    'medium': [
        ('surfaces', {'surfs': True}),
        ('atoms', {'atoms': True}),
        ('edges', {'edges': True}),
        ('logs', {'logs': True}),
        ('vertices', {'verts': True}),
        ('info', {'info': True}),
    ],
    'large': [
        ('atoms', {'atoms': True}),
        ('surfaces', {'surfs': True}),
        ('edges', {'edges': True}),
        ('vertices', {'verts': True}),
        ('info', {'info': True}),
    ],
    'all': [
        ('surfaces', {'surfs': True}),
        ('atoms', {'atoms': True}),
        ('edges', {'edges': True}),
        ('logs', {'logs': True}),
        ('vertices', {'verts': True}),
        ('info', {'info': True}),
    ],
}

GROUP_EXPORT_OPTIONS = {
    'group_atoms': {'atoms': True},
    'atom_surfs': {'atom_surfs': True},
    'atom_edges': {'atom_edges': True},
    'atom_verts': {'atom_verts': True},
    'surfs': {'surfs': True},
    'surfaces': {'surfs': True},
    'sep_surfs': {'sep_surfs': True},
    'shell': {'shell_surfs': True},
    'shell_surfs': {'shell_surfs': True},
    'edges': {'edges': True},
    'sep_edges': {'sep_edges': True},
    'shell_edges': {'shell_edges': True},
    'verts': {'verts': True},
    'vertices': {'verts': True},
    'sep_verts': {'sep_verts': True},
    'shell_verts': {'shell_verts': True},
    'surr_atoms': {'surr_atoms': True},
    'surrounding_atoms': {'surr_atoms': True},
    'surr_resids': {'surr_resids': True},
    'surr_residues': {'surr_resids': True},
    'surrounding_resids': {'surr_resids': True},
    'surrounding_residues': {'surr_resids': True},
    'ext_atoms': {'ext_atoms': True},
    'logs': {'logs': True},
    'lgs': {'logs': True},
}

SYSTEM_EXPORT_OPTIONS = {
    'pdb': {'pdb': True},
    'set_atoms': {'set_atoms': True},
    'mol': {'mol': True},
    'cif': {'cif': True},
    'xyz': {'xyz': True},
    'txt': {'txt': True},
}


class ExportProgress:
    """Report export progress and retain per-operation timing."""

    def __init__(self, total, sys):
        self.total = max(int(total), 1)
        self.current = 0
        self.sys = sys
        self.sys.run_network = None
        self.start = time.perf_counter()
        self.timings = {}
        self.counts = {}
        self.mesh_timings = {}
        self.verbose = bool(getattr(sys, 'verbose', False) or
                            (getattr(sys, 'settings', None) or {}).get('verbose', False) or
                            any((getattr(group, 'settings', None) or {}).get('verbose', False)
                                for group in (getattr(sys, 'groups', None) or [])))

    def show(self, name=None):
        percent = 100.0 * self.current / self.total
        process = 'Export' if name is None else f'Export: {name}'
        self.sys.update_progress(process=process, progress=percent)

    def step(self):
        self.current = min(self.current + 1, self.total)

    def finish(self):
        self.current = self.total
        self.sys.update_progress(process='Export', progress=100.0)
        total_elapsed = time.perf_counter() - self.start
        self.sys.export_timing = self.timings.copy()
        self.sys.export_timing['total'] = total_elapsed
        self.sys.export_mesh_timing = self.mesh_timings.copy()
        debug = str(os.environ.get('VORPY_EXPORT_TIMING', '0')).strip().lower()
        if not self.verbose and debug not in {'1', 'true', 'yes', 'on'}:
            return
        print('\n' + '=' * 70)
        print('EXPORT TIMING')
        print('=' * 70)
        for name, elapsed in sorted(self.timings.items(), key=lambda item: item[1], reverse=True):
            pct = 100.0 * elapsed / total_elapsed if total_elapsed else 0.0
            print(f'{name:<40} {elapsed:10.4f} s  {pct:6.2f} %')
        measured = sum(self.timings.values())
        other = max(total_elapsed - measured, 0.0)
        pct = 100.0 * other / total_elapsed if total_elapsed else 0.0
        print(f'{"Other / export overhead":<40} {other:10.4f} s  {pct:6.2f} %')
        print('-' * 70)
        print(f'{"TOTAL":<40} {total_elapsed:10.4f} s  100.00 %')
        print(f'Export operations: {sum(self.counts.values()):,}')
        print(f'Mesh preparation: {sum(item["prepare"] for item in self.mesh_timings.values()):.4f} s')
        print(f'Mesh file writing: {sum(item["write"] for item in self.mesh_timings.values()):.4f} s')
        print('=' * 70)


def _run_export(progress, name, func, **kwargs):
    progress.show(name)
    start = time.perf_counter()
    owner = getattr(func, '__self__', None)
    net = getattr(owner, 'net', None)
    timing = getattr(net, '_export_mesh_timing', None)
    before = dict(timing) if timing is not None else None
    func(**kwargs)
    elapsed = time.perf_counter() - start
    progress.timings[name] = progress.timings.get(name, 0.0) + elapsed
    progress.counts[name] = progress.counts.get(name, 0) + 1
    details = ''
    if before is not None:
        phases = {phase: timing[phase] - before[phase] for phase in ('prepare', 'write')}
        progress.mesh_timings[name] = phases
        details = f' | mesh preparation={phases["prepare"]:.2f}s; file writing={phases["write"]:.2f}s'
    if progress.verbose:
        print(f'\n[export] {name}: {elapsed:.2f}s{details}', flush=True)
    progress.step()


def _set_group_directory(sys, group):
    root = Path(sys.files['dir'])
    group_indices = {int(index) for index in (getattr(group, 'ball_ndxs', ()) or ())}
    balls = getattr(sys, 'balls', None)
    total_atoms = len(balls) if balls is not None else 0
    if total_atoms and group_indices == set(range(total_atoms)):
        group.dir = str(root)
        root.mkdir(parents=True, exist_ok=True)
        return

    existing = getattr(group, 'dir', None)
    if existing and Path(existing).parent == root:
        group_dir = Path(existing)
    else:
        occupied = [getattr(item, 'dir', None) for item in getattr(sys, 'groups', ()) or ()
                    if getattr(item, 'dir', None)]
        group_dir = unique_output_directory(root, getattr(group, 'name', 'network'), occupied)
    group.dir = str(group_dir)
    group_dir.mkdir(parents=True, exist_ok=True)


def _move_vert_file(sys, group):
    vert_file = group.settings['net_type'] + '_verts.txt'
    source = os.path.join(sys.files['dir'], vert_file)
    destination = os.path.join(group.dir, vert_file)
    if os.path.exists(source) and not os.path.exists(destination):
        shutil.move(source, destination)


def _export_nonpolar_geometry(sys, network, directory, network_id):
    """Export the read-only local geometry analysis tables for one network."""
    destination = os.path.join(directory, 'nonpolar_interface_geometry')
    tables = build_geometry_tables(
        network,
        system_id=getattr(sys, 'name', ''),
        frame_id=getattr(sys, 'frame_index', None),
        network_id=network_id,
    )
    export_geometry_tables(tables, destination)


def _export_dual_geometry_without_launcher(iface):
    """Keep rich dual derivatives, but do not retain preset PyMOL launchers."""
    iface.export(dual=True)
    dual_dir = getattr(iface, 'dir', None)
    if dual_dir is None:
        return
    dual_dir = Path(dual_dir) / 'dual'
    for name in ('apollonius_interface.pml', 'apollonius_interface.py'):
        (dual_dir / name).unlink(missing_ok=True)
    summary = dual_dir / 'interface_dual_summary.json'
    if summary.exists():
        payload = json.loads(summary.read_text(encoding='utf-8'))
        payload['pymol_script'] = None
        summary.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')


@contextmanager
def _group_export_cache(group):
    """Share summaries and color scales only within this export plan."""
    previous = getattr(group, "_export_info_ready", None)
    previous_indices = getattr(group, '_export_topology_indices', None)
    group._export_info_ready = False
    group._export_topology_indices = None
    try:
        with export_color_cache(group.net):
            yield
    finally:
        if previous is None:
            del group._export_info_ready
        else:
            group._export_info_ready = previous
        if previous_indices is None:
            del group._export_topology_indices
        else:
            group._export_topology_indices = previous_indices


def export_preset(sys, preset):
    """Execute one named export plan."""
    groups = [group for group in sys.groups if group.net is not None]
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)
    system_plan = SYSTEM_PRESETS[preset]
    group_plan = GROUP_PRESETS[preset]
    interface_plan = INTERFACE_PRESETS[preset]
    analysis_count = (len(groups) + len(ifaces)) if preset in {'large', 'all'} else 0
    # ``logs`` still includes the legacy sectioned interface CSV.  Whenever it
    # is requested, also emit the canonical Results log: it preserves status,
    # units, identity, and provenance without recalculating geometry.  ``small``
    # has always supplied that canonical log even though its interface plan has
    # no legacy ``logs`` entry.
    write_canonical_interface_logs = (
        preset == 'small'
        or any(kwargs.get('logs', False) for _, kwargs in interface_plan)
    )
    canonical_interface_logs = (
        [iface for iface in ifaces if iface.net is not None]
        if write_canonical_interface_logs else []
    )
    dual_interfaces = [iface for iface in ifaces if iface.net is not None and
                       iface.net.settings.get('net_type', 'aw') == 'aw'] if preset in {'large', 'all'} else []
    archive_enabled = (not getattr(sys, '_export_skip_archive', False)
                       and not getattr(sys, '_export_archive_written', False))
    total = (len(system_plan) + len(group_plan) * len(groups)
             + len(interface_plan) * len(ifaces) + len(canonical_interface_logs)
             + analysis_count + len(dual_interfaces) + int(archive_enabled))
    progress = ExportProgress(total, sys)

    for name, kwargs in system_plan:
        _run_export(progress, f'system {name}', sys.exports, **kwargs)
    for group in groups:
        _set_group_directory(sys, group)
        with _group_export_cache(group):
            for name, kwargs in group_plan:
                _run_export(progress, f'{group.name}: {name}', group.exports, **kwargs)
        # The vertex log is a network-level artifact. Move it once after the
        # group's export plan, rather than checking it from every individual
        # export stage.
        _move_vert_file(sys, group)
        if preset in {'large', 'all'}:
            _run_export(
                progress,
                f'{group.name}: nonpolar interface geometry',
                _export_nonpolar_geometry,
                sys=sys,
                network=group.net,
                directory=group.dir,
                network_id=group.name,
            )
    for iface in ifaces:
        interface_name = getattr(iface, 'name', 'interface')
        with export_color_cache(iface.net):
            for name, kwargs in interface_plan:
                _run_export(progress, f'{interface_name}: {name}', iface.export, **kwargs)
            if iface in canonical_interface_logs:
                _run_export(
                    progress,
                    f'{interface_name}: canonical logs',
                    _export_canonical_interface_log,
                    iface=iface,
                )
            if iface in dual_interfaces:
                _run_export(progress, f'{interface_name}: Apollonius dual',
                            _export_dual_geometry_without_launcher, iface=iface)
            if preset in {'large', 'all'}:
                _run_export(
                    progress,
                    f'{interface_name}: nonpolar interface geometry',
                    _export_nonpolar_geometry,
                    sys=sys,
                    network=iface.net,
                    directory=iface.dir,
                    network_id=interface_name,
                )
    if archive_enabled:
        _run_export(progress, 'network archive', _export_network_archive, sys=sys)
    progress.finish()


def _export_network_archive(sys):
    from vorpy.src.io import save_network
    destination = os.path.join(sys.files['dir'], str(sys.name) + '.vpy')
    save_network(sys, destination)
    if hasattr(sys, '_export_archive_written'):
        sys._export_archive_written = True


def _export_canonical_interface_log(iface):
    """Persist the canonical Results log without rebuilding cached geometry."""
    from vorpy.src.output.visualization import _compact_log

    _compact_log(iface, Path(iface.dir) / 'logs.csv')


def export_micro(sys):
    export_preset(sys, 'micro')


def export_tiny(sys):
    # ``small`` is the historical CLI-small alias; keep its behavior while
    # exposing ``export_preset(..., 'tiny')`` as the metadata-only preset.
    export_preset(sys, 'small')


def export_small(sys):
    export_preset(sys, 'small')


def export_med(sys):
    export_preset(sys, 'medium')


def export_large(sys):
    export_preset(sys, 'large')


def export_all(sys):
    export_preset(sys, 'all')


def other_exports(sys, usr_npt):
    """Run one standalone export without adding preset or system side effects."""
    option = str(usr_npt).strip().lower()
    groups = [group for group in sys.groups if group.net is not None]

    if option in {'none', 'no', 'skip'}:
        return

    if option in {'dual', 'apollonius'}:
        interfaces = [iface for iface in (getattr(sys, 'ifaces', None) or []) if iface.net is not None]
        if not interfaces:
            raise ValueError('-e dual requires a solved AW Interface; use the two-group -i workflow.')
        progress = ExportProgress(len(interfaces), sys)
        for iface in interfaces:
            _run_export(progress, f'{iface.name}: Apollonius dual', iface.export, dual=True)
        progress.finish()
        return

    if option == 'visualize':
        from vorpy.src.output.visualization import export_visualization_bundle
        if getattr(sys, '_compact_interface_workflow', False):
            from vorpy.src.output.visualization import export_compact_visualization_bundle
            export_visualization_bundle = export_compact_visualization_bundle

        progress = ExportProgress(1, sys)
        result = {}

        def write_visualization():
            result['items'] = export_visualization_bundle(sys)

        _run_export(progress, 'visualization bundle', write_visualization)
        progress.finish()
        for item in result.get('items', ()):
            if item.get('launcher'):
                print(f"Visualization launcher written: {item['launcher']}")
                for stage, elapsed in item.get('timings', {}).items():
                    print(f'  {stage}: {elapsed:.2f} s')
                for name, layer in item.get('layers', {}).items():
                    state = layer.get('state', 'NOT_CALCULATED')
                    reason = f" ({layer['reason']})" if layer.get('reason') else ''
                    print(f'  {name}: {state}{reason}')
            else:
                print('Visualization: no solved interface is available; nothing exported.')
        return

    if option in {'water', 'waters', 'interface_water'}:
        interfaces = [iface for iface in (getattr(sys, 'ifaces', None) or [])
                      if getattr(iface, 'net', None) is not None]
        for iface in interfaces:
            state = 'AVAILABLE' if getattr(iface, 'water_topology', None) is not None else 'NOT_CALCULATED'
            print(f'Interface water analysis {iface.name}: {state}')
        return

    if option in {'a', 'atoms', 'atom_cells'}:
        progress = ExportProgress(1, sys)
        _run_export(
            progress, 'atom cells', write_atom_cells,
            net=sys.net, atoms=list(range(len(sys.net.balls))),
            directory=sys.files['dir'], file_type=getattr(sys, 'file_type', 'off'),
        )
        progress.finish()
        return

    if option in SYSTEM_EXPORT_OPTIONS:
        progress = ExportProgress(1, sys)
        _run_export(progress, f'system {option}', sys.exports, **SYSTEM_EXPORT_OPTIONS[option])
        progress.finish()
        return

    kwargs = GROUP_EXPORT_OPTIONS.get(option)
    if kwargs is None:
        raise ValueError(f'Unknown export type: {usr_npt!r}')

    progress = ExportProgress(len(groups), sys)
    for group in groups:
        _set_group_directory(sys, group)
        _run_export(progress, f'{group.name}: {option}', group.exports, **kwargs)
    progress.finish()
