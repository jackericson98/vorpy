"""Data-driven output presets and standalone export commands for VorPy."""

import os
import shutil
import time

from vorpy.src.output.atoms import write_atom_cells


SYSTEM_PRESETS = {
    'micro': [('info', {'info': True})],
    'tiny': [
        ('info', {'info': True}),
        ('PDB', {'pdb': True}),
        ('PyMOL atoms', {'set_atoms': True}),
    ],
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
    'tiny': [
        ('info', {'info': True}),
        ('shell surfaces', {'shell_surfs': True}),
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
    ],
    'large': [
        ('shell vertices', {'shell_verts': True}),
        ('shell edges', {'shell_edges': True}),
        ('shell surfaces', {'shell_surfs': True}),
        ('info', {'info': True}),
        ('edges', {'edges': True}),
        ('vertices', {'verts': True}),
        ('atoms', {'atoms': True}),
        ('surrounding atoms', {'surr_atoms': True}),
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
        ('external atoms', {'ext_atoms': True}),
        ('logs', {'logs': True}),
        ('info', {'info': True}),
        ('atom surfaces', {'atom_surfs': True}),
        ('atom edges', {'atom_edges': True}),
        ('atom vertices', {'atom_verts': True}),
    ],
}

INTERFACE_PRESETS = {
    'micro': [('info', {'info': True})],
    'tiny': [('info', {'info': True})],
    'medium': [
        ('surfaces', {'surfs': True}),
        ('atoms', {'atoms': True}),
        ('edges', {'edges': True}),
        ('logs', {'logs': True}),
        ('vertices', {'verts': True}),
        ('info', {'info': True}),
    ],
    'large': [
        ('balls', {'balls': True}),
        ('surfaces', {'surfs': True}),
        ('edges', {'edges': True}),
        ('vertices', {'verts': True}),
        ('info', {'info': True}),
    ],
    'all': [
        ('balls', {'balls': True}),
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

    def show(self, name=None):
        percent = 100.0 * self.current / self.total
        process = 'Exporting files' if name is None else f'Exporting files: {name}'
        self.sys.update_progress(process=process, progress=percent)

    def step(self):
        self.current = min(self.current + 1, self.total)

    def finish(self):
        self.current = self.total
        self.sys.update_progress(process='Exporting files', progress=100.0)
        total_elapsed = time.perf_counter() - self.start
        self.sys.export_timing = self.timings.copy()
        self.sys.export_timing['total'] = total_elapsed
        debug = str(os.environ.get('VORPY_EXPORT_TIMING', '0')).strip().lower()
        if debug not in {'1', 'true', 'yes', 'on'}:
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
        print('=' * 70)


def _run_export(progress, name, func, **kwargs):
    progress.show(name)
    start = time.perf_counter()
    func(**kwargs)
    elapsed = time.perf_counter() - start
    progress.timings[name] = progress.timings.get(name, 0.0) + elapsed
    progress.counts[name] = progress.counts.get(name, 0) + 1
    progress.step()


def _set_group_directory(sys, group):
    group_dir = os.path.join(sys.files['dir'], group.name)
    group.dir = group_dir
    os.makedirs(group.dir, exist_ok=True)


def _move_vert_file(sys, group):
    vert_file = group.settings['net_type'] + '_verts.txt'
    source = os.path.join(sys.files['dir'], vert_file)
    destination = os.path.join(group.dir, vert_file)
    if os.path.exists(source) and not os.path.exists(destination):
        shutil.move(source, destination)


def export_preset(sys, preset):
    """Execute one named export plan."""
    groups = [group for group in sys.groups if group.net is not None]
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)
    system_plan = SYSTEM_PRESETS[preset]
    group_plan = GROUP_PRESETS[preset]
    interface_plan = INTERFACE_PRESETS[preset]
    total = len(system_plan) + len(group_plan) * len(groups) + len(interface_plan) * len(ifaces)
    progress = ExportProgress(total, sys)

    for name, kwargs in system_plan:
        _run_export(progress, f'system {name}', sys.exports, **kwargs)
    for group in groups:
        _set_group_directory(sys, group)
        for name, kwargs in group_plan:
            _run_export(progress, f'{group.name}: {name}', group.exports, **kwargs)
        if preset in {'large', 'all'}:
            _move_vert_file(sys, group)
    for iface in ifaces:
        interface_name = getattr(iface, 'name', 'interface')
        for name, kwargs in interface_plan:
            _run_export(progress, f'{interface_name}: {name}', iface.export, **kwargs)
    progress.finish()


def export_micro(sys):
    export_preset(sys, 'micro')


def export_tiny(sys):
    export_preset(sys, 'tiny')


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