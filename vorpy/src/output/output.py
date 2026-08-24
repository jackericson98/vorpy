import os
import time
import shutil
from vorpy.src.output import write_atom_cells


def _format_time(seconds):
    """Format elapsed time as H:MM:SS.ss."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = seconds % 60
    return f"{hours}:{minutes:02d}:{seconds:05.2f}"

def _format_time(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = seconds % 60
    return f"{hours}:{minutes}:{seconds:.2f}"


def _get_start_time(sys):
    """Return the earliest network START timestamp for the current run."""
    starts = []

    if getattr(sys, 'net', None) is not None:
        start = sys.net.metrics.get('START')
        if start is not None:
            starts.append(start)

    for group in sys.groups:
        if group.net is None:
            continue

        start = group.net.metrics.get('START')
        if start is not None:
            starts.append(start)

    if starts:
        return min(starts)

    return time.perf_counter()


class ExportProgress:
    """Report export progress through the parent System."""

    def __init__(self, total, sys):
        self.total = max(int(total), 1)
        self.current = 0
        self.sys = sys
        self.sys.run_network = None

    def show(self, name=None):
        percent = 100.0 * self.current / self.total
        process = "Exporting files" if name is None else f"Exporting files: {name}"
        self.sys.update_progress(process=process, progress=percent)

    def step(self):
        self.current = min(self.current + 1, self.total)

    def finish(self):
        self.current = self.total
        self.sys.update_progress(process="Exporting files", progress=100.0)


def _run_export(progress, name, func, **kwargs):
    """Run one export operation and update the progress display."""
    progress.show(name)
    func(**kwargs)
    progress.step()


def _set_group_directory(sys, group):
    """Ensure the group's output directory exists."""
    group_dir = os.path.join(sys.files['dir'], group.name)

    if group.dir is None or not os.path.exists(group_dir):
        group.dir = group_dir
        try:
            os.makedirs(group.dir, exist_ok=True)
        except FileNotFoundError:
            group.dir = os.path.join(sys.files['dir'], 'group')
            os.makedirs(group.dir, exist_ok=True)


def _move_vert_file(sys, group):
    """Move generated vertex text output into the group directory."""
    vert_file = group.settings['net_type'] + '_verts.txt'
    source = os.path.join(sys.files['dir'], vert_file)
    destination = os.path.join(group.dir, vert_file)

    if os.path.exists(source) and not os.path.exists(destination):
        shutil.move(source, destination)


def export_micro(sys):
    """
    Smallest export.

    System:
        - Information

    Groups:
        - Information

    Interfaces:
        - Information
    """

    groups = list(sys.groups)
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)

    progress = ExportProgress(1 + len(groups) + len(ifaces), sys)

    _run_export(progress, "system info", sys.exports, info=True)

    for group in groups:
        _set_group_directory(sys, group)
        _run_export(progress, f"{group.name}: info", group.exports, info=True)

    for iface in ifaces:
        name = getattr(iface, 'name', 'interface')
        _run_export(progress, f"{name}: info", iface.export, info=True)

    progress.finish()


def export_tiny(sys):
    """
    Small export.

    System:
        - Information
        - PDB
        - PyMOL atom script

    Groups:
        - Information
        - Shell surfaces
        - Logs

    Interfaces:
        - Information
    """

    groups = [group for group in sys.groups if group.net is not None]
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)

    # 3 system + 3/group + 1/interface
    progress = ExportProgress(3 + 3 * len(groups) + len(ifaces), sys)

    _run_export(progress, "system info", sys.exports, info=True)
    _run_export(progress, "system PDB", sys.exports, pdb=True)
    _run_export(progress, "system PyMOL atoms", sys.exports, set_atoms=True)

    for group in groups:
        _set_group_directory(sys, group)

        _run_export(progress, f"{group.name}: info", group.exports, info=True)
        _run_export(progress, f"{group.name}: shell surfaces", group.exports, shell_surfs=True)
        _run_export(progress, f"{group.name}: logs", group.exports, logs=True)

    for iface in ifaces:
        name = getattr(iface, 'name', 'interface')
        _run_export(progress, f"{name}: info", iface.export, info=True)

    progress.finish()


def export_med(sys):
    """
    Medium export.

    System:
        - PDB
        - PyMOL atom script
        - Information

    Groups:
        - Shell surfaces
        - Surfaces
        - Shell edges
        - Edges
        - Shell vertices
        - Vertices
        - Logs
        - Atoms
        - Surrounding atoms

    Interfaces:
        - Surfaces
        - Atoms
        - Edges
        - Logs
        - Vertices
        - Information
    """

    groups = [group for group in sys.groups if group.net is not None]
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)

    # 3 system + 9/group + 6/interface
    progress = ExportProgress(3 + 9 * len(groups) + 6 * len(ifaces), sys)

    _run_export(progress, "system PDB", sys.exports, pdb=True)
    _run_export(progress, "system PyMOL atoms", sys.exports, set_atoms=True)
    _run_export(progress, "system info", sys.exports, info=True)

    for group in groups:
        _set_group_directory(sys, group)

        _run_export(progress, f"{group.name}: shell surfaces", group.exports, shell_surfs=True)
        _run_export(progress, f"{group.name}: surfaces", group.exports, surfs=True)
        _run_export(progress, f"{group.name}: shell edges", group.exports, shell_edges=True)
        _run_export(progress, f"{group.name}: edges", group.exports, edges=True)
        _run_export(progress, f"{group.name}: shell vertices", group.exports, shell_verts=True)
        _run_export(progress, f"{group.name}: vertices", group.exports, verts=True)
        _run_export(progress, f"{group.name}: logs", group.exports, logs=True)
        _run_export(progress, f"{group.name}: atoms", group.exports, atoms=True)
        _run_export(progress, f"{group.name}: surrounding atoms", group.exports, surr_atoms=True)

    for iface in ifaces:
        name = getattr(iface, 'name', 'interface')

        _run_export(progress, f"{name}: surfaces", iface.export, surfs=True)
        _run_export(progress, f"{name}: atoms", iface.export, atoms=True)
        _run_export(progress, f"{name}: edges", iface.export, edges=True)
        _run_export(progress, f"{name}: logs", iface.export, logs=True)
        _run_export(progress, f"{name}: vertices", iface.export, verts=True)
        _run_export(progress, f"{name}: info", iface.export, info=True)

    progress.finish()


def export_large(sys):
    """
    Large export.

    System:
        - PDB
        - PyMOL atom script
        - Information

    Groups:
        - Shell vertices
        - Shell edges
        - Shell surfaces
        - Information
        - Edges
        - Vertices
        - Atoms
        - Surrounding atoms
        - Logs
        - Atom surfaces
        - Atom edges
        - Atom vertices

    Interfaces:
        - Balls
        - Surfaces
        - Edges
        - Vertices
        - Information
    """

    groups = [group for group in sys.groups if group.net is not None]
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)

    # 3 system + 12/group + 5/interface
    progress = ExportProgress(3 + 9 * len(groups) + 6 * len(ifaces), sys)

    _run_export(progress, "system PDB", sys.exports, pdb=True)
    _run_export(progress, "system PyMOL atoms", sys.exports, set_atoms=True)
    _run_export(progress, "system info", sys.exports, info=True)

    for group in groups:
        _set_group_directory(sys, group)

        _run_export(progress, f"{group.name}: shell vertices", group.exports, shell_verts=True)
        _run_export(progress, f"{group.name}: shell edges", group.exports, shell_edges=True)
        _run_export(progress, f"{group.name}: shell surfaces", group.exports, shell_surfs=True)
        _run_export(progress, f"{group.name}: info", group.exports, info=True)
        _run_export(progress, f"{group.name}: edges", group.exports, edges=True)
        _run_export(progress, f"{group.name}: vertices", group.exports, verts=True)
        _run_export(progress, f"{group.name}: atoms", group.exports, atoms=True)
        _run_export(progress, f"{group.name}: surrounding atoms", group.exports, surr_atoms=True)
        _run_export(progress, f"{group.name}: logs", group.exports, logs=True)
        _run_export(progress, f"{group.name}: atom surfaces", group.exports, atom_surfs=True)
        _run_export(progress, f"{group.name}: atom edges", group.exports, atom_edges=True)
        _run_export(progress, f"{group.name}: atom vertices", group.exports, atom_verts=True)

        _move_vert_file(sys, group)

    for iface in ifaces:
        name = getattr(iface, 'name', 'interface')

        _run_export(progress, f"{name}: balls", iface.export, balls=True)
        _run_export(progress, f"{name}: surfaces", iface.export, surfs=True)
        _run_export(progress, f"{name}: edges", iface.export, edges=True)
        _run_export(progress, f"{name}: vertices", iface.export, verts=True)
        _run_export(progress, f"{name}: info", iface.export, info=True)

    progress.finish()


def export_all(sys):
    """
    Export all supported outputs.

    Operations are called individually so the user receives progress
    throughout long exports.
    """

    groups = [group for group in sys.groups if group.net is not None]
    ifaces = [] if sys.ifaces is None else list(sys.ifaces)

    group_exports = [
        ("atoms", {"atoms": True}),
        ("shell surfaces", {"shell_surfs": True}),
        ("surfaces", {"surfs": True}),
        ("separate surfaces", {"sep_surfs": True}),
        ("shell edges", {"shell_edges": True}),
        ("edges", {"edges": True}),
        ("separate edges", {"sep_edges": True}),
        ("shell vertices", {"shell_verts": True}),
        ("vertices", {"verts": True}),
        ("separate vertices", {"sep_verts": True}),
        ("surrounding atoms", {"surr_atoms": True}),
        ("external atoms", {"ext_atoms": True}),
        ("logs", {"logs": True}),
        ("info", {"info": True}),
        ("atom surfaces", {"atom_surfs": True}),
        ("atom edges", {"atom_edges": True}),
        ("atom vertices", {"atom_verts": True}),
    ]

    interface_exports = [
        ("balls", {"balls": True}),
        ("surfaces", {"surfs": True}),
        ("atoms", {"atoms": True}),
        ("edges", {"edges": True}),
        ("logs", {"logs": True}),
        ("vertices", {"verts": True}),
        ("info", {"info": True}),
    ]

    total = 3 + len(group_exports) * len(groups) + len(interface_exports) * len(ifaces)
    progress = ExportProgress(total, sys)

    _run_export(progress, "system PDB", sys.exports, pdb=True)
    _run_export(progress, "system PyMOL atoms", sys.exports, set_atoms=True)
    _run_export(progress, "system info", sys.exports, info=True)

    for group in groups:
        _set_group_directory(sys, group)

        for name, kwargs in group_exports:
            _run_export(progress, f"{group.name}: {name}", group.exports, **kwargs)

        _move_vert_file(sys, group)

    for iface in ifaces:
        name = getattr(iface, 'name', 'interface')

        for export_name, kwargs in interface_exports:
            _run_export(progress, f"{name}: {export_name}", iface.export, **kwargs)

    progress.finish()


def other_exports(sys, usr_npt):
    """Run a user-requested standalone export."""

    option = usr_npt.lower()

    if option in {"a", "atoms"}:
        progress = ExportProgress(1, sys)

        _run_export(
            progress,
            "atom cells",
            write_atom_cells,
            net=sys.net,
            atoms=list(range(len(sys.net.balls))),
            directory=sys.files['dir']
        )

        progress.finish()

    elif option in {'logs', 'lgs'}:
        groups = [group for group in sys.groups if group.net is not None]
        progress = ExportProgress(len(groups) + 2, sys)

        for group in groups:
            _set_group_directory(sys, group)
            _run_export(progress, f"{group.name}: logs", group.exports, logs=True)

        _run_export(progress, "system PDB", sys.exports, pdb=True)
        _run_export(progress, "system PyMOL atoms", sys.exports, set_atoms=True)

        progress.finish()

    elif option in {'shell', 'shl'}:
        groups = [group for group in sys.groups if group.net is not None]
        progress = ExportProgress(len(groups), sys)

        for group in groups:
            _set_group_directory(sys, group)
            _run_export(progress, f"{group.name}: shell surfaces", group.exports, shell_surfs=True)

        progress.finish()