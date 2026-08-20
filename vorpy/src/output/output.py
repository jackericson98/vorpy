import os
import time
import shutil
from vorpy.src.output import write_atom_cells


def export_micro(sys):
    """
    Smallest output function. Outputs the information for the system, the groups, and the system's interfaces.
    """
    # Export the information for the system.
    sys.exports(info=True)
    # Loop through the groups in the system
    for group in sys.group:
        # Set up the group directory
        if group.dir is None:
            group.dir = sys.files['dir'] + '/' + group.name
            os.makedirs(group.dir, exist_ok=True)
        # Export the information for the group
        group.exports(info=True)
    # Loop through the interfaces for the groups.
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            # Export the interface information
            iface.export(info=True)


def export_tiny(sys):
    """
    Second smallest of the exports. Outputs are:

    System:
        1. General Information
        2. Set balls script for pymol
        4. The PDB file
        5. The balls file
    Groups:
        1. General Information
        2. Shell for the group
        3. Logs for the group
    Interfaces:
        1.
    """
    sys.exports(info=True, set_atoms=True, pbd=True, balls=True)
    for group in sys.groups:
        group.dir = sys.files['dir'] + '/' + group.name
        os.makedirs(group.dir, exist_ok=True)
        group.export(info=True, shell=True, logs=True)
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(info=True)


def benchmark_exports(group, repeats=3):
    tests = {
        "shell_surfs": {"shell_surfs": True},
        "surfs": {"surfs": True},
        "shell_edges": {"shell_edges": True},
        "edges": {"edges": True},
        "shell_verts": {"shell_verts": True},
        "verts": {"verts": True},
        "logs": {"logs": True},
        "atoms": {"atoms": True},
        "surr_atoms": {"surr_atoms": True},
    }

    net = group.net

    print("\n" + "=" * 80)
    print("EXPORT BENCHMARK")
    print("=" * 80)
    print(f"Atoms:    {len(net.balls):,}")
    print(f"Surfaces: {len(net.surfs):,}")
    print(f"Edges:    {len(net.edges):,}")
    print(f"Vertices: {len(net.verts):,}")
    print("=" * 80)

    results = {}

    for name, kwargs in tests.items():
        times = []

        for n in range(repeats):
            start = time.perf_counter()
            group.exports(**kwargs)
            elapsed = time.perf_counter() - start
            times.append(elapsed)

            print(f"{name:<15} run {n + 1}: {elapsed:10.3f} s")

        avg = sum(times) / len(times)
        results[name] = avg

    print("\n" + "=" * 80)
    print("AVERAGE EXPORT TIMES")
    print("=" * 80)

    for name, elapsed in sorted(results.items(), key=lambda x: x[1], reverse=True):
        print(f"{name:<15} {elapsed:10.3f} s  {elapsed / 60:8.2f} min")


def timed_export(name, func, **kwargs):
    start = time.perf_counter()
    print(f"\nEXPORT: {name}...", flush=True)
    func(**kwargs)
    elapsed = time.perf_counter() - start
    print(f"EXPORT: {name} finished in {elapsed:.3f} s ({elapsed / 60:.2f} min)", flush=True)
    return elapsed


def export_med(sys):
    """Medium export with individual export timing."""
    total_start = time.perf_counter()
    timings = {}

    # System exports
    timings["system_pdb"] = timed_export("System PDB", sys.exports, pdb=True)
    timings["system_set_atoms"] = timed_export("System set atoms", sys.exports, set_atoms=True)
    timings["system_info"] = timed_export("System info", sys.exports, info=True)

    # Group exports
    for group in sys.groups:
        if group.net is None:
            continue

        net = group.net
        print("\nNETWORK EXPORT SIZE")
        print(f"Atoms:    {len(net.balls):,}")
        print(f"Surfaces: {len(net.surfs):,}")
        print(f"Edges:    {len(net.edges):,}")
        print(f"Vertices: {len(net.verts):,}")
        print("\nSURFACE GEOMETRY SIZE")
        surf_points = sum(len(points) for points in net.surfs['points'])
        surf_tris = sum(len(tris) for tris in net.surfs['tris'])
        print(f"Surface points:    {surf_points:,}")
        print(f"Surface triangles: {surf_tris:,}")
        print(f"Points/surface:    {surf_points / len(net.surfs):,.1f}")
        print(f"Tris/surface:      {surf_tris / len(net.surfs):,.1f}")

        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name):
            group.dir = sys.files['dir'] + '/' + group.name
            try:
                os.makedirs(group.dir, exist_ok=True)
            except FileNotFoundError:
                group.dir = sys.files['dir'] + '/group'

        print(f"\n{'=' * 70}\nEXPORTING GROUP: {group.name}\n{'=' * 70}")

        exports = {
            "shell_surfs": {"shell_surfs": True},
            "surfs": {"surfs": True},
            "shell_edges": {"shell_edges": True},
            "edges": {"edges": True},
            "shell_verts": {"shell_verts": True},
            "verts": {"verts": True},
            "logs": {"logs": True},
            "atoms": {"atoms": True},
            "surr_atoms": {"surr_atoms": True},
        }

        for name, kwargs in exports.items():
            key = f"group_{group.name}_{name}"
            timings[key] = timed_export(f"{group.name}: {name}", group.exports, **kwargs)

    # Interface exports
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface_name = getattr(iface, "name", "interface")
            print(f"\n{'=' * 70}\nEXPORTING INTERFACE: {iface_name}\n{'=' * 70}")

            exports = {
                "surfs": {"surfs": True},
                "atoms": {"atoms": True},
                "edges": {"edges": True},
                "logs": {"logs": True},
                "verts": {"verts": True},
                "info": {"info": True},
            }

            for name, kwargs in exports.items():
                key = f"interface_{iface_name}_{name}"
                timings[key] = timed_export(f"{iface_name}: {name}", iface.export, **kwargs)

    total = time.perf_counter() - total_start

    print(f"\n{'=' * 70}")
    print("EXPORT TIMING SUMMARY")
    print(f"{'=' * 70}")

    for name, elapsed in sorted(timings.items(), key=lambda x: x[1], reverse=True):
        print(f"{elapsed:12.3f} s  {elapsed / 60:9.2f} min  {name}")

    print(f"{'-' * 70}")
    print(f"{total:12.3f} s  {total / 60:9.2f} min  TOTAL")


def export_large(sys):
    """
    Large group exports. Exports the basic system files and the shell vertices, the shell surfaces, the information,
    the edges, the vertices, the atosm the surrounding atoms, the logs, the atom surfaces, the atom edges, and the
    atom vertices for each group
    """
    # Export the system exports
    sys.exports(pdb=True, set_atoms=True, info=True)
    # Loop through the groups and export the listed items
    for group in sys.groups:
        # Set and make the group directory
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name):
            group.dir = sys.files['dir'] + '/' + group.name
            os.makedirs(group.dir, exist_ok=True)
        # Export the group exports
        group.exports(shell_verts=True, shell_edges=True, shell_surfs=True, info=True, edges=True, verts=True,
                      atoms=True, surr_atoms=True, logs=True, atom_surfs=True, atom_edges=True, atom_verts=True)
        # Check to see if the verts are in the system directory and if so move them to the group folder
        if os.path.exists(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt'):
            shutil.move(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt',
                        group.dir + '/' + group.settings['net_type'] + '_verts.txt')
    # Export the interfaces
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(balls=True, surfs=True, edges=True, verts=True, info=True)


def export_all(sys):
    """
    Export all. Exports everything there is to export and makes a massive comprehensive set of files that will take a
    lot of space
    """
    # Export the system stuff
    sys.exports(pdb=True, info=True, set_atoms=True)
    # For each group in the system export the
    for group in sys.groups:
        # Set and make the group directory
        if group.dir is None or not os.path.exists(sys.files['dir'] + '/' + group.name):
            group.dir = sys.files['dir'] + '/' + group.name
            os.makedirs(group.dir, exist_ok=True)
        group.dir = sys.files['dir'] + '/' + group.name
        os.makedirs(group.dir, exist_ok=True)
        group.exports(atoms=True, shell_surfs=True, surfs=True, info=True, ext_atoms=True, sep_surfs=True, sep_edges=True,
                      sep_verts=True, verts=True, edges=True, surr_atoms=True, logs=True)

        # Check to see if the verts are in the system directory and if so move them to the group folder
        if os.path.exists(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt'):
            shutil.move(sys.files['dir'] + '/' + group.settings['net_type'] + '_verts.txt',
                        group.dir + '/' + group.settings['net_type'] + '_verts.txt')
    # Make the
    if sys.ifaces is not None:
        for iface in sys.ifaces:
            iface.export(all=True)


def other_exports(sys, usr_npt):
    """

    :param sys:
    :param usr_npt:
    :return:
    """
    # If the first word is atom
    if usr_npt.lower() in {"a", "atoms"}:
        write_atom_cells(sys.net.atoms['num'], sys.files['dir'])
    # If the first word is logs
    elif usr_npt.lower() in {'logs', 'lgs'}:
        for group in sys.groups:
            group.exports(logs=True)
        sys.exports(pdb=True, set_atoms=True)
    # If the first word is shell
    elif usr_npt.lower() in {'shell', 'shl'}:
        for grp in sys.groups:
            grp.exports(shell_surfs=True)
    # If the first word is network
    elif usr_npt.lower() in {'net', 'network'}:
        sys.exports(network=True)

