import os
from vorpy.src.output import write_off_verts
from vorpy.src.output import write_edges
from vorpy.src.output import write_surfs
from vorpy.src.output import write_pdb
from vorpy.src.output import write_interface_logs

def get_interface_atoms(iface):
    return sorted(
        set(iface.group1_indices)
        | set(iface.group2_indices)
    )


def export_info(iface, directory=None):
    if directory is None:
        directory = iface.dir

    os.makedirs(directory, exist_ok=True)

    net = iface.net

    num_verts = 0 if net is None or net.verts is None else len(net.verts)
    num_edges = 0 if net is None or net.edges is None else len(net.edges)
    num_surfs = 0 if net is None or net.surfs is None else len(net.surfs)

    group1_atoms = sorted(set(iface.group1_indices))
    group2_atoms = sorted(set(iface.group2_indices))
    interface_atoms = sorted(set(group1_atoms) | set(group2_atoms))

    file_path = os.path.join(directory, "info.txt")

    with open(file_path, "w", encoding="utf-8") as info:
        info.write(f"{iface.name} - {iface.sys.name}\n\n")

        info.write("Interface groups:\n")
        info.write(f"  Group 1: {iface.group1.name}\n")
        info.write(f"  Group 2: {iface.group2_name}\n\n")

        info.write("Interface atoms:\n")
        info.write(f"  Group 1 atoms: {len(group1_atoms)}\n")
        info.write(f"  Group 2 atoms: {len(group2_atoms)}\n")
        info.write(f"  Unique interface atoms: {len(interface_atoms)}\n\n")

        info.write("Interface network:\n")
        info.write(f"  Vertices: {num_verts}\n")
        info.write(f"  Edges: {num_edges}\n")
        info.write(f"  Surfaces: {num_surfs}\n")


def interface_exports(
        iface,
        all_=False,
        atoms=False,
        surfs=False,
        edges=False,
        verts=False,
        logs=False,
        info=False,
        round_to=3
):
    """
    Export data belonging to an Interface and its dedicated Network.
    """
    if iface.net is None:
        print(f'Interface "{iface.name}" has no network to export.')
        return

    if iface.dir is None:
        iface.dir = os.path.join(
            iface.sys.files["dir"],
            iface.name,
        )

    os.makedirs(iface.dir, exist_ok=True)

    if atoms or all_:
        interface_atoms = get_interface_atoms(iface)

        if iface.sys.files["base_file"][-3:].lower() != "txt":
            write_pdb(
                atoms=interface_atoms,
                file_name="interface_atoms",
                directory=iface.dir,
                sys=iface.sys,
            )

    if info or all_:
        export_info(iface)

    if logs or all_:
        write_interface_logs(
            iface,
            round_to=round_to,
        )

    if verts or all_:
        write_off_verts(
            iface.net,
            list(range(len(iface.net.verts))),
            directory=iface.dir,
            file_name="verts",
            color=iface.net.settings["vert_col"],
        )

    if edges or all_:
        write_edges(
            iface.net,
            list(range(len(iface.net.edges))),
            directory=iface.dir,
            file_name="edges",
            color=iface.net.settings["edge_col"],
        )

    if surfs or all_:
        write_surfs(
            iface.net,
            list(range(len(iface.net.surfs))),
            directory=iface.dir,
            file_name="surfs",
        )