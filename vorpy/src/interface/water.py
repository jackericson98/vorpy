from vorpy.src.group import Group


def get_interface_water_residues(iface):
    """
    Return solvent water residues that contribute at least one atom to a
    retained interface surface.
    """
    if iface.net is None or iface.net.surfs is None:
        return []

    sol = getattr(iface.sys, "sol", None)
    solvent_residues = getattr(sol, "residues", None) or []

    if len(solvent_residues) == 0:
        return []

    surface_balls = set()

    for balls in iface.net.surfs["balls"]:
        surface_balls.update(int(ball) for ball in balls)

    touching_waters = []

    for residue in solvent_residues:
        residue_balls = set(int(ball) for ball in residue.atoms)

        if residue_balls & surface_balls:
            touching_waters.append(residue)

    return touching_waters


def get_water_interface_geometry(iface, residue):
    """
    Extract all retained interface surfaces, edges, and vertices associated
    with one water residue.
    """
    net = iface.net

    if net is None:
        return None

    water_balls = set(int(ball) for ball in residue.atoms)

    direct_surface_indices = []

    for surface_index, row in net.surfs.iterrows():
        defining_balls = set(int(ball) for ball in row["balls"])

        if defining_balls & water_balls:
            direct_surface_indices.append(surface_index)

    edge_indices = set()
    vertex_indices = set()

    for surface_index in direct_surface_indices:
        surface = net.surfs.loc[surface_index]

        edge_indices.update(
            int(edge_index)
            for edge_index in surface["edges"]
        )

        vertex_indices.update(
            int(vertex_index)
            for vertex_index in surface["verts"]
        )

    # Include any additional retained edges directly defined by a water atom.
    for edge_index, row in net.edges.iterrows():
        defining_balls = set(int(ball) for ball in row["balls"])

        if defining_balls & water_balls:
            edge_indices.add(edge_index)

            vertex_indices.update(
                int(vertex_index)
                for vertex_index in row["verts"]
            )

    # Include any retained vertices directly defined by a water atom.
    for vertex_index, row in net.verts.iterrows():
        defining_balls = set(int(ball) for ball in row["balls"])

        if defining_balls & water_balls:
            vertex_indices.add(vertex_index)

    surface_table = net.surfs.loc[
        sorted(direct_surface_indices)
    ].copy()

    edge_table = net.edges.loc[
        sorted(edge_indices)
    ].copy()

    vertex_table = net.verts.loc[
        sorted(vertex_indices)
    ].copy()

    surface_area = (
        float(surface_table["sa"].sum())
        if "sa" in surface_table.columns
        else None
    )

    contact_area = (
        float(surface_table["contact_area"].sum())
        if "contact_area" in surface_table.columns
        else None
    )

    return {
        "residue": residue,
        "water_balls": sorted(water_balls),
        "surface_indices": sorted(direct_surface_indices),
        "edge_indices": sorted(edge_indices),
        "vertex_indices": sorted(vertex_indices),
        "surfaces": surface_table,
        "edges": edge_table,
        "vertices": vertex_table,
        "surface_area": surface_area,
        "contact_area": contact_area,
    }


def analyze_interface_waters(iface):
    """
    Collect interface-associated geometry separately for every water residue
    touching the interface.
    """
    water_geometries = []

    for residue in get_interface_water_residues(iface):
        geometry = get_water_interface_geometry(
            iface=iface,
            residue=residue,
        )

        if geometry is not None:
            water_geometries.append(geometry)

    return water_geometries


def get_water_group_volume(water_group):
    """
    Return the complete volume assigned to the water group.

    This should be finalized against VorPy's authoritative analyzed-volume
    field.
    """
    net = water_group.net

    if (
            hasattr(net, "balls")
            and net.balls is not None
            and "vol" in net.balls.columns
    ):
        return float(
            net.balls["vol"].dropna().sum()
        )

    system_balls = water_group.sys.balls
    water_indices = list(water_group.ball_ndxs)

    if "vol" in system_balls.columns:
        return float(
            system_balls.loc[
                water_indices,
                "vol",
            ].dropna().sum()
        )

    return None


def build_interface_water_groups(
        iface,
        water_geometries,
):
    """
    Fully solve a normal Group for every water residue previously identified
    as touching the interface.
    """
    water_groups = []

    for water_index, water_geometry in enumerate(
            water_geometries
    ):
        residue = water_geometry["residue"]

        water_name = (
            f"{iface.name}_"
            f"{residue.name}_{residue.seq}"
        )

        water_group = Group(
            sys=iface.sys,
            name=water_name,
            residues=[residue],
            settings=iface.group1.settings.copy(),
            make_net=True,
            build_net=False,
        )

        water_group.build()

        water_geometry["full_group_geometry"] = (
            summarize_water_group(
                water_group=water_group,
            )
        )

        water_groups.append(water_group)

    return water_groups


def summarize_water_group(water_group):
    """
    Summarize the complete network built for one water molecule.
    """
    net = water_group.net

    if net is None:
        return None

    surface_area = None

    if (
            net.surfs is not None
            and "sa" in net.surfs.columns
    ):
        surface_area = float(
            net.surfs["sa"].dropna().sum()
        )

    return {
        "group": water_group,
        "group_name": water_group.name,
        "ball_indices": list(water_group.ball_ndxs),

        "vertex_count": (
            len(net.verts)
            if net.verts is not None
            else 0
        ),
        "edge_count": (
            len(net.edges)
            if net.edges is not None
            else 0
        ),
        "surface_count": (
            len(net.surfs)
            if net.surfs is not None
            else 0
        ),

        "surface_area": surface_area,

        # Fill this from the authoritative volume field once confirmed.
        "volume": get_water_group_volume(
            water_group
        ),

        "vertices": net.verts,
        "edges": net.edges,
        "surfaces": net.surfs,
    }