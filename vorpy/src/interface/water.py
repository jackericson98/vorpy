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
    Extract retained interface topology and side-specific surface metrics
    for one water residue.
    """
    net = iface.net

    if net is None:
        return None

    water_balls = set(
        int(ball)
        for ball in residue.atoms
    )

    surface_classes = classify_water_interface_surfaces(
        iface=iface,
        water_balls=water_balls,
    )

    group1_surface_indices = surface_classes["group1"]
    group2_surface_indices = surface_classes["group2"]
    other_surface_indices = surface_classes["other"]

    all_surface_indices = sorted(
        set(group1_surface_indices)
        | set(group2_surface_indices)
        | set(other_surface_indices)
    )

    edge_indices = set()
    vertex_indices = set()

    for surface_index in all_surface_indices:
        surface = net.surfs.loc[surface_index]

        edge_indices.update(
            int(edge_index)
            for edge_index in surface["edges"]
        )

        vertex_indices.update(
            int(vertex_index)
            for vertex_index in surface["verts"]
        )

    # Also capture any retained edge directly defined by a water atom.
    for edge_index, edge in net.edges.iterrows():
        edge_balls = set(
            int(ball)
            for ball in edge["balls"]
        )

        if edge_balls & water_balls:
            edge_indices.add(edge_index)

            vertex_indices.update(
                int(vertex_index)
                for vertex_index in edge["verts"]
            )

    # Capture any retained vertex directly defined by a water atom.
    for vertex_index, vertex in net.verts.iterrows():
        vertex_balls = set(
            int(ball)
            for ball in vertex["balls"]
        )

        if vertex_balls & water_balls:
            vertex_indices.add(vertex_index)

    group1_surfaces = net.surfs.loc[
        group1_surface_indices
    ].copy()

    group2_surfaces = net.surfs.loc[
        group2_surface_indices
    ].copy()

    other_surfaces = net.surfs.loc[
        other_surface_indices
    ].copy()

    all_surfaces = net.surfs.loc[
        all_surface_indices
    ].copy()

    return {
        "residue": residue,
        "water_name": (
            f"{residue.name}_{residue.seq}"
        ),
        "ball_indices": sorted(water_balls),

        "interface_geometry": {
            "surface_indices": all_surface_indices,
            "edge_indices": sorted(edge_indices),
            "vertex_indices": sorted(vertex_indices),

            "group1_surface_indices": (
                group1_surface_indices
            ),
            "group2_surface_indices": (
                group2_surface_indices
            ),
            "other_surface_indices": (
                other_surface_indices
            ),

            "group1": summarize_surfaces(
                group1_surfaces
            ),
            "group2": summarize_surfaces(
                group2_surfaces
            ),
            "other": summarize_surfaces(
                other_surfaces
            ),
            "total": summarize_surfaces(
                all_surfaces
            ),

            "contacts_group1": bool(
                group1_surface_indices
            ),
            "contacts_group2": bool(
                group2_surface_indices
            ),
            "bridging_water": bool(
                group1_surface_indices
                and group2_surface_indices
            ),
        },

        # Filled after the normal water Group has been solved.
        "full_group_geometry": None,
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

    surface_balls = {
        int(ball)
        for balls in iface.net.surfs["balls"]
        for ball in balls
    }

    water_balls = {
        int(ball)
        for residue in getattr(iface.sys.sol, "residues", [])
        for ball in residue.atoms
    }

    print("\n[WATER INDEX VALIDATION]")
    print(
        f"surface ball range: "
        f"{min(surface_balls)} to {max(surface_balls)}"
    )
    print(
        f"water ball range: "
        f"{min(water_balls)} to {max(water_balls)}"
    )
    print(
        f"surface/water overlap: "
        f"{len(surface_balls & water_balls)}"
    )

    return water_geometries

def classify_water_interface_surfaces(
        iface,
        water_balls,
):
    """
    Classify retained interface surfaces involving this water according to
    whether the non-water defining ball belongs to interface Group 1 or
    interface Group 2.
    """
    water_balls = set(int(ball) for ball in water_balls)
    group1_balls = set(int(ball) for ball in iface.group1_indices)
    group2_balls = set(int(ball) for ball in iface.group2_indices)

    group1_surface_indices = []
    group2_surface_indices = []
    other_surface_indices = []

    for surface_index, surface in iface.net.surfs.iterrows():
        surface_balls = set(
            int(ball)
            for ball in surface["balls"]
        )

        water_hits = surface_balls & water_balls

        if not water_hits:
            continue

        non_water_balls = surface_balls - water_balls

        if non_water_balls & group1_balls:
            group1_surface_indices.append(surface_index)

        elif non_water_balls & group2_balls:
            group2_surface_indices.append(surface_index)

        else:
            # Examples include water-water or water-ion surfaces.
            other_surface_indices.append(surface_index)

    return {
        "group1": group1_surface_indices,
        "group2": group2_surface_indices,
        "other": other_surface_indices,
    }

def summarize_surfaces(surface_table):
    """
    Return geometric statistics for a subset of interface surfaces.
    """
    if surface_table is None or len(surface_table) == 0:
        return {
            "surface_count": 0,
            "surface_area": 0.0,
            "contact_area": 0.0,
            "mean_curvature": None,
            "area_weighted_mean_curvature": None,
            "min_mean_curvature": None,
            "max_mean_curvature": None,
            "mean_gaussian_curvature": None,
            "area_weighted_gaussian_curvature": None,
        }

    areas = surface_table["sa"].dropna()
    mean_curvatures = surface_table["avg_mean_curv"].dropna()
    gaussian_curvatures = surface_table["avg_gauss_curv"].dropna()
    contact_areas = surface_table["contact_area"].dropna()

    common_mean_indices = areas.index.intersection(
        mean_curvatures.index
    )
    common_gaussian_indices = areas.index.intersection(
        gaussian_curvatures.index
    )

    weighted_mean_curvature = None

    if (
            len(common_mean_indices) > 0
            and float(areas.loc[common_mean_indices].sum()) > 0
    ):
        weighted_mean_curvature = float(
            (
                mean_curvatures.loc[common_mean_indices]
                * areas.loc[common_mean_indices]
            ).sum()
            / areas.loc[common_mean_indices].sum()
        )

    weighted_gaussian_curvature = None

    if (
            len(common_gaussian_indices) > 0
            and float(areas.loc[common_gaussian_indices].sum()) > 0
    ):
        weighted_gaussian_curvature = float(
            (
                gaussian_curvatures.loc[common_gaussian_indices]
                * areas.loc[common_gaussian_indices]
            ).sum()
            / areas.loc[common_gaussian_indices].sum()
        )

    return {
        "surface_count": len(surface_table),
        "surface_area": (
            float(areas.sum())
            if len(areas) > 0
            else 0.0
        ),
        "contact_area": (
            float(contact_areas.sum())
            if len(contact_areas) > 0
            else 0.0
        ),
        "mean_curvature": (
            float(mean_curvatures.mean())
            if len(mean_curvatures) > 0
            else None
        ),
        "area_weighted_mean_curvature": (
            weighted_mean_curvature
        ),
        "min_mean_curvature": (
            float(mean_curvatures.min())
            if len(mean_curvatures) > 0
            else None
        ),
        "max_mean_curvature": (
            float(mean_curvatures.max())
            if len(mean_curvatures) > 0
            else None
        ),
        "mean_gaussian_curvature": (
            float(gaussian_curvatures.mean())
            if len(gaussian_curvatures) > 0
            else None
        ),
        "area_weighted_gaussian_curvature": (
            weighted_gaussian_curvature
        ),
    }


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


def build_interface_water_groups(
        iface,
        water_geometries,
):
    """
    Fully solve a normal Group for every water residue previously identified
    as touching the interface.
    """
    water_groups = []

    print("\n" + "=" * 70)
    print("BUILDING INTERFACE WATER GROUPS")
    print("=" * 70)
    print(f"touching waters: {len(water_geometries)}")

    for water_index, water_geometry in enumerate(
            water_geometries
    ):
        residue = water_geometry["residue"]

        water_name = (
            f"{iface.name}_"
            f"{residue.name}_{residue.seq}"
        )

        print(
            f"\nWater {water_index + 1}/"
            f"{len(water_geometries)}"
        )
        print(f"  name: {water_name}")
        print(f"  atoms: {list(residue.atoms)}")

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