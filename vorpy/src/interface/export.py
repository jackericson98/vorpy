import os
from vorpy.src.output import write_off_verts
from vorpy.src.output import write_edges
from vorpy.src.output import write_surfs
from vorpy.src.output import write_pdb
from vorpy.src.output import write_interface_logs


def get_interface_atoms(iface):
    """Return the unique atom indices participating in the interface definition."""
    return sorted(
        set(getattr(iface, "group1_indices", []) or [])
        | set(getattr(iface, "group2_indices", []) or [])
    )

def _get_column(dataframe, *candidate_names):
    """
    Return the first matching DataFrame column from candidate_names.

    Matching is case-insensitive and ignores spaces and underscores so this
    works with both internal Network column names and exported-log headings.
    """
    if dataframe is None:
        return None

    normalized_columns = {
        str(column).lower().replace(" ", "").replace("_", ""): column
        for column in dataframe.columns
    }

    for candidate in candidate_names:
        normalized_candidate = (
            str(candidate)
            .lower()
            .replace(" ", "")
            .replace("_", "")
        )

        if normalized_candidate in normalized_columns:
            return normalized_columns[normalized_candidate]

    return None


def _numeric_series(dataframe, *candidate_names):
    """
    Return a numeric Series for the first matching column, or None.
    """
    column = _get_column(dataframe, *candidate_names)

    if column is None:
        return None

    return dataframe[column].dropna()


def _safe_sum(series, empty_value=0.0):
    if series is None:
        return None

    if len(series) == 0:
        return empty_value

    return float(series.sum())


def _safe_mean(series):
    return None if series is None or len(series) == 0 else float(series.mean())


def _safe_min(series):
    return None if series is None or len(series) == 0 else float(series.min())


def _safe_max(series):
    return None if series is None or len(series) == 0 else float(series.max())


def _weighted_mean(values, weights):
    """
    Return the area-weighted mean of values.

    Rows with missing values, missing weights, or non-positive weights are
    excluded.
    """
    if values is None or weights is None:
        return None

    valid_indices = values.index.intersection(weights.index)

    if len(valid_indices) == 0:
        return None

    valid_values = values.loc[valid_indices]
    valid_weights = weights.loc[valid_indices]

    valid_mask = (
        valid_values.notna()
        & valid_weights.notna()
        & (valid_weights > 0)
    )

    valid_values = valid_values[valid_mask]
    valid_weights = valid_weights[valid_mask]

    if len(valid_values) == 0 or float(valid_weights.sum()) == 0:
        return None

    return float(
        (valid_values * valid_weights).sum()
        / valid_weights.sum()
    )


def _format_metric(value, decimals=6):
    if value is None:
        return "Not available"

    return f"{value:.{decimals}f}"


def get_interface_surface_sets(iface):
    """
    Separate interface-network surfaces into:

    direct_surfaces
        One defining ball belongs to Group 1 and the other belongs to Group 2.

    group1_surfaces
        Both defining balls belong to Group 1.

    group2_surfaces
        Both defining balls belong to Group 2.

    support_surfaces
        At least one defining ball is outside the two explicit interface groups.

    The support category may contain geometry required to construct or close
    the dedicated interface network without representing direct inter-group
    contact.
    """
    if iface.net is None or iface.net.surfs is None:
        return None, None, None, None

    surfaces = iface.net.surfs

    ball1_column = _get_column(
        surfaces,
        "Ball 1",
        "ball 1",
        "ball1",
        "b1",
    )
    ball2_column = _get_column(
        surfaces,
        "Ball 2",
        "ball 2",
        "ball2",
        "b2",
    )

    # Some internal surface tables store both balls in one iterable column.
    balls_column = _get_column(
        surfaces,
        "balls",
        "surface balls",
    )

    group1_indices = set(int(index) for index in iface.group1_indices)
    group2_indices = set(int(index) for index in iface.group2_indices)

    if ball1_column is not None and ball2_column is not None:
        ball_pairs = [
            (int(ball1), int(ball2))
            for ball1, ball2 in zip(
                surfaces[ball1_column],
                surfaces[ball2_column],
            )
        ]

    elif balls_column is not None:
        ball_pairs = [
            (int(balls[0]), int(balls[1]))
            for balls in surfaces[balls_column]
        ]

    else:
        # The topology can still be summarized, but surface membership cannot
        # be classified without defining-ball information.
        empty = surfaces.iloc[0:0]
        return empty, empty, empty, surfaces

    direct_mask = []
    group1_mask = []
    group2_mask = []
    support_mask = []

    for ball1, ball2 in ball_pairs:
        ball1_in_group1 = ball1 in group1_indices
        ball2_in_group1 = ball2 in group1_indices
        ball1_in_group2 = ball1 in group2_indices
        ball2_in_group2 = ball2 in group2_indices

        is_direct = (
            (ball1_in_group1 and ball2_in_group2)
            or
            (ball1_in_group2 and ball2_in_group1)
        )

        is_group1 = ball1_in_group1 and ball2_in_group1
        is_group2 = ball1_in_group2 and ball2_in_group2

        is_support = not (
            is_direct
            or is_group1
            or is_group2
        )

        direct_mask.append(is_direct)
        group1_mask.append(is_group1)
        group2_mask.append(is_group2)
        support_mask.append(is_support)

    return (
        surfaces.loc[direct_mask],
        surfaces.loc[group1_mask],
        surfaces.loc[group2_mask],
        surfaces.loc[support_mask],
    )


def write_surface_statistics(info, title, surfaces):
    """
    Write area, curvature, contact, and overlap statistics for a collection
    of surfaces.
    """
    if surfaces is None:
        info.write(f"{title}:\n")
        info.write("  Not available\n\n")
        return

    surface_areas = _numeric_series(
        surfaces,
        "sa",
    )

    mean_curvatures = _numeric_series(
        surfaces,
        "mean_curv",
    )

    average_mean_curvatures = _numeric_series(
        surfaces,
        "avg_mean_curv",
    )

    gaussian_curvatures = _numeric_series(
        surfaces,
        "gauss_curv",
    )

    average_gaussian_curvatures = _numeric_series(
        surfaces,
        "avg_gauss_curv",
    )

    # Integrated curvature geometry
    integrated_mean_curvatures = _numeric_series(
        surfaces,
        "int_mean_curv",
        "Integrated Mean Curvature",
    )

    integrated_mean_curvature_squared = _numeric_series(
        surfaces,
        "int_mean_curv_sq",
        "Integrated Mean Curvature Squared",
    )

    integrated_gaussian_curvatures = _numeric_series(
        surfaces,
        "int_gauss_curv",
        "Integrated Gaussian Curvature",
    )

    contact_areas = _numeric_series(
        surfaces,
        "contact_area",
    )

    contact_areas = _numeric_series(
        surfaces,
        "contact_area",
    )

    overlaps = _numeric_series(
        surfaces,
        "overlap",
    )

    info.write(f"{title}:\n")
    info.write(f"  Number of surfaces: {len(surfaces)}\n")
    info.write(
        f"  Total surface area: "
        f"{_format_metric(_safe_sum(surface_areas))} Å²\n"
    )
    info.write(
        f"  Mean surface area: "
        f"{_format_metric(_safe_mean(surface_areas))} Å²\n"
    )
    info.write(
        f"  Minimum surface area: "
        f"{_format_metric(_safe_min(surface_areas))} Å²\n"
    )
    info.write(
        f"  Maximum surface area: "
        f"{_format_metric(_safe_max(surface_areas))} Å²\n"
    )

    info.write("\n  Mean curvature:\n")
    info.write(
        f"    Mean of surface-average mean curvature: "
        f"{_format_metric(_safe_mean(average_mean_curvatures))} Å⁻¹\n"
    )
    info.write(f"    Area-weighted surface-average mean curvature: {_format_metric(_weighted_mean(average_mean_curvatures, surface_areas))} Å⁻¹\n")
    info.write(
        f"    Minimum surface-average mean curvature: "
        f"{_format_metric(_safe_min(average_mean_curvatures))} Å⁻¹\n"
    )
    info.write(
        f"    Maximum surface-average mean curvature: "
        f"{_format_metric(_safe_max(average_mean_curvatures))} Å⁻¹\n"
    )
    info.write(
        f"    Maximum local mean curvature: "
        f"{_format_metric(_safe_max(mean_curvatures))} Å⁻¹\n"
    )

    info.write("\n  Gaussian curvature:\n")
    info.write(f"    Mean of surface-average Gaussian curvature: "
        f"{_format_metric(_safe_mean(average_gaussian_curvatures))} Å⁻²\n")
    info.write(
        f"    Area-weighted surface-average Gaussian curvature: "
        f"{_format_metric(_weighted_mean(average_gaussian_curvatures, surface_areas))} Å⁻²\n")
    info.write(
        f"    Minimum surface-average Gaussian curvature: "
        f"{_format_metric(_safe_min(average_gaussian_curvatures))} Å⁻²\n"
    )
    info.write(
        f"    Maximum surface-average Gaussian curvature: "
        f"{_format_metric(_safe_max(average_gaussian_curvatures))} Å⁻²\n")
    info.write(
        f"    Maximum local Gaussian curvature: "
        f"{_format_metric(_safe_max(gaussian_curvatures))} Å⁻²\n"
    )

    # ------------------------------------------------------------------
    # Integrated curvature geometry
    # ------------------------------------------------------------------

    info.write("\n  Integrated curvature geometry:\n")

    total_int_mean_curv = _safe_sum(
        integrated_mean_curvatures
    )

    total_int_mean_curv_sq = _safe_sum(
        integrated_mean_curvature_squared
    )

    total_int_gauss_curv = _safe_sum(
        integrated_gaussian_curvatures
    )

    info.write(
        f"    Integrated mean curvature (∫H dA): "
        f"{_format_metric(total_int_mean_curv)} Å\n"
    )

    info.write(
        f"    Integrated squared mean curvature (∫H² dA): "
        f"{_format_metric(total_int_mean_curv_sq)}\n"
    )

    info.write(
        f"    Integrated Gaussian curvature (∫K dA): "
        f"{_format_metric(total_int_gauss_curv)}\n"
    )

    # Derive area-normalized curvature from the integrated quantities
    total_surface_area = _safe_sum(surface_areas)

    if total_surface_area is not None and total_surface_area > 0.0:

        if total_int_mean_curv is not None:
            info.write(
                f"    Area-normalized mean curvature: "
                f"{_format_metric(total_int_mean_curv / total_surface_area)} Å⁻¹\n"
            )

        if total_int_gauss_curv is not None:
            info.write(
                f"    Area-normalized Gaussian curvature: "
                f"{_format_metric(total_int_gauss_curv / total_surface_area)} Å⁻²\n"
            )

    info.write("\n  Contact and overlap:\n")
    info.write(
        f"    Total contact area: "
        f"{_format_metric(_safe_sum(contact_areas))} Å²\n"
    )
    info.write(
        f"    Surfaces with positive contact area: "
        f"{0 if contact_areas is None else int((contact_areas > 0).sum())}\n"
    )
    info.write(
        f"    Total overlap measure: "
        f"{_format_metric(_safe_sum(overlaps))} Å\n"
    )
    info.write(
        f"    Surfaces with positive overlap: "
        f"{0 if overlaps is None else int((overlaps > 0).sum())}\n"
    )

    info.write("\n")


def export_info(iface, directory=None):
    """
    Export a detailed geometric and topological summary of an Interface.
    """
    if directory is None:
        directory = iface.dir

    os.makedirs(directory, exist_ok=True)

    net = iface.net

    num_verts = (
        0
        if net is None or net.verts is None
        else len(net.verts)
    )
    num_edges = (
        0
        if net is None or net.edges is None
        else len(net.edges)
    )
    num_surfs = (
        0
        if net is None or net.surfs is None
        else len(net.surfs)
    )

    group1_atoms = sorted(set(iface.group1_indices))
    group2_atoms = sorted(set(iface.group2_indices))
    interface_atoms = sorted(
        set(group1_atoms) | set(group2_atoms)
    )

    (
        direct_surfaces,
        group1_surfaces,
        group2_surfaces,
        support_surfaces,
    ) = get_interface_surface_sets(iface)

    # One-sided or partially populated interfaces may not have explicit
    # surface-classification sets. Normalize them to empty DataFrames so
    # topology/info export remains valid.
    empty_surfaces = None if net is None or net.surfs is None else net.surfs.iloc[0:0].copy()
    if direct_surfaces is None:
        direct_surfaces = empty_surfaces
    if group1_surfaces is None:
        group1_surfaces = empty_surfaces
    if group2_surfaces is None:
        group2_surfaces = empty_surfaces
    if support_surfaces is None:
        support_surfaces = empty_surfaces

    group1_name = getattr(getattr(iface, "group1", None), "name", "group1")
    group2_name = getattr(iface, "group2_name", None) or getattr(getattr(iface, "group2", None), "name", "surrounding")
    interface_id = getattr(iface, "interface_id", getattr(iface, "name", "interface"))

    file_path = os.path.join(directory, "info.txt")

    with open(file_path, "w", encoding="utf-8") as info:
        info.write(f"{iface.name} - {iface.sys.name}\n\n")

        info.write("Interface definition:\n")
        info.write(f"  Interface ID: {interface_id}\n")
        info.write(f"  Group 1: {group1_name}\n")
        info.write(f"  Group 2: {group2_name}\n\n")

        info.write("Interface atom membership:\n")
        info.write(f"  Group 1 atoms: {len(group1_atoms)}\n")
        info.write(f"  Group 2 atoms: {len(group2_atoms)}\n")
        info.write(
            f"  Unique interface-network atoms: "
            f"{len(interface_atoms)}\n"
        )

        shared_atoms = set(group1_atoms) & set(group2_atoms)

        info.write(
            f"  Atoms shared by both definitions: "
            f"{len(shared_atoms)}\n\n"
        )

        info.write("Interface network topology:\n")
        info.write(f"  Vertices: {num_verts}\n")
        info.write(f"  Edges: {num_edges}\n")
        info.write(f"  Surfaces: {num_surfs}\n")

        if num_verts > 0:
            info.write(
                f"  Edges per vertex: "
                f"{num_edges / num_verts:.6f}\n"
            )
            info.write(
                f"  Surfaces per vertex: "
                f"{num_surfs / num_verts:.6f}\n"
            )

        info.write("\n")

        info.write("Surface classification:\n")
        info.write(
            f"  Direct Group 1–Group 2 surfaces: "
            f"{0 if direct_surfaces is None else len(direct_surfaces)}\n"
        )
        info.write(
            f"  Group 1 internal surfaces: "
            f"{0 if group1_surfaces is None else len(group1_surfaces)}\n"
        )
        info.write(
            f"  Group 2 internal surfaces: "
            f"{0 if group2_surfaces is None else len(group2_surfaces)}\n"
        )
        info.write(
            f"  Supporting/unclassified surfaces: "
            f"{0 if support_surfaces is None else len(support_surfaces)}\n\n"
        )

        write_surface_statistics(
            info,
            "Full interface-network surface statistics",
            net.surfs if net is not None else None,
        )

        write_surface_statistics(
            info,
            "Direct inter-group surface statistics",
            direct_surfaces,
        )

        write_surface_statistics(
            info,
            f"{group1_name} internal surface statistics",
            group1_surfaces,
        )

        write_surface_statistics(
            info,
            f"{group2_name} internal surface statistics",
            group2_surfaces,
        )

        write_surface_statistics(
            info,
            "Supporting surface statistics",
            support_surfaces,
        )

        if net is not None and net.edges is not None:
            edge_lengths = _numeric_series(net.edges, "length")

            info.write("Edge geometry:\n")
            info.write(
                f"  Total edge length: "
                f"{_format_metric(_safe_sum(edge_lengths))} Å\n"
            )
            info.write(
                f"  Average edge length: "
                f"{_format_metric(_safe_mean(edge_lengths))} Å\n"
            )
            info.write(
                f"  Minimum edge length: "
                f"{_format_metric(_safe_min(edge_lengths))} Å\n"
            )
            info.write(
                f"  Maximum edge length: "
                f"{_format_metric(_safe_max(edge_lengths))} Å\n\n"
            )

        if net is not None and net.verts is not None:
            vertex_radii = _numeric_series(net.verts, "rad")

            info.write("Vertex geometry:\n")
            info.write(
                f"  Average vertex radius: "
                f"{_format_metric(_safe_mean(vertex_radii))} Å\n"
            )
            info.write(
                f"  Minimum vertex radius: "
                f"{_format_metric(_safe_min(vertex_radii))} Å\n"
            )
            info.write(
                f"  Maximum vertex radius: "
                f"{_format_metric(_safe_max(vertex_radii))} Å\n"
            )

            negative_radii = (
                0
                if vertex_radii is None
                else int((vertex_radii < 0).sum())
            )

            info.write(
                f"  Vertices with negative radius: "
                f"{negative_radii}\n\n"
            )


def interface_exports(iface, all_=False, atoms=False, surfs=False, edges=False, verts=False, logs=False, info=False,
                      group_info=False, round_to=3):
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

    if group_info or all_:
        exported_group_ids = set()

        for group in (iface.group1, iface.group2):
            if group is None:
                continue

            group_id = str(getattr(group, "group_id", None) or group.name)

            if group_id in exported_group_ids:
                continue

            exported_group_ids.add(group_id)

            group_directory = os.path.join(
                iface.sys.files["dir"],
                group.name,
            )

            export_interface_group_info(
                group=group,
                directory=group_directory,
            )

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
        if iface.net.verts is not None and len(iface.net.verts) > 0:
            write_off_verts(
                iface.net,
                list(range(len(iface.net.verts))),
                directory=iface.dir,
                file_name="verts",
                color=iface.net.settings["vert_col"],
            )

    if edges or all_:
        if iface.net.edges is not None and len(iface.net.edges) > 0:
            write_edges(
                iface.net,
                list(range(len(iface.net.edges))),
                directory=iface.dir,
                file_name="edges",
                color=iface.net.settings["edge_col"],
            )

    if surfs or all_:
        if iface.net.surfs is not None and len(iface.net.surfs) > 0:
            write_surfs(
                iface.net,
                list(range(len(iface.net.surfs))),
                directory=iface.dir,
                file_name="surfs",
            )

def export_interface_group_info(group, directory):
    """
    Export information for a Group that participated in an interface-only
    calculation.

    The group does not need a complete Network. This exports the group's
    selected atoms and interface relationship metadata only.
    """
    os.makedirs(directory, exist_ok=True)

    group_indices = sorted(
        set(int(index) for index in (getattr(group, "ball_ndxs", []) or []))
    )

    info_path = os.path.join(directory, "info.txt")

    with open(info_path, "w", encoding="utf-8") as info:
        info.write(f"{group.name} - {group.sys.name}\n\n")

        info.write("Group definition:\n")
        info.write(f"  Group ID: {getattr(group, 'group_id', group.name)}\n")
        info.write(f"  Selected atoms: {len(group_indices)}\n")
        info.write(f"  Full group network built: {group.net is not None}\n\n")

        info.write("Selections:\n")
        info.write(f"  Atom selections: {len(group.atms or [])}\n")
        info.write(f"  Residue selections: {len(group.rsds or [])}\n")
        info.write(f"  Chain selections: {len(group.chns or [])}\n")
        info.write(f"  Molecule selections: {len(group.mols or [])}\n\n")

        metadata = getattr(group, "interface_metadata", {}) or {}

        info.write("Interfaces:\n")

        if not metadata:
            info.write("  None\n")
        else:
            for interface_id, record in metadata.items():
                info.write(f"  Interface ID: {interface_id}\n")
                info.write(
                    f"    Name: "
                    f"{record.get('interface_name', interface_id)}\n"
                )
                info.write(
                    f"    Other group: "
                    f"{record.get('other_group_name', 'surrounding')}\n"
                )
                info.write(
                    f"    Side: {record.get('side')}\n"
                )
                info.write(
                    f"    Status: {record.get('status', 'defined')}\n"
                )
                info.write(
                    f"    Built: {record.get('built', False)}\n"
                )

    base_file = group.sys.files.get("base_file")

    if (
            base_file is not None
            and not str(base_file).lower().endswith(".txt")
    ):
        write_pdb(
            atoms=group_indices,
            file_name="group_atoms",
            directory=directory,
            sys=group.sys,
        )