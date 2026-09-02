import os
import pandas as pd
import csv
import numpy as np
from datetime import datetime
from vorpy.src.calculations import round_func
from vorpy.src.version import __version__


def write_logs(group, net_name=None, round_to=None):
    """
    Exports a comprehensive log file containing detailed information about the network analysis.

    This function generates a CSV log file with multiple sections:

    Build Information:
    - Network name, location, and completion date
    - Network type and key parameters (surface resolution, box size, max vertices)
    - Performance metrics (total time, vertex processing time, connection time, etc.)
    - Maximum vertex radius found

    Group Information:
    - Basic properties (name, volume, surface area, mass, density)
    - Center of mass (both standard and VDW)
    - Moment of inertia tensors (standard and spatial)

    Atoms:
    - Basic atom properties (index, name, residue info, chain, mass)
    - Spatial information (coordinates, radius, volumes)
    - Curvature metrics (mean and Gaussian curvatures)
    - Topological properties (sphericity, isometric quotient)
    - Neighbor analysis (count, distances, overlaps)
    - Contact areas and volumes
    - Center of mass and moment of inertia

    Args:
        group: Group object containing the network and system information
        net_name (str, optional): Additional identifier for the log file name
        round_to (int, optional): Number of decimal places to round numerical values to. Defaults to 3.
    """

    if round_to is None:
        round_to = group.sys.round_to

    net = group.net
    r = round_func(round_to)

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    if "system_num" in net.balls.columns:
        sys_nums = [int(_) for _ in net.balls["system_num"].tolist()]
    else:
        sys_nums = [int(_) for _ in net.balls["num"].tolist()]

    sys_balls = group.sys.balls.iloc[sys_nums].to_dict(orient='records')

    # Pull surface topology out of pandas once.
    surf_balls = net.surfs['balls'].tolist()

    # ------------------------------------------------------------------
    # Write log
    # ------------------------------------------------------------------

    with open(group.settings['net_type'] + "_logs.csv", 'w', newline='', buffering=1024 * 1024) as log_file:
        lg_fl = csv.writer(log_file, lineterminator='\n')

        # ==============================================================
        # Build + group information
        # ==============================================================

        lg_fl.writerow(["build informaiton"])
        lg_fl.writerow([
            "Name", "Location", "Completion Date", "Network Type", "Surface Resolution", "Box Size",
            "Maximum Allowable Vertex", "Total Time", "Vertex Time", "Connect Time",
            "Surface Building Time", "Analysis time", "Maximum Found Vertex", "vorPy version"
        ])

        lg_fl.writerow([
            group.sys.name,
            group.sys.files['base_file'],
            datetime.now(),
            net.settings['net_type'],
            net.settings['surf_res'],
            net.settings['box_size'],
            net.settings['max_vert'],
            r(net.metrics['tot']),
            r(net.metrics['vert']),
            r(net.metrics['con']),
            r(net.metrics['surf']),
            r(net.metrics['anal']),
            r(max(net.verts['rad'])),
            __version__
        ])

        lg_fl.writerow(["group information"])
        lg_fl.writerow([
            "Name", "Volume", "Surface Area", "Mass", "Density", "Center of Mass",
            "VDW Volume", "VDW Center of Mass", "Moment of Inertia", "Spatial Moment of Inertia"
        ])

        group.get_info()

        log_name = group.sys.name
        if getattr(net, "rebuilt_from_logs", False):
            log_name = f"{log_name}_rebuilt_from_logs"

        lg_fl.writerow([
            log_name,
            r(group.vol),
            r(group.sa),
            float(group.mass),
            r(group.density),
            [float(r(_)) for _ in group.com],
            r(group.vdw_vol),
            [float(r(_)) for _ in group.vdw_com],
            [[float(r(__)) for __ in _] for _ in group.moi],
            [[float(r(__)) for __ in _] for _ in group.spatial_moment]
        ])

        # ==============================================================
        # Atoms
        # ==============================================================

        lg_fl.writerow(["Atoms"])
        lg_fl.writerow(["Index", "Name", "Residue", "Residue Sequence", "Chain", "Mass", "X", "Y", "Z", "Radius",
                        "Volume", "Van Der Waals Volume", "Surface Area", "Complete Cell?", "Maximum Mean Curvature",
                        "Average Mean Surface Curvature", "Maximum Gaussian Curvature",
                        "Average Gaussian Surface Curvature", "Integrated Mean Curvature",
                        "Integrated Mean Curvature Squared", "Integrated Gaussian Curvature",
                        "Representative Surface Energy", "Sphericity", "Isometric Quotient", "Inner Ball?",
                        "Number of Neighbors", "Closest Neighbor", "Closest Neighbor Distance",
                        "Layer Distance Average", "Layer Distance RMSD", "Minimum Point Distance",
                        "Maximum Point Distance", "Number of Overlaps", "Contact Area", "Non-Overlap Volume",
                        "Overlap Volume", "Center of Mass", "Moment of Inertia Tensor", "Bounding Box", "neighbors"])

        atom_rows = net.balls.itertuples(index=True, name='AtomRow')

        for atom in atom_rows:
            i = atom.Index
            sys_ball = sys_balls[i]

            if atom.sa == 0:
                continue

            if not atom.complete:
                continue

            atom_num = atom.num

            nbrs = []
            for surf_ndx in atom.surfs:
                balls = surf_balls[surf_ndx]
                nbrs.append(balls[0] if balls[0] != atom_num else balls[1])

            lg_fl.writerow([
                i,
                sys_ball['name'],
                sys_ball['res_name'],
                sys_ball['res_seq'],
                sys_ball['chain_name'],
                sys_ball['mass'],
                atom.loc[0],
                atom.loc[1],
                atom.loc[2],
                atom.rad,
                r(atom.vol),
                r(atom.vdw_vol),
                r(atom.sa),
                atom.complete,
                r(atom.max_mean_curv),
                r(atom.avg_mean_surf_curv),
                r(atom.max_gauss_curv),
                r(atom.avg_gauss_surf_curv),
                r(atom.int_mean_curv),
                r(atom.int_mean_curv_sq),
                r(atom.int_gauss_curv),
                r(getattr(atom, 'surf_energy', 2.0 * atom.int_mean_curv_sq)),
                r(atom.sphericity),
                r(atom.isometric_quotient),
                atom.ball_inside,
                atom.number_of_neighbors,
                atom.nearest_neighbor,
                atom.nearest_neighbor_distance,
                [float(_) for _ in r(atom.neighbor_distance_average)],
                [float(_) for _ in r(atom.neighbor_distance_rmsd)],
                r(atom.min_spike),
                r(atom.max_spike),
                atom.number_of_olaps,
                r(atom.contact_area),
                r(atom.vdw_vol),
                r(atom.olap_vol),
                [float(r(_)) for _ in atom.com],
                [[float(r(__)) for __ in _] for _ in atom.moi],
                [
                    [float(r(_)) for _ in atom.bounding_box[0]],
                    [float(r(_)) for _ in atom.bounding_box[1]]
                ],
                nbrs
            ])

        # ==============================================================
        # Surfaces
        # ==============================================================

        lg_fl.writerow(["Surfaces"])
        lg_fl.writerow(["Index", "Ball 1", "Ball 2", "Surface Area", "Mean Curvature", "Average Mean Curvature",
                        "Gaussian Curvature", "Average Gaussian Curvature", "Integrated Mean Curvature",
                        "Integrated Mean Curvature Squared", "Integrated Gaussian Curvature",
                        "Representative Surface Energy", "Ball 1 Volume Contribution", "Ball 2 Volume Contribution",
                        "Contact Area", "Overlap"])

        for surf in net.surfs.itertuples(index=True, name='SurfRow'):
            ball1, ball2 = surf.balls
            vols = surf.vols

            lg_fl.writerow([surf.Index, ball1, ball2, r(surf.sa), r(surf.mean_curv), r(surf.avg_mean_curv),
                            r(surf.gauss_curv), r(surf.avg_gauss_curv), r(surf.int_mean_curv), r(surf.int_mean_curv_sq),
                            r(surf.int_gauss_curv), r(getattr(surf, 'surf_energy', 2.0 * surf.int_mean_curv_sq)),
                            r(vols[ball1]), r(vols[ball2]), r(surf.contact_area), r(surf.overlap)])

        # ==============================================================
        # Edges
        # ==============================================================

        lg_fl.writerow(["Edges"])
        lg_fl.writerow(["Index", "Ball 1", "Ball 2", "Ball 3", "Length"])

        for edge in net.edges.itertuples(index=True, name='EdgeRow'):
            balls = edge.balls

            lg_fl.writerow([
                edge.Index,
                balls[0],
                balls[1],
                balls[2],
                r(edge.length)
            ])

        # ==============================================================
        # Vertices
        # ==============================================================

        lg_fl.writerow(["Vertices"])
        lg_fl.writerow([
            "Index", "Ball 1", "Ball 2", "Ball 3", "Ball 4",
            "x", "y", "z", "r"
        ])

        for vert in net.verts.itertuples(index=True, name='VertRow'):
            balls = vert.balls
            loc = r(vert.loc)

            lg_fl.writerow([
                vert.Index,
                balls[0],
                balls[1],
                balls[2],
                balls[3],
                loc[0],
                loc[1],
                loc[2],
                r(vert.rad)
            ])


def write_interface_logs(iface, net_name=None, round_to=None):
    _debug_interface_logs = os.environ.get("VORPY_INTERFACE_LOG_DEBUG", "0") == "1"
    import time
    _log_t_total = time.perf_counter()
    """
    Export an Interface network using the standard vorPy log format.

    The generated CSV preserves the same row ordering, section labels, column
    ordering, and topology schemas used by the standard Group log writer. This
    allows the existing read_logs() function to parse interface logs without
    modification.

    Parameters
    ----------
    iface : Interface
        Interface containing the network to export.

    net_name : str, optional
        Optional suffix added to the output filename.

    round_to : int, optional
        Number of decimal places used for numerical output. Defaults to the
        system's round_to setting.

    Returns
    -------
    str or None
        Path to the exported log file, or None if no network is available.
    """

    # ------------------------------------------------------------------
    # Basic validation
    # ------------------------------------------------------------------

    if iface is None:
        print("Cannot export interface logs: interface is None.")
        return None

    net = getattr(iface, "net", None)

    if net is None:
        print(
            f'Cannot export interface logs for '
            f'"{getattr(iface, "name", "unnamed interface")}": '
            f"the interface has no network."
        )
        return None

    sys = getattr(iface, "sys", None)

    if sys is None:
        print(
            f'Cannot export interface logs for '
            f'"{getattr(iface, "name", "unnamed interface")}": '
            f"the interface has no system."
        )
        return None

    if round_to is None:
        round_to = getattr(sys, "round_to", 3)

    r = round_func(round_to)

    # ------------------------------------------------------------------
    # Small defensive helpers
    # ------------------------------------------------------------------

    def safe_round(value, default=0.0):
        """
        Round a value while safely handling None, NaN, and incompatible types.
        """
        if value is None:
            value = default

        try:
            return r(value)
        except (TypeError, ValueError, KeyError):
            try:
                return r(default)
            except (TypeError, ValueError):
                return default

    def safe_float(value, default=0.0):
        """
        Convert a value to float without allowing missing data to stop export.
        """
        if value is None:
            return float(default)

        try:
            return float(value)
        except (TypeError, ValueError):
            return float(default)

    def safe_list(value, default=None):
        """
        Return a normal Python list suitable for CSV serialization.
        """
        if default is None:
            default = []

        if value is None:
            return default

        if isinstance(value, np.ndarray):
            return value.tolist()

        if isinstance(value, tuple):
            return list(value)

        if isinstance(value, list):
            return value

        try:
            return list(value)
        except TypeError:
            return default

    def safe_nested_list(value, default=None):
        """
        Convert nested array-like values to ordinary Python lists.
        """
        if default is None:
            default = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        if value is None:
            return default

        if isinstance(value, np.ndarray):
            return value.tolist()

        try:
            return [safe_list(row, default=[0.0, 0.0, 0.0]) for row in value]
        except TypeError:
            return default

    def get_topology_index(atom_row, dataframe_index):
        """
        Return the atom identifier used by surfaces, edges, and vertices.

        net.balls['num'] should normally use the same index space as the
        topology tables. The DataFrame index is retained as a fallback.
        """
        if "num" in atom_row and atom_row["num"] is not None:
            try:
                return int(atom_row["num"])
            except (TypeError, ValueError):
                pass

        return int(dataframe_index)

    def get_system_index(atom_row, topology_index):
        """
        Return the corresponding index into sys.balls.

        system_num is used for metadata lookup when available. It is not
        automatically used as the topology identifier.
        """
        if "system_num" in atom_row and atom_row["system_num"] is not None:
            try:
                return int(atom_row["system_num"])
            except (TypeError, ValueError):
                pass

        return int(topology_index)

    def get_metric(name, default=0.0):
        metrics = getattr(net, "metrics", None)

        if metrics is None:
            return default

        try:
            return metrics[name]
        except (KeyError, TypeError):
            return default

    def get_setting(name, default=None):
        settings = getattr(net, "settings", None)

        if settings is None:
            return default

        try:
            return settings[name]
        except (KeyError, TypeError):
            return default

    # ------------------------------------------------------------------
    # Output directory and filename
    # ------------------------------------------------------------------

    interface_name = getattr(iface, "name", None)

    if interface_name is None:
        group1_name = getattr(getattr(iface, "group1", None), "name", "group1")
        group2_name = getattr(getattr(iface, "group2", None), "name", "group2")
        interface_name = f"{group1_name}_{group2_name}_interface"

    interface_directory = getattr(iface, "dir", None)

    if interface_directory is None:
        system_directory = sys.files["dir"]
        interface_directory = os.path.join(system_directory, interface_name)
        iface.dir = interface_directory

    os.makedirs(interface_directory, exist_ok=True)

    network_type = get_setting("net_type", "network")

    filename_suffix = ""
    if net_name is not None:
        filename_suffix = f"_{net_name}"

    log_filename = f"{network_type}{filename_suffix}_logs.csv"
    log_path = os.path.join(interface_directory, log_filename)

    # ------------------------------------------------------------------
    # Determine network contents
    # ------------------------------------------------------------------

    balls = getattr(net, "balls", None)
    surfs = getattr(net, "surfs", None)
    edges = getattr(net, "edges", None)
    verts = getattr(net, "verts", None)

    num_balls = 0 if balls is None else len(balls)
    num_surfs = 0 if surfs is None else len(surfs)
    num_edges = 0 if edges is None else len(edges)
    num_verts = 0 if verts is None else len(verts)

    if num_verts > 0 and "rad" in verts.columns:
        max_found_vertex = safe_round(verts["rad"].max())
    else:
        max_found_vertex = 0.0

    # ------------------------------------------------------------------
    # Gather interface atom mappings and aggregate information
    # ------------------------------------------------------------------

    _log_t_mapping = time.perf_counter()

    interface_atom_records = []
    system_atom_cache = {}

    if balls is not None and num_balls > 0:
        dataframe_indices = list(balls.index)

        if "num" in balls.columns:
            topology_values = balls["num"].to_numpy(copy=False)
        else:
            topology_values = np.asarray(dataframe_indices)

        if "system_num" in balls.columns:
            raw_system_values = balls["system_num"].to_numpy(copy=False)
        else:
            raw_system_values = topology_values

        # Resolve system indices once without materializing every network-ball
        # row as a Python dictionary.
        system_indices = np.empty(num_balls, dtype=int)
        for position, value in enumerate(raw_system_values):
            try:
                system_indices[position] = int(value)
            except (TypeError, ValueError):
                try:
                    system_indices[position] = int(topology_values[position])
                except (TypeError, ValueError):
                    system_indices[position] = int(dataframe_indices[position])

        # Interface mass/COM are vector operations across all network balls.
        # This replaces the old all-row DataFrame->dict conversion.
        try:
            locations = np.asarray(balls["loc"].tolist(), dtype=float)
            parent_masses = np.asarray(sys.balls["mass"].to_numpy(), dtype=float)
            masses = parent_masses[system_indices]
        except (KeyError, IndexError, TypeError, ValueError):
            locations = np.asarray(balls["loc"].tolist(), dtype=float)
            try:
                masses = np.asarray(balls["mass"].to_numpy(), dtype=float)
            except (KeyError, TypeError, ValueError):
                masses = np.zeros(num_balls, dtype=float)

        finite_mass = np.isfinite(masses)
        masses = np.where(finite_mass, masses, 0.0)
        total_mass = float(masses.sum())
        if total_mass > 0.0 and len(locations) == len(masses):
            interface_com = (
                (locations * masses[:, None]).sum(axis=0) / total_mass
            ).tolist()
        else:
            interface_com = [0.0, 0.0, 0.0]

        # Only atoms with nonzero interface surface area are emitted in the
        # atom section. Convert only those rows, and only columns the writer
        # actually consumes.
        if "sa" in balls.columns:
            sa_values = pd.to_numeric(balls["sa"], errors="coerce").fillna(0.0).to_numpy()
            active_positions = np.flatnonzero(sa_values != 0.0)
        else:
            active_positions = np.arange(num_balls, dtype=int)

        atom_columns = [
            "loc", "rad", "mass", "sa", "complete", "bounding_box",
            "com", "moi", "vol", "vdw_vol", "max_mean_curv",
            "avg_mean_surf_curv", "max_gauss_curv", "avg_gauss_surf_curv",
            "int_mean_curv", "int_mean_curv_sq", "int_gauss_curv",
            "surf_energy", "sphericity", "isometric_quotient", "ball_inside",
            "number_of_neighbors", "nearest_neighbor",
            "nearest_neighbor_distance", "neighbor_distance_average",
            "neighbor_distance_rmsd", "min_spike", "max_spike",
            "number_of_olaps", "contact_area", "olap_vol",
        ]
        atom_columns = [column for column in atom_columns if column in balls.columns]

        if len(active_positions) > 0:
            active_records = balls.iloc[active_positions][atom_columns].to_dict(orient="records")
            active_system_indices = sorted(set(int(system_indices[pos]) for pos in active_positions))

            try:
                selected_system_rows = sys.balls.iloc[active_system_indices]
                selected_columns = [
                    column for column in
                    ("name", "res_name", "res_seq", "chain_name", "mass")
                    if column in selected_system_rows.columns
                ]
                selected_records = selected_system_rows[selected_columns].to_dict(orient="records")
                system_atom_cache = {
                    int(system_index): record
                    for system_index, record in zip(active_system_indices, selected_records)
                }
            except (IndexError, KeyError, TypeError, ValueError):
                system_atom_cache = {}

            for position, atom in zip(active_positions, active_records):
                dataframe_index = dataframe_indices[position]
                try:
                    topology_index = int(topology_values[position])
                except (TypeError, ValueError):
                    topology_index = int(dataframe_index)

                interface_atom_records.append({
                    "dataframe_index": dataframe_index,
                    "topology_index": topology_index,
                    "system_index": int(system_indices[position]),
                    "atom": atom,
                })
    else:
        total_mass = 0.0
        interface_com = [0.0, 0.0, 0.0]

    if num_surfs > 0 and "sa" in surfs.columns:
        interface_surface_area = safe_float(surfs["sa"].sum())
    else:
        interface_surface_area = 0.0

    _log_mapping_seconds = time.perf_counter() - _log_t_mapping

    # An interface-only network does not necessarily define a closed volume.
    undefined_scalar = float("nan")

    zero_tensor = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]

    # ------------------------------------------------------------------
    # Precompute lookup tables used by the atom section.
    #
    # The legacy implementation performed a pandas surface lookup for every
    # atom/surface relationship.  On large interface networks that dominates
    # export time.  Build the topology neighbor map once instead.
    # ------------------------------------------------------------------

    _log_t_prepare = time.perf_counter()

    neighbor_map = {}
    if surfs is not None and num_surfs > 0:
        # Pull the small topology column out of pandas once and build the
        # adjacency map using native Python containers.
        for surface_balls in surfs["balls"].tolist():
            surface_balls = safe_list(surface_balls, default=[])
            if len(surface_balls) < 2:
                continue
            try:
                ball1 = int(surface_balls[0])
                ball2 = int(surface_balls[1])
            except (TypeError, ValueError):
                continue
            neighbor_map.setdefault(ball1, set()).add(ball2)
            neighbor_map.setdefault(ball2, set()).add(ball1)

    neighbor_map = {
        atom_index: sorted(neighbors)
        for atom_index, neighbors in neighbor_map.items()
    }

    _log_prepare_seconds = time.perf_counter() - _log_t_prepare

    # ------------------------------------------------------------------
    # Write the log
    # ------------------------------------------------------------------

    with open(log_path, "w", newline="", encoding="utf-8") as log_file:

        lg_fl = csv.writer(log_file, lineterminator="\n")

        # ==============================================================
        # Build information
        #
        # These must remain the first three rows because read_logs()
        # accesses the build-data row by its absolute row number.
        # ==============================================================

        lg_fl.writerow(["build information"])

        lg_fl.writerow(["Name", "Location", "Completion Date", "Network Type", "Surface Resolution", "Box Size",
                        "Maximum Allowable Vertex", "Total Time", "Vertex Time", "Connect Time",
                        "Surface Building Time", "Analysis time", "Maximum Found Vertex", "vorPy version"])

        lg_fl.writerow([interface_name, sys.files.get("base_file", ""), datetime.now(), network_type,
                        get_setting("surf_res", 0.0), get_setting("box_size", 0.0),
                        get_setting("max_vert", 0.0), safe_round(get_metric("tot", 0.0)),
                        safe_round(get_metric("vert", 0.0)), safe_round(get_metric("con", 0.0)),
                        safe_round(get_metric("surf", 0.0)), safe_round(get_metric("anal", 0.0)),
                        max_found_vertex, __version__])

        # ==============================================================
        # Group information compatibility section
        #
        # The section must retain this name and contain ten values because
        # read_logs() parses row 5 positionally as group information.
        # ==============================================================

        lg_fl.writerow(["group information"])

        lg_fl.writerow(["Name", "Volume", "Surface Area", "Mass", "Density", "Center of Mass", "VDW Volume",
                        "VDW Center of Mass", "Moment of Inertia", "Spatial Moment of Inertia"])

        lg_fl.writerow([f"INTERFACE:{interface_name}", undefined_scalar, safe_round(interface_surface_area),
                        safe_float(total_mass), undefined_scalar, [safe_float(value) for value in interface_com],
                        undefined_scalar, [safe_float(value) for value in interface_com], zero_tensor, zero_tensor])

        # ==============================================================
        # Atoms
        # ==============================================================
        # Write the atom header
        lg_fl.writerow(["Atoms"])
        # Write the column labels
        lg_fl.writerow(["Index", "Name", "Residue", "Residue Sequence", "Chain", "Mass", "X", "Y", "Z", "Radius",
                        "Volume", "Van Der Waals Volume", "Surface Area", "Complete Cell?",
                        "Maximum Mean Curvature", "Average Mean Surface Curvature", "Maximum Gaussian Curvature",
                        "Average Gaussian Surface Curvature", "Integrated Mean Curvature",
                        "Integrated Mean Curvature Squared", "Integrated Gaussian Curvature",
                        "Representative Surface Energy", "Sphericity", "Isometric Quotient", "Inner Ball?",
                        "Number of Neighbors", "Closest Neighbor",
                        "Closest Neighbor Distance", "Layer Distance Average", "Layer Distance RMSD",
                        "Minimum Point Distance", "Maximum Point Distance", "Number of Overlaps", "Contact Area",
                        "Non-Overlap Volume", "Overlap Volume", "Center of Mass", "Moment of Inertia Tensor",
                        "Bounding Box", "neighbors"])

        _log_t_atoms = time.perf_counter()
        _log_atom_rows = 0

        for record in interface_atom_records:
            atom = record["atom"]
            topology_index = record["topology_index"]
            system_index = record["system_index"]

            # Skip non-participating atoms before any metadata or topology work.
            atom_surface_area = safe_float(atom.get("sa", 0.0))
            if atom_surface_area == 0:
                continue

            sys_atom = system_atom_cache.get(int(system_index))

            if sys_atom is None:
                atom_name = str(atom.get("name", ""))
                residue_name = str(atom.get("res_name", ""))
                residue_sequence = int(atom.get("res_seq", 0))
                chain_name = str(atom.get("chain_name", ""))
                atom_mass = safe_float(atom.get("mass", 0.0))
            else:
                atom_name = str(sys_atom.get("name", ""))
                residue_name = str(sys_atom.get("res_name", ""))
                residue_sequence = int(sys_atom.get("res_seq", 0))
                chain_name = str(sys_atom.get("chain_name", ""))
                atom_mass = safe_float(sys_atom.get("mass", 0.0))

            location = safe_list(
                atom.get("loc", [0.0, 0.0, 0.0]),
                default=[0.0, 0.0, 0.0],
            )

            if len(location) < 3:
                location += [0.0] * (3 - len(location))

            # Neighbor relationships were precomputed from surfaces once.
            neighbors = neighbor_map.get(int(topology_index), [])

            complete = bool(atom.get("complete", False))

            bounding_box = safe_nested_list(atom.get("bounding_box"), default=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

            atom_com = safe_list(atom.get("com", location[:3]), default=location[:3])

            atom_moi = safe_nested_list(atom.get("moi"), default=zero_tensor)

            lg_fl.writerow([
                topology_index,
                atom_name,
                residue_name,
                residue_sequence,
                chain_name,
                atom_mass,
                safe_float(location[0]),
                safe_float(location[1]),
                safe_float(location[2]),
                safe_float(atom.get("rad", 0.0)),
                safe_round(atom.get("vol", 0.0)),
                safe_round(atom.get("vdw_vol", 0.0)),
                safe_round(atom_surface_area),
                complete,
                safe_round(atom.get("max_mean_curv", 0.0)),
                safe_round(atom.get("avg_mean_surf_curv", 0.0)),
                safe_round(atom.get("max_gauss_curv", 0.0)),
                safe_round(atom.get("avg_gauss_surf_curv", 0.0)),
                safe_round(atom.get("int_mean_curv", 0.0)),
                safe_round(atom.get("int_mean_curv_sq", 0.0)),
                safe_round(atom.get("int_gauss_curv", 0.0)),
                safe_round(atom.get(
                    "surf_energy",
                    2.0 * safe_float(atom.get("int_mean_curv_sq", 0.0))
                )),
                safe_round(atom.get("sphericity", 0.0)),
                safe_round(atom.get("isometric_quotient", 0.0)),
                bool(atom.get("ball_inside", False)),
                int(atom.get("number_of_neighbors", len(neighbors))),
                int(atom.get("nearest_neighbor", -1)),
                safe_round(
                    atom.get("nearest_neighbor_distance", 0.0)
                ),
                safe_list(
                    atom.get("neighbor_distance_average", []),
                    default=[],
                ),
                safe_list(
                    atom.get("neighbor_distance_rmsd", []),
                    default=[],
                ),
                safe_round(atom.get("min_spike", 0.0)),
                safe_round(atom.get("max_spike", 0.0)),
                int(atom.get("number_of_olaps", 0)),
                safe_round(atom.get("contact_area", 0.0)),

                # These two positions deliberately preserve the current
                # Group log writer's output ordering.
                safe_round(atom.get("olap_vol", 0.0)),
                safe_round(atom.get("vdw_vol", 0.0)),

                [safe_round(value) for value in atom_com],
                [
                    [safe_round(value) for value in row]
                    for row in atom_moi
                ],
                [
                    [safe_round(value) for value in row]
                    for row in bounding_box
                ],
                neighbors,
            ]
            )
            _log_atom_rows += 1

        _log_atom_seconds = time.perf_counter() - _log_t_atoms
        if _debug_interface_logs:
            print(f"[INTERFACE LOG] atoms: {_log_atom_rows} rows in {_log_atom_seconds:.3f} s")

        # ==============================================================
        # Surfaces
        # ==============================================================

        lg_fl.writerow(["Surfaces"])

        lg_fl.writerow(
            [
                "Index",
                "Ball 1",
                "Ball 2",
                "Surface Area",
                "Mean Curvature",
                "Average Mean Curvature",
                "Gaussian Curvature",
                "Average Gaussian Curvature",
                "Integrated Mean Curvature",
                "Integrated Mean Curvature Squared",
                "Integrated Gaussian Curvature",
                "Representative Surface Energy",
                "Ball 1 Volume Contribution",
                "Ball 2 Volume Contribution",
                "Contact Area",
                "Overlap",
            ]
        )

        _log_t_surfaces = time.perf_counter()
        _log_surface_rows = 0
        if surfs is not None:
            for surface_index, surf in surfs.iterrows():
                surface_balls = safe_list(
                    surf.get("balls", []),
                    default=[],
                )

                if len(surface_balls) < 2:
                    continue

                ball1 = int(surface_balls[0])
                ball2 = int(surface_balls[1])

                surface_volumes = surf.get("vols", {})

                try:
                    ball1_volume = surface_volumes[ball1]
                except (KeyError, TypeError, IndexError):
                    try:
                        ball1_volume = surface_volumes[str(ball1)]
                    except (KeyError, TypeError, IndexError):
                        ball1_volume = 0.0

                try:
                    ball2_volume = surface_volumes[ball2]
                except (KeyError, TypeError, IndexError):
                    try:
                        ball2_volume = surface_volumes[str(ball2)]
                    except (KeyError, TypeError, IndexError):
                        ball2_volume = 0.0

                lg_fl.writerow(
                    [
                        int(surface_index),
                        ball1,
                        ball2,
                        safe_round(surf.get("sa", 0.0)),
                        safe_round(surf.get("mean_curv", 0.0)),
                        safe_round(surf.get("avg_mean_curv", 0.0)),
                        safe_round(surf.get("gauss_curv", 0.0)),
                        safe_round(surf.get("avg_gauss_curv", 0.0)),
                        safe_round(surf.get("int_mean_curv", 0.0)),
                        safe_round(surf.get("int_mean_curv_sq", 0.0)),
                        safe_round(surf.get("int_gauss_curv", 0.0)),
                        safe_round(surf.get(
                            "surf_energy",
                            2.0 * safe_float(surf.get("int_mean_curv_sq", 0.0))
                        )),
                        safe_round(ball1_volume),
                        safe_round(ball2_volume),
                        safe_round(surf.get("contact_area", 0.0)),
                        safe_round(surf.get("overlap", 0.0)),
                    ]
                )
                _log_surface_rows += 1

        _log_surface_seconds = time.perf_counter() - _log_t_surfaces
        if _debug_interface_logs:
            print(f"[INTERFACE LOG] surfaces: {_log_surface_rows} rows in {_log_surface_seconds:.3f} s")

        # ==============================================================
        # Edges
        # ==============================================================

        lg_fl.writerow(["Edges"])

        lg_fl.writerow(
            [
                "Index",
                "Ball 1",
                "Ball 2",
                "Ball 3",
                "Length",
            ]
        )

        _log_t_edges = time.perf_counter()
        _log_edge_rows = 0
        if edges is not None:
            for edge_index, edge in edges.iterrows():
                edge_balls = safe_list(
                    edge.get("balls", []),
                    default=[],
                )

                if len(edge_balls) < 3:
                    continue

                lg_fl.writerow(
                    [
                        int(edge_index),
                        int(edge_balls[0]),
                        int(edge_balls[1]),
                        int(edge_balls[2]),
                        safe_round(edge.get("length", 0.0)),
                    ]
                )
                _log_edge_rows += 1

        _log_edge_seconds = time.perf_counter() - _log_t_edges
        if _debug_interface_logs:
            print(f"[INTERFACE LOG] edges: {_log_edge_rows} rows in {_log_edge_seconds:.3f} s")

        # ==============================================================
        # Vertices
        # ==============================================================

        lg_fl.writerow(["Vertices"])

        lg_fl.writerow(
            [
                "Index",
                "Ball 1",
                "Ball 2",
                "Ball 3",
                "Ball 4",
                "x",
                "y",
                "z",
                "r",
            ]
        )

        _log_t_vertices = time.perf_counter()
        _log_vertex_rows = 0
        if verts is not None:
            for vertex_index, vert in verts.iterrows():
                vertex_balls = safe_list(
                    vert.get("balls", []),
                    default=[],
                )

                vertex_location = safe_list(
                    vert.get("loc", [0.0, 0.0, 0.0]),
                    default=[0.0, 0.0, 0.0],
                )

                if len(vertex_balls) < 4:
                    continue

                if len(vertex_location) < 3:
                    vertex_location += [0.0] * (
                            3 - len(vertex_location)
                    )

                lg_fl.writerow(
                    [
                        int(vertex_index),
                        int(vertex_balls[0]),
                        int(vertex_balls[1]),
                        int(vertex_balls[2]),
                        int(vertex_balls[3]),
                        safe_round(vertex_location[0]),
                        safe_round(vertex_location[1]),
                        safe_round(vertex_location[2]),
                        safe_round(vert.get("rad", 0.0)),
                    ]
                )
                _log_vertex_rows += 1

        _log_vertex_seconds = time.perf_counter() - _log_t_vertices
        if _debug_interface_logs:
            print(f"[INTERFACE LOG] vertices: {_log_vertex_rows} rows in {_log_vertex_seconds:.3f} s")

    _log_total_seconds = time.perf_counter() - _log_t_total
    if _debug_interface_logs:
        print("\n" + "=" * 72)
        print("INTERFACE LOG EXPORT TIMING")
        print("=" * 72)
        print(f"Atom mapping / metadata     {_log_mapping_seconds:10.3f} s")
        print(f"Preparation / caches        {_log_prepare_seconds:10.3f} s")
        print(f"Atoms                       {_log_atom_seconds:10.3f} s")
        print(f"Surfaces                    {_log_surface_seconds:10.3f} s")
        print(f"Edges                       {_log_edge_seconds:10.3f} s")
        print(f"Vertices                    {_log_vertex_seconds:10.3f} s")
        print("-" * 72)
        print(f"TOTAL                       {_log_total_seconds:10.3f} s")
        print("=" * 72 + "\n")

    return log_path
