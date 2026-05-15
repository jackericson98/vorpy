import os.path
import os
import time
import pandas as pd
from vorpy.src.inputs.logs import read_logs
import numpy as np
from vorpy.src.calculations.vert import calc_vert, calc_flat_vert
from vorpy.src.calculations.calcs import get_time


def _rename_log_columns(df):
    rename_map = {
        "Index": "num",
        "Name": "name",
        "Mass": "mass",
        "Radius": "rad",
        "Volume": "vol",
        "Surface Area": "sa",
        "Complete Cell?": "complete",
        "Balls": "balls",
        "Length": "length",
        "Ball Volumes": "ball_vols",
        "Mean Curvature": "mean_curv",
        "Average Mean Curvature": "avg_mean_curv",
        "Gauss Curvature": "gauss_curv",
        "Average Gauss Curvature": "avg_gauss_curv",
        "Contact Area": "contact_area",
        "Overlap": "overlap",
    }

    return df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})


def _standardize_log_geometry_columns(balls, verts, edges, surfs):
    balls = balls.copy()
    verts = verts.copy()
    edges = edges.copy()
    surfs = surfs.copy()

    if "loc" not in balls.columns and {"X", "Y", "Z"}.issubset(balls.columns):
        balls["loc"] = balls[["X", "Y", "Z"]].values.tolist()

    if "rad" not in balls.columns and "Radius" in balls.columns:
        balls["rad"] = balls["Radius"]

    if "balls" not in verts.columns:
        vert_ball_cols = ["Ball 1", "Ball 2", "Ball 3", "Ball 4"]

        if all(col in verts.columns for col in vert_ball_cols):
            verts["balls"] = verts[vert_ball_cols].astype(int).values.tolist()

    if "loc" not in verts.columns and {"x", "y", "z"}.issubset(verts.columns):
        verts["loc"] = verts[["x", "y", "z"]].values.tolist()

    if "rad" not in verts.columns and "r" in verts.columns:
        verts["rad"] = verts["r"]

    if "balls" not in edges.columns:
        edge_ball_cols = ["Ball 1", "Ball 2", "Ball 3"]

        if all(col in edges.columns for col in edge_ball_cols):
            edges["balls"] = edges[edge_ball_cols].astype(int).values.tolist()

    if "balls" not in surfs.columns:
        surf_ball_cols = ["Ball 1", "Ball 2"]

        if all(col in surfs.columns for col in surf_ball_cols):
            surfs["balls"] = surfs[surf_ball_cols].astype(int).values.tolist()

    balls = balls.loc[:, ~balls.columns.duplicated()].copy()
    verts = verts.loc[:, ~verts.columns.duplicated()].copy()
    edges = edges.loc[:, ~edges.columns.duplicated()].copy()
    surfs = surfs.loc[:, ~surfs.columns.duplicated()].copy()

    return balls, verts, edges, surfs


def _solution_is_sane(new_loc, new_rad, max_abs_rad=10.0, allow_negative_rad=True):
    # print("\n=== SANITY CHECK DEBUG ===")
    # print(f"new_loc = {new_loc}")
    # print(f"new_rad = {new_rad}")
    # print(f"max_abs_rad = {max_abs_rad}")
    # print(f"allow_negative_rad = {allow_negative_rad}")

    if new_loc is None or new_rad is None:
        print("FAIL: loc or rad is None")
        return False

    new_loc = np.array(new_loc, dtype=float)
    new_rad = float(new_rad)

    # print(f"finite loc = {np.all(np.isfinite(new_loc))}")
    # print(f"finite rad = {np.isfinite(new_rad)}")
    # print(f"abs rad = {abs(new_rad)}")

    if not np.all(np.isfinite(new_loc)):
        # print("FAIL: loc not finite")
        return False

    if not np.isfinite(new_rad):
        # print("FAIL: rad not finite")
        return False

    if allow_negative_rad:
        if abs(new_rad) > max_abs_rad:
            # print("FAIL: abs(rad) > max_abs_rad")
            return False

    else:
        if new_rad < 0:
            # print("FAIL: negative rad not allowed")
            return False

        if new_rad > max_abs_rad:
            # print("FAIL: rad > max_abs_rad")
            return False

    # print("PASS")
    return True


def _recalculate_loaded_vertices(net, verts):
    """
    Recalculate loaded vertex loc/rad from saved 4-ball definitions.

    Important:
    - Uses current system ball coordinates/radii.
    - Preserves doublet pairs by using both calc_vert solutions.
    - Rejects recalculated solutions that are far from saved log coordinates.
    """
    if verts is None or len(verts) == 0:
        return verts

    if "balls" not in verts.columns:
        return verts

    verts = verts.copy()

    if "dub" not in verts.columns:
        from vorpy.src.network.net_logs_connect import get_dubs
        verts["dub"] = get_dubs(verts)

    net_type = net.settings.get("net_type", "aw") if net.settings else "aw"
    allow_negative_rad = net_type in {"aw", "pow"}
    use_flat = net_type in {"prm", "pow"}
    use_power = net_type == "pow"

    # print("\n=== VERTEX MODE DEBUG ===")
    # print(f"net_type = {net_type}")
    # print(f"use_flat = {use_flat}")
    # print(f"use_power = {use_power}")
    # print(f"allow_negative_rad = {allow_negative_rad}")

    start_time = time.perf_counter()
    total_verts = len(verts)
    processed_verts = 0

    new_locs = list(verts["loc"])
    new_rads = list(verts["rad"])

    grouped = {}

    # print(f"\n\n\n\nUse Flat = {use_flat}\n\n\n\n")

    for vert_i, vert in verts.iterrows():
        key = tuple(int(_) for _ in vert["balls"])
        grouped.setdefault(key, []).append(vert_i)

    failed = []
    failed_details = []

    for vballs_key, vert_indices in grouped.items():
        vballs = list(vballs_key)

        if len(vballs) != 4:
            failed.extend(vert_indices)

            for vi in vert_indices:
                failed_details.append((vi, "invalid 4-ball vertex definition", vballs))

            processed_verts += len(vert_indices)
            continue

        locs = [
            np.array(net.balls.iloc[ball_i]["loc"], dtype=float)
            for ball_i in vballs
        ]

        rads = [
            float(net.balls.iloc[ball_i]["rad"])
            for ball_i in vballs
        ]

        try:
            if use_flat:
                loc, rad = calc_flat_vert(locs, rads, power=use_power)
                loc2, rad2 = None, None
            else:
                loc, rad, loc2, rad2 = calc_vert(locs, rads)

            if loc is None or rad is None:
                failed.extend(vert_indices)
                processed_verts += len(vert_indices)
                continue

            # Normal / single vertex case
            if len(vert_indices) == 1:
                vi = vert_indices[0]

                saved_loc = np.array(verts.iloc[vi]["loc"], dtype=float)
                saved_rad = float(verts.iloc[vi]["rad"])

                candidates = []

                if loc is not None and rad is not None:
                    candidates.append((np.array(loc, dtype=float), float(rad)))

                if loc2 is not None and rad2 is not None:
                    candidates.append((np.array(loc2, dtype=float), float(rad2)))

                if not candidates:
                    new_locs[vi] = saved_loc
                    new_rads[vi] = saved_rad
                    failed.append(vi)
                    processed_verts += len(vert_indices)
                    continue

                best_loc, best_rad = min(
                    candidates,
                    key=lambda x: np.linalg.norm(saved_loc - x[0]) + abs(saved_rad - x[1])
                )

                if _solution_is_sane(best_loc, best_rad, allow_negative_rad=allow_negative_rad):
                    new_locs[vi] = best_loc
                    new_rads[vi] = best_rad
                else:
                    new_locs[vi] = saved_loc
                    new_rads[vi] = saved_rad
                    failed.append(vi)
                    failed_details.append((vi, "best candidate failed sanity check", vballs, candidates))

                processed_verts += len(vert_indices)

            # Doublet pair case
            elif len(vert_indices) == 2:
                vi0, vi1 = vert_indices

                saved0 = np.array(verts.iloc[vi0]["loc"], dtype=float)
                saved1 = np.array(verts.iloc[vi1]["loc"], dtype=float)

                sol0 = np.array(loc, dtype=float)
                sol1 = np.array(loc2, dtype=float)

                assignment_a = (
                    np.linalg.norm(saved0 - sol0) +
                    np.linalg.norm(saved1 - sol1)
                )

                assignment_b = (
                    np.linalg.norm(saved0 - sol1) +
                    np.linalg.norm(saved1 - sol0)
                )

                if assignment_a <= assignment_b:
                    proposed = [
                        (vi0, sol0, float(rad)),
                        (vi1, sol1, float(rad2)),
                    ]
                else:
                    proposed = [
                        (vi0, sol1, float(rad2)),
                        (vi1, sol0, float(rad)),
                    ]

                for vi, proposed_loc, proposed_rad in proposed:
                    saved_loc = verts.iloc[vi]["loc"]
                    saved_rad = verts.iloc[vi]["rad"]

                    if _solution_is_sane(proposed_loc, proposed_rad, allow_negative_rad=allow_negative_rad):
                        new_locs[vi] = np.array(proposed_loc, dtype=float)
                        new_rads[vi] = float(proposed_rad)
                    else:
                        new_locs[vi] = np.array(saved_loc, dtype=float)
                        new_rads[vi] = float(saved_rad)
                        failed.append(vi)
                        failed_details.append(
                            (
                                vi,
                                "doublet candidate failed sanity check",
                                vballs,
                                [(loc, rad), (loc2, rad2)]
                            )
                        )

                processed_verts += len(vert_indices)

            # Unexpected duplicate count
            else:
                failed.extend(vert_indices)
                processed_verts += len(vert_indices)

        except Exception:
            failed.extend(vert_indices)
            processed_verts += len(vert_indices)
            if len(failed) < 10:
                print("\n=== VERT RECALC FAILURE ===")
                print(f"vert_indices = {vert_indices}")
                print(f"vballs = {vballs}")
                print(f"saved loc/rad:")
                for vi in vert_indices:
                    print(
                        f"  vert {vi}: loc={verts.iloc[vi]['loc']} rad={verts.iloc[vi]['rad']} balls={verts.iloc[vi]['balls']}")

                print("system balls used:")
                for bi in vballs:
                    b = net.balls.iloc[bi]
                    print(
                        f"  ball {bi}: name={b.get('name', b.get('Name', None))} "
                        f"res={b.get('res_name', b.get('Residue', None))} "
                        f"res_seq={b.get('res_seq', b.get('Residue Sequence', None))} "
                        f"chain={b.get('chain', b.get('Chain', None))} "
                        f"loc={b.get('loc', None)} rad={b.get('rad', b.get('Radius', None))}"
                    )

                print("===========================\n")

        my_time = time.perf_counter() - start_time
        h, m, s = get_time(my_time)
        percentage = 100.0 * processed_verts / total_verts

        print(
            "\rRun Time = {}:{:02d}:{:2.2f} - Process: recalculating vertices: {} verts - {:.2f} %"
            .format(int(h), int(m), round(s, 2), processed_verts, percentage),
            end=""
        )

    verts["loc"] = new_locs
    verts["rad"] = new_rads

    if failed:
        failed_unique = sorted(set(failed))

        print(
            f"WARNING: {len(failed_unique)} loaded vertices could not be safely recalculated "
            f"and kept saved log coordinates."
        )

        print("\n=== FAILED VERTEX INDICES ===")
        print(failed_unique)
        print("=============================\n")

        for detail in failed_details[:50]:
            vi = detail[0]
            reason = detail[1]
            vballs = detail[2] if len(detail) > 2 else None
            candidates = detail[3] if len(detail) > 3 else None

            _print_bad_loaded_vertex(
                net=net,
                verts=verts,
                vi=vi,
                reason=reason,
                candidates=candidates,
                vballs=vballs
            )

    return verts


def _print_bad_loaded_vertex(net, verts, vi, reason, saved_loc=None, saved_rad=None,
                             candidates=None, vballs=None):
    print("\n=== BAD LOADED VERTEX ===")
    print(f"reason = {reason}")
    print(f"vert index = {vi}")

    if vballs is None:
        try:
            vballs = [int(_) for _ in verts.iloc[vi]["balls"]]
        except Exception:
            vballs = None

    print(f"balls = {vballs}")

    if "dub" in verts.columns:
        print(f"dub = {verts.iloc[vi].get('dub', None)}")

    if saved_loc is None:
        saved_loc = verts.iloc[vi].get("loc", None)

    if saved_rad is None:
        saved_rad = verts.iloc[vi].get("rad", None)

    print(f"saved loc = {saved_loc}")
    print(f"saved rad = {saved_rad}")

    if candidates:
        print("candidate solutions:")
        for cand_i, cand in enumerate(candidates):
            cand_loc, cand_rad = cand
            dist = np.linalg.norm(np.array(saved_loc, dtype=float) - np.array(cand_loc, dtype=float))
            rdiff = abs(float(saved_rad) - float(cand_rad))

            print(
                f"  candidate {cand_i}: "
                f"loc={np.array(cand_loc).tolist()} "
                f"rad={cand_rad} "
                f"loc_delta={dist:.6f} "
                f"rad_delta={rdiff:.6f}"
            )

    if vballs is not None:
        print("system balls:")
        for bi in vballs:
            if not (0 <= bi < len(net.balls)):
                print(f"  ball {bi}: OUTSIDE SYSTEM RANGE")
                continue

            b = net.balls.iloc[bi]
            print(
                f"  ball {bi}: "
                f"name={b.get('name', b.get('Name', None))} "
                f"res={b.get('res_name', b.get('Residue', None))} "
                f"res_seq={b.get('res_seq', b.get('Residue Sequence', None))} "
                f"chain={b.get('chain', b.get('Chain', None))} "
                f"loc={b.get('loc', None)} "
                f"rad={b.get('rad', b.get('Radius', None))}"
            )

    print("=========================\n")


def read_net(group, net, file_name, rebuild_edges=True, rebuild_surfs=True, analyze=True,
             store_points=True, index_map=None, index_offset=None):
    """
    Load a vorpy logs file into an existing Network object.

    Important model:
    - net.balls is the full system balls table.
    - net.group is the logged/network subset.
    - topology ball refs are global system indices unless proven otherwise.
    - log_balls only stores complete-cell/network ball metadata.
    """
    if file_name is None or not os.path.exists(file_name):
        raise FileNotFoundError(f"Network log file not found: {file_name}")

    def _get_system_balls(group, net):
        if getattr(group, "sys", None) is not None and getattr(group.sys, "balls", None) is not None:
            return group.sys.balls.copy().reset_index(drop=True)

        if getattr(net, "balls", None) is not None:
            return net.balls.copy().reset_index(drop=True)

        raise ValueError(
            "Cannot load logs because full system balls are unavailable. "
            "Load the PDB/system before calling read_net()."
        )

    def _collect_topology_refs(verts, edges, surfs):
        refs = set()

        for df in (verts, edges, surfs):
            if df is None or "balls" not in df.columns:
                continue

            for ball_list in df["balls"]:
                refs.update(int(_) for _ in ball_list)

        return refs

    def _validate_topology_refs(topology_refs, n_balls):
        bad_refs = sorted(_ for _ in topology_refs if _ < 0 or _ >= n_balls)

        if bad_refs:
            raise ValueError(
                "Loaded topology references balls outside the full system ball table.\n"
                f"System ball count = {n_balls}\n"
                f"Bad ref count = {len(bad_refs)}\n"
                f"First bad refs = {bad_refs[:50]}"
            )

    def _remap_ball_list(ball_list, index_map=None, n_balls=None, index_offset=None):
        remapped = []

        for ball_i in ball_list:
            ball_i = int(ball_i)

            if index_map is not None and ball_i in index_map:
                mapped_i = int(index_map[ball_i])

            elif index_offset is not None:
                mapped_i = ball_i + int(index_offset)

            else:
                mapped_i = ball_i

            if n_balls is not None and not (0 <= mapped_i < n_balls):
                raise KeyError(
                    f"Topology ball index {ball_i} maps to {mapped_i}, "
                    f"outside system range 0-{n_balls - 1}."
                )

            remapped.append(mapped_i)

        return remapped

    log_data = read_logs(file_name, all_=True)

    build_data = log_data["data"]
    group_data = log_data["group data"]

    log_balls = _rename_log_columns(log_data["atoms"])
    verts = _rename_log_columns(log_data["verts"])
    edges = _rename_log_columns(log_data["edges"])
    surfs = _rename_log_columns(log_data["surfs"])

    log_balls, verts, edges, surfs = _standardize_log_geometry_columns(
        log_balls,
        verts,
        edges,
        surfs
    )

    # Critical:
    # Loaded networks must mimic normal-built networks.
    # net.balls must be the full system balls table, not only logged/network balls.
    net.balls = _get_system_balls(group, net)

    n_balls = len(net.balls)

    topology_refs = _collect_topology_refs(verts, edges, surfs)
    logged_refs = set(int(_) for _ in log_balls["num"]) if "num" in log_balls.columns else set()
    missing_topology_refs = sorted(topology_refs - logged_refs)

    # print("\n=== TOPOLOGY COVERAGE DEBUG ===")
    # print(f"system ball count = {n_balls}")
    # print(f"logged/network atom count = {len(logged_refs)}")
    # print(f"topology referenced atoms = {len(topology_refs)}")
    # print(f"topology refs not in logged atoms = {len(missing_topology_refs)}")

    # if missing_topology_refs:
    #     print(f"first 50 topology-only refs = {missing_topology_refs[:50]}")
    #
    # print("================================\n")

    remap_topology = False

    try:
        _validate_topology_refs(topology_refs, n_balls)

        # print(
        #     "\nTopology references are valid against full system balls.\n"
        #     "Treating topology as global system indices.\n"
        #     "Using index_map only for logged/network group assignment.\n"
        # )

    except ValueError:
        if index_map is None and index_offset is None:
            raise

        remap_topology = True

    if remap_topology:
        # print(
        #     "\nTopology references do not fit the system directly.\n"
        #     "Attempting topology remap using index_map/index_offset.\n"
        # )

        verts["balls"] = verts["balls"].apply(
            lambda x: _remap_ball_list(x, index_map, n_balls, index_offset)
        )
        edges["balls"] = edges["balls"].apply(
            lambda x: _remap_ball_list(x, index_map, n_balls, index_offset)
        )
        surfs["balls"] = surfs["balls"].apply(
            lambda x: _remap_ball_list(x, index_map, n_balls, index_offset)
        )

        topology_refs = _collect_topology_refs(verts, edges, surfs)
        _validate_topology_refs(topology_refs, n_balls)

    logged_group_indices = sorted(int(_) for _ in log_balls["num"].dropna().tolist())

    if index_map is not None:
        net.group = sorted(set(int(index_map[int(_)]) for _ in logged_group_indices))
    else:
        net.group = logged_group_indices

    group.ball_ndxs = net.group
    group.loaded_log_ball_indices = logged_group_indices
    net.loaded_log_ball_indices = logged_group_indices

    # Existing system columns initialized as integer zeros may later receive
    # computed float log metrics. These are true numeric metrics, so cast only
    # known continuous columns to float, not object.
    float_metric_cols = {
        "vol",
        "sa",
        "contact_area",
        "mean_curv",
        "avg_mean_curv",
        "gauss_curv",
        "avg_gauss_curv",
    }

    for col in float_metric_cols:
        if col in net.balls.columns:
            net.balls[col] = net.balls[col].astype(float)

    skip_log_ball_cols = {
        "Index",
        "num",
        "loc",
        "rad",
        "name",
        "mass",
        "X",
        "Y",
        "Z",
    }

    for _, log_ball in log_balls.iterrows():
        old_ball_num = int(log_ball["num"])
        new_ball_num = int(index_map[old_ball_num]) if index_map is not None else old_ball_num

        if not (0 <= new_ball_num < n_balls):
            raise ValueError(
                f"Logged ball {old_ball_num} maps to system ball {new_ball_num}, "
                f"outside range 0-{n_balls - 1}."
            )

        for col in log_balls.columns:
            if col in skip_log_ball_cols:
                continue

            if col not in net.balls.columns:
                net.balls[col] = None

            net.balls.at[new_ball_num, col] = log_ball[col]

    if net.settings is None:
        net.default_settings()

    net.settings["net_type"] = build_data.get("network_type", net.settings.get("net_type", "aw"))
    net.settings["surf_res"] = build_data.get("surface_resolution", net.settings.get("surf_res", 0.2))
    net.settings["box_size"] = build_data.get("box_size", net.settings.get("box_size", 1.5))
    net.settings["max_vert"] = build_data.get("max_vert", net.settings.get("max_vert", 40))

    verts = _recalculate_loaded_vertices(net, verts)

    net.verts = verts
    net.edges = edges
    net.surfs = surfs

    net.log_data = build_data
    net.group_data = group_data

    net.original_metrics = {
        "tot": build_data.get("Total_Time", 0),
        "vert": build_data.get("vert_time", 0),
        "con": build_data.get("connect_time", 0),
        "surf": build_data.get("surf_time", 0),
        "anal": build_data.get("analysis_time", 0),
    }

    net.metrics.setdefault("start", time.perf_counter())
    net.metrics["vert"] = 0
    net.metrics["con"] = 0
    net.metrics["surf"] = 0
    net.metrics["anal"] = 0
    net.metrics["tot"] = 0

    net.loaded_from_logs = True
    net.rebuilt_from_logs = True

    net.logs_connect()

    if rebuild_edges:
        net.build_edges()

    if rebuild_surfs:
        net.build_surfaces(store_points=store_points)

    if analyze:
        net.analyze()

    net.loaded_from_logs = True

    return net


# Input index function. Takes in an index file and loads it into the list of indices
def read_ndx(sys, file=None):
    """
    Read and process an index file into a system object.

    This function parses index files, which contain atom group definitions and metadata.
    The function converts the index data into a standardized format for further processing.

    The index format includes:
    - Group names in square brackets
    - Atom indices for each group
    - Multiple groups can be defined

    Parameters:
    -----------
    sys : System
        The system object to populate with index data
    file : str, optional
        Path to the index file. If None, uses sys.ndx_file

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom groups from index data
        - Storing group names in sys.ndx_names
        - Storing atom indices in sys.ndxs
    """
    # If no file is provided, check the system
    if file is None:
        file = sys.ndx_file
    # Get the file information and make sure to close the file when done
    try:
        with open(file, 'r') as f:
            my_file = f.readlines()
    except FileNotFoundError:
        return
    # Set up the indices lists and the current index
    curr_ndx = -1
    indices = []
    names = []
    # Go through the lines in the file
    for line in my_file:
        # Split the line into
        line = line.split()
        # Add the
        if line[0] == "[":
            curr_ndx += 1
            names.append([line[1]])
        else:
            for i in range(len(line)):
                indices[curr_ndx].append(line[i])
    # Set the systems indices
    sys.ndx_names = names
    # Set the systems indices
    sys.ndxs = [[sys.atoms[ndx] for ndx in indices[i]] for i in range(len(indices))]
