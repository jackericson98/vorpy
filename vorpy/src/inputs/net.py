import os.path
import os
import time
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


def _recalculate_loaded_vertices(net, verts):
    """
    Recalculate loaded vertex loc/rad from saved 4-ball definitions.

    Important:
    - Normal vertex: use loc/rad.
    - Doublet vertex pair: use both loc/rad and loc2/rad2.
    """
    start_time = net.metrics['start']
    if verts is None or len(verts) == 0:
        return verts

    if "balls" not in verts.columns:
        return verts

    verts = verts.copy()
    total_verts = len(verts)
    processed_verts = 0

    if "dub" not in verts.columns:
        from vorpy.src.network.net_logs_connect import get_dubs
        verts["dub"] = get_dubs(verts)

    net_type = net.settings.get("net_type", "aw") if net.settings else "aw"
    use_flat = net_type in {"prm", "pow"}
    use_power = net_type == "pow"

    new_locs = list(verts["loc"])
    new_rads = list(verts["rad"])

    grouped = {}

    for vert_i, vert in verts.iterrows():
        key = tuple(int(_) for _ in vert["balls"])
        grouped.setdefault(key, []).append(vert_i)

    failed = []

    for vballs_key, vert_indices in grouped.items():

        processed_verts += len(vert_indices)

        my_time = time.perf_counter() - start_time
        h, m, s = get_time(my_time)
        percentage = 100.0 * processed_verts / total_verts

        print(
            "\rRun Time = {}:{:02d}:{:2.2f} - Process: recalculating vertices: {} verts - {:.2f} %"
            .format(int(h), int(m), round(s, 2), processed_verts, percentage),
            end=""
        )

        vballs = list(vballs_key)

        if len(vballs) != 4:
            failed.extend(vert_indices)
            continue

        locs = [np.array(net.balls.iloc[ball_i]["loc"], dtype=float) for ball_i in vballs]
        rads = [float(net.balls.iloc[ball_i]["rad"]) for ball_i in vballs]

        try:
            if use_flat:
                loc, rad = calc_flat_vert(locs, rads, power=use_power)
                loc2, rad2 = None, None
            else:
                loc, rad, loc2, rad2 = calc_vert(locs, rads)

            if loc is None or rad is None:
                failed.extend(vert_indices)
                continue

            # Normal vertex
            if len(vert_indices) == 1 or loc2 is None or rad2 is None:
                vi = vert_indices[0]
                new_locs[vi] = np.array(loc, dtype=float)
                new_rads[vi] = float(rad)
                continue

            # Doublet vertex pair
            if len(vert_indices) == 2:
                vi0, vi1 = vert_indices

                saved0 = np.array(verts.iloc[vi0]["loc"], dtype=float)
                saved1 = np.array(verts.iloc[vi1]["loc"], dtype=float)

                sol0 = np.array(loc, dtype=float)
                sol1 = np.array(loc2, dtype=float)

                # Prefer assigning by nearest saved rounded coordinates.
                assignment_a = (
                    np.linalg.norm(saved0 - sol0) +
                    np.linalg.norm(saved1 - sol1)
                )

                assignment_b = (
                    np.linalg.norm(saved0 - sol1) +
                    np.linalg.norm(saved1 - sol0)
                )

                if assignment_a <= assignment_b:
                    new_locs[vi0] = sol0
                    new_rads[vi0] = float(rad)
                    new_locs[vi1] = sol1
                    new_rads[vi1] = float(rad2)
                else:
                    new_locs[vi0] = sol1
                    new_rads[vi0] = float(rad2)
                    new_locs[vi1] = sol0
                    new_rads[vi1] = float(rad)

                continue

            # More than two identical 4-ball rows is unexpected.
            failed.extend(vert_indices)

        except Exception:
            failed.extend(vert_indices)

    verts["loc"] = new_locs
    verts["rad"] = new_rads

    if failed:
        print(f"WARNING: {len(failed)} loaded vertices could not be recalculated and kept saved log coordinates.")

    return verts


def read_net(net, file_name, rebuild_edges=True, rebuild_surfs=True, analyze=True, store_points=True):
    """
    Load a vorpy logs file into an existing Network object.

    This reads the saved ball / vertex / edge / surface topology from aw_logs.csv-style
    logs, attaches it to `net`, and then optionally rebuilds missing edge points,
    surface points, triangles, and derived analysis fields.
    """
    if file_name is None or not os.path.exists(file_name):
        raise FileNotFoundError(f"Network log file not found: {file_name}")

    log_data = read_logs(file_name, all_=True)

    build_data = log_data["data"]
    group_data = log_data["group data"]

    balls = _rename_log_columns(log_data["atoms"])
    verts = _rename_log_columns(log_data["verts"])
    edges = _rename_log_columns(log_data["edges"])
    surfs = _rename_log_columns(log_data["surfs"])

    if "loc" not in balls.columns and {"X", "Y", "Z"}.issubset(balls.columns):
        balls["loc"] = balls[["X", "Y", "Z"]].values.tolist()

    if "rad" not in balls.columns and "Radius" in balls.columns:
        balls["rad"] = balls["Radius"]

    if "balls" not in verts.columns and "Balls" in verts.columns:
        verts["balls"] = verts["Balls"]

    if "balls" not in edges.columns and "Balls" in edges.columns:
        edges["balls"] = edges["Balls"]

    if "balls" not in surfs.columns and "Balls" in surfs.columns:
        surfs["balls"] = surfs["Balls"]

    # Logs use system/global ball indices.
    # Keep the full system ball table so vertices/edges/surfs can reference surrounding balls.
    logged_group_indices = sorted([int(_) for _ in balls["num"].dropna().tolist()])
    net.group = logged_group_indices

    if net.balls is None:
        net.balls = balls.copy().reset_index(drop=True)
    else:
        for _, log_ball in balls.iterrows():
            ball_num = int(log_ball["num"])

            for col in balls.columns:
                if col in {"Index", "num", "loc", "rad", "name", "mass"}:
                    continue

                if col not in net.balls.columns:
                    net.balls[col] = None

                net.balls.at[ball_num, col] = log_ball[col]

    net.loaded_log_ball_indices = logged_group_indices

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
    net.edges = edges
    net.surfs = surfs

    if net.settings is None:
        net.default_settings()

    net.settings["net_type"] = build_data.get("network_type", net.settings.get("net_type", "aw"))
    net.settings["surf_res"] = build_data.get("surface_resolution", net.settings.get("surf_res", 0.2))
    net.settings["box_size"] = build_data.get("box_size", net.settings.get("box_size", 1.5))
    net.settings["max_vert"] = build_data.get("max_vert", net.settings.get("max_vert", 40))

    net.log_data = build_data
    net.group_data = group_data

    # Preserve original build metrics from logs.
    net.original_metrics = {
        "tot": build_data.get("Total_Time", 0),
        "vert": build_data.get("vert_time", 0),
        "con": build_data.get("connect_time", 0),
        "surf": build_data.get("surf_time", 0),
        "anal": build_data.get("analysis_time", 0),
    }

    # Metrics for this reload/rebuild session.
    net.metrics.setdefault("start", 0)
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
