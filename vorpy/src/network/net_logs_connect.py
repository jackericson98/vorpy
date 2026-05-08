import time
import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import combinations


def _empty_index_lists(n):
    return [[] for _ in range(n)]


def _as_int_list(value):
    if value is None:
        return []

    if isinstance(value, list):
        return [int(_) for _ in value]

    if isinstance(value, tuple):
        return [int(_) for _ in value]

    if pd.isna(value):
        return []

    return [int(value)]


def _dist(p0, p1):
    return float(np.linalg.norm(np.array(p0, dtype=float) - np.array(p1, dtype=float)))


def _build_subset_lookup(df, subset_size):
    """
    Map every subset of size subset_size from each row's balls to row indices.

    Example:
        verts have 4 balls.
        For edges, we need 3-ball subsets.
        For surfs, we need 2-ball subsets.
    """
    lookup = defaultdict(list)

    for row_i, row in df.iterrows():
        balls = tuple(sorted(int(_) for _ in row["balls"]))

        for subset in combinations(balls, subset_size):
            lookup[subset].append(row_i)

    return lookup


def _split_doublet_duplicate_edges(net, raw_edge_verts):
    """
    Fix duplicate log edges caused by doublets.

    For normal edges:
        one edge row with 2 verts -> unchanged

    For doublet split edges:
        two edge rows with same 3 balls may initially each see 4 candidate verts.
        We split them into two 2-vert edges using the original doublify logic:
        assign each outer vertex to the closer member of the doublet pair.
    """
    edge_groups = defaultdict(list)

    for edge_i, edge in net.edges.iterrows():
        key = tuple(sorted([int(_) for _ in edge["balls"]]))
        edge_groups[key].append(edge_i)

    fixed_edge_verts = list(raw_edge_verts)

    for edge_key, edge_indices in edge_groups.items():
        if len(edge_indices) == 1:
            continue

        candidate_verts = sorted(set(
            vert_i
            for edge_i in edge_indices
            for vert_i in raw_edge_verts[edge_i]
        ))

        if len(candidate_verts) != 4:
            continue

        doublet_verts = [
            vert_i
            for vert_i in candidate_verts
            if int(net.verts.iloc[vert_i].get("dub", 0)) in {1, 2}
        ]

        outer_verts = [
            vert_i
            for vert_i in candidate_verts
            if vert_i not in doublet_verts
        ]

        if len(doublet_verts) != 2 or len(outer_verts) != 2:
            continue

        # Preserve original convention: dub=2 is primary, dub=1 is secondary.
        doublet_verts = sorted(
            doublet_verts,
            key=lambda v: int(net.verts.iloc[v].get("dub", 0)),
            reverse=True
        )

        primary_dub = doublet_verts[0]
        secondary_dub = doublet_verts[1]

        primary_loc = net.verts.iloc[primary_dub]["loc"]
        secondary_loc = net.verts.iloc[secondary_dub]["loc"]

        primary_outer = []
        secondary_outer = []

        for outer_vert in outer_verts:
            outer_loc = net.verts.iloc[outer_vert]["loc"]

            if _dist(outer_loc, primary_loc) < _dist(outer_loc, secondary_loc):
                primary_outer.append(outer_vert)
            else:
                secondary_outer.append(outer_vert)

        if len(primary_outer) != 1 or len(secondary_outer) != 1:
            continue

        split_pairs = [
            sorted([primary_dub, primary_outer[0]]),
            sorted([secondary_dub, secondary_outer[0]]),
        ]

        # Assign split pairs to duplicate edge rows by matching saved edge length.
        available_pairs = list(split_pairs)

        for edge_i in edge_indices:
            saved_length = float(net.edges.iloc[edge_i].get("length", np.nan))

            if np.isnan(saved_length):
                fixed_edge_verts[edge_i] = available_pairs.pop(0)
                continue

            best_pair_i = None
            best_diff = float("inf")

            for pair_i, pair in enumerate(available_pairs):
                v0 = np.array(net.verts.iloc[pair[0]]["loc"], dtype=float)
                v1 = np.array(net.verts.iloc[pair[1]]["loc"], dtype=float)

                pair_dist = _dist(v0, v1)
                diff = abs(pair_dist - saved_length)

                if diff < best_diff:
                    best_diff = diff
                    best_pair_i = pair_i

            fixed_edge_verts[edge_i] = available_pairs.pop(best_pair_i)

    return fixed_edge_verts


def get_dubs(verts):
    """
    Reconstruct the dub column for vertices loaded from logs.

    Consecutive vertices with the same 4 defining balls are treated as a doublet.
    The first gets dub=2 and the second gets dub=1 to match the original
    find_net_verts assignment.

    Parameters
    ----------
    verts : pandas.DataFrame
        Vertex dataframe with a 'balls' column.

    Returns
    -------
    list[int]
        dub values aligned to verts.index.
    """
    dubs = [0 for _ in range(len(verts))]

    i = 0
    while i < len(verts) - 1:
        balls_i = list(verts.iloc[i]["balls"])
        balls_j = list(verts.iloc[i + 1]["balls"])

        if balls_i == balls_j:
            dubs[i] = 2
            dubs[i + 1] = 1
            i += 2
        else:
            i += 1

    return dubs


def net_logs_connect(net):
    """
    Connect a network loaded from logs.

    This does not solve vertices or rebuild topology. It only fills:
        balls -> verts, edges, surfs
        verts -> edges, surfs, dub
        edges -> verts, surfs
        surfs -> verts, edges
    """

    if net.balls is None or net.verts is None or net.edges is None or net.surfs is None:
        raise ValueError("Cannot logs-connect network. balls, verts, edges, and surfs must all exist.")

    net.balls = net.balls.reset_index(drop=True)
    net.verts = net.verts.reset_index(drop=True)
    net.edges = net.edges.reset_index(drop=True)
    net.surfs = net.surfs.reset_index(drop=True)

    for df_name in ["verts", "edges", "surfs"]:
        df = getattr(net, df_name)

        if "Balls" in df.columns and "balls" not in df.columns:
            df["balls"] = df["Balls"]

        if "balls" not in df.columns:
            raise KeyError(
                f"{df_name} dataframe has no 'balls' or 'Balls' column. "
                f"Columns found: {list(df.columns)}"
            )

        df["balls"] = df["balls"].apply(_as_int_list)

    n_balls = len(net.balls)
    n_verts = len(net.verts)
    n_edges = len(net.edges)
    n_surfs = len(net.surfs)

    ball_verts = _empty_index_lists(n_balls)
    ball_edges = _empty_index_lists(n_balls)
    ball_surfs = _empty_index_lists(n_balls)

    vert_edges = _empty_index_lists(n_verts)
    vert_surfs = _empty_index_lists(n_verts)

    edge_verts = _empty_index_lists(n_edges)
    edge_surfs = _empty_index_lists(n_edges)

    surf_verts = _empty_index_lists(n_surfs)
    surf_edges = _empty_index_lists(n_surfs)

    # balls <-> verts
    for vert_i, vert in net.verts.iterrows():
        for ball_i in vert["balls"]:
            if 0 <= ball_i < n_balls:
                ball_verts[ball_i].append(vert_i)

    vert_by_edge_balls = _build_subset_lookup(net.verts, 3)

    for edge_i, edge in net.edges.iterrows():
        edge_key = tuple(sorted(int(_) for _ in edge["balls"]))

        for ball_i in edge_key:
            if 0 <= ball_i < n_balls:
                ball_edges[ball_i].append(edge_i)

        matched_verts = vert_by_edge_balls.get(edge_key, [])

        edge_verts[edge_i] = list(matched_verts)

        for vert_i in matched_verts:
            vert_edges[vert_i].append(edge_i)

    vert_by_surf_balls = _build_subset_lookup(net.verts, 2)
    edge_by_surf_balls = _build_subset_lookup(net.edges, 2)

    # balls <-> surfs, surfs <-> verts, surfs <-> edges
    for surf_i, surf in net.surfs.iterrows():
        surf_key = tuple(sorted(int(_) for _ in surf["balls"]))

        for ball_i in surf_key:
            if 0 <= ball_i < n_balls:
                ball_surfs[ball_i].append(surf_i)

        matched_verts = vert_by_surf_balls.get(surf_key, [])
        matched_edges = edge_by_surf_balls.get(surf_key, [])

        surf_verts[surf_i] = list(matched_verts)
        surf_edges[surf_i] = list(matched_edges)

        for vert_i in matched_verts:
            vert_surfs[vert_i].append(surf_i)

        for edge_i in matched_edges:
            edge_surfs[edge_i].append(surf_i)

    net.balls["verts"] = ball_verts
    net.balls["edges"] = ball_edges
    net.balls["surfs"] = ball_surfs

    net.verts["edges"] = vert_edges
    net.verts["surfs"] = vert_surfs

    # Logs do not save dub. Use 0 as safe default unless rebuilt later.
    if "dub" not in net.verts.columns:
        net.verts["dub"] = get_dubs(net.verts)

    edge_verts = _split_doublet_duplicate_edges(net, edge_verts)

    net.edges["verts"] = edge_verts
    net.edges["surfs"] = edge_surfs
    vert_edges = [[] for _ in range(len(net.verts))]

    for edge_i, verts in enumerate(edge_verts):
        for vert_i in verts:
            vert_edges[vert_i].append(edge_i)

    net.verts["edges"] = vert_edges

    net.surfs["verts"] = surf_verts
    net.surfs["edges"] = surf_edges

    return net
