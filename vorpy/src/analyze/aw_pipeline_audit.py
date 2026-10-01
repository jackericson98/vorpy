"""Forensic comparison of full and interface-filtered AW networks.

This module is diagnostic only. It does not change AW construction or the
eligibility decisions in :mod:`aw_interface_curvature`.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


def stable_atom_identity(network, generator_id):
    """Return a source-record identity independent of network row numbering."""
    # Network.balls is a geometry-focused table and may omit molecular labels.
    # Its generator index is the system row index, so resolve the source atom
    # from Network.sys.balls whenever that table is available.
    atom_table = getattr(getattr(network, "sys", None), "balls", None)
    if atom_table is None or int(generator_id) not in atom_table.index:
        atom_table = network.balls
    row = atom_table.loc[int(generator_id)]
    xyz = tuple(round(float(value), 3) for value in row["loc"])
    return "|".join((
        str(row.get("chain_name", "")), str(row.get("res_name", "")),
        str(row.get("res_seq", "")), str(row.get("pdb_ins_code", "")),
        str(row.get("name", "")), str(row.get("element", "")),
        ",".join(map(str, xyz)),
    ))


def _ids(value):
    try:
        return tuple(int(item) for item in value)
    except (TypeError, ValueError):
        return ()


def audit_interface_surfaces(network, group_a, group_b, *, full_reference=None):
    """Describe candidate interface surfaces and test their own boundary.

    A surface passes the feature-level topology check when it has a finite
    triangulated patch and finite positive area, and its listed AW boundary
    edges form disjoint closed cycles with resolved endpoints and reciprocal
    surface/vertex incidence. This deliberately does not inspect participating
    cells' global ``complete`` flags.
    """
    group_a, group_b = set(map(int, group_a)), set(map(int, group_b))
    reference = {}
    if full_reference is not None:
        for sid, row in full_reference.surfs.iterrows():
            ids = _ids(row.get("balls", ()))
            if len(ids) == 2:
                reference[tuple(sorted(stable_atom_identity(full_reference, i) for i in ids))] = (int(sid), row)

    edge_geometry_cache = {}
    rows = []
    for sid, surf in network.surfs.iterrows():
        ids = _ids(surf.get("balls", ()))
        if len(ids) != 2 or not ((ids[0] in group_a and ids[1] in group_b)
                                  or (ids[1] in group_a and ids[0] in group_b)):
            continue
        edge_ids, vert_ids = _ids(surf.get("edges", ())), _ids(surf.get("verts", ()))
        points = np.asarray(surf.get("points", ()), dtype=float)
        tris = np.asarray(surf.get("tris", ()), dtype=int)
        area = _float(surf.get("sa"))
        errors = []
        endpoint_degree = {vid: 0 for vid in vert_ids}
        all_edges_supported = bool(edge_ids)
        edge_types, edge_lengths = [], []
        for eid in edge_ids:
            if eid not in network.edges.index:
                all_edges_supported = False
                errors.append(f"missing_edge:{eid}")
                continue
            edge = network.edges.loc[eid]
            endpoints = _ids(edge.get("verts", ()))
            if len(endpoints) != 2 or any(v not in vert_ids or v not in network.verts.index for v in endpoints):
                all_edges_supported = False
                errors.append(f"unresolved_edge_endpoints:{eid}")
                continue
            for vid in endpoints:
                endpoint_degree[vid] = endpoint_degree.get(vid, 0) + 1
            if int(sid) not in _ids(edge.get("surfs", ())):
                all_edges_supported = False
                errors.append(f"missing_reciprocal_edge_surface:{eid}")
            if any(int(sid) not in _ids(network.verts.loc[v].get("surfs", ())) for v in endpoints):
                all_edges_supported = False
                errors.append(f"missing_reciprocal_vertex_surface:{eid}")
            try:
                if eid not in edge_geometry_cache:
                    from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge
                    edge_geometry_cache[eid] = resolve_aw_network_edge(network, eid)
                resolved = edge_geometry_cache[eid]
                geometry = resolved.geometry
                finite_interval = (math.isfinite(float(geometry.t_min)) and
                                   math.isfinite(float(geometry.t_max)) and
                                   float(geometry.t_max) > float(geometry.t_min))
                if not finite_interval:
                    all_edges_supported = False
                    errors.append(f"nonfinite_edge_interval:{eid}")
                edge_types.append(type(geometry).__name__)
                edge_lengths.append(_float(getattr(resolved, "length", None)))
            except Exception as exc:
                all_edges_supported = False
                edge_types.append("unresolved")
                errors.append(f"edge_geometry:{eid}:{type(exc).__name__}")
        closed = bool(endpoint_degree) and all(degree == 2 for degree in endpoint_degree.values())
        if not closed:
            errors.append("boundary_not_disjoint_closed_cycles")
        finite_mesh = (points.ndim == 2 and points.shape[1:] == (3,) and len(points) > 0
                       and np.isfinite(points).all() and tris.ndim == 2 and tris.shape[1:] == (3,)
                       and len(tris) > 0 and (tris >= 0).all() and (tris < len(points)).all())
        if not finite_mesh:
            errors.append("missing_or_invalid_surface_mesh")
        if area is None or area <= 0:
            errors.append("nonpositive_or_nonfinite_area")
        cells_complete = all(bool(network.balls.at[i, "complete"]) for i in ids)
        analytic_surface = surf.get("func") is not None
        atom_pair = tuple(sorted(stable_atom_identity(network, i) for i in ids))
        ref = reference.get(atom_pair)
        full_cell_flags = ()
        if ref is not None and full_reference is not None:
            full_ids = _ids(ref[1].get("balls", ()))
            full_cell_flags = tuple(bool(full_reference.balls.at[i, "complete"]) for i in full_ids)
        rows.append({
            "surface_id": int(sid), "generator_ids": ";".join(map(str, ids)),
            "stable_generator_ids": " || ".join(atom_pair),
            "atom_numbers": ";".join(str(int(network.balls.at[i, "num"])) for i in ids),
            "residue_identities": " || ".join(_residue_identity(network, i) for i in ids),
            "cell_a_complete": bool(network.balls.at[ids[0], "complete"]),
            "cell_b_complete": bool(network.balls.at[ids[1], "complete"]),
            "both_cells_complete": cells_complete,
            "surface_has_vertices": bool(vert_ids), "vertex_count": len(vert_ids),
            "surface_has_edges": bool(edge_ids), "edge_count": len(edge_ids),
            "boundary_closed_cycles": closed, "analytic_aw_surface_available": analytic_surface,
            "all_boundary_edges_supported": all_edges_supported,
            "boundary_edge_geometry_types": ";".join(sorted(set(edge_types))),
            "surface_area": area, "mesh_point_count": len(points) if points.ndim else 0,
            "mesh_triangle_count": len(tris) if tris.ndim else 0,
            "feature_complete": bool(finite_mesh and area is not None and area > 0 and closed and all_edges_supported),
            "feature_errors": ";".join(errors),
            "full_network_surface_id": "" if ref is None else ref[0],
            "full_network_area": "" if ref is None else _float(ref[1].get("sa")),
            "full_network_cell_flags": ";".join(map(str, full_cell_flags)),
            "full_network_cells_complete": all(full_cell_flags) if full_cell_flags else "",
            "area_difference_vs_full": "" if ref is None else area - _float(ref[1].get("sa")),
        })
    return rows


def compare_surface_tables(filtered_rows, full_network, group_a, group_b, previous_csv=None):
    """Match candidate surfaces to a full network by stable atom identities."""
    full = {}
    for sid, surf in full_network.surfs.iterrows():
        ids = _ids(surf.get("balls", ()))
        if len(ids) == 2 and ((ids[0] in group_a and ids[1] in group_b)
                              or (ids[1] in group_a and ids[0] in group_b)):
            key = tuple(sorted(stable_atom_identity(full_network, i) for i in ids))
            full[key] = (int(sid), surf)
    old = {}
    if previous_csv and Path(previous_csv).exists():
        import pandas as pd
        for _, row in pd.read_csv(previous_csv).iterrows():
            try:
                ids = tuple(int(value) for value in str(row.generator_ids).split(";"))
            except ValueError:
                continue
            key = tuple(sorted(stable_atom_identity(full_network, i) for i in ids))
            old[key] = row
    result = []
    for item in filtered_rows:
        key = tuple(sorted(item["stable_generator_ids"].split(" || ")))
        full_item, prior = full.get(key), old.get(key)
        result.append({
            "stable_generator_ids": " || ".join(key),
            "filtered_surface_id": item["surface_id"],
            "full_surface_id": "" if full_item is None else full_item[0],
            "previous_standalone_present": prior is not None,
            "filtered_area": item["surface_area"],
            "full_area": "" if full_item is None else _float(full_item[1].get("sa")),
            "previous_standalone_area": "" if prior is None else float(prior.area),
            "filtered_cells_complete": item["both_cells_complete"],
            "full_cells_complete": "" if full_item is None else all(
                bool(full_network.balls.at[i, "complete"])
                for i in _ids(full_item[1].get("balls", ()))
            ),
            "filtered_feature_complete": item["feature_complete"],
            "full_geometry_present": full_item is not None,
            "previous_standalone_included": "" if prior is None else bool(prior.included),
            "difference_class": "present_in_both" if full_item else "filtered_only",
        })
    filtered_keys = {tuple(sorted(x["stable_generator_ids"].split(" || "))) for x in filtered_rows}
    for key, (sid, surf) in full.items():
        if key not in filtered_keys:
            prior = old.get(key)
            result.append({
                "stable_generator_ids": " || ".join(key), "filtered_surface_id": "",
                "full_surface_id": sid, "previous_standalone_present": prior is not None,
                "filtered_area": "", "full_area": _float(surf.get("sa")),
                "previous_standalone_area": "" if prior is None else float(prior.area),
                "filtered_cells_complete": "", "full_cells_complete": all(
                    bool(full_network.balls.at[i, "complete"]) for i in _ids(surf.get("balls", ()))
                ), "filtered_feature_complete": "", "full_geometry_present": True,
                "previous_standalone_included": "" if prior is None else bool(prior.included),
                "difference_class": "full_network_only",
            })
    return result


def write_rows(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if fields:
            writer.writeheader()
            writer.writerows(rows)


def _float(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _residue_identity(network, index):
    atom_table = getattr(getattr(network, "sys", None), "balls", None)
    if atom_table is None or int(index) not in atom_table.index:
        atom_table = network.balls
    row = atom_table.loc[int(index)]
    return f"{row.get('chain_name', '')}:{row.get('res_name', '')}{row.get('res_seq', '')}{row.get('pdb_ins_code', '')}:{row.get('name', '')}"
