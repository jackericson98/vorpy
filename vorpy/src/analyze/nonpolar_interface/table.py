"""Build local geometry analysis tables without mutating VorPy networks."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

from .classification import (
    NONPOLAR_RULE_VERSION,
    classify_atom,
    classify_surface_pair,
)


def _as_int(value):
    return int(value)


def _record(table, index, fallback=None):
    """Return a row by network index, then by its ``num`` field."""
    if table is None:
        return {}
    try:
        row = table.loc[index]
        return row.to_dict() if hasattr(row, "to_dict") else dict(row)
    except (KeyError, TypeError, IndexError):
        pass
    if fallback is not None and "num" in table.columns:
        matches = table[table["num"].map(lambda value: int(value) == int(index))]
        if len(matches):
            return matches.iloc[0].to_dict()
    return {}


def _float(row, *names, default=np.nan):
    for name in names:
        value = row.get(name)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                continue
    return default


def _interface_sides(network, iface_grps):
    groups = iface_grps if iface_grps is not None else getattr(network, "iface_grps", None)
    if groups is None:
        return None
    if len(groups) != 2:
        raise ValueError("Interface groups must contain exactly two collections.")
    return frozenset(int(value) for value in groups[0]), frozenset(int(value) for value in groups[1])


def _surface_interface(balls, sides):
    if sides is None:
        return None, None, None
    side_a, side_b = sides
    members = set(int(value) for value in balls)
    in_a = bool(members & side_a)
    in_b = bool(members & side_b)
    if in_a and in_b:
        return 1, "cross", "cross_side"
    if members & side_a:
        return 0, "side_1", "same_side"
    if members & side_b:
        return 0, "side_2", "same_side"
    return 0, "outside", "outside_interface_groups"


def _identity(row, system_row=None):
    merged = dict(system_row or {})
    merged.update(row)
    return merged


def build_geometry_tables(
        network,
        *,
        system_id="",
        frame_id=None,
        network_id=None,
        iface_grps=None,
        nonpolar_rule=NONPOLAR_RULE_VERSION,
        complete_only=True,
):
    """Return ``{"surfaces": DataFrame, "atoms": DataFrame, "metadata": dict}``.

    The function only reads ``network``. Interface labels are based on the
    network's two explicit interface groups and surface ownership, never on a
    new distance cutoff. A missing interface-group definition produces missing
    labels (``None``), rather than falsely labeling all surfaces non-interface.
    """
    if network is None or getattr(network, "balls", None) is None:
        raise ValueError("A solved network with a balls table is required.")
    if getattr(network, "surfs", None) is None:
        raise ValueError("A solved network with a surfaces table is required.")

    balls = network.balls
    system = getattr(network, "sys", None)
    system_balls = getattr(system, "balls", None)
    sides = _interface_sides(network, iface_grps)
    resolved_network_id = network_id if network_id is not None else getattr(network, "group_name", "")

    atom_records = {}
    for index in balls.index:
        index = _as_int(index)
        local = _record(balls, index)
        parent = _record(system_balls, index, fallback=balls)
        atom_records[index] = _identity(local, parent)

    surfaces = []
    atom_surface_area = defaultdict(float)
    atom_interface_area = defaultdict(float)
    atom_interface_count = defaultdict(int)
    for surface_index, surface in network.surfs.iterrows():
        row = surface.to_dict()
        raw_balls = row.get("balls", row.get("Balls", ()))
        if raw_balls is None or len(raw_balls) < 2:
            continue
        surface_balls = tuple(_as_int(value) for value in raw_balls)
        if complete_only and len(surface_balls) != 2:
            continue
        area = _float(row, "sa", "Surface Area", default=0.0)
        if not np.isfinite(area) or area <= 0:
            continue
        first, second = (atom_records.get(index, {}) for index in surface_balls[:2])
        interface_label, interface_side, interface_kind = _surface_interface(surface_balls, sides)
        chemistry = classify_surface_pair(first, second, nonpolar_rule)
        record = {
            "system_id": system_id,
            "frame_id": frame_id,
            "network_id": resolved_network_id,
            "surface_id": _as_int(surface_index),
            "ball_1": surface_balls[0],
            "ball_2": surface_balls[1],
            "surface_area": area,
            "mean_curvature": _float(row, "mean_curv", "Mean Curvature"),
            "average_mean_curvature": _float(row, "avg_mean_curv", "Average Mean Curvature"),
            "gaussian_curvature": _float(row, "gauss_curv", "Gaussian Curvature"),
            "average_gaussian_curvature": _float(row, "avg_gauss_curv", "Average Gaussian Curvature"),
            "integrated_mean_curvature": _float(row, "int_mean_curv", "Integrated Mean Curvature"),
            "integrated_mean_curvature_squared": _float(row, "int_mean_curv_sq", "Integrated Mean Curvature Squared"),
            "integrated_gaussian_curvature": _float(row, "int_gauss_curv", "Integrated Gaussian Curvature"),
            "representative_surface_energy": _float(row, "surf_energy", "Representative Surface Energy"),
            "contact_area": _float(row, "contact_area", "Contact Area", default=0.0),
            "overlap": _float(row, "overlap", "Overlap", default=0.0),
            "interface_label": interface_label,
            "interface_side": interface_side,
            "interface_kind": interface_kind,
            "nonpolar_surface_label": chemistry,
            "nonpolar_rule": nonpolar_rule,
        }
        for position, atom_index in enumerate(surface_balls[:2], start=1):
            atom = atom_records.get(atom_index, {})
            record[f"ball_{position}_nonpolar_class"] = classify_atom(atom, nonpolar_rule)
            record[f"ball_{position}_atom_name"] = atom.get("name", "")
            record[f"ball_{position}_element"] = atom.get("element", "")
            record[f"ball_{position}_residue_name"] = atom.get("res_name", "")
            record[f"ball_{position}_residue_number"] = atom.get("res_seq", "")
            record[f"ball_{position}_chain"] = atom.get("chain_name", atom.get("chain", ""))
            atom_surface_area[atom_index] += area
        if interface_label == 1:
            for atom_index in surface_balls[:2]:
                atom_interface_area[atom_index] += area
                atom_interface_count[atom_index] += 1
        surfaces.append(record)

    atom_rows = []
    for atom_index, atom in atom_records.items():
        exposed_area = _float(atom, "sa", "Surface Area", default=atom_surface_area[atom_index])
        if complete_only and "complete" in atom and not bool(atom["complete"]):
            continue
        if not np.isfinite(exposed_area) or exposed_area <= 0:
            continue
        interface_area = atom_interface_area[atom_index]
        atom_rows.append({
            "system_id": system_id,
            "frame_id": frame_id,
            "network_id": resolved_network_id,
            "atom_index": atom_index,
            "chain": atom.get("chain_name", atom.get("chain", "")),
            "residue_name": atom.get("res_name", ""),
            "residue_number": atom.get("res_seq", ""),
            "atom_name": atom.get("name", ""),
            "element": atom.get("element", ""),
            "radius": _float(atom, "rad", "Radius"),
            "volume": _float(atom, "vol", "Volume"),
            "surface_area": exposed_area,
            "interface_area": interface_area,
            "interface_fraction": interface_area / exposed_area,
            "interface_surface_count": atom_interface_count[atom_index],
            "nonpolar_class": classify_atom(atom, nonpolar_rule),
            "interface_label": None if sides is None else int(interface_area > 0),
            "max_mean_curvature": _float(atom, "max_mean_curv", "Maximum Mean Curvature"),
            "average_mean_curvature": _float(atom, "avg_mean_surf_curv", "Average Mean Surface Curvature"),
            "max_gaussian_curvature": _float(atom, "max_gauss_curv", "Maximum Gaussian Curvature"),
            "average_gaussian_curvature": _float(atom, "avg_gauss_surf_curv", "Average Gaussian Surface Curvature"),
            "integrated_mean_curvature": _float(atom, "int_mean_curv", "Integrated Mean Curvature"),
            "integrated_mean_curvature_squared": _float(atom, "int_mean_curv_sq", "Integrated Mean Curvature Squared"),
            "integrated_gaussian_curvature": _float(atom, "int_gauss_curv", "Integrated Gaussian Curvature"),
            "neighbor_count": atom.get("number_of_neighbors", atom.get("neighbor_count", np.nan)),
        })

    metadata = {
        "system_id": system_id,
        "frame_id": frame_id,
        "network_id": resolved_network_id,
        "nonpolar_rule": nonpolar_rule,
        "complete_only": complete_only,
        "interface_groups_supplied": sides is not None,
        "surface_count": len(surfaces),
        "atom_count": len(atom_rows),
        "units": {"surface_area": "A^2", "volume": "A^3", "mean_curvature": "A^-1", "gaussian_curvature": "A^-2"},
    }
    return {
        "surfaces": pd.DataFrame(surfaces),
        "atoms": pd.DataFrame(atom_rows),
        "metadata": metadata,
    }
