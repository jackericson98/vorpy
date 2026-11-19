"""Export detailed atlas tables from an existing network archive.

This module is deliberately read-only with respect to the solved network: it
calls :func:`load_network`, copies no executable solver state, and never calls
the network builder or an interface producer.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np

from vorpy.src.io.network_archive import load_network


_SURFACE_CURVATURE = (
    ("mean_curv", "mean_curvature_A^-1"),
    ("avg_mean_curv", "average_mean_curvature_A^-1"),
    ("gauss_curv", "gaussian_curvature_A^-2"),
    ("avg_gauss_curv", "average_gaussian_curvature_A^-2"),
    ("int_mean_curv", "integrated_mean_curvature_A"),
    ("int_mean_curv_sq", "integrated_mean_curvature_squared_A2"),
    ("int_gauss_curv", "integrated_gaussian_curvature"),
)


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _ids(value):
    try:
        return tuple(int(item) for item in value)
    except (TypeError, ValueError):
        return ()


def _number(value):
    try:
        value = np.asarray(value, dtype=float)
        if value.size != 1:
            return None
        value = float(value.reshape(-1)[0])
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def _array(value, shape_last=None, dtype=float):
    try:
        array = np.asarray(value, dtype=dtype)
    except (TypeError, ValueError):
        return np.empty((0,), dtype=dtype)
    if shape_last is not None and (array.ndim != 2 or array.shape[1] != shape_last):
        return np.empty((0, shape_last), dtype=dtype)
    return array


def _text(row, *names):
    for name in names:
        value = row.get(name) if hasattr(row, "get") else None
        if value is not None and str(value).strip() not in {"", "nan", "None"}:
            return str(value).strip()
    return ""


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_atom(net, generator_id, namespace):
    """Return archive-backed atom identity for a local network generator."""
    local = int(generator_id)
    network_row = net.balls.loc[local]
    system_id = _number(network_row.get("system_num"))
    source = None
    if system_id is not None and getattr(net, "sys", None) is not None:
        try:
            source = net.sys.balls.loc[int(system_id)]
        except (KeyError, TypeError, ValueError):
            source = None
    row = source if source is not None else network_row
    loc = _array(row.get("loc"), dtype=float).reshape(-1)
    xyz = tuple(round(float(value), 3) for value in loc[:3]) if len(loc) >= 3 else ()
    chain = _text(row, "chain_name", "chain", "auth_asym_id")
    res_name = _text(row, "res_name", "auth_comp_id", "residue_name")
    res_seq = _text(row, "res_seq", "auth_seq_id", "residue_number")
    insertion = _text(row, "pdb_ins_code", "insertion_code", "icode")
    atom_name = _text(row, "name", "auth_atom_id", "atom_name")
    element = _text(row, "element").upper()
    stable_generator = (
        (getattr(net, "generator_stable_ids", None) or ())[local]
        if getattr(net, "generator_stable_ids", None) is not None
        and local < len(net.generator_stable_ids)
        else f"generator:{local}"
    )
    identity = "|".join((chain, res_name, res_seq, insertion, atom_name, element, str(xyz)))
    atom_id = "atom-v1:" + hashlib.sha256(identity.encode("utf-8")).hexdigest()[:24]
    archive_solvent_flag = bool(network_row.get("is_solvent_generator", False))
    residue_group = _group_name(chain, res_name, element, archive_solvent_flag)
    return {
        "generator_id": local,
        "generator_stable_id": str(stable_generator),
        "source_atom_id": atom_id,
        "source_identity": identity,
        "atom_name": atom_name,
        "residue_name": res_name,
        "residue_number": res_seq,
        "insertion_code": insertion,
        "chain": chain,
        "molecule": residue_group.split("_", 1)[0] if "_" in residue_group else residue_group,
        "group": residue_group,
        "element": element,
        "xyz": xyz,
        "radius_A": _number(row.get("rad")),
        "system_index": None if system_id is None else int(system_id),
        "is_solvent": archive_solvent_flag,
        "group_assignment_status": (
            "CHEMICAL_IDENTITY_OVERRIDES_ARCHIVE_SOLVENT_FLAG"
            if archive_solvent_flag and residue_group == "Zn"
            else "ARCHIVE_SOLVENT_FLAG" if archive_solvent_flag else "SOURCE_METADATA"
        ),
    }


def _group_name(chain, residue, element, is_solvent):
    if element.upper() == "ZN" or residue.upper() == "ZN":
        return "Zn"
    if is_solvent or residue.upper() in {"HOH", "WAT", "SOL", "TIP3", "TIP4"}:
        return "solvent"
    if chain.upper() == "A":
        return "protein_A"
    if chain.upper() in {"B", "C"}:
        return f"DNA_{chain.upper()}"
    return chain or "unclassified"


def _generator_records(net, ids, namespace):
    return [_source_atom(net, generator_id, namespace) for generator_id in ids]


def _groups(records):
    return tuple(dict.fromkeys(record["group"] for record in records if record["group"]))


def _classification(groups):
    groups = set(groups)
    if not groups:
        return "unclassified"
    if "solvent" in groups and "Zn" in groups:
        return "solvent_ion_contact"
    if "solvent" in groups:
        return "solvent_contact"
    if "Zn" in groups:
        return "ion_contact"
    return "cross_group" if len(groups) > 1 else "intra_group"


def _write_csv(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _json(value) if isinstance(value, (tuple, list, dict)) else value
                             for key, value in row.items()})


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
                          encoding="utf-8")


def _feature_status(mesh_ok, refs_ok, cells_complete):
    if not mesh_ok or not refs_ok:
        return "UNRESOLVED"
    return "CERTIFIED" if cells_complete else "PARTIAL"


def _write_surface_mesh(path, parts):
    points, triangles = [], []
    feature_map = []
    for record_id, point_array, triangle_array in parts:
        point_array = _array(point_array, 3)
        triangle_array = _array(triangle_array, 3, dtype=int)
        if not len(point_array) or not len(triangle_array):
            feature_map.append({"record_id": record_id, "point_count": 0, "triangle_count": 0})
            continue
        offset = len(points)
        points.extend(point_array.tolist())
        triangles.extend((triangle_array + offset).tolist())
        feature_map.append({"record_id": record_id, "point_offset": offset,
                            "point_count": len(point_array), "triangle_count": len(triangle_array)})
    with Path(path).open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"OFF\n{len(points)} {len(triangles)} 0\n")
        for point in points:
            stream.write("%.6f %.6f %.6f\n" % tuple(point))
        for triangle in triangles:
            stream.write("3 %d %d %d\n" % tuple(triangle))
    return feature_map


def _write_edges(path, edge_geometry):
    vertex_count = 0
    mapping = []
    for record_id, points in edge_geometry:
        points = _array(points, 3)
        if len(points) < 2:
            continue
        start = vertex_count + 1
        mapping.append({"record_id": record_id, "point_offset": start - 1,
                        "point_count": len(points)})
        vertex_count += len(points)
    with Path(path).open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# VorPy Molecular Interface Atlas edge layer\n")
        offset = 0
        for record_id, points in edge_geometry:
            points = _array(points, 3)
            for point in points:
                stream.write("v %.6f %.6f %.6f\n" % tuple(point))
            if len(points) >= 2:
                stream.write("l %s\n" % " ".join(str(offset + i + 1) for i in range(len(points))))
            offset += len(points)
    return mapping


def _pdb_atom(serial, atom):
    x, y, z = (atom["xyz"] + (0.0, 0.0, 0.0))[:3]
    chain = (atom["chain"] or "X")[:1]
    res = (atom["residue_name"] or "UNK")[:3].rjust(3)
    name = (atom["atom_name"] or "X")[:4].rjust(4)
    element = (atom["element"] or "X")[:2].rjust(2)
    try:
        seq = int(float(atom["residue_number"]))
    except (TypeError, ValueError):
        seq = 1
    return (f"ATOM  {serial:5d} {name} {res} {chain}{seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element}\n")


def _write_pdb(path, atoms, record_ids=None):
    mapping = []
    with Path(path).open("w", encoding="utf-8", newline="\n") as stream:
        for serial, atom in enumerate(atoms, 1):
            stream.write(_pdb_atom(serial, atom))
            if record_ids is not None and serial <= len(record_ids):
                mapping.append({"record_id": record_ids[serial - 1], "serial": serial})
        stream.write("END\n")
    return mapping


def _write_pml(path, surface_files, selected=False):
    lines = [
        "# Cache-only VorPy Molecular Interface Atlas demonstration",
        "load source_atoms.pdb, source_atoms",
        "load geometry/edges.obj, atlas_edges",
        "load geometry/vertices.pdb, atlas_vertices",
        "show sticks, source_atoms",
        "show lines, atlas_edges",
        "show spheres, atlas_vertices",
        "select protein_A, source_atoms and chain A",
        "select DNA_B, source_atoms and chain B",
        "select DNA_C, source_atoms and chain C",
        "select zinc, source_atoms and resn ZN",
        "color cyan, protein_A",
        "color yellow, DNA_B or DNA_C",
        "color orange, zinc",
    ]
    for name in sorted(surface_files):
        lines.extend([f"load geometry/{name}, {Path(name).stem}",
                      f"show surface, {Path(name).stem}"])
    if selected:
        lines.extend([
            "load geometry/selected_surface.off, selected_surface",
            "load geometry/selected_edges.obj, selected_incident_edges",
            "load geometry/selected_vertices.pdb, selected_incident_vertices",
            "color red, selected_surface or selected_incident_edges or selected_incident_vertices",
            "show surface, selected_surface",
            "show lines, selected_incident_edges",
            "show spheres, selected_incident_vertices",
        ])
    lines.extend([
        "# Surface objects are independent by contact class.",
        "# Table record IDs and feature-to-geometry offsets are in geometry_index.json.",
        "# Pairwise contacts are adjacency classifications, not independent interface solves.",
    ])
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _file_manifest(root):
    output = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = str(path.relative_to(root)).replace("\\", "/")
        output[relative] = {"bytes": path.stat().st_size, "sha256": _sha256(path)}
    return output


def export_atlas_archive(archive_path, output_dir, *, source_path=None):
    """Create the Phase 2 tables and visualization layers from one archive.

    ``output_dir`` must not exist. This is the exporter equivalent of
    ``--fail-if-exists`` and protects prior atlas evidence from replacement.
    """
    archive_path = Path(archive_path).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite atlas output: {output_dir}")
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)

    archive_hash = _sha256(archive_path)
    net = load_network(archive_path)
    output_dir.mkdir(parents=True)
    geometry_dir = output_dir / "geometry"
    geometry_dir.mkdir()
    namespace = f"archive:{getattr(net, 'archive_id', 'unknown')}"
    surface_map = {int(index): f"{namespace}:surface:{int(index)}" for index in net.surfs.index}
    edge_map = {int(index): f"{namespace}:edge:{int(index)}" for index in net.edges.index}
    vertex_map = {int(index): f"{namespace}:vertex:{int(index)}" for index in net.verts.index}
    atoms = {int(index): _source_atom(net, int(index), namespace) for index in net.balls.index}

    def atom_records(ids):
        records = []
        for generator_id in ids:
            if int(generator_id) not in atoms:
                atoms[int(generator_id)] = _source_atom(net, int(generator_id), namespace)
            records.append(atoms[int(generator_id)])
        return records

    surface_rows, surface_geometry = [], {name: [] for name in
                                          ("intra_group", "cross_group", "ion_contact",
                                           "solvent_contact", "solvent_ion_contact", "unclassified")}
    for index, row in net.surfs.iterrows():
        index = int(index)
        generators = _ids(row.get("balls"))
        generator_records = atom_records(generators)
        groups = _groups(generator_records)
        edge_ids, vertex_ids = _ids(row.get("edges")), _ids(row.get("verts"))
        mesh_points, mesh_triangles = _array(row.get("points"), 3), _array(row.get("tris"), 3, dtype=int)
        refs_ok = all(edge_id in edge_map for edge_id in edge_ids) and all(vertex_id in vertex_map for vertex_id in vertex_ids)
        cells_complete = all(bool(net.balls.loc[generator_id].get("complete", False))
                             for generator_id in generators if generator_id in net.balls.index)
        topology_ok = all(
            len(_ids(net.edges.loc[edge_id].get("surfs"))) == 2
            for edge_id in edge_ids if edge_id in net.edges.index
        )
        mesh_ok = bool(len(mesh_points) and len(mesh_triangles) and np.isfinite(mesh_points).all())
        status = _feature_status(mesh_ok, refs_ok, cells_complete and topology_ok)
        classification = _classification(groups)
        if classification not in surface_geometry:
            classification = "unclassified"
        row_data = {
            "record_id": f"{namespace}:surface:{index}", "feature_index": index,
            "generator_ids": generators,
            "generator_stable_ids": tuple(record["generator_stable_id"] for record in generator_records),
            "source_atom_ids": tuple(record["source_atom_id"] for record in generator_records),
            "source_atoms": tuple({key: value for key, value in record.items() if key != "xyz"}
                                   for record in generator_records),
            "generator_ordered_groups": groups, "group_pair": tuple(sorted(groups)),
            "scheme": str((net.settings or {}).get("net_type", "unknown")),
            "representation": "analytic_plus_mesh" if row.get("func") is not None else "mesh_only",
            "surface_type": "cross_group" if len(set(groups)) > 1 else "intra_group" if groups else "unclassified",
            "contact_classification": classification, "area_A2": _number(row.get("sa")),
            "volume_A3": _number(row.get("vols")), "incident_edge_ids": tuple(edge_map[i] for i in edge_ids if i in edge_map),
            "incident_vertex_ids": tuple(vertex_map[i] for i in vertex_ids if i in vertex_map),
            "point_count": len(mesh_points), "triangle_count": len(mesh_triangles),
            "bounded": bool(edge_ids) and all(len(_ids(net.edges.loc[i].get("verts"))) == 2
                                               for i in edge_ids if i in net.edges.index),
            "bounded_status": "INFERRED_FROM_EDGE_ENDPOINTS",
            "complete": status == "CERTIFIED", "scientific_status": status,
            "unsupported_incidence": not refs_ok,
            "curvature_status": "PARTIAL_RAW_NETWORK_FIELDS" if any(_number(row.get(key)) is not None for key, _ in _SURFACE_CURVATURE) else "NOT_CALCULATED",
            "source_hash": archive_hash, "archive_identity": getattr(net, "archive_id", None),
            "geometry_lineage": "RAW_ARCHIVE_NETWORK_SURFACE;NO_RECALCULATION",
        }
        for source_name, export_name in _SURFACE_CURVATURE:
            row_data[export_name] = _number(row.get(source_name))
            row_data[f"{export_name}_unit"] = "1" if source_name == "int_gauss_curv" else "Å" if source_name.startswith("int_") else "Å^-2" if "gauss" in source_name else "Å^-1"
        surface_rows.append(row_data)
        surface_geometry[classification].append((row_data["record_id"], row.get("points"), row.get("tris")))

    edge_rows, edge_geometry = [], []
    for index, row in net.edges.iterrows():
        index = int(index)
        generators = _ids(row.get("balls")); generator_records = atom_records(generators)
        vertices, surfaces = _ids(row.get("verts")), _ids(row.get("surfs"))
        points = _array(row.get("points"), 3)
        refs_ok = all(vertex in vertex_map for vertex in vertices) and all(surface in surface_map for surface in surfaces)
        incidence = "complete" if len(surfaces) == 2 else "boundary" if len(surfaces) == 1 else "nonmanifold" if len(surfaces) > 2 else "unresolved"
        edge_status = ("CERTIFIED" if refs_ok and len(vertices) == 2 and incidence == "complete"
                       else "PARTIAL" if refs_ok else "UNRESOLVED")
        record_id = edge_map[index]
        edge_rows.append({
            "record_id": record_id, "feature_index": index, "generator_ids": generators,
            "generator_stable_ids": tuple(record["generator_stable_id"] for record in generator_records),
            "source_atom_ids": tuple(record["source_atom_id"] for record in generator_records),
            "endpoint_vertex_ids": tuple(vertex_map[v] for v in vertices if v in vertex_map),
            "adjacent_surface_ids": tuple(surface_map[s] for s in surfaces if s in surface_map),
            "analytic_type": "bounded_segment_or_polyline" if len(points) >= 2 else "UNRESOLVED",
            "length_A": _number(row.get("length")), "point_count": len(points),
            "bounded": len(vertices) == 2, "incidence_status": incidence,
            "unsupported_incidence": incidence in {"nonmanifold", "unresolved"},
            "angle_contribution": None, "curvature_contribution": None,
            "curvature_status": "NOT_CALCULATED_FROM_ARCHIVE_EDGE_FIELDS",
            "scientific_status": edge_status, "source_hash": archive_hash,
            "archive_identity": getattr(net, "archive_id", None),
        })
        edge_geometry.append((record_id, points))

    vertex_rows = []
    for index, row in net.verts.iterrows():
        index = int(index); generators = _ids(row.get("balls")); generator_records = atom_records(generators)
        edges, surfaces = _ids(row.get("edges")), _ids(row.get("surfs"))
        xyz = _array(row.get("loc"), dtype=float).reshape(-1)
        xyz = tuple(float(value) for value in xyz[:3]) if len(xyz) >= 3 and np.isfinite(xyz[:3]).all() else None
        refs_ok = all(edge in edge_map for edge in edges) and all(surface in surface_map for surface in surfaces)
        vertex_rows.append({
            "record_id": vertex_map[index], "feature_index": index, "xyz_A": xyz,
            "generator_ids": generators,
            "generator_stable_ids": tuple(record["generator_stable_id"] for record in generator_records),
            "source_atom_ids": tuple(record["source_atom_id"] for record in generator_records),
            "incident_edge_ids": tuple(edge_map[e] for e in edges if e in edge_map),
            "incident_surface_ids": tuple(surface_map[s] for s in surfaces if s in surface_map),
            "angle_defect": None, "angle_defect_unit": "rad", "curvature_status": "NOT_CALCULATED_FROM_ARCHIVE_VERTEX_FIELDS",
            "bounded": bool(edges) and all(len(_ids(net.edges.loc[e].get("verts"))) == 2 for e in edges if e in net.edges.index),
            "unsupported_incidence": not refs_ok, "scientific_status": "CERTIFIED" if xyz is not None and refs_ok else "UNRESOLVED",
            "source_hash": archive_hash, "archive_identity": getattr(net, "archive_id", None),
        })

    contacts = []
    for row in surface_rows:
        groups = tuple(row["generator_ordered_groups"])
        contacts.append({
            "contact_id": row["record_id"] + ":contact", "surface_id": row["record_id"],
            "ordered_group_a": groups[0] if groups else None, "ordered_group_b": groups[1] if len(groups) > 1 else None,
            "groups": groups, "group_pair": row["group_pair"], "area_A2": row["area_A2"],
            "contact_classification": row["contact_classification"],
            "derivation_status": "DERIVED_FROM_FULL_NETWORK_PRIMAL_SURFACES",
            "ordered_interface_status": "NOT_AN_INDEPENDENT_INTERFACE_SOLUTION",
            "scientific_status": row["scientific_status"],
        })

    summary_keys = set()
    for row in contacts:
        for pair in combinations(sorted(set(row["groups"])), 2):
            summary_keys.add(pair)
    interface_rows = []
    for pair in sorted(summary_keys):
        matching = [row for row in contacts if set(pair).issubset(set(row["groups"]))]
        interface_rows.append({
            "interface": " ↔ ".join(pair), "group_a": pair[0], "group_b": pair[1],
            "surface_count": len(matching), "area_A2_nonexclusive_sum": sum(row["area_A2"] or 0 for row in matching),
            "surface_membership": "NONEXCLUSIVE_FOR_MULTI_GENERATOR_SURFACES",
            "ordered_interface_solved": False, "derivation_status": "FULL_NETWORK_ADJACENCY_CLASSIFICATION",
            "scientific_status": "PARTIAL_CONTACT_CLASSIFICATION",
        })

    _write_csv(output_dir / "surfaces.csv", surface_rows)
    _write_csv(output_dir / "edges.csv", edge_rows)
    _write_csv(output_dir / "vertices.csv", vertex_rows)
    _write_csv(output_dir / "contacts.csv", contacts)
    _write_csv(output_dir / "interface_summary.csv", interface_rows)

    geometry_index = {"surfaces": {}, "edges": {}, "vertices": {}}
    surface_files = []
    for classification, parts in surface_geometry.items():
        filename = f"surfaces_{classification}.off"
        mapping = _write_surface_mesh(geometry_dir / filename, parts)
        surface_files.append(filename)
        geometry_index["surfaces"][filename] = mapping
    geometry_index["edges"] = _write_edges(geometry_dir / "edges.obj", edge_geometry)
    vertex_atoms = [{"xyz": tuple(row["xyz_A"] or (0, 0, 0)), "atom_name": "V",
                     "residue_name": "VTX", "residue_number": row["feature_index"],
                     "chain": "V", "element": "V"} for row in vertex_rows]
    geometry_index["vertices"] = _write_pdb(
        geometry_dir / "vertices.pdb", vertex_atoms,
        [row["record_id"] for row in vertex_rows],
    )
    _write_pdb(output_dir / "source_atoms.pdb", list(atoms.values()))
    vertex_atoms_by_index = {
        row["feature_index"]: atom for row, atom in zip(vertex_rows, vertex_atoms)
    }
    selected = next((row for row in surface_rows if row["contact_classification"] == "cross_group"), None)
    if selected is not None:
        selected_index = selected["feature_index"]
        selected_surface_raw = net.surfs.loc[selected_index]
        _write_surface_mesh(geometry_dir / "selected_surface.off", [
            (selected["record_id"], selected_surface_raw.get("points"), selected_surface_raw.get("tris"))
        ])
        selected_edges = [
            (edge_map[edge_id], net.edges.loc[edge_id].get("points"))
            for edge_id in _ids(selected_surface_raw.get("edges")) if edge_id in edge_map
        ]
        selected_vertices = [
            vertex_atoms_by_index[vertex_id] for vertex_id in _ids(selected_surface_raw.get("verts"))
            if vertex_id in vertex_atoms_by_index
        ]
        _write_edges(geometry_dir / "selected_edges.obj", selected_edges)
        _write_pdb(geometry_dir / "selected_vertices.pdb", selected_vertices)
        geometry_index["selected_surface"] = {
            "surface_id": selected["record_id"],
            "edge_ids": [item[0] for item in selected_edges],
            "vertex_ids": [vertex_map[item] for item in _ids(selected_surface_raw.get("verts")) if item in vertex_map],
        }
    _write_pml(output_dir / "visualization.pml", surface_files, selected is not None)
    _write_json(geometry_dir / "geometry_index.json", geometry_index)

    raw_area = sum(row["area_A2"] or 0.0 for row in surface_rows)
    group = (getattr(net.sys, "groups", None) or [None])[0]
    group_area = _number(getattr(group, "sa", None)) if group is not None else None
    group_volume = _number(getattr(group, "vol", None)) if group is not None else None
    group_gauss = _number(getattr(group, "int_gauss_curv", None)) if group is not None else None
    missing = {
        "edge_angle_or_curvature_contributions": len(edge_rows),
        "vertex_angle_defect": len(vertex_rows),
        "typed_results_snapshot": 1 if not getattr(net, "typed_results", None) else 0,
        "geometry_lineage_snapshot": 1 if not getattr(net, "geometry_lineage", None) else 0,
    }
    validation = {
        "schema": "molecular_interface_atlas.v1",
        "archive": {"path": str(archive_path), "sha256": archive_hash,
                    "archive_id": getattr(net, "archive_id", None), "source_network_id": getattr(net, "source_network_id", None)},
        "geometry_recalculated": False, "loader": "vorpy.src.io.network_archive.load_network",
        "settings": getattr(net, "settings", None),
        "counts": {"surfaces": len(surface_rows), "edges": len(edge_rows), "vertices": len(vertex_rows), "generators": len(atoms), "contacts": len(contacts)},
        "incidence": {
            "surface_edge_missing": sum(1 for row in surface_rows if row["unsupported_incidence"]),
            "edge_vertex_or_surface_missing": sum(1 for row in edge_rows if row["unsupported_incidence"]),
            "vertex_incidence_missing": sum(1 for row in vertex_rows if row["unsupported_incidence"]),
            "edge_status_counts": {status: sum(row["incidence_status"] == status for row in edge_rows) for status in ("complete", "boundary", "nonmanifold", "unresolved")},
        },
        "area_reconciliation": {
            "raw_surface_table_sum_A2": raw_area, "group_summary_sa_A2": group_area,
            "difference_A2": None if group_area is None else raw_area - group_area,
            "status": "NOT_RECONCILED_DIFFERENT_AGGREGATES" if group_area is not None and not math.isclose(raw_area, group_area, rel_tol=1e-9, abs_tol=1e-9) else "RECONCILED",
            "interpretation": "group.sa is a group-level aggregate; raw surfs.sa contains all archived surface rows and is not certified as the same physical boundary.",
        },
        "scientific_audit": {
            "group_volume_A3": group_volume, "group_surface_area_A2": group_area,
            "integrated_gaussian_curvature": group_gauss, "integrated_gaussian_curvature_unit": "1",
            "gaussian_convention": "intrinsic surface integral; dimensionless under the standard curvature convention",
            "gaussian_status": "RAW_GROUP_VALUE_PARTIAL_ARCHIVE_SUPPORT",
            "solvent_cells_complete": getattr(net, "solvent_cells_complete", None),
            "archive_solvent_identity_conflicts": sum(
                record["group_assignment_status"] == "CHEMICAL_IDENTITY_OVERRIDES_ARCHIVE_SOLVENT_FLAG"
                for record in atoms.values()
            ),
        },
        "missing_or_unsupported": missing,
        "interface_feasibility": {
            "directly_derivable": ["pairwise contact classification from full-network surface generators"],
            "not_derivable_as_solved_interface": ["protein_A↔DNA_B", "protein_A↔DNA_C", "protein_A↔DNA_B+C", "DNA_B↔DNA_C", "protein_or_DNA↔Zn"],
            "reason": "The archive is a full ordinary POW network, not an ordered two-group interface solve; pairwise rows are nonexclusive adjacency classifications.",
        },
    }
    _write_json(output_dir / "validation_report.json", validation)

    manifest = {
        "schema": "molecular_interface_atlas.v1", "exporter": "cache_export",
        "generated_from": str(archive_path), "archive_sha256": archive_hash,
        "archive_id": getattr(net, "archive_id", None), "source_path": None if source_path is None else str(Path(source_path).resolve()),
        "source_sha256": None if source_path is None else _sha256(source_path),
        "geometry_recalculated": False, "output_fail_if_exists": True,
        "settings": getattr(net, "settings", None), "counts": validation["counts"],
        "files": _file_manifest(output_dir),
        "limitations": validation["missing_or_unsupported"],
    }
    _write_json(output_dir / "atlas_manifest.json", manifest)
    return output_dir

