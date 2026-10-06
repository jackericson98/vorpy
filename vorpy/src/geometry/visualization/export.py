"""Visualization-only exports derived from solved geometry and alpha records.

The dual uses generator-center coordinates as an abstract geometric
realization. Physical surfaces, edges, and vertices are copied from the solved
network through VorPy's existing surface mesh and edge tube builders.
"""

from __future__ import annotations

import csv
import itertools
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from vorpy.src.output.draw import DEFAULT_EDGE_RADIUS, draw_joint, draw_line
from vorpy.src.output.mesh import MeshData, combine_mesh_parts, write_mesh
from vorpy.src.output.surfs import prepare_surfs

_COLORS = {
    "dual_full_edges": (0.35, 0.45, 0.55),
    "dual_full_triangles": (0.55, 0.62, 0.70),
    "dual_full_tetrahedra": (0.72, 0.75, 0.78),
    "alpha_edges": (1.0, 0.55, 0.05),
    "alpha_triangles": (1.0, 0.72, 0.25),
    "alpha_tetrahedra": (1.0, 0.82, 0.45),
    "interface_full_surfaces": (0.25, 0.75, 0.8),
    "interface_alpha_surfaces": (1.0, 0.15, 0.25),
    "interface_alpha_physical_edges": (0.9, 0.25, 0.8),
    "dual_mapping_edges": (1.0, 0.9, 0.1),
}


def export_dual_visualization(network, dual, output_dir, *, filtration=None,
                              alpha=None, max_dimension=None, interface_full=None,
                              interface_alpha=None, structure_path=None,
                              edge_radius=DEFAULT_EDGE_RADIUS):
    """Write independently toggleable full-dual, alpha, and interface layers.

    ``filtration.simplices_at(alpha)`` is the sole source of alpha membership.
    No birth or geometry is recalculated here. ``structure_path`` is optional;
    absent it, a pseudoatom PDB records the exact generator centers.
    """
    if dual.network is not network:
        raise ValueError("DualComplex must reference the supplied solved network.")
    if (filtration is None) != (alpha is None):
        raise ValueError("filtration and alpha must be supplied together.")
    if filtration is not None and filtration.network is not network:
        raise ValueError("AlphaFiltration must reference the supplied solved network.")
    if max_dimension is not None and not 0 <= int(max_dimension) <= 3:
        raise ValueError("max_dimension must be between 0 and 3")

    root = Path(output_dir).expanduser().resolve()
    scheme_root = root / dual.scheme
    dual_dir = scheme_root / "dual"
    dual_dir.mkdir(parents=True, exist_ok=True)
    coordinate_map = {int(row["generator_id"]): np.asarray(row["coordinates"], dtype=float)
                      for row in dual.generators}
    full_records = [simplex for dimension in range(4)
                    for simplex in dual.simplices[dimension].values()]
    _write_simplex_layers(dual_dir, "dual_full", full_records, coordinate_map, edge_radius)
    _write_simplex_tables(dual_dir, "dual_full", full_records, coordinate_map)
    generator_count = _write_generators(network, dual, dual_dir, coordinate_map)

    alpha_counts = {dimension: 0 for dimension in range(4)}
    blocked_count = 0
    alpha_dir = None
    if filtration is not None:
        alpha_subcomplex = filtration.simplices_at(float(alpha))
        records = [row for row in alpha_subcomplex.records
                   if max_dimension is None or row.dimension <= int(max_dimension)]
        blocked_count = sum(
            row.dimension >= 2 and (max_dimension is None or row.dimension <= int(max_dimension))
            for row in alpha_subcomplex.blocked
        )
        alpha_counts = {dimension: sum(row.dimension == dimension for row in records)
                        for dimension in range(4)}
        alpha_dir = scheme_root / f"alpha_{_alpha_slug(alpha)}"
        alpha_dir.mkdir(parents=True, exist_ok=True)
        _write_simplex_layers(alpha_dir, "alpha", records, coordinate_map, edge_radius)
        _write_simplex_tables(alpha_dir, "alpha", records, coordinate_map)
        _write_alpha_vertices(alpha_dir, records, coordinate_map)

    mapping_dir = scheme_root / "mapping"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    mapping_rows = _dual_mapping_rows(dual)
    _write_csv(mapping_dir / "dual_voronoi_mapping.csv", mapping_rows)

    interface_layers = {}
    if interface_full is not None:
        if interface_full.scheme != dual.scheme:
            raise ValueError("Full interface scheme must match the dual.")
        interface_layers.update(_write_interface_layers(
            network, interface_full, scheme_root / "interfaces" / "full",
            "interface_full", coordinate_map, edge_radius,
        ))
    if interface_alpha is not None:
        if interface_alpha.scheme != dual.scheme or interface_alpha.selection_mode != "alpha":
            raise ValueError("Alpha interface must match the dual scheme and use selection_mode='alpha'.")
        interface_layers.update(_write_interface_layers(
            network, interface_alpha, scheme_root / f"alpha_{_alpha_slug(alpha if alpha is not None else interface_alpha.alpha)}",
            "interface_alpha", coordinate_map, edge_radius,
        ))
        _write_csv(mapping_dir / "alpha_interface_surface_mapping.csv",
                   interface_alpha.selection_mappings)

    pml = scheme_root / "visualize.pml"
    _write_pymol_script(pml, structure_path, dual_dir, alpha_dir, interface_layers)
    metadata = {
        "scheme": dual.scheme,
        "scheme_name": {"prm": "ordinary Voronoi / Delaunay", "pow": "Power-Laguerre / regular", "aw": "AW Apollonius"}.get(dual.scheme, dual.scheme),
        "dual_geometry": "straight generator-center realization of abstract dual simplices",
        "aw_warning": "AW dual simplices are visualization only; physical AW geometry is exported from the solved network." if dual.scheme == "aw" else None,
        "alpha_native": None if alpha is None else float(alpha),
        "alpha_units": None if filtration is None else filtration.metadata.get("native_alpha_units"),
        "alpha_convention": None if filtration is None else filtration.metadata.get("filtration_kind"),
        "alpha_experimental": bool(filtration and filtration.metadata.get("experimental", False)),
        "full_dual_counts": {str(k): len(dual.simplices[k]) for k in range(4)},
        "alpha_counts": {str(k): alpha_counts[k] for k in range(4)} if filtration is not None else None,
        "alpha_blocked_records_not_rendered": blocked_count,
        "alpha_higher_dimensions_not_calculated": (
            int(filtration.metadata.get("birth_dimensions_calculated", 3))
            if filtration is not None else None
        ),
        "alpha_uncomputed_dual_triangles_not_rendered": (
            len(dual.simplices[2]) if filtration is not None and
            int(filtration.metadata.get("birth_dimensions_calculated", 3)) < 2 else 0
        ),
        "alpha_uncomputed_dual_tetrahedra_not_rendered": (
            len(dual.simplices[3]) if filtration is not None and
            int(filtration.metadata.get("birth_dimensions_calculated", 3)) < 3 else 0
        ),
        "generator_count": generator_count,
        "rendered_vertices": {"full": generator_count,
                              "alpha": alpha_counts[0] if filtration is not None else None},
        "mapped_physical_features": len(mapping_rows),
        "dual_audit": "not run by visualization exporter; use DualComplex.audit() explicitly when required",
        "unresolved_dual_features": len(dual.unresolved_features),
        "interface_layers": interface_layers,
        "pymol_script": str(pml),
    }
    for directory in (dual_dir, alpha_dir, mapping_dir):
        if directory is not None:
            (directory / "visualization_metadata.json").write_text(
                json.dumps(metadata, indent=2, sort_keys=True, default=str) + "\n",
                encoding="utf-8",
            )
    return metadata


def export_alpha_complex_visualization(network, filtration, alpha, output_dir, *,
                                       dual=None, max_dimension=None, **kwargs):
    """Export an alpha query and its full dual context as independent layers."""
    if dual is None:
        from vorpy.src.geometry.duals import build_dual
        dual = build_dual(network)
    return export_dual_visualization(
        network, dual, output_dir, filtration=filtration, alpha=alpha,
        max_dimension=max_dimension, **kwargs,
    )


def export_dual_voronoi_mapping(network, dual, interface, output_dir, *,
                                alpha_filtration=None, edge_radius=DEFAULT_EDGE_RADIUS):
    """Export physical interface patches and dual-to-primal feature records."""
    if interface.scheme != dual.scheme or dual.network is not network:
        raise ValueError("Network, dual, and interface must belong to the same solve/scheme.")
    directory = Path(output_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    rows = _dual_mapping_rows(dual)
    _write_csv(directory / "dual_voronoi_mapping.csv", rows)
    physical_layers = _write_interface_layers(
        network, interface, directory,
        "interface_alpha" if interface.selection_mode == "alpha" else "interface_full",
        {int(row["generator_id"]): np.asarray(row["coordinates"], dtype=float)
         for row in dual.generators}, edge_radius,
    )
    if interface.selection_mode == "alpha":
        _write_csv(directory / "alpha_interface_surface_mapping.csv", interface.selection_mappings)
    else:
        _write_csv(directory / "full_interface_surface_mapping.csv",
                   _interface_mapping_rows(interface))
    result = {
        "scheme": dual.scheme,
        "selection": interface.selection_mode,
        "alpha_native": interface.alpha,
        "alpha_units": (alpha_filtration.metadata.get("native_alpha_units")
                        if alpha_filtration is not None else None),
        "dual_simplex_count": len(rows),
        "mapped_physical_features": sum(
            bool(row.get("mapped", row.get("physical_feature_id") is not None))
            for row in rows
        ),
        "unresolved_features": sum(not bool(row["supported"]) for row in rows),
        "interface_pairs": len(interface.selected_surfaces),
        "interface_area_A2": interface.area_A2,
        "layers": physical_layers,
        "aw_geometry_note": "Straight dual edges show center connectivity; exported surface/edge meshes are solved AW geometry." if dual.scheme == "aw" else None,
    }
    (directory / "visualization_metadata.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return result


def visualize_simplex_mapping(network, dual, generator_tuple, output_dir, *,
                              edge_radius=DEFAULT_EDGE_RADIUS):
    """Export one abstract simplex and its linked solved physical feature(s)."""
    ids = tuple(sorted(map(int, generator_tuple)))
    if not 1 <= len(ids) <= 4:
        raise ValueError("A dual simplex must contain one through four generators.")
    simplex = dual.simplex(len(ids) - 1, ids)
    if simplex is None:
        raise KeyError(f"No dual simplex exists for generator tuple {ids}.")
    directory = Path(output_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    coords = {int(row["generator_id"]): np.asarray(row["coordinates"], dtype=float)
              for row in dual.generators}
    layers = _write_simplex_layers(directory, "selected_simplex", [simplex], coords, edge_radius)
    features = []
    for feature in simplex.primal_features:
        features.append({"kind": feature.kind, "feature_id": feature.feature_id,
                         "bounded": feature.bounded, "complete": feature.complete,
                         "supported": bool(simplex.supported and feature.bounded
                                           and feature.complete
                                           and feature.generator_cells_complete),
                         "generator_tuple": ids})
        if feature.kind == "surf":
            position = _position_by_feature(network.surfs, feature.feature_id)
            if position is not None:
                mesh = _surface_mesh(network, [position])
                if mesh is not None:
                    write_mesh(mesh, "physical_surface.off", directory=directory)
        elif feature.kind == "edge":
            position = _position_by_feature(network.edges, feature.feature_id)
            if position is not None:
                mesh = _physical_edge_mesh(network, [position], edge_radius)
                if mesh is not None:
                    write_mesh(mesh, "physical_edge.off", directory=directory)
        elif feature.kind == "vert":
            row = network.verts.iloc[_position_by_feature(network.verts, feature.feature_id)]
            xyz = np.asarray(row.get("loc", row.get("point", ())), dtype=float)
            if xyz.shape == (3,):
                pts, tris = draw_joint(xyz, radius=edge_radius * 2)
                write_mesh(MeshData(np.asarray(pts), np.asarray(tris)),
                           "physical_vertex.off", directory=directory)
    result = {"scheme": dual.scheme, "generator_tuple": ids,
              "simplex_id": simplex.simplex_id, "dimension": simplex.dimension,
              "dual_center_coordinates": simplex.generator_coordinates,
              "physical_features": features, "layers": layers,
              "aw_note": "Center-connected simplex is abstract visualization; physical feature remains solved AW geometry." if dual.scheme == "aw" else None}
    (directory / "selected_simplex_mapping.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return result


def _write_simplex_layers(directory, prefix, records, coordinates, edge_radius):
    by_dim = {dimension: [row for row in records if row.dimension == dimension]
              for dimension in range(4)}
    layer_counts = {str(dim): len(by_dim[dim]) for dim in range(4)}
    all_pairs = set()
    for record in records:
        all_pairs.update(itertools.combinations(_record_ids(record), 2))
    edge_mesh = _edge_mesh(sorted(all_pairs), coordinates, edge_radius)
    edge_name = "dual_full_edges" if prefix == "dual_full" else "alpha_edges" if prefix == "alpha" else f"{prefix}_edges"
    if edge_mesh is not None:
        write_mesh(edge_mesh, f"{edge_name}.off", directory=directory)
    triangle_mesh = _triangle_mesh(by_dim[2], coordinates)
    triangle_name = "dual_full_triangles" if prefix == "dual_full" else "alpha_triangles" if prefix == "alpha" else f"{prefix}_triangles"
    if triangle_mesh is not None:
        write_mesh(triangle_mesh, f"{triangle_name}.off", directory=directory)
    tetra_edges = set()
    tetra_faces = []
    for record in by_dim[3]:
        ids = _record_ids(record)
        tetra_edges.update(itertools.combinations(ids, 2))
        tetra_faces.extend(itertools.combinations(ids, 3))
    tetra_meshes = []
    tetra_edge_mesh = _edge_mesh(sorted(tetra_edges), coordinates, edge_radius)
    tetra_face_mesh = _triangle_mesh_from_tuples(tetra_faces, coordinates)
    if tetra_edge_mesh is not None:
        write_mesh(tetra_edge_mesh, f"{prefix}_tetrahedra_edges.off", directory=directory)
        tetra_meshes.append(f"{prefix}_tetrahedra_edges.off")
    if tetra_face_mesh is not None:
        write_mesh(tetra_face_mesh, f"{prefix}_tetrahedra_faces.off", directory=directory)
        tetra_meshes.append(f"{prefix}_tetrahedra_faces.off")
    return {"counts": layer_counts,
            "edges_file": None if edge_mesh is None else f"{edge_name}.off",
            "triangles_file": None if triangle_mesh is None else f"{triangle_name}.off",
            "tetrahedra_files": tetra_meshes}


def _write_simplex_tables(directory, prefix, records, coordinates):
    for dimension, name in enumerate(("vertices", "edges", "triangles", "tetrahedra")):
        rows = []
        for record in records:
            if record.dimension != dimension:
                continue
            ids = tuple(record.generator_ids if hasattr(record, "generator_ids")
                        else record.generator_tuple)
            row = {"dimension": dimension,
                   "simplex_id": getattr(record, "simplex_id", f"{dimension}:" + ",".join(map(str, ids))),
                   "generator_tuple": ids,
                   "visualization_kind": "straight generator-center realization"}
            for i, generator_id in enumerate(ids):
                xyz = coordinates[int(generator_id)]
                row.update({f"generator_{i}": int(generator_id), f"x{i}_A": float(xyz[0]),
                            f"y{i}_A": float(xyz[1]), f"z{i}_A": float(xyz[2])})
            rows.append(row)
        _write_csv(Path(directory) / f"{prefix}_{name}.csv", rows)


def _write_alpha_vertices(directory, records, coordinates):
    vertices = [tuple(row.generator_tuple) for row in records if row.dimension == 0]
    path = Path(directory) / "alpha_vertices.pdb"
    with path.open("w", encoding="ascii", newline="\n") as stream:
        for serial, (generator_id,) in enumerate(vertices, 1):
            xyz = coordinates[int(generator_id)]
            stream.write(f"HETATM{serial % 100000:5d}  C   ALP G{generator_id % 9999:4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          C \n")
        stream.write("END\n")
    return len(vertices)


def _edge_mesh(pairs, coordinates, radius):
    parts_points, parts_tris = [], []
    for pair in pairs:
        points, triangles = draw_line([coordinates[int(pair[0])], coordinates[int(pair[1])]], radius=radius)
        parts_points.append(points)
        parts_tris.append(triangles)
    if not parts_points:
        return None
    return combine_mesh_parts(parts_points, parts_tris)


def _triangle_mesh(records, coordinates):
    return _triangle_mesh_from_tuples([_record_ids(row) for row in records], coordinates)


def _record_ids(record):
    return tuple(record.generator_ids if hasattr(record, "generator_ids")
                 else record.generator_tuple)


def _triangle_mesh_from_tuples(tuples, coordinates):
    points, triangles = [], []
    for ids in tuples:
        if len(ids) != 3:
            continue
        start = len(points)
        points.extend(coordinates[int(index)] for index in ids)
        triangles.append((start, start + 1, start + 2))
    if not triangles:
        return None
    return MeshData(np.asarray(points), np.asarray(triangles))


def _write_generators(network, dual, directory, coordinates):
    rows = []
    balls = network.balls
    with (directory / "generators.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["generator_id", "x_A", "y_A", "z_A", "radius_A", "weight_A2"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for item in dual.generators:
            idx = int(item["generator_id"]); xyz = coordinates[idx]
            writer.writerow({"generator_id": idx, "x_A": xyz[0], "y_A": xyz[1], "z_A": xyz[2],
                             "radius_A": item["radius"], "weight_A2": item["weight"]})
            rows.append((idx, xyz, item["radius"]))
    with (directory / "generators.pdb").open("w", encoding="ascii", newline="\n") as stream:
        for serial, (idx, xyz, _radius) in enumerate(rows, 1):
            atom_name = "C"
            if "name" in balls.columns:
                atom_name = str(balls.loc[idx, "name"]).strip()[:4] or "C"
            resname = "GEN"
            chain = "G"
            if "res" in balls.columns:
                resname = str(balls.loc[idx, "res"]).strip()[:3] or "GEN"
            if "chn" in balls.columns:
                chain = str(balls.loc[idx, "chn"]).strip()[:1] or "G"
            stream.write(f"HETATM{serial % 100000:5d} {atom_name:>4s} {resname:>3s} {chain}{(idx % 9999):4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          C \n")
        stream.write("END\n")
    return len(rows)


def _dual_mapping_rows(dual):
    rows = []
    for dimension in range(4):
        for ids, simplex in sorted(dual.simplices[dimension].items()):
            for feature in simplex.primal_features:
                rows.append({"scheme": dual.scheme, "dimension": dimension,
                             "generator_tuple": ids, "simplex_id": simplex.simplex_id,
                             "physical_feature_kind": feature.kind,
                             "physical_feature_id": feature.feature_id,
                             "bounded": feature.bounded, "complete": feature.complete,
                             "generator_cells_complete": feature.generator_cells_complete,
                             "supported": simplex.supported,
                             "status": simplex.status,
                             "dual_visualization": "straight generator-center realization",
                             "physical_geometry": "solved VorPy Voronoi network"})
    return rows


def _write_interface_layers(network, interface, directory, prefix, coordinates, edge_radius):
    directory.mkdir(parents=True, exist_ok=True)
    selected = interface.selected_surfaces
    surf_positions = [_position_by_feature(network.surfs, row.feature_id) for row in selected]
    surf_positions = [index for index in surf_positions if index is not None]
    mesh = _surface_mesh(network, surf_positions)
    surface_file = f"{prefix}_surfaces.off"
    if mesh is not None:
        write_mesh(mesh, surface_file, directory=directory)
    edge_ids = {edge.feature_id for edge in interface.interface_edges}
    edge_positions = [_position_by_feature(network.edges, edge_id) for edge_id in sorted(edge_ids)]
    edge_positions = [index for index in edge_positions if index is not None]
    edge_mesh = _physical_edge_mesh(network, edge_positions, edge_radius)
    edge_file = f"{prefix}_physical_edges.off"
    if edge_mesh is not None:
        write_mesh(edge_mesh, edge_file, directory=directory)
    center_pairs = sorted({tuple(sorted(surface.generator_ids)) for surface in selected})
    mapping_mesh = _edge_mesh(center_pairs, coordinates, edge_radius)
    mapping_edge_file = "dual_mapping_edges.off"
    if mapping_mesh is not None:
        write_mesh(mapping_mesh, mapping_edge_file, directory=directory)
    vertex_rows = [row for row in interface.interface_vertices if row.xyz is not None]
    vertex_parts, vertex_faces = [], []
    for row in vertex_rows:
        pts, tris = draw_joint(row.xyz, radius=edge_radius * 2)
        vertex_parts.append(pts); vertex_faces.append(tris)
    vertex_mesh = combine_mesh_parts(vertex_parts, vertex_faces) if vertex_parts else None
    vertex_file = f"{prefix}_physical_vertices.off"
    if vertex_mesh is not None:
        write_mesh(vertex_mesh, vertex_file, directory=directory)
    mapping_file = directory / f"{prefix}_dual_edges.csv"
    _write_csv(mapping_file, (interface.selection_mappings
                              if interface.selection_mappings
                              else _interface_mapping_rows(interface)))
    for name, records in (
        ("interface_surfaces.csv", interface.surfaces),
        ("interface_edges.csv", interface.edges),
        ("interface_vertices.csv", interface.vertices),
    ):
        _write_csv(directory / name, [asdict(record) for record in records])
    _write_csv(directory / "alpha_incidence_audit.csv", interface.alpha_incidence_audit)
    return {prefix: {"directory": str(directory), "surfaces": surface_file if mesh is not None else None,
                     "physical_edges": edge_file if edge_mesh is not None else None,
                     "dual_mapping_edges": mapping_edge_file if mapping_mesh is not None else None,
                     "physical_vertices": vertex_file if vertex_mesh is not None else None,
                     "surface_count": len(selected), "edge_count": len(edge_ids),
                     "mapping_csv": mapping_file.name}}


def _surface_mesh(network, positions):
    if not positions:
        return None
    # prepare_surfs is VorPy's shared physical-surface mesh builder. Give it a
    # shallow network clone with a copied surface table so export-side color
    # caches cannot mutate the solved network.
    import copy
    clone = copy.copy(network)
    clone.surfs = network.surfs.copy(deep=True)
    if "tri_colors" in clone.surfs.columns:
        clone.surfs["tri_colors"] = [[] for _ in range(len(clone.surfs))]
    return prepare_surfs(clone, positions, color=(0.25, 0.75, 0.8),
                         concave_colors=False, ref_surfs=[], universal_max=True,
                         color_scheme="none")


def _physical_edge_mesh(network, positions, radius):
    parts_points, parts_tris = [], []
    for index in positions:
        points = network.edges.iloc[index].get("points", ())
        if len(points) >= 2:
            vertices, faces = draw_line(points, radius=radius)
            parts_points.append(vertices); parts_tris.append(faces)
    return combine_mesh_parts(parts_points, parts_tris) if parts_points else None


def _position_by_feature(table, feature_id):
    try:
        return int(table.index.get_loc(feature_id))
    except (KeyError, TypeError, ValueError):
        return None


def _interface_mapping_rows(interface):
    return [{
        "scheme": interface.scheme,
        "selection": interface.selection_mode,
        "alpha_native": interface.alpha,
        "generator_i": surface.generator_ids[0],
        "generator_j": surface.generator_ids[1],
        "dual_edge_id": surface.dual_simplex_id,
        "surface_id": surface.feature_id,
        "mapped": True,
        "surface_complete": surface.complete,
        "surface_supported": surface.supported,
        "exclusion_reason": "",
    } for surface in interface.selected_surfaces]


def _write_csv(path, rows):
    rows = list(rows)
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        if not keys:
            return
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows({key: json.dumps(value, separators=(",", ":"))
                          if isinstance(value, (tuple, list, dict, set, frozenset)) else value
                          for key, value in row.items()} for row in rows)


def _alpha_slug(alpha):
    value = f"{float(alpha):g}".replace("-", "m").replace(".", "p")
    return f"{value}"


def _write_pymol_script(path, structure_path, dual_dir, alpha_dir, interface_layers):
    structure = Path(structure_path).resolve() if structure_path else dual_dir / "generators.pdb"
    lines = [
        "from pymol import cmd",
        "from pymol.cgo import BEGIN, END, TRIANGLES, COLOR, VERTEX",
        "import os",
        f"cmd.load({str(structure)!r}, 'generators')",
        "cmd.hide('everything', 'generators')",
        "cmd.show('spheres', 'generators')",
        "cmd.set('sphere_scale', 0.22, 'generators')",
        "cmd.color('gray70', 'generators')",
        "def _load_off_layer(path, name, color):",
        "    if not path or not os.path.exists(path): return",
        "    with open(path, 'r', encoding='utf-8') as f: lines=f.readlines()",
        "    if not lines or lines[0].strip() != 'OFF': return",
        "    nv,nf,_=map(int, lines[1].split()); pts=[tuple(map(float, x.split())) for x in lines[2:2+nv]]",
        "    obj=[]; obj.extend([BEGIN, TRIANGLES, COLOR, *color])",
        "    for line in lines[2+nv:2+nv+nf]:",
        "        fields=line.split(); n=int(fields[0]); ids=list(map(int, fields[1:1+n]))",
        "        for k in range(1, n-1): obj.extend([VERTEX, *pts[ids[0]], VERTEX, *pts[ids[k]], VERTEX, *pts[ids[k+1]]])",
        "    obj.append(END); cmd.load_cgo(obj, name)",
    ]
    if alpha_dir is not None and (alpha_dir / "alpha_vertices.pdb").exists():
        lines.extend([
            f"cmd.load({str(alpha_dir / 'alpha_vertices.pdb')!r}, 'alpha_vertices')",
            "cmd.show('spheres', 'alpha_vertices')",
            "cmd.set('sphere_scale', 0.35, 'alpha_vertices')",
            "cmd.color('orange', 'alpha_vertices')",
        ])
    scheme = dual_dir.parent.name
    for filename, name, color in _mesh_layers(dual_dir, "dual_full"):
        name = f"{scheme}_{name}"
        lines.append(f"_load_off_layer({str(filename)!r}, {name!r}, {color!r})")
    if alpha_dir is not None:
        for filename, name, color in _mesh_layers(alpha_dir, "alpha"):
            name = f"{scheme}_{name}"
            lines.append(f"_load_off_layer({str(filename)!r}, {name!r}, {color!r})")
    for layer_key, layer_info in interface_layers.items():
        directory = Path(layer_info["directory"])
        for key in ("surfaces", "physical_edges", "dual_mapping_edges", "physical_vertices"):
            filename = layer_info.get(key)
            if filename:
                name = f"{scheme}_{layer_key}_{key}"
                if key == "dual_mapping_edges":
                    color = _COLORS["dual_mapping_edges"]
                else:
                    color = _COLORS["interface_alpha_surfaces" if layer_key == "interface_alpha"
                                    else "interface_full_surfaces"]
                lines.append(f"_load_off_layer({str(directory / filename)!r}, {name!r}, {color!r})")
    lines.append("cmd.zoom('generators')")
    path = Path(path)
    python_path = path.with_suffix(".py")
    python_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.write_text(f'run "{python_path.resolve().as_posix()}"\n', encoding="utf-8")
    _write_pymol_views(path, lines, scheme, alpha_dir, interface_layers)


def _write_pymol_views(master_path, master_lines, scheme, alpha_dir, interface_layers):
    """Provide ready-to-open full, alpha, and mapped-interface PyMOL views."""
    directory = Path(master_path).parent / "views"
    directory.mkdir(parents=True, exist_ok=True)
    full_objects = [f"{scheme}_dual_full_{name}" for name in
                    ("edges", "triangles", "tetrahedra_edges", "tetrahedra_faces")]
    alpha_objects = ([f"{scheme}_alpha_{name}" for name in
                      ("edges", "triangles", "tetrahedra_edges", "tetrahedra_faces")]
                     if alpha_dir is not None else [])
    if alpha_dir is not None and (Path(alpha_dir) / "alpha_vertices.pdb").exists():
        alpha_objects.append("alpha_vertices")
    mapping_objects = []
    for layer_key, layer_info in interface_layers.items():
        for key in ("surfaces", "physical_edges", "dual_mapping_edges", "physical_vertices"):
            if layer_info.get(key):
                mapping_objects.append(f"{scheme}_{layer_key}_{key}")
    views = {
        "view_A_full_dual.pml": full_objects,
        "view_B_alpha_complex.pml": alpha_objects,
        "view_C_dual_to_physical_mapping.pml": mapping_objects,
    }
    for filename, visible in views.items():
        names = ["generators", *visible]
        lines = [*master_lines, "cmd.disable('all')"]
        lines.extend(f"cmd.enable({name!r})" for name in names)
        pml_path = directory / filename
        python_path = pml_path.with_suffix(".py")
        python_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        pml_path.write_text(f'run "{python_path.resolve().as_posix()}"\n', encoding="utf-8")


def _mesh_layers(directory, prefix):
    names = (["dual_full_edges.off", "dual_full_triangles.off", "dual_full_tetrahedra_edges.off", "dual_full_tetrahedra_faces.off"]
             if prefix == "dual_full" else ["alpha_edges.off", "alpha_triangles.off", "alpha_tetrahedra_edges.off", "alpha_tetrahedra_faces.off"])
    return [(Path(directory) / name, Path(name).stem, _COLORS.get(Path(name).stem, _COLORS.get(prefix + "_edges", (0.8, 0.5, 0.1))))
            for name in names if (Path(directory) / name).exists()]
