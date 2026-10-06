"""Controlled 2KAI matched-probe Power/AW restriction audit.

Analysis-only: reads the saved H-tier AW network, uses its loaded radii as the
shared base model, and delegates births and curvature to existing analyzers.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

NETWORK = Path("output/comparative_study/cases/2KAI/H/full_aw.vpy")
RADIUS_AUDIT = Path("output/comparative_study/cases/2KAI/H/radii_audit.csv")
OUTPUT = Path("output/comparative_study/matched_probe_control/2KAI")
PDB = Path("vorpy/data/2KAI.pdb")
PROBE = 1.4
AW_ALPHA = 1.4


def _csv(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(",", ":"), allow_nan=False)
                             if isinstance(v, (list, tuple, dict, set)) else v
                             for k, v in row.items()})


def _atom_maps(path):
    aw, power = {}, {}
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            ident = tuple(json.loads(row["atom_id"]))
            item = {"atom_id": ident, "stable": "|".join(map(str, ident)),
                    "chain": ident[0], "resnum": ident[1], "icode": ident[2],
                    "resname": ident[3], "atom": ident[4],
                    "aw_base": float(row["aw_radius_A"]),
                    "power_base": float(row["power_base_radius_A"])}
            i, j = int(row["aw_generator_id"]), int(row["power_generator_id"])
            if i in aw or j in power:
                raise ValueError(f"Duplicate identity in radius audit: {ident}")
            aw[i], power[j] = item, item
    return aw, power


def _pdb_elements(path):
    from vorpy.src.analyze.comparative_geometry_study import canonical_identity, heavy_atoms, read_pdb
    atoms, _headers = read_pdb(path)
    return {canonical_identity(atom, {}): atom.element for atom in heavy_atoms(atoms)}


def validate_matched_inputs(atom_ids_power, atom_ids_aw, coordinates_power,
                            coordinates_aw, base_radii_power, base_radii_aw,
                            probe=PROBE, aw_alpha=AW_ALPHA):
    """Fail before geometry unless stable IDs, coordinates, and radii match."""
    if list(atom_ids_power) != list(atom_ids_aw):
        raise AssertionError("Power/AW stable atom identities differ.")
    if not np.array_equal(np.asarray(coordinates_power), np.asarray(coordinates_aw)):
        raise AssertionError("Power/AW coordinates differ.")
    rp, ra = np.asarray(base_radii_power, dtype=float), np.asarray(base_radii_aw, dtype=float)
    if not np.array_equal(rp, ra):
        raise AssertionError("Power/AW base radii differ.")
    if not np.array_equal(rp + float(probe), ra + float(aw_alpha)):
        raise AssertionError("Power and AW effective physical radii differ.")
    return True


def pair_partition(power_pairs, aw_pairs):
    """Return the exact intersection/disjoint decomposition used in reports."""
    power_pairs = {tuple(sorted(map(int, p))) for p in power_pairs}
    aw_pairs = {tuple(sorted(map(int, p))) for p in aw_pairs}
    shared, power_only, aw_only = power_pairs & aw_pairs, power_pairs - aw_pairs, aw_pairs - power_pairs
    if shared | power_only | aw_only != power_pairs | aw_pairs:
        raise AssertionError("Pair partition is incomplete.")
    if power_only & aw_pairs or aw_only & power_pairs:
        raise AssertionError("Pair partition is not disjoint.")
    return shared, power_only, aw_only


def bicolor_pairs(pairs, group_a, group_b):
    a, b = set(map(int, group_a)), set(map(int, group_b))
    return {tuple(sorted(map(int, pair))) for pair in pairs
            if (int(pair[0]) in a and int(pair[1]) in b)
            or (int(pair[1]) in a and int(pair[0]) in b)}


def _aw_surfaces(network):
    """Return every solved pair surface with a reference to its true geometry."""
    from vorpy.src.analyze.apollonius import PrimalFeatureRef
    result = {}
    cells = {int(i): bool(row.get("complete", True)) for i, row in network.balls.iterrows()}
    for sid, row in network.surfs.iterrows():
        pair = tuple(sorted(int(x) for x in row.get("balls", ())))
        if len(pair) != 2 or pair[0] == pair[1]:
            continue
        bounded = bool(row.get("edges", ())) and bool(row.get("verts", ()))
        complete = bounded
        cell_complete = all(cells.get(i, False) for i in pair)
        ref = PrimalFeatureRef("surf", int(sid), pair, bounded, complete,
                               "supported" if cell_complete else "incomplete_network_geometry",
                               generator_cells_complete=cell_complete)
        result.setdefault(pair, []).append((ref, row))
    return result


def _power_areas(result, pairs):
    from itertools import combinations
    from collections import defaultdict
    from vorpy.src.analyze.cazals_validation import power_facet_polygon_area
    edge_triangles, triangle_tetrahedra = {}, defaultdict(list)
    for tri in result.full_simplices[2]:
        for edge in combinations(tri, 2):
            edge_triangles.setdefault(tuple(sorted(edge)), []).append(tuple(sorted(tri)))
    for tet in result.full_simplices[3]:
        for tri in combinations(tet, 3):
            triangle_tetrahedra[tuple(sorted(tri))].append(tuple(sorted(tet)))
    areas = {}
    for pair in pairs:
        facet = SimpleNamespace(generator_ids=pair,
                                incident_regular_triangles=tuple(sorted(edge_triangles.get(pair, ()))))
        areas[pair] = power_facet_polygon_area(result, facet,
                                                triangle_tetrahedra=triangle_tetrahedra)
    return areas


def _dist_stats(rows):
    x = [float(r["center_distance_A"]) for r in rows]
    if not x:
        return {"count": 0, "minimum_A": None, "median_A": None, "mean_A": None,
                "p95_A": None, "maximum_A": None}
    return {"count": len(x), "minimum_A": min(x), "median_A": float(np.median(x)),
            "mean_A": float(np.mean(x)), "p95_A": float(np.percentile(x, 95)),
            "maximum_A": max(x)}


def _aw_metrics(aw_result, chosen, area):
    selected_surfs = {int(s.surface_id): s for s in aw_result.surfaces
                      if tuple(sorted(s.generator_ids)) in chosen}
    smooth_unknown = sum((not s.included) or s.integrated_mean_curvature is None
                         for s in selected_surfs.values())
    smooth = None if smooth_unknown else sum(float(s.integrated_mean_curvature)
                                              for s in selected_surfs.values())
    selected_ids = set(selected_surfs)
    signed = unsigned = positive = negative = boundary_length = 0.0
    interior = boundary = unresolved = 0
    for edge in aw_result.edges:
        n = len(set(edge.surface_ids) & selected_ids)
        if n == 1:
            boundary += 1
            if edge.length is not None:
                boundary_length += float(edge.length)
        elif n == 2:
            if edge.status != "included" or edge.signed_contribution is None:
                unresolved += 1
            else:
                interior += 1
                signed += float(edge.signed_contribution)
                unsigned += float(edge.unsigned_contribution or 0.0)
                positive += float(edge.positive_contribution or 0.0)
                negative += float(edge.negative_contribution or 0.0)
    edge_h = .5 * signed
    combined = None if smooth is None else smooth + edge_h
    return {"scheme": "aw", "area_A2": area, "smooth_H_A": smooth,
            "raw_signed_edge_A_rad": signed, "raw_unsigned_edge_A_rad": unsigned,
            "positive_edge_A_rad": positive, "negative_edge_A_rad": negative,
            "conventional_edge_H_A": edge_h, "combined_H_A": combined,
            "H_per_area_inv_A": combined / area if combined is not None and area else None,
            "interior_edges": interior, "boundary_edges": boundary,
            "boundary_length_A": boundary_length, "unresolved_edges": unresolved,
            "unresolved_smooth_surfaces": smooth_unknown}


def _power_metrics(result, chosen, areas):
    finite = {p for p in chosen if areas.get(p, {}).get("status") == "finite_bounded"}
    area = sum(float(areas[p]["area_A2"]) for p in finite)
    signed = unsigned = positive = negative = boundary_length = 0.0
    interior = boundary = unresolved = 0
    for edge in result.edges:
        n = len(set(edge.bicolor_facets) & finite)
        if n == 1:
            boundary += 1
            if edge.length is not None:
                boundary_length += float(edge.length)
        elif n == 2:
            if edge.geometry_status != "finite" or edge.length_beta is None:
                unresolved += 1
            else:
                interior += 1
                value = float(edge.length_beta)
                signed += value; unsigned += float(edge.length_abs_beta or 0.0)
                positive += max(0.0, value); negative += min(0.0, value)
    edge_h = .5 * signed
    return {"scheme": "power", "area_A2": area, "smooth_H_A": 0.0,
            "raw_signed_edge_A_rad": signed, "raw_unsigned_edge_A_rad": unsigned,
            "positive_edge_A_rad": positive, "negative_edge_A_rad": negative,
            "conventional_edge_H_A": edge_h, "combined_H_A": edge_h,
            "H_per_area_inv_A": edge_h / area if area else None,
            "interior_edges": interior, "boundary_edges": boundary,
            "boundary_length_A": boundary_length, "unresolved_edges": unresolved,
            "unresolved_smooth_surfaces": 0}


def _export_visualization(network, power_rows, aw_rows, shared, power_only, aw_only, out, pdb_path):
    """Use VorPy's existing edge and physical-surface mesh builders."""
    from vorpy.src.geometry.visualization.export import _edge_mesh, _surface_mesh
    from vorpy.src.output.mesh import MeshData, combine_mesh_parts, write_mesh
    root = Path(out) / "visualization"; root.mkdir(parents=True, exist_ok=True)
    xyz = {int(i): np.asarray(row["loc"], dtype=float) for i, row in network.balls.iterrows()}
    pfull = {tuple(sorted((r["generator_i"], r["generator_j"]))) for r in power_rows}
    psel = {tuple(sorted((r["generator_i"], r["generator_j"]))) for r in power_rows if r["restricted"]}
    afull = {tuple(sorted((r["generator_i"], r["generator_j"]))) for r in aw_rows}
    asel = {tuple(sorted((r["generator_i"], r["generator_j"]))) for r in aw_rows if r["restriction_status"] == "retained"}
    layers = {"full_Power_dual": pfull, "restricted_Power_dual": psel,
              "full_AW_dual": afull, "restricted_AW_dual": asel,
              "shared_pairs": shared, "Power_only_pairs": power_only, "AW_only_pairs": aw_only}
    for name, pairs in layers.items():
        mesh = _edge_mesh(sorted(pairs), xyz, .07)
        if mesh is not None:
            write_mesh(mesh, f"{name}.off", directory=root)
    pair_surfaces = {}
    for sid, row in network.surfs.iterrows():
        pair = tuple(sorted(int(x) for x in row.get("balls", ())))
        if len(pair) == 2:
            pair_surfaces.setdefault(pair, []).append(int(sid))
    positions = [network.surfs.index.get_loc(sid) for p in sorted(shared)
                 for sid in pair_surfaces.get(p, ())]
    mesh = _surface_mesh(network, positions)
    if mesh is not None:
        write_mesh(mesh, "shared_AW_physical_surfaces.off", directory=root)
    # The Power polygon coordinates come from the same regular-complex power
    # vertices used by the validated facet-area routine.
    power_parts = []
    for row in power_rows:
        pair = tuple(sorted((int(row["generator_i"]), int(row["generator_j"]))))
        if pair not in shared or not row.get("polygon_vertices"):
            continue
        vertices = np.asarray(row["polygon_vertices"], dtype=float)
        if len(vertices) >= 3:
            faces = np.asarray([(0, k, k + 1) for k in range(1, len(vertices)-1)], dtype=int)
            power_parts.append((vertices, faces))
    if power_parts:
        write_mesh(combine_mesh_parts([p[0] for p in power_parts], [p[1] for p in power_parts]),
                   "shared_Power_physical_surfaces.off", directory=root)
    # Standalone toggleable PyMOL loader using VorPy's existing OFF mesh layer
    # convention. It preserves the source molecule as a separate object.
    pymol_lines = [
        "from pymol import cmd", "from pymol.cgo import BEGIN, TRIANGLES, COLOR, VERTEX, END",
        "import os", f"cmd.load({str(Path(pdb_path).resolve())!r}, '2KAI_atoms')",
        "def _load_off(path,name,color):", "    if not os.path.exists(path): return",
        "    rows=open(path,encoding='utf-8').read().splitlines()",
        "    nv,nf,_=map(int,rows[1].split())", "    pts=[tuple(map(float,x.split())) for x in rows[2:2+nv]]",
        "    obj=[BEGIN,TRIANGLES,COLOR,*color]", "    for line in rows[2+nv:2+nv+nf]:",
        "        f=list(map(int,line.split())); ids=f[1:]",
        "        for k in range(1,len(ids)-1): obj += [VERTEX,*pts[ids[0]],VERTEX,*pts[ids[k]],VERTEX,*pts[ids[k+1]]]",
        "    obj.append(END); cmd.load_cgo(obj,name)"]
    colors = {"full_Power_dual": (0.48,0.50,0.58), "restricted_Power_dual": (0.2,0.45,1.0),
              "full_AW_dual": (0.58,0.58,0.60), "restricted_AW_dual": (0.0,0.85,0.85),
              "shared_pairs": (0.2,1.0,0.2), "Power_only_pairs": (1.0,0.35,0.0),
              "AW_only_pairs": (0.8,0.1,0.85), "shared_AW_physical_surfaces": (0.1,0.8,0.8),
              "shared_Power_physical_surfaces": (0.1,0.35,1.0)}
    for name, color in colors.items():
        mesh_path = root / f"{name}.off"
        if mesh_path.exists():
            pymol_lines.append(f"_load_off({str(mesh_path.resolve())!r},{name!r},{color!r})")
    pymol_lines += ["cmd.hide('everything','2KAI_atoms')", "cmd.show('spheres','2KAI_atoms and polymer.protein')",
                    "cmd.set('sphere_scale',0.22,'2KAI_atoms')", "cmd.zoom('2KAI_atoms')"]
    (root / "matched_probe_layers.py").write_text("\n".join(pymol_lines)+"\n", encoding="utf-8")
    (root / "matched_probe_layers.pml").write_text(
        f'run "{(root / "matched_probe_layers.py").resolve().as_posix()}"\n', encoding="utf-8")
    metadata = {"layers": {name: len(pairs) for name, pairs in layers.items()},
                "AW_physical_surface_mesh": "shared_AW_physical_surfaces.off uses saved solved AW surfaces",
                "Power_physical_surface_mesh": "shared_Power_physical_surfaces.off uses bounded polygons from weighted regular-complex power vertices",
                "Power_physical_surface_mesh_available": bool(power_parts)}
    (root / "visualization_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return str(root)


def run(network_path=NETWORK, radius_audit=RADIUS_AUDIT, output_dir=OUTPUT,
        *, pdb_path=PDB, aw_alpha=AW_ALPHA, calculate_curvature=True):
    from vorpy.src.io import load_network
    from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    from vorpy.src.analyze.aw_interface_curvature import _oriented_surface_curvatures
    from vorpy.src.geometry.aw_alpha.births import (
        surface_birth, _pair_contact_candidate, _clearances,
    )

    network_path, radius_audit, out, pdb_path = [Path(p).resolve() for p in
                                                  (network_path, radius_audit, output_dir, pdb_path)]
    out.mkdir(parents=True, exist_ok=True)
    network = load_network(network_path)
    aw_atoms, power_atoms = _atom_maps(radius_audit)
    ids = [int(i) for i in network.balls.index]
    if ids != list(range(len(ids))):
        raise AssertionError("This saved 2KAI solve is expected to use contiguous generator IDs.")
    if set(ids) != set(aw_atoms) or set(ids) != set(power_atoms):
        raise AssertionError("Generator/identity sets differ; stopping before geometry.")
    coords = np.asarray(network.balls.loc[ids, "loc"].tolist(), dtype=float)
    base = np.asarray(network.balls.loc[ids, "rad"], dtype=float)
    expanded = base + PROBE
    aw_effective = base + aw_alpha
    if not np.array_equal(base, np.asarray([aw_atoms[i]["aw_base"] for i in ids])):
        raise AssertionError("Saved AW network radii differ from the auditable AW base-radius column.")
    elements = _pdb_elements(pdb_path)
    if set(a["atom_id"] for a in aw_atoms.values()) != set(elements):
        raise AssertionError("H-tier atom identities do not exactly match protein-heavy atoms in the source PDB.")
    validate_matched_inputs([power_atoms[i]["atom_id"] for i in ids],
                            [aw_atoms[i]["atom_id"] for i in ids],
                            coords, coords, base, base, PROBE, aw_alpha)
    group_a = {i for i in ids if aw_atoms[i]["chain"] in {"A", "B"}}
    group_b = {i for i in ids if aw_atoms[i]["chain"] == "I"}
    if (len(ids), len(group_a), len(group_b)) != (2236, 1799, 437):
        raise AssertionError("Unexpected 2KAI H-tier atom/group counts.")
    input_rows = []
    for k, i in enumerate(ids):
        a = aw_atoms[i]
        input_rows.append({"generator_id": i, "stable_atom_id": a["stable"], "chain": a["chain"],
            "residue_number": a["resnum"], "insertion_code": a["icode"], "residue_name": a["resname"],
            "atom_name": a["atom"], "element": elements[a["atom_id"]],
            "x_A": coords[k, 0], "y_A": coords[k, 1], "z_A": coords[k, 2],
            "base_radius_A": base[k], "expanded_radius_A": expanded[k],
            "radius_source": "saved AW-loaded H-tier base radius"})
    _csv(out / "input_identity_audit.csv", input_rows)

    # Persist the expensive weighted regular result so AW/export fixes can be
    # rerun without solving Power geometry again. The cache is input-keyed.
    cache_key = hashlib.sha256(coords.tobytes() + expanded.tobytes() +
                               repr((sorted(group_a), sorted(group_b), 1e300)).encode()).hexdigest()
    cache_path = out / "power_geometry_cache.pkl"
    power = None
    if cache_path.exists():
        try:
            with cache_path.open("rb") as stream:
                cached_key, power = pickle.load(stream)
            if cached_key != cache_key:
                power = None
        except (OSError, EOFError, pickle.PickleError, ValueError):
            power = None
    if power is None:
        print("Building matched weighted regular complex and Power edge records…", flush=True)
        power = analyze_power_interface_curvature(coords, expanded, group_a, group_b,
                                                   condition_beta_m=1e300)
        with cache_path.open("wb") as stream:
            pickle.dump((cache_key, power), stream, protocol=pickle.HIGHEST_PROTOCOL)
    if not np.array_equal(power.points, coords) or not np.array_equal(power.expanded_radii, expanded):
        raise AssertionError("Power analyzer changed the controlled input arrays.")
    pfull = {tuple(sorted(p)) for p in power.full_simplices[1]}
    pres = {tuple(sorted(p)) for p in power.alpha_simplices[1]}
    if not pres <= pfull:
        raise AssertionError("Power restricted pair set is not a subset of its own full dual.")
    p_bicolor_full = {tuple(sorted(pair)) for pair in power.alpha_simplices[1]
                      if (pair[0] in group_a and pair[1] in group_b)
                      or (pair[1] in group_a and pair[0] in group_b)}
    parea = _power_areas(power, p_bicolor_full)
    power_rows = []
    for pair in sorted(pfull):
        i, j = pair; ai, aj = power_atoms[i], power_atoms[j]
        pa = parea.get(pair, {"status": ("not_calculated_rejected_bicolor" if
            ((pair[0] in group_a and pair[1] in group_b) or
             (pair[1] in group_a and pair[0] in group_b)) else "not_calculated_non_bicolor"),
            "area_A2": None})
        birth = float(power.filtration[pair])
        power_rows.append({"scheme": "power", "generator_i": i, "generator_j": j,
            "stable_atom_i": ai["stable"], "stable_atom_j": aj["stable"],
            "chain_i": ai["chain"], "chain_j": aj["chain"],
            "center_distance_A": float(np.linalg.norm(coords[i]-coords[j])),
            "base_radius_i_A": base[i], "base_radius_j_A": base[j],
            "expanded_radius_i_A": expanded[i], "expanded_radius_j_A": expanded[j],
            "full_dual": True, "restricted": pair in pres, "alpha_native": birth,
            "alpha_units": "A^2", "threshold": 0.0, "margin": -birth,
            "bicolor": (i in group_a and j in group_b) or (i in group_b and j in group_a),
            "physical_surface_id": f"power_pair:{i}:{j}",
            "surface_status": pa["status"], "surface_area_A2": pa["area_A2"],
            "bounded": pa["status"] == "finite_bounded", "complete": pa["status"] == "finite_bounded",
            "supported": pa["status"] == "finite_bounded", "restriction_status": "retained" if pair in pres else "rejected",
            "polygon_vertices": pa.get("polygon_vertices")})

    # Exact AW clipped-surface pair birth. No Power alpha values or guessed births.
    surfaces = _aw_surfaces(network)
    aw_rows = []
    edge_birth_cache, vertex_birth_cache = {}, {}
    all_aw_locations = np.asarray(network.balls["loc"].tolist(), dtype=float)
    all_aw_radii = np.asarray(network.balls["rad"], dtype=float)
    for n, (pair, features) in enumerate(sorted(surfaces.items()), 1):
        refs = [x[0] for x in features]
        simplex = SimpleNamespace(generator_ids=pair, primal_features=refs)
        i, j = pair
        is_bicolor = (i in group_a and j in group_b) or (i in group_b and j in group_a)
        candidate, candidate_diag = _pair_contact_candidate(network, pair, 1e-6)
        # The unconstrained pairwise AW minimum is a rigorous lower bound for
        # the minimum on the clipped surface. It can certify rejection only;
        # bicolor contacts and possible retained contacts always get the exact
        # clipped-surface minimization.
        bounded_complete = all(ref.bounded and ref.complete for ref in refs)
        cells_complete = all(ref.generator_cells_complete for ref in refs)
        birth_geometry_supported = bounded_complete and cells_complete
        certified_reject = (not is_bicolor and candidate is not None and
                            float(candidate["alpha"]) > aw_alpha + 1e-6)
        if certified_reject:
            birth, diag = None, {"birth_source": "pair_global_lower_bound",
                                 "reason": "Unconstrained pair minimum exceeds the AW threshold; clipped minimum cannot be lower."}
            lower_bound = float(candidate["alpha"])
        elif candidate is not None:
            candidate_clearances = _clearances(candidate["point"], all_aw_locations, all_aw_radii)
            candidate_alpha = float(candidate["alpha"])
            tie_count = int(np.sum(np.abs(candidate_clearances-candidate_alpha) <= 1e-6))
            if float(np.min(candidate_clearances)) >= candidate_alpha-1e-6 and tie_count == 2:
                birth, diag = candidate_alpha, {"birth_source":"interior_analytic",
                    "reason":"Analytic pair minimum is in the actual AW cell intersection.",
                    "residual_A":candidate["residual_A"]}
                lower_bound = candidate_alpha
            elif not birth_geometry_supported:
                birth, diag = None, {"birth_source":"unresolved",
                    "reason":"Pair minimum is clipped by other generators, but pair surface/cell support is incomplete."}
                lower_bound = candidate_alpha
            else:
                birth, diag = surface_birth(network, simplex, 1e-6,
                    edge_birth_cache=edge_birth_cache, vertex_birth_cache=vertex_birth_cache)
                lower_bound = candidate_alpha
        else:
            if not birth_geometry_supported:
                birth, diag = None, {"birth_source":"unresolved",
                    "reason":"No analytic pair contact and pair surface/cell support is incomplete."}
            else:
                birth, diag = surface_birth(network, simplex, 1e-6,
                    edge_birth_cache=edge_birth_cache, vertex_birth_cache=vertex_birth_cache)
            lower_bound = None
        filtered = None if birth is None else max(float(birth),
            -float(network.balls.loc[pair[0], "rad"]), -float(network.balls.loc[pair[1], "rad"]))
        selected = (not certified_reject and filtered is not None and filtered <= aw_alpha + 1e-6)
        i, j = pair; ai, aj = aw_atoms[i], aw_atoms[j]
        area = sum(float(row["sa"]) for _, row in features
                   if row.get("sa") is not None and math.isfinite(float(row["sa"])))
        bounded = any(ref.bounded for ref in refs); complete = any(ref.complete for ref in refs)
        supported = any(ref.bounded and ref.complete and ref.generator_cells_complete
                        for ref in refs)
        aw_rows.append({"scheme": "aw", "generator_i": i, "generator_j": j,
            "stable_atom_i": ai["stable"], "stable_atom_j": aj["stable"],
            "chain_i": ai["chain"], "chain_j": aj["chain"],
            "center_distance_A": float(np.linalg.norm(coords[i]-coords[j])),
            "base_radius_i_A": float(base[i]), "base_radius_j_A": float(base[j]),
            "expanded_radius_i_A": float(base[i]+aw_alpha), "expanded_radius_j_A": float(base[j]+aw_alpha),
            "full_dual": True, "restricted": bool(selected), "surface_birth_A": birth,
            "birth_lower_bound_A": lower_bound,
            "alpha_native": filtered, "alpha_units": "A", "threshold": aw_alpha,
            "margin": (aw_alpha-lower_bound if certified_reject else
                       None if filtered is None else aw_alpha-filtered),
            "bicolor": is_bicolor,
            "physical_surface_id": ";".join(str(ref.feature_id) for ref in refs),
            "surface_status": "supported" if supported else "incomplete_or_unsupported",
            "surface_area_A2": area, "bounded": bounded, "complete": complete, "supported": supported,
            "birth_source": diag.get("birth_source", "unresolved"), "birth_reason": diag.get("reason", ""),
            "restriction_status": "rejected" if certified_reject else
                "unresolved" if filtered is None else "retained" if selected else "rejected",
            "restriction_certificate": "global_pair_minimum_lower_bound" if certified_reject else "exact_clipped_surface_birth",
            "boundary_edge_count": sum(len(row.get("edges", ())) for _, row in features)})
        if n % 500 == 0:
            print(f"AW exact pair births {n}/{len(surfaces)}", flush=True)
    _csv(out / "full_dual_pairs.csv", power_rows + aw_rows)
    _csv(out / "restricted_pairs.csv", power_rows + aw_rows)

    p_bi = {tuple(sorted((r["generator_i"], r["generator_j"]))) for r in power_rows
            if r["restricted"] and r["bicolor"]}
    w_bi = {tuple(sorted((r["generator_i"], r["generator_j"]))) for r in aw_rows
            if r["restricted"] and r["bicolor"]}
    shared, p_only, a_only = pair_partition(p_bi, w_bi)
    pm = {tuple((r["generator_i"], r["generator_j"])): r for r in power_rows}
    am = {tuple((r["generator_i"], r["generator_j"])): r for r in aw_rows}
    comp_rows = []
    for pair in sorted(p_bi | w_bi):
        pr, ar = pm.get(pair), am.get(pair)
        comp_rows.append({"generator_i": pair[0], "generator_j": pair[1],
            "stable_atom_i": aw_atoms[pair[0]]["stable"], "stable_atom_j": aw_atoms[pair[1]]["stable"],
            "center_distance_A": float(np.linalg.norm(coords[pair[0]]-coords[pair[1]])),
            "in_Power": pair in p_bi, "in_AW": pair in w_bi,
            "classification": "shared" if pair in shared else "Power_only" if pair in p_only else "AW_only",
            "Power_birth_A2": pr["alpha_native"] if pr else None,
            "Power_margin_A2": pr["margin"] if pr else None,
            "Power_surface_area_A2": pr["surface_area_A2"] if pr else None,
            "Power_surface_status": pr["surface_status"] if pr else "not_in_full_dual",
            "AW_birth_A": ar["alpha_native"] if ar else None,
            "AW_surface_birth_A": ar["surface_birth_A"] if ar else None,
            "AW_threshold_A": aw_alpha, "AW_margin_A": ar["margin"] if ar else None,
            "AW_surface_area_A2": ar["surface_area_A2"] if ar else None,
            "AW_surface_status": ar["surface_status"] if ar else "not_in_full_dual",
            "AW_birth_source": ar["birth_source"] if ar else None,
            "AW_birth_reason": ar["birth_reason"] if ar else None})
    _csv(out / "bicolor_pair_comparison.csv", comp_rows)
    _csv(out / "disagreement_pairs.csv", [r for r in comp_rows if r["classification"] != "shared"])
    shared_rows = []
    for pair in sorted(shared):
        pr, ar = pm[pair], am[pair]; pa, aa = pr["surface_area_A2"], ar["surface_area_A2"]
        shared_rows.append({"generator_i": pair[0], "generator_j": pair[1],
            "stable_atom_i": aw_atoms[pair[0]]["stable"], "stable_atom_j": aw_atoms[pair[1]]["stable"],
            "Power_area_A2": pa, "AW_area_A2": aa,
            "delta_area_A2": None if pa is None or aa is None else aa-pa,
            "relative_delta_area": None if pa in (None, 0) or aa is None else (aa-pa)/pa,
            "Power_surface_id": pr["physical_surface_id"], "AW_surface_id": ar["physical_surface_id"],
            "Power_status": pr["surface_status"], "AW_status": ar["surface_status"]})
    _csv(out / "shared_pair_surfaces.csv", shared_rows)

    if calculate_curvature:
        aw_curve = analyze_aw_interface_curvature(network, group_a, group_b, selection_mode="supported")
        power_shared = _power_metrics(power, shared, parea)
        power_total = _power_metrics(power, p_bi, parea)
        aw_shared_area = sum(float(am[p]["surface_area_A2"] or 0) for p in shared)
        aw_total_area = sum(float(am[p]["surface_area_A2"] or 0) for p in w_bi)
        aw_shared = _aw_metrics(aw_curve, shared, aw_shared_area)
        aw_total = _aw_metrics(aw_curve, w_bi, aw_total_area)
        curvature_rows = [{"comparison": "shared_pair_geometry", **power_shared},
                          {"comparison": "shared_pair_geometry", **aw_shared},
                          {"comparison": "total_restricted_interface", **power_total},
                          {"comparison": "total_restricted_interface", **aw_total}]
        _csv(out / "shared_pair_curvature.csv", curvature_rows[:2])
        _csv(out / "total_interface_metrics.csv", curvature_rows)
    else:
        curvature_rows = []

    distance_rows = []
    for scheme, rows, selected in (("power", power_rows, {tuple((r["generator_i"],r["generator_j"])) for r in power_rows if r["restricted"]}),
                                   ("aw", aw_rows, {tuple((r["generator_i"],r["generator_j"])) for r in aw_rows if r["restricted"]})):
        bicolor_set = p_bi if scheme == "power" else w_bi
        for stage, pairs in (("full", set(tuple((r["generator_i"],r["generator_j"])) for r in rows)),
                             ("molecular_restricted", selected), ("interface", bicolor_set)):
            for pair in sorted(pairs):
                r = (pm if scheme == "power" else am)[pair]
                distance_rows.append({"scheme": scheme, "stage": stage, "pair_id": f"{pair[0]}:{pair[1]}",
                    "bicolor": r["bicolor"], "center_distance_A": r["center_distance_A"], "selected": True})
    _csv(out / "connection_distances.csv", distance_rows)

    distributions = {}
    for scheme, rows in (("power", power_rows), ("aw", aw_rows)):
        distributions[scheme] = {}
        for stage in ("full", "rejected", "restricted", "unresolved"):
            part = [r for r in rows if stage == "full" or
                    (stage == "restricted" and r["restriction_status"] == "retained") or
                    (stage == "rejected" and r["restriction_status"] == "rejected") or
                    (stage == "unresolved" and r["restriction_status"] == "unresolved")]
            distributions[scheme][stage] = _dist_stats(part)
    total_pair_power_area = sum(float(pm[p]["surface_area_A2"] or 0) for p in p_bi
                                if pm[p]["surface_area_A2"] is not None)
    total_pair_aw_area = sum(float(am[p]["surface_area_A2"] or 0) for p in w_bi)
    shared_p_area = sum(float(pm[p]["surface_area_A2"] or 0) for p in shared
                        if pm[p]["surface_area_A2"] is not None)
    shared_a_area = sum(float(am[p]["surface_area_A2"] or 0) for p in shared)
    summary = {"system": "2KAI", "atom_count": len(ids), "group_A_atoms": len(group_a), "group_B_atoms": len(group_b),
        "input_control": {"stable_atom_ids_equal": True, "coordinates_exactly_equal": True,
            "base_radii_exactly_equal": True, "expanded_radii_exactly_equal": bool(np.array_equal(expanded, aw_effective)),
            "base_radius_source": "saved H-tier AW-loaded radius", "Power_expansion": "R=r+1.4; weight=R^2; alpha=0 A^2",
            "AW_expansion": "base r; actual restricted surface birth <=1.4 A; no pre-expansion",
            "frozen_Cazals_reference_separate": {"full_bicolor_regular":690,"alpha0":362,"M5_final":362,
                "radius_model":"literature-reference Cazals/Chothia + 1.4 A"}},
        "full_dual_counts": {"power":len(pfull),"aw":len(surfaces)},
        "restricted_counts": {"power":{"retained":len(pres),"rejected":len(pfull-pres),"unresolved":0},
            "aw":{"retained":sum(r["restriction_status"]=="retained" for r in aw_rows),
                 "rejected":sum(r["restriction_status"]=="rejected" for r in aw_rows),
                 "unresolved":sum(r["restriction_status"]=="unresolved" for r in aw_rows)}},
        "distance_distributions_A":distributions,
        "bicolor":{"Power":len(p_bi),"AW":len(w_bi),"shared":len(shared),"Power_only":len(p_only),
            "AW_only":len(a_only),"union":len(p_bi|w_bi),"Jaccard":len(shared)/len(p_bi|w_bi) if p_bi|w_bi else None,
            "fraction_Power_shared":len(shared)/len(p_bi) if p_bi else None,"fraction_AW_shared":len(shared)/len(w_bi) if w_bi else None},
        "shared_pair_area_A2":{"Power":shared_p_area,"AW":shared_a_area,"delta_AW_minus_Power":shared_a_area-shared_p_area},
        "total_interface_area_A2":{"Power":total_pair_power_area,"AW":total_pair_aw_area,"delta_AW_minus_Power":total_pair_aw_area-total_pair_power_area},
        "curvature":curvature_rows,
        "longest_retained_connections":{"power":sorted([r for r in power_rows if r["restricted"]],key=lambda r:r["center_distance_A"],reverse=True)[:10],
            "aw":sorted([r for r in aw_rows if r["restriction_status"]=="retained"],key=lambda r:r["center_distance_A"],reverse=True)[:10]},
        "boundary_policy":"No boundary curvature term; integrate only seams shared by two selected physical surfaces.",
        "visualization":_export_visualization(network,power_rows,aw_rows,shared,p_only,a_only,out,pdb_path)}
    (out/"comparison_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    (out/"comparison_report.txt").write_text(_report(summary),encoding="utf-8")
    return summary


def _report(s):
    b=s["bicolor"]; p=s["shared_pair_area_A2"]; t=s["total_interface_area_A2"]
    return ("2KAI matched-probe Power-vs-AW molecular restriction comparison\n"
        "Shared base radii are the saved H-tier AW-loaded radii; Power uses R=r+1.4 and AW uses base r with alpha=1.4 A.\n"
        "Power weighted-alpha units are A^2; AW additive-offset units are A. The frozen Cazals/Chothia result is context only.\n"
        f"Full pairs: Power {s['full_dual_counts']['power']}, AW {s['full_dual_counts']['aw']}.\n"
        f"Restricted pairs: Power {s['restricted_counts']['power']}; AW {s['restricted_counts']['aw']}.\n"
        f"Bicolor: Power {b['Power']}, AW {b['AW']}, shared {b['shared']}, Power-only {b['Power_only']}, AW-only {b['AW_only']}, Jaccard {b['Jaccard']}.\n"
        f"Shared-pair areas A^2: Power {p['Power']:.8f}, AW {p['AW']:.8f}, delta {p['delta_AW_minus_Power']:.8f}.\n"
        f"Total restricted interface areas A^2: Power {t['Power']:.8f}, AW {t['AW']:.8f}, delta {t['delta_AW_minus_Power']:.8f}.\n"
        "Boundary curvature is not assigned; unresolved geometry remains reported rather than replaced by a heuristic.\n"
        f"Visualization layers: {s['visualization']}\n")


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network",default=str(NETWORK)); parser.add_argument("--radius-audit",default=str(RADIUS_AUDIT))
    parser.add_argument("--output-dir",default=str(OUTPUT)); parser.add_argument("--pdb",default=str(PDB))
    parser.add_argument("--aw-alpha",type=float,default=AW_ALPHA)
    parser.add_argument("--skip-curvature",action="store_true")
    a=parser.parse_args(argv)
    result=run(a.network,a.radius_audit,a.output_dir,pdb_path=a.pdb,aw_alpha=a.aw_alpha,calculate_curvature=not a.skip_curvature)
    print(json.dumps({"bicolor":result["bicolor"],"full":result["full_dual_counts"],"restricted":result["restricted_counts"],"output":str(Path(a.output_dir).resolve())},indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
