"""Independent stage-by-stage 2KAI Cazals reproduction audit.

This diagnostic does not edit or replace any production implementation. It
parses PDB records itself, assigns the already-published VorPy Cazals radius
table values, constructs a GUDHI weighted alpha complex directly, and
independently reconstructs power vertices, interface-edge lengths, and
unsigned dihedrals. At each stage it compares to the frozen analyzer and
stops the comparison at the first disagreement.

Run from the repository root:
    python -m vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_2kai_independent_audit

Exports under output/cazals_2kai_independent_audit are designed so later
stages can be reproduced from their CSV inputs without importing VorPy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from vorpy.src.chemistry.Cazals import (
    OTHER_RADIUS,
    PROBE_RADIUS,
    element_radii,
    special_radii,
)

GROUP_A_CHAINS = frozenset({"A", "B"})
GROUP_B_CHAINS = frozenset({"I"})
ALPHA_TOLERANCE = 1e-12
M5_THRESHOLD = 5.0


def _infer_element(line):
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    name = line[12:16].strip().upper().lstrip("0123456789")
    return name[:1] if name and name[0] in "CNO SHP".replace(" ", "") else name[:2]


def _radius(resname, atomname, element):
    """Local radius lookup using the frozen tabulated data, not its helper."""
    residue = resname.upper()
    atom = atomname.upper()
    if residue in special_radii and atom in special_radii[residue]:
        value = float(special_radii[residue][atom])
        if element == "C":
            kind = "aliphatic_C" if np.isclose(value, 1.87) else "trigonal_C" if np.isclose(value, 1.76) else "explicit_residue_atom"
        elif element == "N":
            kind = "neutral_N" if np.isclose(value, 1.65) else "charged_N" if np.isclose(value, 1.50) else "explicit_residue_atom"
        elif element == "O":
            kind = "oxygen"
        elif element == "S":
            kind = "sulfur"
        else:
            kind = "explicit_residue_atom"
        return value, kind, True
    value = float(element_radii.get(element, OTHER_RADIUS))
    kind = {"O": "oxygen_element_fallback", "S": "sulfur_element_fallback",
            "C": "ambiguous_carbon_fallback", "N": "ambiguous_nitrogen_fallback"}.get(element, "other_fallback")
    return value, kind, False


def parse_pdb_independently(pdb_path):
    """Parse all coordinate records, then select the heavy protein model."""
    path = Path(pdb_path)
    raw = []
    inventory = Counter()
    by_chain_record = Counter()
    for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        record = line[:6].strip().upper()
        if record not in {"ATOM", "HETATM"}:
            continue
        inventory[record] += 1
        chain = line[21:22].strip()
        resname = line[17:20].strip().upper()
        by_chain_record[(chain, record, resname)] += 1
        if record != "ATOM" or chain not in GROUP_A_CHAINS | GROUP_B_CHAINS:
            continue
        element = _infer_element(line)
        atom_name = line[12:16].strip().upper()
        if element == "H" or atom_name.lstrip("0123456789").startswith("H"):
            continue
        try:
            raw.append({
                "line_number": line_number,
                "pdb_serial": int(line[6:11]),
                "record_type": record,
                "atom_name": atom_name,
                "altloc": line[16:17].strip().upper(),
                "residue_name": resname,
                "chain": chain,
                "residue_number": int(line[22:26]),
                "insertion_code": line[26:27].strip(),
                "element": element,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
        except (ValueError, IndexError) as error:
            raise ValueError(f"Invalid atom record at {path}:{line_number}") from error

    # Independent implementation of the documented altloc preference.
    chosen = {}
    priority = {"": 0, "A": 1}
    for atom in raw:
        key = (atom["chain"], atom["residue_number"], atom["insertion_code"],
               atom["residue_name"], atom["atom_name"])
        old = chosen.get(key)
        if old is None or priority.get(atom["altloc"], 2) < priority.get(old["altloc"], 2):
            chosen[key] = atom
    atoms = sorted(chosen.values(), key=lambda row: row["pdb_serial"])
    altloc_removed = len(raw) - len(atoms)
    for index, atom in enumerate(atoms):
        base, radius_class, explicit = _radius(atom["residue_name"], atom["atom_name"], atom["element"])
        group = "A" if atom["chain"] in GROUP_A_CHAINS else "B"
        atom.update({
            "atom_id": index,
            "group": group,
            "radius_class": radius_class,
            "explicit_assignment": explicit,
            "base_radius_A": base,
            "probe_radius_A": float(PROBE_RADIUS),
            "expanded_radius_A": base + float(PROBE_RADIUS),
            "power_weight_A2": (base + float(PROBE_RADIUS)) ** 2,
        })
    census = {
        "filename": path.name,
        "total_coordinate_records": sum(inventory.values()),
        "ATOM_records": inventory["ATOM"],
        "HETATM_records": inventory["HETATM"],
        "protein_heavy_atoms_all_chains": sum(1 for row in raw),
        "selected_heavy_protein_atoms_after_altloc": len(atoms),
        "selected_group_A_atoms": sum(row["group"] == "A" for row in atoms),
        "selected_group_B_atoms": sum(row["group"] == "B" for row in atoms),
        "altloc_duplicates_removed": altloc_removed,
        "chain_record_residue_counts": {"|".join(key): value for key, value in sorted(by_chain_record.items())},
    }
    return atoms, census


def _write_csv(path, rows, fields=None):
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def _csv_value(value):
    if isinstance(value, (tuple, list, set, frozenset, dict)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    if isinstance(value, np.ndarray):
        return json.dumps(value.tolist(), separators=(",", ":"))
    return value


def _atom_stage_compare(atoms, production_prepared, tolerance=0.0):
    if len(atoms) != len(production_prepared.atoms):
        return {"match": False, "reason": f"atom count {len(atoms)} != {len(production_prepared.atoms)}"}
    for index, (left, right) in enumerate(zip(atoms, production_prepared.atoms)):
        identity = (left["pdb_serial"], left["chain"], left["residue_name"], left["residue_number"], left["insertion_code"], left["atom_name"])
        production_identity = (right.pdb_serial, right.chain, right.residue_name, right.residue_number, right.insertion_code, right.atom_name)
        if identity != production_identity:
            return {"match": False, "reason": f"atom identity mismatch at selected index {index}: {identity} != {production_identity}"}
        for key, value in (("x", right.x), ("y", right.y), ("z", right.z),
                           ("base_radius_A", right.base_radius),
                           ("expanded_radius_A", right.expanded_radius), ("power_weight_A2", right.weight)):
            if not math.isclose(float(left[key]), float(value), rel_tol=0.0, abs_tol=tolerance):
                return {"match": False, "reason": f"{key} mismatch for atom {identity}: {left[key]} != {value}"}
        if left["group"] != right.group:
            return {"match": False, "reason": f"group mismatch for atom {identity}"}
    return {"match": True, "reason": "all selected heavy-atom identities, coordinates, groups, radii, and weights agree"}


def _direct_gudhi(points, weights):
    """Construct the weighted alpha/regular complex without VorPy helpers."""
    import gudhi

    ac = gudhi.AlphaComplex(points=np.asarray(points, dtype=float).tolist(),
                            weights=np.asarray(weights, dtype=float).tolist(), precision="safe")
    tree = ac.create_simplex_tree()
    simplices = {dimension: {} for dimension in range(4)}
    for simplex, filtration in tree.get_filtration():
        simplex = tuple(sorted(map(int, simplex)))
        dimension = len(simplex) - 1
        if dimension > 3:
            raise RuntimeError(f"Unexpected regular simplex dimension {dimension}: {simplex}")
        simplices[dimension][simplex] = float(filtration)
    return tree, simplices


def _power_vertex_independent(tetrahedron, points, weights):
    """Solve equal-power equations using a separate least-squares formulation."""
    ids = tuple(sorted(map(int, tetrahedron)))
    p = np.asarray(points, dtype=float)[list(ids)]
    w = np.asarray(weights, dtype=float)[list(ids)]
    lhs = 2.0 * (p[1:] - p[0])
    rhs = np.sum(p[1:] * p[1:], axis=1) - w[1:] - (float(p[0] @ p[0]) - w[0])
    position, _, rank, singular = np.linalg.lstsq(lhs, rhs, rcond=None)
    if rank != 3:
        raise np.linalg.LinAlgError(f"rank deficient tetrahedron {ids}: rank={rank}")
    residuals = np.sum((position[None, :] - p) ** 2, axis=1) - w
    return position, float(np.ptp(residuals)), float(np.linalg.cond(lhs)), tuple(map(float, singular))


def _incident_tetrahedra(simplices):
    tri_to_tets = defaultdict(list)
    for tet in simplices[3]:
        for omit in range(4):
            tri_to_tets[tuple(tet[j] for j in range(4) if j != omit)].append(tet)
    edge_to_triangles = defaultdict(list)
    for tri in simplices[2]:
        for omit in range(3):
            edge_to_triangles[tuple(tri[j] for j in range(3) if j != omit)].append(tri)
    return tri_to_tets, edge_to_triangles


def _is_bicolor(pair, group_a, group_b):
    i, j = pair
    return (i in group_a and j in group_b) or (i in group_b and j in group_a)


def _pair_pattern(simplex, group_a, group_b):
    return "".join("A" if i in group_a else "B" if i in group_b else "?" for i in simplex)


def _facet_m5_independent(pair, filtration, simplices, tetra_q, tri_to_tets,
                          smaller_radius, alpha_tolerance=ALPHA_TOLERANCE):
    alpha_value = filtration.get(pair, math.inf)
    alpha_selected = alpha_value <= alpha_tolerance
    if not alpha_selected:
        return {"pair": pair, "alpha_filtration": alpha_value,
                "alpha0_selected": False, "boundary_triangles": [],
                "candidate_q": [], "m_A": None, "m_over_r": None,
                "M5_accepted": None, "M5_status": "not_evaluated_not_alpha0",
                "reason": "facet is absent at alpha=0"}
    triangles = [tri for tri in simplices[2] if set(pair).issubset(tri)]
    candidate_q = []
    reason = ""
    for tri in triangles:
        incident = [tet for tet in tri_to_tets.get(tri, ())]
        if len(incident) != 2:
            reason = f"facet boundary triangle {tri} has {len(incident)} incident tetrahedra"
            continue
        if any(tet not in tetra_q for tet in incident):
            reason = f"unresolved tetrahedron in facet star at triangle {tri}"
            continue
        candidate_q.extend((tet, tetra_q[tet]) for tet in incident)
    if not triangles:
        reason = "no regular triangles in facet star"
    if reason:
        m = ratio = None
        accepted = False if alpha_selected else None
        status = "unresolved"
    elif not candidate_q:
        m = ratio = None
        accepted = False if alpha_selected else None
        status = "unresolved"
        reason = "facet star has no orthogonal-ball candidates"
    else:
        max_tet, max_q = max(candidate_q, key=lambda item: item[1])
        m = math.sqrt(max(0.0, max_q))
        ratio = m / smaller_radius
        accepted = bool(ratio <= M5_THRESHOLD)
        status = "accepted" if accepted else "rejected_m_over_r"
        reason = f"max-q tetrahedron {max_tet}"
    return {"pair": pair, "alpha_filtration": alpha_value,
            "alpha0_selected": alpha_selected, "boundary_triangles": triangles,
            "candidate_q": candidate_q, "m_A": m, "m_over_r": ratio,
            "M5_accepted": accepted, "M5_status": status, "reason": reason}


def _interface_facets_for_triangle(triangle, bicolor_pairs):
    tri = tuple(triangle)
    sides = [tuple(sorted((tri[i], tri[j]))) for i, j in ((0, 1), (0, 2), (1, 2))]
    return sorted(pair for pair in sides if pair in bicolor_pairs)


def _unsigned_angle_from_normals(pair_1, pair_2, points, group_a, group_b):
    def normal(pair):
        i, j = pair
        if i in group_a and j in group_b:
            p, q = points[i], points[j]
        elif j in group_a and i in group_b:
            p, q = points[j], points[i]
        else:
            raise ValueError(f"Not a bicolor pair: {pair}")
        return (q - p) / np.linalg.norm(q - p)
    n1, n2 = normal(pair_1), normal(pair_2)
    return float(np.arccos(np.clip(float(n1 @ n2), -1.0, 1.0))), n1, n2


def _unsigned_triangle_angle(triangle, points, group_a, group_b):
    a = [index for index in triangle if index in group_a]
    b = [index for index in triangle if index in group_b]
    if len(a) == 1:
        singleton, others = a[0], b
    elif len(b) == 1:
        singleton, others = b[0], a
    else:
        raise ValueError(f"Not a 1:2 bicolor triangle: {triangle}")
    u = points[others[0]] - points[singleton]
    v = points[others[1]] - points[singleton]
    return float(np.arccos(np.clip(float(u @ v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1.0, 1.0))), singleton


def run_independent_2kai_audit(pdb_path=None, output_dir="output/cazals_2kai_independent_audit"):
    from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import (
        prepare_cazals_pdb,
    )
    from vorpy.src.analyze.power_interface_curvature import (
        analyze_power_interface_curvature,
    )

    pdb_path = Path(pdb_path or Path(__file__).resolve().parents[4] / "data" / "2KAI.pdb")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stages = []

    atoms, census = parse_pdb_independently(pdb_path)
    points = np.asarray([[row[axis] for axis in ("x", "y", "z")] for row in atoms], dtype=float)
    radii = np.asarray([row["expanded_radius_A"] for row in atoms], dtype=float)
    weights = np.asarray([row["power_weight_A2"] for row in atoms], dtype=float)
    group_a = {i for i, row in enumerate(atoms) if row["group"] == "A"}
    group_b = {i for i, row in enumerate(atoms) if row["group"] == "B"}
    atom_rows = [{**row, "selected_atom_id": i} for i, row in enumerate(atoms)]
    _write_csv(output / "stage_01_atoms_and_radii.csv", atom_rows)
    production_prepared = prepare_cazals_pdb(
        pdb_path, GROUP_A_CHAINS, GROUP_B_CHAINS,
        probe_radius=PROBE_RADIUS, include_hydrogens=False,
        include_hetatm=False, strict_heavy_fallbacks=False, verbose=False,
    )
    prepared_comparison = _atom_stage_compare(atoms, production_prepared)
    census["group_A_chain_definition"] = sorted(GROUP_A_CHAINS)
    census["group_B_chain_definition"] = sorted(GROUP_B_CHAINS)
    census["radius_class_counts"] = dict(Counter(row["radius_class"] for row in atoms))
    census["explicit_radius_assignments"] = sum(row["explicit_assignment"] for row in atoms)
    census["fallback_radius_assignments"] = sum(not row["explicit_assignment"] for row in atoms)
    stage = {"stage": "01_pdb_atoms_radii", "independent_count": len(atoms),
             "production_count": len(production_prepared.atoms), **prepared_comparison}
    stages.append(stage)
    _write_csv(output / "stage_01_census.csv", [census])
    if not prepared_comparison["match"]:
        return _finish(output, stages, "01_pdb_atoms_radii", prepared_comparison["reason"], census)

    # Production result is computed only after independent PDB/radius equality.
    production = analyze_power_interface_curvature(
        production_prepared.points, production_prepared.expanded_radii,
        production_prepared.group_a, production_prepared.group_b,
        alpha=0.0, condition_beta_m=M5_THRESHOLD,
    )
    _tree, independent_simplices = _direct_gudhi(points, weights)
    filtration = {simplex: value for dimension in independent_simplices
                  for simplex, value in independent_simplices[dimension].items()}
    alpha_simplices = {dimension: {simplex for simplex, value in independent_simplices[dimension].items()
                                   if value <= ALPHA_TOLERANCE}
                       for dimension in range(4)}
    simplex_rows = []
    for dimension in range(4):
        for simplex, value in sorted(independent_simplices[dimension].items()):
            simplex_rows.append({"dimension": dimension, "generator_ids": simplex,
                                 "gudhi_filtration_native": value,
                                 "alpha0_threshold": 0.0,
                                 "alpha0_tolerance": ALPHA_TOLERANCE,
                                 "alpha0_member": simplex in alpha_simplices[dimension]})
    _write_csv(output / "stage_02_gudhi_regular_and_alpha_simplices.csv", simplex_rows)
    full_sets_match = all(independent_simplices[d].keys() == production.full_simplices[d] for d in range(4))
    alpha_sets_match = all(alpha_simplices[d] == production.alpha_simplices[d] for d in range(4))
    filtration_max_difference = max((abs(filtration[key] - production.filtration[key])
                                     for key in filtration.keys() & production.filtration.keys()), default=0.0)
    gudhi_match = full_sets_match and alpha_sets_match and filtration.keys() == production.filtration.keys() and filtration_max_difference <= 1e-12
    stages.append({"stage": "02_gudhi_regular_alpha", "independent_counts": {d: len(independent_simplices[d]) for d in range(4)},
                   "production_counts": production.full_counts,
                   "independent_alpha0_counts": {d: len(alpha_simplices[d]) for d in range(4)},
                   "production_alpha0_counts": production.alpha_counts,
                   "full_membership_match": full_sets_match, "alpha0_membership_match": alpha_sets_match,
                   "filtration_key_match": filtration.keys() == production.filtration.keys(),
                   "max_filtration_difference": filtration_max_difference, "match": gudhi_match})
    if not gudhi_match:
        return _finish(output, stages, "02_gudhi_regular_alpha", "GUDHI simplex membership or filtration differs from production", census)

    tri_to_tets, _edge_to_triangles = _incident_tetrahedra(independent_simplices)
    tetra_vertices = {}
    tetra_rows = []
    for tetrahedron in sorted(independent_simplices[3]):
        try:
            xyz, residual, condition, singular = _power_vertex_independent(tetrahedron, points, weights)
            qvalues = [float(np.sum((xyz - points[i]) ** 2) - weights[i]) for i in tetrahedron]
            q = float(np.mean(qvalues))
            tetra_vertices[tetrahedron] = (xyz, q)
            tetra_rows.append({"tetrahedron": tetrahedron, "x": xyz[0], "y": xyz[1], "z": xyz[2],
                               "q_power_A2": q, "equal_power_spread_A2": max(qvalues) - min(qvalues),
                               "linear_residual": residual, "condition_number": condition,
                               "singular_values": singular, "status": "resolved"})
        except np.linalg.LinAlgError as error:
            tetra_rows.append({"tetrahedron": tetrahedron, "status": "unresolved", "reason": str(error)})
    _write_csv(output / "stage_03_independent_power_vertices.csv", tetra_rows)
    tetra_q = {tet: item[1] for tet, item in tetra_vertices.items()}
    production_power_vertices = production.power_vertices
    independent_tetrahedra = set(tetra_vertices)
    production_tetrahedra = set(production_power_vertices)
    power_vertex_position_difference = max((
        float(np.linalg.norm(tetra_vertices[tet][0] - np.asarray(production_power_vertices[tet].position)))
        for tet in independent_tetrahedra & production_tetrahedra
    ), default=0.0)
    power_vertex_population_match = independent_tetrahedra == production_tetrahedra
    power_vertex_match = power_vertex_population_match and power_vertex_position_difference <= 1e-8
    stages.append({
        "stage": "03_power_vertex_solutions",
        "independent_regular_tetrahedra": len(independent_simplices[3]),
        "production_regular_tetrahedra": production.full_counts[3],
        "independent_resolved_power_vertices": len(independent_tetrahedra),
        "production_resolved_power_vertices": len(production_tetrahedra),
        "resolved_tetrahedron_membership_match": power_vertex_population_match,
        "maximum_endpoint_coordinate_difference_A": power_vertex_position_difference,
        "match": power_vertex_match,
    })
    if not power_vertex_match:
        return _finish(output, stages, "03_power_vertex_solutions",
                       "independently solved regular-tetrahedron power vertices disagree", census)

    bicolor_pairs = {pair for pair in independent_simplices[1] if _is_bicolor(pair, group_a, group_b)}
    facet_rows = []
    facet_tet_rows = []
    prod_facets = {row.generator_ids: row for row in production.facets}
    for pair in sorted(bicolor_pairs):
        independent = _facet_m5_independent(
            pair, filtration, independent_simplices, tetra_q, tri_to_tets,
            min(radii[list(pair)]),
        )
        independent["r_smaller_A"] = min(radii[list(pair)])
        production_facet = prod_facets[pair]
        independent.update({
            "production_alpha0_selected": production_facet.alpha_zero_selected,
            "production_M5_status": production_facet.condition_beta_status,
            "production_m_A": production_facet.largest_orthogonal_ball_radius,
            "production_m_over_r": production_facet.m_over_r,
            "production_M5_accepted": production_facet.condition_beta_accepted,
            "production_final_selected": production_facet.final_selected,
            "independent_final_selected": bool(independent["alpha0_selected"] and independent["M5_accepted"] is True),
        })
        facet_rows.append(independent)
        for tet, q in independent["candidate_q"]:
            facet_tet_rows.append({"pair": pair, "tetrahedron": tet, "q_power_A2": q,
                                   "m_candidate_A": math.sqrt(max(0.0, q)),
                                   "candidate_is_maximum": independent["reason"].endswith(str(tet))})
    _write_csv(output / "stage_04_bicolor_m5_facets.csv", facet_rows)
    _write_csv(output / "stage_04_facet_orthogonal_ball_candidates.csv", facet_tet_rows)
    independent_final = {tuple(row["pair"]) for row in facet_rows if row["independent_final_selected"]}
    production_final = {tuple(row.generator_ids) for row in production.final_facets}
    facet_stage_match = (independent_final == production_final and all(
        bool(row["alpha0_selected"]) == bool(row["production_alpha0_selected"])
        and (not row["alpha0_selected"] or row["M5_accepted"] == row["production_M5_accepted"])
        and (not row["alpha0_selected"] or row["m_A"] is None or
             (row["production_m_A"] is not None and
              math.isclose(row["m_A"], row["production_m_A"], rel_tol=0.0, abs_tol=1e-8)))
        for row in facet_rows
    ))
    stages.append({"stage": "04_alpha0_bicolor_m5_selection", "regular_bicolor_pairs": len(bicolor_pairs),
                   "independent_alpha0_bicolor_pairs": sum(row["alpha0_selected"] for row in facet_rows),
                   "production_alpha0_bicolor_pairs": len(production.alpha_zero_facets),
                   "independent_final_pairs": len(independent_final), "production_final_pairs": len(production_final),
                   "independent_M5_rejected": sum(row["M5_accepted"] is False for row in facet_rows if row["alpha0_selected"]),
                   "production_M5_rejected": sum(row.alpha_zero_selected and row.condition_beta_accepted is False for row in production.facets),
                   "match": facet_stage_match})
    if not facet_stage_match:
        return _finish(output, stages, "04_alpha0_bicolor_m5_selection", "alpha-zero or M=5 facet membership differs", census)

    points_by_id = points
    independent_edges = []
    for triangle in sorted(independent_simplices[2]):
        interface_pairs = _interface_facets_for_triangle(triangle, independent_final)
        if len(interface_pairs) != 2:
            continue
        incident_tets = tuple(sorted(tri_to_tets.get(triangle, ())))
        if len(incident_tets) != 2 or any(tet not in tetra_vertices for tet in incident_tets):
            continue
        p1, q1 = tetra_vertices[incident_tets[0]]
        p2, q2 = tetra_vertices[incident_tets[1]]
        length = float(np.linalg.norm(p2 - p1))
        beta_planes, normal1, normal2 = _unsigned_angle_from_normals(
            interface_pairs[0], interface_pairs[1], points_by_id, group_a, group_b
        )
        beta_triangle, singleton = _unsigned_triangle_angle(triangle, points_by_id, group_a, group_b)
        pattern = _pair_pattern(triangle, group_a, group_b)
        sign = -1 if pattern == "AAB" else 1 if pattern == "ABB" else 0
        signed_beta = sign * beta_triangle
        independent_edges.append({
            "triangle_generators": triangle,
            "group_pattern": pattern,
            "singleton_generator": singleton,
            "incident_bicolor_facets": interface_pairs,
            "incident_tetrahedra": incident_tets,
            "endpoint_1_x": p1[0], "endpoint_1_y": p1[1], "endpoint_1_z": p1[2],
            "endpoint_2_x": p2[0], "endpoint_2_y": p2[1], "endpoint_2_z": p2[2],
            "power_q_endpoint_1_A2": q1, "power_q_endpoint_2_A2": q2,
            "edge_length_independent_A": length,
            "normal_1_xyz": tuple(map(float, normal1)), "normal_2_xyz": tuple(map(float, normal2)),
            "beta_unsigned_power_plane_normals_rad": beta_planes,
            "beta_unsigned_power_plane_normals_deg": math.degrees(beta_planes),
            "beta_unsigned_dual_triangle_rad": beta_triangle,
            "beta_unsigned_dual_triangle_deg": math.degrees(beta_triangle),
            "unsigned_angle_difference_rad": abs(beta_planes - beta_triangle),
            "Cazals_sign": sign,
            "beta_signed_Cazals_rad": signed_beta,
            "signed_beta_l_independent_A_rad": signed_beta * length,
            "unsigned_beta_l_independent_A_rad": beta_triangle * length,
        })

    _write_csv(output / "stage_05_independent_interface_edges.csv", independent_edges)
    prod_included = {tuple(row.regular_triangle): row for row in production.included_edges}
    independent_edge_ids = {tuple(row["triangle_generators"]) for row in independent_edges}
    production_edge_ids = set(prod_included)
    if independent_edge_ids != production_edge_ids:
        stages.append({"stage": "05_interface_edge_population", "independent_edge_count": len(independent_edge_ids),
                       "production_edge_count": len(production_edge_ids),
                       "independent_only": sorted(independent_edge_ids - production_edge_ids),
                       "production_only": sorted(production_edge_ids - independent_edge_ids), "match": False})
        return _finish(output, stages, "05_interface_edge_population", "finite selected edge populations differ", census)
    edge_rows_by_id = {tuple(row["triangle_generators"]): row for row in independent_edges}
    length_differences = []
    plane_angle_differences = []
    triangle_angle_differences = []
    for tri, prod in prod_included.items():
        row = edge_rows_by_id[tri]
        length_differences.append(abs(row["edge_length_independent_A"] - prod.length))
        plane_angle_differences.append(abs(row["beta_unsigned_power_plane_normals_rad"] - abs(prod.beta_radians)))
        triangle_angle_differences.append(abs(row["beta_unsigned_dual_triangle_rad"] - abs(prod.beta_radians)))
    geometry_match = max(length_differences, default=0.0) <= 1e-8
    stages.append({"stage": "05_independent_edge_lengths", "edge_count": len(independent_edges),
                   "independent_total_length_A": sum(row["edge_length_independent_A"] for row in independent_edges),
                   "production_total_length_A": production.edge_totals["length"],
                   "maximum_edge_length_difference_A": max(length_differences, default=0.0), "match": geometry_match})
    if not geometry_match:
        return _finish(output, stages, "05_independent_edge_lengths", "one or more independently solved dual-edge lengths disagree", census)

    plane_match = max(plane_angle_differences, default=0.0) <= 1e-10
    triangle_match = max(triangle_angle_differences, default=0.0) <= 1e-10
    stages.append({"stage": "06_unsigned_dihedral_angles", "edge_count": len(independent_edges),
                   "max_power_plane_normal_difference_rad": max(plane_angle_differences, default=0.0),
                   "max_dual_triangle_angle_difference_rad": max(triangle_angle_differences, default=0.0),
                   "power_plane_vs_triangle_max_difference_rad": max((row["unsigned_angle_difference_rad"] for row in independent_edges), default=0.0),
                   "match_production": plane_match and triangle_match})
    if not (plane_match and triangle_match):
        return _finish(output, stages, "06_unsigned_dihedral_angles", "independent unsigned angle disagrees with production beta magnitude", census)

    production_signed_differences = [abs(edge_rows_by_id[tri]["beta_signed_Cazals_rad"] - row.beta_radians)
                                     for tri, row in prod_included.items()]
    signed_match = max(production_signed_differences, default=0.0) <= 1e-10
    signed_numerator = sum(row["signed_beta_l_independent_A_rad"] for row in independent_edges)
    unsigned_numerator = sum(row["unsigned_beta_l_independent_A_rad"] for row in independent_edges)
    denominator = sum(row["edge_length_independent_A"] for row in independent_edges)
    mean_rad = signed_numerator / denominator if denominator else float("nan")
    signed_summary = {
        "edge_count": len(independent_edges),
        "total_edge_length_A": denominator,
        "signed_numerator_A_rad": signed_numerator,
        "unsigned_numerator_A_rad": unsigned_numerator,
        "s_H_signed_rad": mean_rad,
        "s_H_signed_deg": math.degrees(mean_rad),
        "production_signed_numerator_A_rad": production.edge_totals["signed"],
        "production_unsigned_numerator_A_rad": production.edge_totals["unsigned"],
        "production_total_edge_length_A": production.edge_totals["length"],
        "production_s_H_signed_deg": production.edge_totals["signed_mean_deg"],
        "maximum_signed_beta_difference_rad": max(production_signed_differences, default=0.0),
        "sign_classification_match": signed_match,
    }
    (output / "stage_07_final_sH.json").write_text(json.dumps(signed_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stages.append({"stage": "07_cazals_sign_and_final_sH", **signed_summary, "match": signed_match and
                   math.isclose(denominator, production.edge_totals["length"], rel_tol=0.0, abs_tol=1e-8) and
                   math.isclose(signed_numerator, production.edge_totals["signed"], rel_tol=0.0, abs_tol=1e-8)})
    if not stages[-1]["match"]:
        return _finish(output, stages, "07_cazals_sign_and_final_sH", "signed classification or final s_H differs", census)
    return _finish(output, stages, None, "all independently audited stages agree with production", census)


def _finish(output, stages, stop_stage, reason, census):
    _write_csv(output / "stage_counts_and_comparisons.csv", stages)
    metadata = {
        "purpose": "Independent diagnostic reproduction; no production implementation is changed.",
        "input_census": census,
        "atom_selection": {
            "records": "ATOM only; HETATM excluded",
            "hydrogens": "excluded by element or atom-name prefix",
            "altloc_preference": "blank, then A, then first encountered alternate",
            "group_A_chains": sorted(GROUP_A_CHAINS),
            "group_B_chains": sorted(GROUP_B_CHAINS),
            "stable_row_mapping": "selected_atom_id in stage_01_atoms_and_radii.csv; PDB serial also retained",
        },
        "radii": {
            "base_radius_table": "vorpy.src.chemistry.Cazals.special_radii and element_radii values copied into the atom export",
            "probe_radius_A": float(PROBE_RADIUS),
            "expanded_radius": "R_i = r_i + probe_radius",
            "power_weight": "w_i = R_i^2",
        },
        "regular_and_alpha_complex": {
            "builder": "direct gudhi.AlphaComplex(points, weights, precision='safe') call",
            "alpha0_native_threshold": 0.0,
            "alpha_tolerance": ALPHA_TOLERANCE,
            "alpha0_membership": "GUDHI filtration <= alpha0 + tolerance",
            "simplex_export": "stage_02_gudhi_regular_and_alpha_simplices.csv; one row per regular simplex with native filtration",
        },
        "m5_selection": {
            "facet_membership": "bicolor regular 1-simplex with alpha0 membership",
            "orthogonal_ball_squared_radius": "q = ||x-p_i||^2 - w_i at a regular tetrahedron power vertex",
            "orthogonal_ball_radius": "m = sqrt(max(0, max incident q))",
            "reference_radius": "r = min expanded radius of the pair",
            "acceptance": "m/r <= M",
            "M": M5_THRESHOLD,
            "candidate_export": "stage_04_facet_orthogonal_ball_candidates.csv",
        },
        "interface_edge_and_curvature": {
            "edge_population": "regular triangles incident to exactly two final selected bicolor facets and exactly two resolved regular tetrahedra",
            "power_vertex_equations": "2(p_j-p_i) dot x = (||p_j||^2-w_j) - (||p_i||^2-w_i); solved independently by least squares",
            "edge_length": "Euclidean distance between independently solved power vertices",
            "plane_normal": "unit vector from group-A generator center to group-B generator center",
            "unsigned_dihedral": "acos(clamp(n_1 dot n_2, -1, 1))",
            "dual_triangle_unsigned_angle": "angle at the unique singleton-group generator in the regular triangle",
            "Cazals_sign": "AAB=-1; ABB=+1",
            "final_statistic": "s_H = sum(sign * beta_abs * length) / sum(length)",
        },
        "comparison_policy": "Stop stage progression at the first mismatch; stage outputs before that point are sufficient to reproduce the next audited stage offline.",
    }
    (output / "reproduction_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
    )
    result = {"pdb_census": census, "stop_stage": stop_stage, "reason": reason,
              "stages": stages}
    (output / "audit_result.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    lines = ["Independent 2KAI Cazals reproduction audit", "",
             f"Input census: {census}", f"First disagreement: {stop_stage or 'none'}", f"Result: {reason}", "",
             "Stage-by-stage results:"]
    for row in stages:
        lines.append(json.dumps(row, sort_keys=True, default=str))
    (output / "audit_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=Path("output/cazals_2kai_independent_audit"))
    args = parser.parse_args(argv)
    run_independent_2kai_audit(args.pdb, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
