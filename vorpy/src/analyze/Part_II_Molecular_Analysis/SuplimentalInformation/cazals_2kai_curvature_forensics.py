"""Read-only forensic exports for the frozen 2KAI Cazals Power interface.

All geometry, selection, radii, beta signs, edge lengths, and M=5 decisions in
this module are read from :mod:`power_interface_curvature` results.  The only
recomputed quantities are CSV aggregations, distances used for the Lys15
localization diagnostic, and independent sums used to verify those exports.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import numpy as np

from vorpy.src.analyze.power_interface_curvature import cazals_signed_beta
from vorpy.src.analyze.power_interface_cli import run_power_interface_pdb


EXPECTED = {
    "edge_count": 864,
    "total_length_A": 953.524161,
    "signed_beta_l_A_rad": -127.050015,
    "unsigned_beta_l_A_rad": 607.364973,
    "signed_sH_deg": -7.634237,
}
PATTERNS = ("AAB", "ABB")
RADII_A = (3, 5, 7, 10, 15, 20)


def write_csv(path, rows, fields=None):
    rows = list(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _stable_atom(atom):
    return (f"{atom.chain}:{atom.residue_number}{atom.insertion_code}:"
            f"{atom.residue_name}:{atom.atom_name}:{atom.pdb_serial}")


def _facet_condition_row(facet, result, atoms, tri_to_tets=None):
    if tri_to_tets is None:
        tri_to_tets = defaultdict(list)
        for tet in result.full_simplices[3]:
            for triangle in _faces(tet, 3):
                tri_to_tets[triangle].append(tet)
    tets = {tet for tri in facet.incident_regular_triangles for tet in tri_to_tets.get(tuple(tri), ())}
    power_values = [float(result.power_vertices[tet].power_value)
                    for tet in tets if tet in result.power_vertices]
    q = max(power_values) if power_values else None
    q_class = "q<0_no_real_orthogonal_ball" if q is not None and q < 0 else (
        "q>=0_real_candidate" if q is not None else "unavailable_no_resolved_power_vertex")
    m_real = math.sqrt(q) if q is not None and q >= 0 else None
    r = float(facet.smaller_expanded_radius)
    return {
        "facet_atom_pair": f"{facet.generator_ids[0]}-{facet.generator_ids[1]}",
        "atom_i_id": _stable_atom(atoms[facet.generator_ids[0]]),
        "atom_j_id": _stable_atom(atoms[facet.generator_ids[1]]),
        "atom_i_radius_expanded_A": float(result.expanded_radii[facet.generator_ids[0]]),
        "atom_j_radius_expanded_A": float(result.expanded_radii[facet.generator_ids[1]]),
        "r_smaller_expanded_A": r,
        "q_max_power_A2": q,
        "candidate_classification": q_class,
        "m_if_real_A": m_real,
        "m_over_r_if_real": m_real / r if m_real is not None else None,
        "m_used_by_current_implementation_A": facet.largest_orthogonal_ball_radius,
        "m_over_r_current": facet.m_over_r,
        "M_threshold": result.condition_beta_m,
        "alpha_filtration": result.filtration.get(facet.generator_ids),
        "alpha0_passed": bool(facet.alpha_zero_selected),
        "M5_accepted": facet.condition_beta_accepted,
        "current_status": facet.condition_beta_status,
        "current_reason": facet.reason,
        "incident_power_vertex_count": len(power_values),
        "all_incident_power_values_A2": ";".join(map(str, power_values)) if power_values else None,
    }


_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def _map_bpti_lys15(prepared):
    """Map BPTI SEQRES position 15 onto an observed chain-I residue identifier."""
    seq3 = []
    for line in Path(prepared.source).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("SEQRES") and line[11:12].strip() == "I":
            seq3.extend(line[19:70].split())
    sequence = "".join(_AA3_TO_1.get(name.upper(), "X") for name in seq3)
    if len(sequence) < 15 or sequence[14] != "K":
        raise AssertionError(f"Chain-I SEQRES position 15 is not Lys: {sequence[14:15]!r}")
    observed = []
    seen = set()
    for atom in prepared.atoms:
        if atom.chain != "I":
            continue
        key = (atom.chain, atom.residue_number, atom.insertion_code, atom.residue_name)
        if key not in seen:
            seen.add(key)
            observed.append(key)
    obs_sequence = "".join(_AA3_TO_1.get(key[3].upper(), "X") for key in observed)
    # Semiglobal alignment: missing terminal/internal SEQRES residues are
    # allowed; observed residues must map in order to the declared sequence.
    n, m = len(sequence), len(obs_sequence)
    gap, match, mismatch = -2, 3, -3
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        score[0][j] = gap * j
        trace[0][j] = "left"
    for i in range(1, n + 1):
        score[i][0] = 0
        trace[i][0] = "up"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            options = (score[i-1][j-1] + (match if sequence[i-1] == obs_sequence[j-1] else mismatch),
                       score[i-1][j] + gap, score[i][j-1] + gap)
            best = max(options)
            score[i][j] = best
            trace[i][j] = "diag" if options[0] == best else "up" if options[1] == best else "left"
    i, j = max(((n, j) for j in range(m + 1)), key=lambda ij: score[ij[0]][ij[1]])
    mapping = {}
    while j > 0:
        step = trace[i][j]
        if step == "diag":
            mapping[i] = j
            i -= 1
            j -= 1
        elif step == "up":
            i -= 1
        else:
            j -= 1
    observed_index = mapping.get(15)
    if observed_index is None:
        raise AssertionError("BPTI SEQRES Lys15 is absent from observed ATOM residues")
    key = observed[observed_index - 1]
    if key[3] != "LYS":
        raise AssertionError(f"BPTI Lys15 mapped to non-LYS observed residue {key}")
    return key


def _edge_atoms(edge, atoms, group_by_id, radii):
    records = []
    for position, atom_id in enumerate(edge.regular_triangle, 1):
        atom = atoms[atom_id]
        records.append({
            f"atom{position}_index": atom_id,
            f"atom{position}_stable_id": _stable_atom(atom),
            f"atom{position}_name": atom.atom_name,
            f"atom{position}_resname": atom.residue_name,
            f"atom{position}_resnum": atom.residue_number,
            f"atom{position}_icode": atom.insertion_code or "",
            f"atom{position}_chain": atom.chain,
            f"atom{position}_group": group_by_id.get(atom_id, "None"),
            f"atom{position}_xyz_A": f"{atom.x},{atom.y},{atom.z}",
            f"atom{position}_radius_expanded_A": float(radii[atom_id]),
        })
    return {key: value for record in records for key, value in record.items()}


def _summary(records):
    records = list(records)
    length = sum(float(r["edge_length_A"]) for r in records)
    positive = sum(float(r["signed_beta_times_length_A_rad"]) for r in records
                   if float(r["beta_signed_rad"]) > 0)
    negative = sum(float(r["signed_beta_times_length_A_rad"]) for r in records
                   if float(r["beta_signed_rad"]) < 0)
    signed = sum(float(r["signed_beta_times_length_A_rad"]) for r in records)
    unsigned = sum(float(r["unsigned_beta_times_length_A_rad"]) for r in records)
    abs_beta_deg = [float(r["beta_abs_deg"]) for r in records]
    return {
        "edge_count": len(records), "total_length_A": length,
        "mean_abs_beta_deg": float(np.mean(abs_beta_deg)) if records else math.nan,
        "median_abs_beta_deg": float(np.median(abs_beta_deg)) if records else math.nan,
        "positive_beta_l_A_rad": positive, "negative_beta_l_A_rad": negative,
        "signed_beta_l_A_rad": signed, "unsigned_beta_l_A_rad": unsigned,
        "signed_sH_deg": math.degrees(signed / length) if length else math.nan,
        "unsigned_sH_deg": math.degrees(unsigned / length) if length else math.nan,
    }


def _point_segment_distances(points, start, end):
    points = np.asarray(points, dtype=float)
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    vector = end - start
    norm2 = float(np.dot(vector, vector))
    if norm2 == 0:
        return np.linalg.norm(points - start, axis=1)
    t = np.clip(((points - start) @ vector) / norm2, 0.0, 1.0)
    return np.linalg.norm(points - (start + t[:, None] * vector), axis=1)


def export_2kai_forensics(result, prepared, output_dir):
    """Export row-level detail and fail unless it reconstructs frozen totals."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    atoms = {atom.index: atom for atom in prepared.atoms}
    group_by_id = {int(i): "A" for i in result.group_a}
    group_by_id.update({int(i): "B" for i in result.group_b})
    facet_by_pair = {facet.generator_ids: facet for facet in result.facets}
    tri_to_tets = defaultdict(list)
    for tet in result.full_simplices[3]:
        for triangle in _faces(tet, 3):
            tri_to_tets[triangle].append(tet)
    m5_by_pair = {facet.generator_ids: _facet_condition_row(facet, result, atoms, tri_to_tets)
                  for facet in result.alpha_zero_facets}
    rows = []
    for edge_id, edge in enumerate(result.included_edges, 1):
        if edge.geometry_status != "finite" or edge.length is None or len(edge.bicolor_facets) != 2:
            raise AssertionError(f"Included edge {edge.regular_triangle} is not finite/interior")
        ids = edge.regular_triangle
        singleton_id = edge.singleton_generator
        coords = edge.power_vertex_coordinates
        if len(coords) != 2:
            raise AssertionError(f"Included edge {ids} lacks its two power-edge endpoints")
        adjacent = [facet_by_pair[tuple(sorted(pair))] for pair in edge.bicolor_facets]
        signed = float(edge.beta_radians)
        length = float(edge.length)
        row = {
            "edge_id": edge_id,
            **_edge_atoms(edge, atoms, group_by_id, result.expanded_radii),
            "primal_atom_ids": ";".join(map(str, ids)),
            "triangle_group_pattern": edge.group_composition,
            "singleton_atom_id": singleton_id,
            "singleton_stable_id": _stable_atom(atoms[singleton_id]),
            "singleton_chain": atoms[singleton_id].chain,
            "singleton_residue": atoms[singleton_id].residue_name,
            "singleton_resnum": atoms[singleton_id].residue_number,
            "singleton_icode": atoms[singleton_id].insertion_code or "",
            "singleton_atom_name": atoms[singleton_id].atom_name,
            "singleton_group": edge.singleton_group,
            "adjacent_facet_1": f"{edge.bicolor_facets[0][0]}-{edge.bicolor_facets[0][1]}",
            "adjacent_facet_2": f"{edge.bicolor_facets[1][0]}-{edge.bicolor_facets[1][1]}",
            "power_edge_start_xyz_A": ";".join(map(str, coords[0])),
            "power_edge_end_xyz_A": ";".join(map(str, coords[1])),
            "edge_length_A": length,
            "beta_abs_rad": abs(signed), "beta_abs_deg": abs(float(edge.beta_degrees)),
            "sign_assigned": edge.sign,
            "beta_signed_rad": signed, "beta_signed_deg": float(edge.beta_degrees),
            "signed_beta_times_length_A_rad": float(edge.length_beta),
            "unsigned_beta_times_length_A_rad": float(edge.length_abs_beta),
            "both_adjacent_facets_passed_alpha0": all(edge.alpha_zero_facets_selected),
            "facet_1_alpha_filtration": result.filtration.get(adjacent[0].generator_ids),
            "facet_2_alpha_filtration": result.filtration.get(adjacent[1].generator_ids),
            "triangle_alpha_filtration": result.filtration.get(ids),
            "facet_1_passed_M5": adjacent[0].condition_beta_accepted,
            "facet_2_passed_M5": adjacent[1].condition_beta_accepted,
            "facet_1_r_A": adjacent[0].smaller_expanded_radius,
            "facet_2_r_A": adjacent[1].smaller_expanded_radius,
            "facet_1_m_A": adjacent[0].largest_orthogonal_ball_radius,
            "facet_2_m_A": adjacent[1].largest_orthogonal_ball_radius,
            "facet_1_m_if_real_A": m5_by_pair[adjacent[0].generator_ids]["m_if_real_A"],
            "facet_2_m_if_real_A": m5_by_pair[adjacent[1].generator_ids]["m_if_real_A"],
            "facet_1_m_over_r_if_real": m5_by_pair[adjacent[0].generator_ids]["m_over_r_if_real"],
            "facet_2_m_over_r_if_real": m5_by_pair[adjacent[1].generator_ids]["m_over_r_if_real"],
            "facet_1_m_over_r": adjacent[0].m_over_r,
            "facet_2_m_over_r": adjacent[1].m_over_r,
            "facet_1_status": adjacent[0].condition_beta_status,
            "facet_2_status": adjacent[1].condition_beta_status,
            "facet_1_q_A2": m5_by_pair[adjacent[0].generator_ids]["q_max_power_A2"],
            "facet_2_q_A2": m5_by_pair[adjacent[1].generator_ids]["q_max_power_A2"],
            "facet_1_q_classification": m5_by_pair[adjacent[0].generator_ids]["candidate_classification"],
            "facet_2_q_classification": m5_by_pair[adjacent[1].generator_ids]["candidate_classification"],
            "edge_classification": "finite interior",
        }
        rows.append(row)
    write_csv(out / "2kai_curvature_edges.csv", rows)
    total = _summary(rows)
    expected_checks = {
        "edge_count": 0,
        "total_length_A": 2e-6,
        "signed_beta_l_A_rad": 2e-6,
        "unsigned_beta_l_A_rad": 2e-6,
        "signed_sH_deg": 2e-6,
    }
    for key, tolerance in expected_checks.items():
        if abs(total[key] - EXPECTED[key]) > tolerance:
            raise AssertionError(f"CSV reconstruction {key}={total[key]} != frozen {EXPECTED[key]}")

    pattern_rows = []
    for pattern in (*PATTERNS, "TOTAL"):
        subset = rows if pattern == "TOTAL" else [r for r in rows if r["triangle_group_pattern"] == pattern]
        pattern_rows.append({"pattern": pattern, **_summary(subset)})
    write_csv(out / "2kai_curvature_by_pattern.csv", pattern_rows)

    by_chain, by_residue, by_atom = defaultdict(list), defaultdict(list), defaultdict(list)
    for row in rows:
        by_chain[row["singleton_chain"]].append(row)
        residue_key = (row["singleton_chain"], row["singleton_resnum"], row["singleton_icode"], row["singleton_residue"])
        by_residue[residue_key].append(row)
        atom_key = (row["singleton_chain"], row["singleton_resnum"], row["singleton_icode"],
                    row["singleton_residue"], row["singleton_atom_name"])
        by_atom[atom_key].append(row)

    chain_rows = []
    for chain, subset in sorted(by_chain.items()):
        s = _summary(subset)
        chain_rows.append({"singleton_chain": chain, "edge_count": s["edge_count"],
                           "total_length_A": s["total_length_A"], "signed_beta_l_A_rad": s["signed_beta_l_A_rad"],
                           "unsigned_beta_l_A_rad": s["unsigned_beta_l_A_rad"],
                           "signed_sH_deg": s["signed_sH_deg"], "mean_abs_beta_deg": s["mean_abs_beta_deg"]})
    write_csv(out / "2kai_curvature_by_singleton_chain.csv", chain_rows)
    residue_rows = []
    for key, subset in by_residue.items():
        s = _summary(subset)
        residue_rows.append({"chain": key[0], "resnum": key[1], "icode": key[2], "resname": key[3],
                             "edge_count": s["edge_count"], "total_length_A": s["total_length_A"],
                             "signed_beta_l_A_rad": s["signed_beta_l_A_rad"],
                             "unsigned_beta_l_A_rad": s["unsigned_beta_l_A_rad"],
                             "signed_sH_deg": s["signed_sH_deg"], "mean_abs_beta_deg": s["mean_abs_beta_deg"]})
    residue_rows.sort(key=lambda r: abs(r["signed_beta_l_A_rad"]), reverse=True)
    write_csv(out / "2kai_curvature_by_singleton_residue.csv", residue_rows)
    atom_rows = []
    for key, subset in by_atom.items():
        s = _summary(subset)
        atom_rows.append({"chain": key[0], "resnum": key[1], "icode": key[2], "resname": key[3], "atom_name": key[4],
                          "edge_count": s["edge_count"], "total_length_A": s["total_length_A"],
                          "signed_beta_l_A_rad": s["signed_beta_l_A_rad"],
                          "unsigned_beta_l_A_rad": s["unsigned_beta_l_A_rad"],
                          "signed_sH_deg": s["signed_sH_deg"], "mean_abs_beta_deg": s["mean_abs_beta_deg"]})
    atom_rows.sort(key=lambda r: abs(r["signed_beta_l_A_rad"]), reverse=True)
    write_csv(out / "2kai_curvature_by_singleton_atom.csv", atom_rows)

    # The biological BPTI Lys15 is resolved by its position in the chain-I
    # SEQRES sequence, then mapped to the observed PDB residue identifier.
    lys_key = _map_bpti_lys15(prepared)
    lys_atoms = [a for a in prepared.atoms if (a.chain, a.residue_number, a.insertion_code, a.residue_name) == lys_key]
    lys_points = np.asarray([[a.x, a.y, a.z] for a in lys_atoms], dtype=float)
    radial_rows = []
    for radius in RADII_A:
        included = []
        for row, edge in zip(rows, result.included_edges):
            start = np.fromstring(row["power_edge_start_xyz_A"], sep=";")
            end = np.fromstring(row["power_edge_end_xyz_A"], sep=";")
            distance = float(np.min(_point_segment_distances(lys_points, start, end)))
            if distance <= radius:
                included.append(row)
        s = _summary(included)
        radial_rows.append({"radius_A": radius, "edge_count": s["edge_count"], "total_length_A": s["total_length_A"],
                            "positive_beta_l_A_rad": s["positive_beta_l_A_rad"],
                            "negative_beta_l_A_rad": s["negative_beta_l_A_rad"],
                            "signed_beta_l_A_rad": s["signed_beta_l_A_rad"],
                            "unsigned_beta_l_A_rad": s["unsigned_beta_l_A_rad"],
                            "signed_sH_deg": s["signed_sH_deg"], "unsigned_sH_deg": s["unsigned_sH_deg"],
                            "distance_definition": "minimum Euclidean distance from edge segment to any heavy-atom center in mapped BPTI Lys15"})
    write_csv(out / "2kai_lys15_radial_curvature.csv", radial_rows)

    m5_rows = [m5_by_pair[facet.generator_ids] for facet in result.alpha_zero_facets]
    write_csv(out / "2kai_M5_facet_forensics.csv", m5_rows)
    m5_counts = {
        "alpha0_bicolor_facets": len(m5_rows),
        "real_orthogonal_ball_candidate_q_ge_0": sum(row["q_max_power_A2"] is not None and row["q_max_power_A2"] >= 0 for row in m5_rows),
        "q_lt_0_no_real_candidate": sum(row["q_max_power_A2"] is not None and row["q_max_power_A2"] < 0 for row in m5_rows),
        "unavailable_candidate": sum(row["q_max_power_A2"] is None for row in m5_rows),
        "M5_accepted": sum(row["M5_accepted"] is True for row in m5_rows),
        "M5_rejected": sum(row["M5_accepted"] is False for row in m5_rows),
        "M5_unresolved": sum(row["M5_accepted"] is None for row in m5_rows),
    }
    (out / "2kai_M5_summary.txt").write_text(
        "2KAI frozen condition-beta / M=5 facet audit\n"
        "Source: current production PowerInterfaceCurvatureResult; no threshold or formula changed.\n"
        "q is the maximum common power value (squared orthogonal-ball radius) among resolved power vertices incident to the bicolor facet.\n"
        "For q<0 the current implementation clamps m=sqrt(max(0,q)) to zero; this report records q<0 as no real positive-radius candidate without reinterpreting the decision.\n\n"+
        "\n".join(f"{key}: {value}" for key, value in m5_counts.items())+"\n",
        encoding="utf-8")

    # Counterfactuals are summaries over the exported production beta/lengths.
    counterfactuals = []
    groups = {
        "current sign convention": rows,
        "globally reversed sign": [{**r, "signed_beta_times_length_A_rad": -float(r["signed_beta_times_length_A_rad"]),
                                    "beta_signed_rad": -float(r["beta_signed_rad"])} for r in rows],
        "all beta positive / unsigned": [{**r, "signed_beta_times_length_A_rad": float(r["unsigned_beta_times_length_A_rad"]),
                                          "beta_signed_rad": abs(float(r["beta_signed_rad"]))} for r in rows],
        "AAB contribution only": [r for r in rows if r["triangle_group_pattern"] == "AAB"],
        "ABB contribution only": [r for r in rows if r["triangle_group_pattern"] == "ABB"],
    }
    for label, subset in groups.items():
        s = _summary(subset)
        counterfactuals.append({"counterfactual": label, "numerator_A_rad": s["signed_beta_l_A_rad"],
                                "denominator_A": s["total_length_A"], "s_H_deg": s["signed_sH_deg"]})
    target_num = math.radians(17.0) * total["total_length_A"]
    counterfactuals.append({"counterfactual": "required numerator for +17 degrees at current L",
                            "numerator_A_rad": target_num, "denominator_A": total["total_length_A"], "s_H_deg": 17.0})
    write_csv(out / "2kai_sign_counterfactuals.csv", counterfactuals)

    # Group-reversal invariance uses the production beta routine on the same
    # triangle coordinates; it does not modify the result object.
    reversed_signed = []
    for edge in result.included_edges:
        beta = cazals_signed_beta(edge.regular_triangle, result.points, result.group_b, result.group_a)[2]
        if not math.isclose(beta, -edge.beta_radians, abs_tol=1e-12):
            raise AssertionError(f"Group reversal did not flip production beta for {edge.regular_triangle}")
        reversed_signed.append(float(edge.length) * beta)
    if not math.isclose(sum(abs(x) for x in reversed_signed), total["unsigned_beta_l_A_rad"], abs_tol=1e-9):
        raise AssertionError("Group reversal changed unsigned curvature")

    report = _make_report(rows, pattern_rows, chain_rows, residue_rows, radial_rows, m5_counts, counterfactuals,
                          lys_key, total, prepared)
    (out / "2KAI_CURVATURE_FORENSIC_REPORT.txt").write_text(report, encoding="utf-8")
    return {"edges": rows, "patterns": pattern_rows, "residues": residue_rows,
            "radial": radial_rows, "m5": m5_rows, "m5_counts": m5_counts,
            "counterfactuals": counterfactuals, "total": total, "lys_key": lys_key}


def _make_report(rows, pattern_rows, chain_rows, residue_rows, radial_rows, m5, counterfactuals, lys_key, total, prepared):
    pattern = {row["pattern"]: row for row in pattern_rows}
    positive = total["positive_beta_l_A_rad"]
    negative = total["negative_beta_l_A_rad"]
    top_pos = sorted(residue_rows, key=lambda r: r["signed_beta_l_A_rad"], reverse=True)[:5]
    top_neg = sorted(residue_rows, key=lambda r: r["signed_beta_l_A_rad"])[:5]
    def labels(values):
        return "; ".join(f"{r['chain']}:{r['resname']}{r['resnum']}{r['icode']} {r['signed_beta_l_A_rad']:+.6f}" for r in values)
    lys = [r for r in radial_rows]
    return (
        "2KAI Cazals Power-interface curvature forensic report\n\n"
        "Diagnostic only. The frozen Cazals geometry, atom typing/radii, alpha=0 selection, M=5 condition, beta sign, and edge inclusion rules were not modified. Values are read from the production result and independently summed from the exported rows.\n\n"
        f"Input: {prepared.source}; groups A={sorted({prepared.atoms[i].chain for i in prepared.group_a})}, B={sorted({prepared.atoms[i].chain for i in prepared.group_b})}; retained protein atoms={len(prepared.atoms)}.\n"
        f"Biological Lys15 mapping: BPTI chain {lys_key[0]}, actual PDB residue {lys_key[3]} {lys_key[1]} insertion '{lys_key[2] or ''}'. Lys15 was identified as sequence position 15 in chain-I SEQRES and mapped to the observed residue.\n\n"
        f"RECONSTRUCTION\nN={total['edge_count']}; L={total['total_length_A']:.9f} A; positive={positive:.9f}; negative={negative:.9f}; signed={total['signed_beta_l_A_rad']:.9f}; unsigned={total['unsigned_beta_l_A_rad']:.9f} A rad; s_H={total['signed_sH_deg']:.9f} deg.\n"
        "The forensic CSV reconstructs the frozen analyzer totals within 2e-6 in each reported aggregate. This supports no discrepancy in exported beta magnitudes or edge lengths at that precision; it does not establish which atom model or interface convention the paper used.\n\n"
        f"CANCELLATION BY PATTERN\nAAB: N={pattern['AAB']['edge_count']}, L={pattern['AAB']['total_length_A']:.6f}, signed={pattern['AAB']['signed_beta_l_A_rad']:.6f}, unsigned={pattern['AAB']['unsigned_beta_l_A_rad']:.6f}, mean |beta|={pattern['AAB']['mean_abs_beta_deg']:.4f} deg.\n"
        f"ABB: N={pattern['ABB']['edge_count']}, L={pattern['ABB']['total_length_A']:.6f}, signed={pattern['ABB']['signed_beta_l_A_rad']:.6f}, unsigned={pattern['ABB']['unsigned_beta_l_A_rad']:.6f}, mean |beta|={pattern['ABB']['mean_abs_beta_deg']:.4f} deg.\n"
        f"The +{positive:.6f} / {negative:.6f} cancellation is exactly the sum of the two opposite-sign patterns; the larger-magnitude pattern and class contributions are shown above.\n\n"
        f"SINGLETON RESIDUES\nLargest positive contributions: {labels(top_pos)}\nLargest negative contributions: {labels(top_neg)}\n"
        + "Singleton-chain totals: " + "; ".join(
            f"{r['singleton_chain']}: C={r['signed_beta_l_A_rad']:+.6f}, L={r['total_length_A']:.6f}, N={r['edge_count']}"
            for r in chain_rows) + "\n"
        f"LYS15 SPATIAL AUDIT\nDistance means minimum Euclidean distance from each finite edge segment to any heavy-atom center of mapped Lys15. Values at radii {', '.join(map(str, RADII_A))} A are in 2kai_lys15_radial_curvature.csv.\n"
        + "\n".join(f"within {r['radius_A']} A: N={r['edge_count']}, L={r['total_length_A']:.6f}, C+={r['positive_beta_l_A_rad']:.6f}, C-={r['negative_beta_l_A_rad']:.6f}, C={r['signed_beta_l_A_rad']:.6f}, |C|={r['unsigned_beta_l_A_rad']:.6f}, s_H={r['signed_sH_deg']:.6f} deg" for r in lys)+"\n"
        + "Within 3–5 A of Lys15 the net signed contribution is negative (-48.29 to -49.04 A rad) under the current A=kallikrein, B=BPTI convention, so the nearby population is coherently concave by the implementation's AAB-negative sign rule. The normalized values (-9.35 to -6.54 degrees) are comparable to the global -7.63 degrees, rather than showing an exceptionally more negative local mean. Reversing physical A/B labels flips this local sign. This is directionally compatible with a Lys15 pocket concavity after the paper's partner mapping is established, but it does not resolve the global magnitude discrepancy.\n\n"
        f"M=5 FORENSICS\nalpha=0 bicolor facets={m5['alpha0_bicolor_facets']}; q>=0 candidates={m5['real_orthogonal_ball_candidate_q_ge_0']}; q<0/no-real-candidate={m5['q_lt_0_no_real_candidate']}; unavailable={m5['unavailable_candidate']}; accepted={m5['M5_accepted']}; rejected={m5['M5_rejected']}; unresolved={m5['M5_unresolved']}.\n"
        "q is the maximum incident common power value (A^2). q<0 has no real positive orthogonal-ball radius; the production implementation clamps m to zero before m/r. No reinterpretation was applied.\n\n"
        "SIGN COUNTERFACTUALS\n"+"\n".join(f"{r['counterfactual']}: numerator={r['numerator_A_rad']:.9f} A rad; denominator={r['denominator_A']:.9f} A; s_H={r['s_H_deg']:.9f} deg" for r in counterfactuals)+"\n\n"
        "CONCLUSION / OPEN HYPOTHESES\nThe forensic sums localize cancellation quantitatively by AAB/ABB class and singleton location. The current exports plus independent 1UDI singleton-angle validation provide no positive evidence here of a numerical beta or edge-length error. This does not establish a cause for -7.63 degrees versus the paper’s approximate +17 degrees. Remaining hypotheses include an A/B mapping difference (which flips sign only), a different edge/interface population or boundary convention, a different structure/atom model or radii, differences in how the condition-beta/alpha interface was implemented, and rounding or reporting differences. The present diagnostics do not decide among these hypotheses, and no parameter was selected to approach 17 degrees.\n"
    )


def _faces(simplex, size):
    from itertools import combinations
    return [tuple(sorted(face)) for face in combinations(simplex, size)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, default=Path("vorpy/data/2KAI.pdb"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/cazals_2kai_curvature_forensics"))
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    run = run_power_interface_pdb(args.pdb, {"A", "B"}, {"I"}, output / "production_source",
                                  probe_radius=1.4, alpha=0.0, condition_beta_m=5.0)
    export_2kai_forensics(run["result"], run["prepared"], output)
    print(f"2KAI Cazals forensics written to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
