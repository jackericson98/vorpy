"""Forensic, non-production audits of Cazals power-interface edge sets.

The production curvature analyzer is deliberately not changed here.  This
module compares its regular-triangle incidence against edges recovered from
the actual finite polygons dual to the selected power facets.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
import csv
import math
from pathlib import Path

import numpy as np


def classify_interface_edge(selected_facet_count, *, geometry_status="finite"):
    if geometry_status not in {"finite", ""}:
        if "unbounded" in geometry_status or "boundary_power_edge" in geometry_status:
            return "UNBOUNDED"
        if "degenerate" in geometry_status or "unresolved" in geometry_status or "nonmanifold" in geometry_status:
            return "DEGENERATE_UNRESOLVED"
    if selected_facet_count == 1:
        return "BOUNDARY"
    if selected_facet_count == 2:
        return "INTERIOR"
    if selected_facet_count > 2:
        return "NONMANIFOLD"
    return "NOT_ON_SELECTED_INTERFACE"


def summarize_edge_records(records):
    records = list(records)
    length = sum(float(r["length"]) for r in records)
    plus = sum(float(r["length_beta"]) for r in records if float(r["beta"]) > 0)
    minus = sum(float(r["length_beta"]) for r in records if float(r["beta"]) < 0)
    signed, unsigned = plus + minus, sum(abs(float(r["length_beta"])) for r in records)
    return {"N_edges": len(records), "L": length, "positive_edges": sum(float(r["beta"]) > 0 for r in records),
            "negative_edges": sum(float(r["beta"]) < 0 for r in records), "C_plus": plus,
            "C_minus": minus, "C_signed": signed, "C_unsigned": unsigned,
            "s_H_signed_rad": signed / length if length else math.nan,
            "s_H_signed_deg": math.degrees(signed / length) if length else math.nan,
            "s_H_unsigned_deg": math.degrees(unsigned / length) if length else math.nan}


def classify_triangle_subset(records):
    return {tuple(r["triangle"]): bool(r["triangle_alpha_zero"]) for r in records}


def select_interior_finite_edges(records, *, require_triangle_alpha=False):
    """Filter already classified rows without manufacturing dual edges."""
    return [r for r in records
            if r.get("classification") == "INTERIOR"
            and r.get("geometry_status") == "finite"
            and (not require_triangle_alpha or bool(r.get("triangle_alpha_zero")))]


def condition_beta_audit_rows(result, atom_labels):
    """Expose m/r for every alpha-selected bicolor facet, including failures."""
    tri_to_tets=defaultdict(list)
    for tet in result.full_simplices[3]:
        for tri in _faces(tet,3):
            tri_to_tets[tri].append(tet)
    rows=[]
    for facet in result.facets:
        if not facet.alpha_zero_selected:
            continue
        i,j=facet.generator_ids
        incident={tet for tri in facet.incident_regular_triangles for tet in tri_to_tets.get(tuple(tri),())}
        power_values=[result.power_vertices[tet].power_value for tet in incident if tet in result.power_vertices]
        raw_max_power=max(power_values) if power_values else None
        rows.append({"atom_i":i,"atom_j":j,"atom_i_label":atom_labels[i],"atom_j_label":atom_labels[j],
                     "r_smaller_expanded_A":facet.smaller_expanded_radius,
                     "m_largest_orthogonal_ball_A":facet.largest_orthogonal_ball_radius,
                     "max_power_radius_squared_A2":raw_max_power,
                     "m_clamped_to_zero":bool(raw_max_power is not None and raw_max_power <= 0),
                     "m_over_r":facet.m_over_r,"threshold_M":result.condition_beta_m,
                     "accepted":facet.condition_beta_accepted,"status":facet.condition_beta_status})
    return rows


def _faces(simplex, n):
    return [tuple(sorted(f)) for f in combinations(simplex, n)]


def _coordinate_key(point, tolerance):
    return tuple(int(round(float(x) / tolerance)) for x in point)


def _ordered_facet_polygon(result, facet, tri_to_tets, tolerance):
    tets = set()
    for tri in facet.incident_regular_triangles:
        incident = tri_to_tets.get(tuple(tri), ())
        if len(incident) != 2:
            return None, f"facet boundary triangle has {len(incident)} incident tetrahedra"
        tets.update(incident)
    if any(tet not in result.power_vertices for tet in tets):
        return None, "facet contains unresolved power vertex"
    unique = {}
    for tet in tets:
        p = np.asarray(result.power_vertices[tet].position, dtype=float)
        unique.setdefault(_coordinate_key(p, tolerance), (p, tet))
    if len(unique) < 3:
        return None, "fewer than three distinct power vertices"
    vertices = list(unique.values())
    ids = facet.generator_ids
    normal = result.points[ids[1]] - result.points[ids[0]]
    norm = np.linalg.norm(normal)
    if norm <= tolerance:
        return None, "coincident generator centers"
    normal = normal / norm
    center = np.mean([v[0] for v in vertices], axis=0)
    axis = np.eye(3)[int(np.argmin(np.abs(normal)))]
    u = np.cross(normal, axis); u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    vertices.sort(key=lambda item: math.atan2(float(np.dot(item[0]-center, v)), float(np.dot(item[0]-center, u))))
    return vertices, ""


def audit_final_interface_polygons(result, *, coordinate_tolerance=2e-6):
    """Build the final selected interface as a polygonal 2-complex.

    Returns each geometric segment with incident selected facets, its dual
    regular triangle, and the segment length.  No distance-based adjacency is
    used; equal endpoint keys require coordinate agreement within tolerance.
    """
    tri_to_tets = defaultdict(list)
    for tet in result.full_simplices[3]:
        for tri in _faces(tet, 3):
            tri_to_tets[tri].append(tet)
    segments = defaultdict(list)
    polygon_issues = []
    for facet in result.final_facets:
        polygon, reason = _ordered_facet_polygon(result, facet, tri_to_tets, coordinate_tolerance)
        if polygon is None:
            polygon_issues.append({"facet": facet.generator_ids, "reason": reason})
            continue
        for index, (p, tet1) in enumerate(polygon):
            q, tet2 = polygon[(index + 1) % len(polygon)]
            common = set(tet1) & set(tet2)
            if len(common) != 3:
                polygon_issues.append({"facet": facet.generator_ids, "reason": f"consecutive dual vertices share {len(common)} generators"})
                continue
            triangle = tuple(sorted(common))
            k1, k2 = _coordinate_key(p, coordinate_tolerance), _coordinate_key(q, coordinate_tolerance)
            key = tuple(sorted((k1, k2)))
            segments[key].append({"facet": facet.generator_ids, "triangle": triangle,
                                  "length": float(np.linalg.norm(q-p)),
                                  "p": tuple(map(float, p)), "q": tuple(map(float, q))})
    rows = []
    for key, incidences in sorted(segments.items()):
        triangles = sorted({x["triangle"] for x in incidences})
        rows.append({"segment_key": key, "incident_facets": tuple(sorted(x["facet"] for x in incidences)),
                     "incident_selected_facet_count": len(incidences),
                     "triangle": triangles[0] if len(triangles) == 1 else None,
                     "triangle_candidates": tuple(triangles),
                     "length": float(np.mean([x["length"] for x in incidences])),
                     "length_spread": float(max(x["length"] for x in incidences)-min(x["length"] for x in incidences)),
                     "coordinates": (incidences[0]["p"], incidences[0]["q"]),
                     "classification": classify_interface_edge(len(incidences))})
    return rows, polygon_issues


def _atom_label(atom):
    return f"{atom.chain}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}:{atom.atom_name}[{atom.pdb_serial}]"


def run_2kai_forensics(result, prepared, output_dir, *, significant_facets=None):
    """Write the requested detailed 2KAI edge-set and polygon-incidence audit."""
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    atoms = {a.index: a for a in prepared.atoms}
    geom_rows, polygon_issues = audit_final_interface_polygons(result)
    geom_by_triangle = defaultdict(list)
    geom_by_facet_pair = defaultdict(list)
    for row in geom_rows:
        if row["triangle"]:
            geom_by_triangle[row["triangle"]].append(row)
        if len(row["incident_facets"]) == 2:
            geom_by_facet_pair[tuple(sorted(row["incident_facets"]))].append(row)
    edge_by_tri = {e.regular_triangle: e for e in result.edges}
    final_facet_ids = {f.generator_ids for f in result.final_facets}
    component_ids = {tuple(f): 1 for f in final_facet_ids}  # 2KAI has one measured component
    rows = []
    for edge in result.edges:
        tri = edge.regular_triangle
        facets = tuple(edge.bicolor_facets)
        selected = [f for f in facets if f in final_facet_ids]
        geom = geom_by_triangle.get(tri, [])
        count = len(selected)
        classification = classify_interface_edge(count, geometry_status=edge.geometry_status)
        if edge.geometry_status != "finite":
            classification = classify_interface_edge(count, geometry_status=edge.geometry_status)
        ids = tri
        rows.append({
            "triangle": tri, "generator_ids": ";".join(map(str, tri)),
            "atom_labels": ";".join(_atom_label(atoms[i]) for i in tri),
            "group_pattern": edge.group_composition,
            "singleton_id": edge.singleton_generator,
            "singleton_atom": _atom_label(atoms[edge.singleton_generator]),
            "incident_selected_facets": ";".join(f"{i}-{j}" for i,j in selected),
            "incident_selected_facet_count": count,
            "incident_full_power_facets": 3 if tri in result.full_simplices[2] else 0,
            "incident_geometric_segment_count": len(geom),
            "geometry_status": edge.geometry_status,
            "triangle_alpha_filtration": result.filtration.get(tri),
            "triangle_alpha_zero": edge.triangle_alpha_zero_selected,
            "facet_alpha_flags": repr(edge.alpha_zero_facets_selected),
            "condition_b_flags": repr(edge.condition_beta_facets_accepted),
            "component_id": 1 if count else "",
            "significant_component": bool(count),
            "length_A": edge.length,
            "beta_rad": edge.beta_radians, "beta_deg": edge.beta_degrees,
            "l_beta": edge.length_beta,
            "classification": classification,
            "polygon_geometric_support": bool(geom),
            "polygon_length": geom[0]["length"] if geom else None,
        })
    fields = list(rows[0]) if rows else []
    _write_csv(out / "2KAI_interface_edge_forensics.csv", rows, fields)

    selected_final = select_interior_finite_edges(rows)
    tri_in_alpha = select_interior_finite_edges(rows, require_triangle_alpha=True)
    edge_obj = {e.regular_triangle: e for e in result.edges}
    def summary_for(row_set, reverse=False):
        data=[]
        for r in row_set:
            e=edge_obj[r["triangle"]]
            data.append({"length": e.length, "beta": -e.beta_radians if reverse else e.beta_radians,
                         "length_beta": -(e.length_beta or 0.) if reverse else (e.length_beta or 0.)})
        return summarize_edge_records(data)
    variants = []
    for name, subset, justification in [
        ("A_CURRENT_864", selected_final, "Current analyzer: finite full dual edges shared by two final alpha=0/condition-b facets"),
        ("B_TRIANGLE_IN_K0", tri_in_alpha, "Additional triangle alpha=0 restriction; not stated for Equation 6"),
    ]:
        for reversed_labels in (False, True):
            s=summary_for(subset, reversed_labels)
            variants.append({"method":name + ("_GROUP_REVERSED" if reversed_labels else ""),
                             "edge_definition":"finite power dual edge with two selected incident bicolor facets" + ("; sign reversed for diagnostic" if reversed_labels else ""), **s,
                             "interpretation_status":"DIAGNOSTIC ONLY" if reversed_labels or name.startswith("B_") else "DIRECT CONSEQUENCE OF PUBLISHED DEFINITION",
                             "source_justification":"sign diagnostic only" if reversed_labels else justification})
    poly_interior_tri = {r["triangle"] for r in geom_rows if r["classification"] == "INTERIOR" and r["triangle"]}
    polygon_rows = [r for r in rows if r["triangle"] in poly_interior_tri and r["geometry_status"] == "finite"
                    and all(f in final_facet_ids for f in edge_obj[r["triangle"]].bicolor_facets)]
    # De-duplicate if multiple coincident triangulated segments represent a higher-order feature.
    polygon_rows = list({r["triangle"]: r for r in polygon_rows}.values())
    variants.append({"method":"C_POLYGONAL_TOPOLOGICAL_INTERIOR",
                     "edge_definition":"polygon edge with exactly two selected facet polygons, matched to finite regular-triangle dual",
                     **summary_for(polygon_rows),
                     "interpretation_status":"DIRECT CONSEQUENCE OF PUBLISHED DEFINITION",
                     "source_justification":"Direct polygon-edge incidence after final facet selection"})
    variants.append({"method":"D_SIGNIFICANT_COMPONENT_ONLY",
                     "edge_definition":"current population restricted to significant connected components",
                     **summary_for(selected_final),
                     "interpretation_status":"DIAGNOSTIC ONLY",
                     "source_justification":"Diagnostic; the sole 2KAI component is significant"})
    _write_csv(out / "2KAI_edge_set_comparison.csv", variants)

    # Attach alpha-triangle discrepancy details.
    discrepant = [r for r in selected_final if not r["triangle_alpha_zero"]]
    discrepancy_fields = ["generator_ids", "group_pattern", "singleton_atom", "atom_labels", "triangle_alpha_filtration",
                          "length_A", "beta_rad", "beta_deg", "l_beta"]
    _write_csv(out / "2KAI_13_alpha_triangle_discrepant_edges.csv",
               [{k:r[k] for k in discrepancy_fields} for r in discrepant], discrepancy_fields)

    condition_rows=condition_beta_audit_rows(result, {i:_atom_label(a) for i,a in atoms.items()})
    _write_csv(out / "2KAI_condition_b_audit.csv", condition_rows)

    # Incidence audit compares each combinatorial selected edge with matching polygon segments.
    inc=[]
    for r in rows:
        if r["incident_selected_facet_count"] < 1:
            continue
        edge=edge_obj[r["triangle"]]
        geom=geom_by_triangle.get(r["triangle"],[])
        matching=[g for g in geom if set(g["incident_facets"]) == set(edge.bicolor_facets)]
        inc.append({"triangle":r["generator_ids"],"group_pattern":r["group_pattern"],
                    "combinatorial_selected_facets":";".join(f"{i}-{j}" for i,j in edge.bicolor_facets),
                    "polygon_segment_matches":len(matching),"expected_shared":r["incident_selected_facet_count"]==2,
                    "geometrically_shared":bool(matching),"regular_triangle_in_K0":r["triangle_alpha_zero"],
                    "power_edge_status":edge.geometry_status,"combinatorial_length_A":edge.length,
                    "polygon_length_A":matching[0]["length"] if matching else None,
                    "length_difference_A":abs(edge.length-matching[0]["length"]) if matching and edge.length is not None else None,
                    "facet_pair_match":";".join(f"{i}-{j}" for i,j in edge.bicolor_facets) if matching else ""})
    _write_csv(out / "2KAI_incidence_audit.csv", inc)

    # Group-pattern curvature balance.
    pattern_rows=[]
    for pattern in ("AAB","ABB"):
        subset=[r for r in selected_final if r["group_pattern"]==pattern]
        s=summary_for(subset)
        s["mean_edge_length_A"] = s["L"]/s["N_edges"] if s["N_edges"] else math.nan
        s["mean_abs_beta_deg"] = float(np.mean([abs(edge_obj[r["triangle"]].beta_degrees) for r in subset])) if subset else math.nan
        all_s=summary_for(selected_final)
        s["fraction_total_L"] = s["L"]/all_s["L"] if all_s["L"] else math.nan
        s["fraction_total_abs_C"] = s["C_unsigned"]/all_s["C_unsigned"] if all_s["C_unsigned"] else math.nan
        pattern_rows.append({"pattern":pattern,**s})
    _write_csv(out / "2KAI_pattern_curvature.csv", pattern_rows)

    geo_shared_tri = {g["triangle"] for g in geom_rows if g["classification"] == "INTERIOR" and g["triangle"]}
    comb_shared_tri = {tuple(r["triangle"]) for r in selected_final}
    geometric_only = geo_shared_tri - comb_shared_tri
    combinatorial_only = comb_shared_tri - geo_shared_tri
    report = _report(result, rows, selected_final, tri_in_alpha, discrepant, variants, condition_rows,
                     geom_rows, polygon_issues, inc, pattern_rows, geometric_only, combinatorial_only)
    (out / "2KAI_curvature_forensic_report.txt").write_text(report, encoding="utf-8")
    (out / "2KAI_sign_convention_audit.txt").write_text(
        "Cazals 2006 states: ABB (one A, two B) positive; AAB negative. The supplied 2KAI chain/group assignment is A=kallikrein chains A+B and B=BPTI chain I. The paper names partners as A/B generically but does not map those labels to kallikrein/BPTI for 2KAI: NOT SPECIFIED. Reversing labels reverses beta and s_H sign only: -7.634237 deg becomes +7.634237 deg, while magnitude remains 7.634237 deg. This does not explain the +17 deg magnitude.\n",
        encoding="utf-8")
    return {"rows":rows,"geometric_edges":geom_rows,"discrepant":discrepant,"variants":variants,
            "condition_rows":condition_rows,"polygon_issues":polygon_issues,"incidence":inc}


def _write_csv(path, rows, fields=None):
    rows=list(rows); path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else []
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction="ignore"); writer.writeheader()
        for row in rows:
            writer.writerow({k:(repr(v) if isinstance(v,(tuple,list,dict)) else v) for k,v in row.items()})


def _report(result, rows, current, tri_alpha, discrepant, variants, condition, geom, issues, incidence, patterns,
            geometric_only=(), combinatorial_only=()):
    classes={}
    for row in rows: classes[row["classification"]]=classes.get(row["classification"],0)+1
    ratios=[x["m_over_r"] for x in condition if x["m_over_r"] is not None]
    closest_condition=min(condition, key=lambda row: abs(row["m_over_r"]-result.condition_beta_m)
                          if row["m_over_r"] is not None else float("inf")) if condition else None
    diffs=[x["length_difference_A"] for x in incidence if x["length_difference_A"] is not None]
    geo_shared=sum(x["classification"]=="INTERIOR" for x in geom)
    comb_tri={r["triangle"] for r in current}
    return ("2KAI Cazals power-interface edge forensics\n\n"
        "Primary source: Cazals et al. (2006), Protein Science 15:2082-2092, doi:10.1110/ps.062245906.\n\n"
        "SOURCE STATEMENTS\n"
        "EXPLICITLY STATED: alpha=0 selects bicolor facets dual to bicolor regular edges in K(0); condition beta rejects facets using m/r and M=5; facets form the interface. The beta magnitude is the singleton angle; ABB is positive, AAB negative. Equation 6 sums beta*l over all 'interior edges', divided by total edge length.\n"
        "STRONGLY IMPLIED: an interior edge is a finite power edge shared by two facets of the selected interface; a boundary edge with only one selected facet is excluded. The selected facets determine the surface first, then edge curvature is accumulated on that final interface.\n"
        "NOT SPECIFIED: a requirement that the regular triangle itself lie in K(0); mapping physical 2KAI proteins to A/B; restriction of Equation 6 to significant components.\n\n"
        f"EDGE POPULATION\nall AAB/ABB triangles: {sum(r['group_pattern'] in ('AAB','ABB') for r in rows)}; finite regular duals: {sum(r['geometry_status']=='finite' for r in rows)}; unbounded dual edges: {sum(r['geometry_status']=='unbounded_boundary_power_edge' for r in rows)}.\n"
        f"all 1713 regular-triangle records classified: {classes}\n"
        f"Selected polygon complex: {sum(x['classification']=='INTERIOR' for x in geom)} interior segments, {sum(x['classification']=='BOUNDARY' for x in geom)} boundary segments, {sum(x['classification']=='NONMANIFOLD' for x in geom)} nonmanifold segments.\n"
        f"current shared finite final-facet edges: {len(current)}; triangle also in K(0): {len(tri_alpha)}; discrepant triangles: {len(discrepant)}.\n"
        f"polygon segments: {len(geom)}; polygonal interior segments: {geo_shared}; polygon-construction issues: {len(issues)}.\n"
        f"polygon incidence length audit rows: {len(incidence)}; max |length difference|={max(diffs,default=float('nan')):.9g} A; RMS={float(np.sqrt(np.mean(np.square(diffs)))) if diffs else float('nan'):.9g} A.\n"
        f"Expected shared facet pairs={sum(1 for r in incidence if r['expected_shared'])}; geometrically confirmed={sum(1 for r in incidence if r['expected_shared'] and r['geometrically_shared'])}; combinatorial-only={len(combinatorial_only)}; geometric-only={len(geometric_only)}; mismatched lengths >1e-9 A={sum(x>1e-9 for x in diffs)}.\n\n"
        f"ALPHA TRIANGLE DISCREPANCY\nN={len(discrepant)}; see 2KAI_13_alpha_triangle_discrepant_edges.csv. L={sum(float(r['length_A']) for r in discrepant):.9g} A; positive/negative edge counts={sum(float(r['beta_rad'])>0 for r in discrepant)}/{sum(float(r['beta_rad'])<0 for r in discrepant)}; C+={sum(float(r['l_beta']) for r in discrepant if float(r['beta_rad'])>0):.9g}; C-={sum(float(r['l_beta']) for r in discrepant if float(r['beta_rad'])<0):.9g}; C_signed={sum(float(r['l_beta']) for r in discrepant):.9g}; C_unsigned={sum(abs(float(r['l_beta'])) for r in discrepant):.9g} A rad. Keeping these 13 yields -7.634237 deg; dropping them yields -6.561559 deg, a +1.072678 degree shift toward zero (unsigned 36.495614 -> 35.807273 deg). These edges have both selected bicolor facets, are finite, but the regular triangle filtration exceeds zero. Requiring triangle K(0) membership is not stated by Cazals and would impose extra clipping/selection.\n\n"
        f"CONDITION B\nalpha-selected bicolor facets audited: {len(condition)}; rejected={sum(not x['accepted'] for x in condition)}; m/r min={min(ratios,default=float('nan')):.9g}; max={max(ratios,default=float('nan')):.9g}; closest to M={closest_condition['m_over_r'] if closest_condition else float('nan'):.9g} for {closest_condition['atom_i_label'] if closest_condition else 'n/a'}--{closest_condition['atom_j_label'] if closest_condition else 'n/a'}; M={result.condition_beta_m:g}. Raw maximum squared power is negative for {sum(x['max_power_radius_squared_A2'] is not None and x['max_power_radius_squared_A2'] < 0 for x in condition)} facets; the current implementation reports m=0 for those nonpositive squared radii. Zero rejections are supported by the exported per-facet ratios, not inferred from summary alone.\n\n"
        "EDGE-SET VARIANTS\n"+"\n".join(f"{v['method']}: N={v['N_edges']}, L={v['L']:.9g}, C_signed={v['C_signed']:.9g}, sH={v['s_H_signed_deg']:.9g} deg [{v['interpretation_status']}]; {v['edge_definition']}. {v['source_justification']}" for v in variants)+"\n\n"
        "GROUP-PATTERN CURVATURE\n"+"\n".join(f"{p['pattern']}: N={p['N_edges']}, L={p['L']:.9g}, mean length={p['mean_edge_length_A']:.9g}, C_signed={p['C_signed']:.9g}, C_unsigned={p['C_unsigned']:.9g}, mean |beta|={p['mean_abs_beta_deg']:.9g} deg, L fraction={p['fraction_total_L']:.6%}, |C| fraction={p['fraction_total_abs_C']:.6%}" for p in patterns)+"\n\n"
        "CONCAVITY CHECK\nWith the current labeling, 522/864 edges (60.42%) are AAB and 342/864 (39.58%) ABB; the negative AAB contribution carries 60.46% of |C| and 57.88% of length. Since AAB has the singleton in group B (the inhibitor here), this pattern is consistent with inhibitor-side concavity only if singleton-side geometry is taken as that concavity orientation. The 2006 text reports concavity toward inhibitor but does not explicitly equate that phrase to the algebraic beta sign; treat the linkage as an inference.\n\n"
        "LABEL FORENSICS\nCazals' sign rule gives current -7.634237 deg for A=protease/B=inhibitor; reversing labels gives +7.634237 deg, changing sign only. The paper's reported +17 deg remains a magnitude discrepancy. The supplied original 2KAI COMPND records assign chains A+B to kallikrein and I to BPTI; the 2006 article uses generic partners A and B and does not name their mapping for this complex.\n\n"
        "HISTORICAL IMPLEMENTATION SEARCH\nThe original 2006 paper identifies CGAL Alpha_shape_3 and points to Intervor; it does not publish source code or spell out a separate K(0) requirement for dual triangles. The 2009 Cazals-Loriot Intervor technical report (INRIA RR-7069 / HAL inria-00426041) describes interface tiles and explicitly records boundary loops before/after component merging, but the accessible report/searchable material does not expose the Equation 6 accumulation code. No historical source-level rule resolving the 864-vs-851 choice or physical A/B mapping was verified.\n\n"
        "INTERPRETATION\nThe current 864 records are the analyzer's finite regular edges whose two bicolor incident facets passed alpha and condition beta. The independent polygon audit confirms all 864 by exact dual-vertex endpoint matching; the polygonal interface also has 111 boundary segments. No triangle K(0) condition is added. The remaining magnitude discrepancy is not resolved here. A different curvature population/normalization remains a hypothesis, not an established explanation.\n")
