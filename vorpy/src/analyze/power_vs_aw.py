"""Transparent side-by-side reporting for Power and AW interface results."""
from __future__ import annotations

import csv
from collections import defaultdict, deque
from pathlib import Path
import numpy as np


def compare_power_and_aw(power_run, aw_result, output_dir, *, aw_lys15_ids=(), input_notes=(), group_labels=None):
    """Write comparable interface metrics while keeping curvature terms separate.

    ``power_run`` is the mapping returned by ``run_power_interface_pdb``.
    No combined Power/AW curvature value is emitted: the Power quantity is
    discrete turning on planar facets, while AW has curved surfaces and its
    analyzer intentionally withholds a combined signed curvature measure.
    """
    output_dir=Path(output_dir); output_dir.mkdir(parents=True,exist_ok=True)
    power=power_run["result"]; prepared=power_run["prepared"]
    p_atoms_a,p_atoms_b=_power_interface_atoms(power)
    power_components=_power_components(power_run)
    p_coords=prepared.points[sorted(p_atoms_a|p_atoms_b)] if p_atoms_a or p_atoms_b else np.empty((0,3))
    p_centroid=_centroid(p_coords)

    aw_surfaces=aw_result.selected_surfaces
    aw_atom_ids={i for s in aw_surfaces for i in s.generator_ids}
    aw_surface_candidates = list(getattr(aw_result, "surfaces", ()))
    aw_excluded = [surface for surface in aw_surface_candidates if not surface.included]
    aw_excluded_reasons = defaultdict(int)
    for surface in aw_excluded:
        aw_excluded_reasons[surface.reason or "unspecified"] += 1
    aw_components=_aw_surface_components(aw_result)
    network=aw_result.network
    try:
        aw_coords=np.asarray(network.balls.loc[sorted(aw_atom_ids),"loc"].tolist(),dtype=float)
    except Exception:
        aw_coords=np.empty((0,3))
    aw_centroid=_centroid(aw_coords)
    edge=aw_result.edge_curvature
    p_edge=power.edge_totals
    aw_supported=bool(aw_surfaces)
    rows=[]
    def add(model,metric,value,unit,notes=""):
        rows.append({"model":model,"metric":metric,"value":value,"unit":unit,"notes":notes})
    add("Power","retained_atoms",len(prepared.atoms),"atoms","Cazals PDB preparation")
    add("AW","retained_atoms",len(network.balls),"atoms","AW network input generators")
    if group_labels is not None:
        add("Power","group_A_selection",",".join(sorted(group_labels[0])),"chains")
        add("Power","group_B_selection",",".join(sorted(group_labels[1])),"chains")
        add("AW","group_A_selection",",".join(sorted(group_labels[0])),"chains")
        add("AW","group_B_selection",",".join(sorted(group_labels[1])),"chains")
    add("Power","interface_atoms_A",len(p_atoms_a),"atoms")
    add("Power","interface_atoms_B",len(p_atoms_b),"atoms")
    add("Power","interface_atoms_total",len(p_atoms_a|p_atoms_b),"atoms")
    add("AW","interface_atoms_total",len(aw_atom_ids),"atoms","Unique generators on selected supported AW surfaces")
    power_area=_power_area(power_run)
    add("Power","interface_area",power_area,"A^2","Sum of bounded selected power-facet polygon areas")
    add("AW","interface_area",aw_result.interface_area,"A^2","Selected supported AW surfaces")
    add("AW","interface_surface_candidates",len(aw_surface_candidates),"surfaces")
    add("AW","excluded_interface_surfaces",len(aw_excluded),"surfaces",repr(dict(aw_excluded_reasons)))
    add("Power","connected_components",len(power_components),"components")
    add("AW","connected_components",len(aw_components),"components","Surface adjacency through selected supported AW edges")
    add("Power","significant_components",sum(c["significant"] for c in power_components),"components","7.5% interface-area threshold")
    add("AW","significant_components",sum(c["significant"] for c in aw_components),"components","7.5% interface-area threshold")
    add("Power","interface_centroid",p_centroid,"A","Centroid of unique interface atom centers")
    add("AW","interface_centroid",aw_centroid,"A","Centroid of unique generators on selected AW surfaces")
    add("Power","planar_facets",len(power.final_facets),"facets")
    add("Power","C_edge_raw_signed",p_edge["signed"],"A rad","Cazals raw sum(l beta)")
    add("Power","C_edge_H_signed",p_edge["signed"]/2.,"A","Conventional half-scaled edge measure")
    add("AW","curved_surfaces",len(aw_surfaces),"surfaces")
    add("AW","C_edge_raw_signed",edge["signed"] if aw_supported else None,"A",
        "AW integral of beta ds" if aw_supported else "Unavailable: no selected supported AW surfaces")
    add("AW","C_edge_H_signed",edge["signed"]/2. if aw_supported else None,"A",
        "Half-scaled edge contribution for H=(k1+k2)/2" if aw_supported else "Unavailable: no selected supported AW surfaces")
    add("AW","C_surface_H",aw_result.surface_curvature if aw_supported else None,"A",
        "Existing AW integral of H dA; H=(k1+k2)/2" if aw_supported else "Unavailable: no selected supported AW surfaces")
    add("AW","C_total_H",None,"A","Withheld: edge orientation and surface-normal sign conventions have not been reconciled")
    lys_ids=set(aw_lys15_ids)
    lys_power = sum(any(prepared.atoms[i].chain=="I" and prepared.atoms[i].residue_name=="LYS" and prepared.atoms[i].residue_number==15 for i in f.generator_ids) for f in power.final_facets)
    lys_aw_surfaces = [s for s in aw_surfaces if lys_ids.intersection(s.generator_ids)]
    lys_aw = len(lys_aw_surfaces)
    lys_aw_area = sum(s.area or 0.0 for s in lys_aw_surfaces)
    lys_aw_surface_h = sum(s.integrated_mean_curvature or 0.0 for s in lys_aw_surfaces)
    lys_aw_edges = [e for e in aw_result.selected_edges if lys_ids.intersection(e.generator_ids)]
    lys_power_measurements = [m for m in getattr(power, "measurements", ())
                              if lys_ids.intersection(m.triangle)]
    lys_power_facets = [r for r in power_run.get("facet_records", ())
                        if "LYS15" in str(r.get("residue_i", ""))
                        or "LYS15" in str(r.get("residue_j", ""))]
    lys_power_area = sum(float(r.get("area_A2") or 0.0) for r in lys_power_facets)
    add("Power","Lys15_facets",lys_power,"facets","BPTI chain I Lys15 pocket localization")
    add("Power","Lys15_interface_area",lys_power_area,"A^2","Final selected power facets incident to a BPTI chain I Lys15 atom")
    add("Power","Lys15_curvature_edges",len(lys_power_measurements),"edges","Selected finite power interface edges incident to a Lys15 atom")
    add("Power","Lys15_C_signed",sum(m.length_times_beta for m in lys_power_measurements),"A rad")
    add("Power","Lys15_C_unsigned",sum(m.length_times_abs_beta for m in lys_power_measurements),"A rad")
    add("AW","Lys15_surfaces",lys_aw,"surfaces","Selected supported AW surfaces incident to BPTI chain I Lys15")
    add("AW","Lys15_interface_area",lys_aw_area,"A^2")
    add("AW","Lys15_curvature_edges",len(lys_aw_edges),"edges","Included AW analytic curvature edges incident to a Lys15 atom")
    add("AW","Lys15_C_edge_signed",sum(e.signed_contribution or 0.0 for e in lys_aw_edges),"A")
    add("AW","Lys15_C_edge_unsigned",sum(e.unsigned_contribution or 0.0 for e in lys_aw_edges),"A")
    add("AW","Lys15_C_surface_H",lys_aw_surface_h,"A","Existing integral of H dA over selected Lys15-associated AW surfaces")
    for note in input_notes:
        add("Input","note",note,"text")
    _write_csv(output_dir/"power_vs_aw_metrics.csv",rows)
    input_note_text = "Input note: " + " ".join(input_notes) + "\n\n" if input_notes else ""
    groups_text = (f"Selections: group A chains={sorted(group_labels[0])}; group B chains={sorted(group_labels[1])}.\n"
                   if group_labels is not None else "")
    summary=("Power vs AW Interface Comparison\n"
        f"Power: {len(prepared.atoms)} Cazals/Chothia generators; groups A/B={len(prepared.group_a)}/{len(prepared.group_b)}; probe={prepared.probe_radius:.2f} A; alpha={power.alpha:g}; M={power.condition_beta_m:g}.\n"
        f"AW: {len(network.balls)} generators; selection={aw_result.selection_mode}; groups A/B={len(aw_result.group_a)}/{len(aw_result.group_b)}.\n\n"
        + groups_text
        + input_note_text
        + f"Power interface atoms={len(p_atoms_a)}/{len(p_atoms_b)}; facets={len(power.final_facets)}; area={power_area:.9g} A^2; components={len(power_components)}; significant={sum(c['significant'] for c in power_components)}.\n"
        f"AW interface atoms={len(aw_atom_ids)}; selected surfaces={len(aw_surfaces)}/{len(aw_surface_candidates)} candidates; excluded={len(aw_excluded)} {dict(aw_excluded_reasons)}; area={aw_result.interface_area:.9g} A^2; components={len(aw_components)}; significant={sum(c['significant'] for c in aw_components)}.\n"
        f"Interface-center centroid distance={_format_distance(_distance(p_centroid,aw_centroid))} A (localization diagnostic only).\n\n"
        f"Power curvature: C_edge_raw={p_edge['signed']:.9g} A rad; C_edge_H=C_raw/2={p_edge['signed']/2.:.9g} A.\n"
        + (f"AW curvature: C_edge_raw={edge['signed']:.9g} A; C_edge_H=C_raw/2={edge['signed']/2.:.9g} A; C_surface_H={aw_result.surface_curvature:.9g} A.\n"
           if aw_supported else "AW curvature: unavailable; no selected supported AW surfaces.\n")
        +
        "C_total_H is withheld because the edge orientation and surface-normal signs are not yet reconciled. Power s_H and AW edge-normalized beta statistics are not treated as interchangeable.\n"
        f"Lys15 pocket diagnostic: Power facets={lys_power}, area={lys_power_area:.9g} A^2, curvature edges={len(lys_power_measurements)}, C_signed={sum(m.length_times_beta for m in lys_power_measurements):.9g} A rad; "
        f"AW selected surfaces={lys_aw}, area={lys_aw_area:.9g} A^2, curvature edges={len(lys_aw_edges)}, "
        f"C_edge_signed={sum(e.signed_contribution or 0.0 for e in lys_aw_edges):.9g} A, "
        f"C_surface_H={lys_aw_surface_h:.9g} A.\n")
    (output_dir/"power_vs_aw_summary.txt").write_text(summary,encoding="utf-8")
    return {"summary":summary,"metrics":rows}


def _power_interface_atoms(result):
    a=set();b=set()
    for facet in result.final_facets:
        for i in facet.generator_ids:
            (a if i in result.group_a else b).add(i)
    return a,b


def _power_components(power_run):
    from vorpy.src.analyze.cazals_validation import interface_components
    return interface_components(power_run["result"],power_run.get("facet_area_by_facet",{}))


def _power_area(power_run):
    return float(sum(float(row["area_A2"] or 0.) for row in power_run.get("facet_records",())))


def _aw_surface_components(result):
    surfaces={s.surface_id:s for s in result.selected_surfaces}
    adjacency={sid:set() for sid in surfaces}
    for edge in result.selected_edges:
        u,v=edge.surface_ids
        if u in surfaces and v in surfaces:
            adjacency[u].add(v);adjacency[v].add(u)
    unseen=set(surfaces); out=[]
    while unseen:
        seed=min(unseen);unseen.remove(seed);queue=deque([seed]);ids=set()
        while queue:
            node=queue.popleft();ids.add(node)
            for nxt in adjacency[node]&unseen:unseen.remove(nxt);queue.append(nxt)
        area=sum(surfaces[i].area or 0. for i in ids);out.append({"area":area,"surfaces":ids})
    total=sum(c["area"] for c in out)
    for c in out:c["fraction"]=c["area"]/total if total else 0.;c["significant"]=c["fraction"]>=.075
    return out


def _centroid(points):
    return tuple(map(float,np.mean(points,axis=0))) if len(points) else None


def _distance(first,second):
    if first is None or second is None:return None
    return float(np.linalg.norm(np.asarray(first)-np.asarray(second)))


def _format_distance(value):
    return "n/a" if value is None else f"{value:.9g}"


def _write_csv(path,rows):
    rows=list(rows);fields=list(rows[0]) if rows else []
    with Path(path).open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)
