"""Read-only audit of molecular filtering stages on dual connections.

Power membership and curvature are supplied by the frozen Cazals analyzer.
This module only labels/exports its records and computes center distances for
diagnostics; it never applies a geometric distance cutoff.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
from vorpy.src.analyze.power_interface_cli import run_power_interface_pdb


def _write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def classify_power_edges(result, atoms):
    """Return one exhaustive row per full regular edge, sourced from result."""
    from vorpy.src.analyze.cazals_validation import power_facet_polygon_area

    groups_a, groups_b = set(result.group_a), set(result.group_b)
    facets = {tuple(record.generator_ids): record for record in result.facets}
    rows = []
    for pair in sorted(result.full_simplices[1]):
        i, j = pair
        ai, aj = atoms[i], atoms[j]
        bicolor = (i in groups_a and j in groups_b) or (j in groups_a and i in groups_b)
        facet = facets.get(pair)
        polygon = power_facet_polygon_area(result, facet) if facet is not None else {}
        filtration = result.filtration.get(pair)
        alpha0 = pair in result.alpha_simplices[1]
        stage_path = ["full_regular"]
        if alpha0:
            stage_path.append("alpha0")
        if alpha0 and bicolor:
            stage_path.append("bicolor_alpha0")
        if not alpha0:
            stage = "rejected_alpha0"
            stage_path.append("rejected_alpha0")
            reason = "regular edge filtration exceeds alpha threshold"
        elif not bicolor:
            stage = "rejected_not_bicolor"
            stage_path.append("rejected_not_bicolor")
            reason = "alpha0 regular edge does not connect the requested groups"
        elif facet is None:
            raise RuntimeError(f"Alpha0 bicolor regular edge {pair} lacks production facet record")
        elif facet.condition_beta_accepted is False:
            stage = "rejected_M5"
            stage_path.append("rejected_M5")
            reason = facet.reason or facet.condition_beta_status
        elif facet.condition_beta_accepted is None:
            stage = "unresolved_M5"
            stage_path.append("unresolved_M5")
            reason = facet.reason or facet.condition_beta_status
        else:
            stage = "cazals_final"
            stage_path.extend(("condition_b_pass", "cazals_final"))
            reason = "selected by frozen alpha0 + bicolor + condition-beta implementation"
        delta = result.points[i] - result.points[j]
        rows.append({
            "generator_i": i, "generator_j": j,
            "atom_i_stable_id": _stable_id(ai), "atom_j_stable_id": _stable_id(aj),
            "atom_i_chain": ai.chain, "atom_i_residue": ai.residue_name,
            "atom_i_resnum": ai.residue_number, "atom_i_icode": ai.insertion_code,
            "atom_i_name": ai.atom_name, "atom_j_chain": aj.chain,
            "atom_j_residue": aj.residue_name, "atom_j_resnum": aj.residue_number,
            "atom_j_icode": aj.insertion_code, "atom_j_name": aj.atom_name,
            "group_i": ai.group, "group_j": aj.group,
            "center_distance_A": float(np.linalg.norm(delta)),
            "base_radius_i_A": ai.base_radius, "base_radius_j_A": aj.base_radius,
            "expanded_radius_i_A": ai.expanded_radius,
            "expanded_radius_j_A": aj.expanded_radius,
            "power_weight_i_A2": ai.weight, "power_weight_j_A2": aj.weight,
            "full_regular": True,
            "stage_path": ";".join(stage_path),
            "alpha0_membership": bool(alpha0), "alpha_filtration_A2": filtration,
            "bicolor": bool(bicolor),
            "condition_b_applicable": bool(bicolor and alpha0),
            "condition_b_status": facet.condition_beta_status if facet else "not_evaluated",
            "condition_b_r_A": facet.smaller_expanded_radius if facet else None,
            "condition_b_m_A": facet.largest_orthogonal_ball_radius if facet else None,
            "condition_b_m_over_r": facet.m_over_r if facet else None,
            "condition_b_accepted": facet.condition_beta_accepted if facet else None,
            "final_cazals_membership": bool(facet.final_selected) if facet else False,
            "facet_area_A2": polygon.get("area_A2"),
            "facet_area_status": polygon.get("status", "not_a_bicolor_facet"),
            "facet_bounded_status": (facet.condition_beta_status if facet else "not_evaluated"),
            "stage": stage, "rejection_reason": reason,
        })
    if len(rows) != len(result.full_simplices[1]):
        raise AssertionError("A full-regular edge disappeared from the diagnostic export")
    return rows


def _stable_id(atom):
    return f"{atom.chain}:{atom.residue_number}{atom.insertion_code}:{atom.residue_name}:{atom.atom_name}:{atom.pdb_serial}"


def summarize_distances(rows):
    stages = {
        "full_regular": rows,
        "rejected_alpha0": [r for r in rows if r["stage"] == "rejected_alpha0"],
        "retained_alpha0": [r for r in rows if r["alpha0_membership"]],
        "rejected_M5": [r for r in rows if r["stage"] == "rejected_M5"],
        "final_cazals": [r for r in rows if r["final_cazals_membership"]],
    }
    summaries = []
    for name, selected in stages.items():
        values = np.asarray([row["center_distance_A"] for row in selected], dtype=float)
        summaries.append({
            "stage": name, "count": len(values),
            "min_A": float(np.min(values)) if len(values) else None,
            "median_A": float(np.median(values)) if len(values) else None,
            "mean_A": float(np.mean(values)) if len(values) else None,
            "p95_A": float(np.percentile(values, 95)) if len(values) else None,
            "max_A": float(np.max(values)) if len(values) else None,
        })
    return summaries


def build_probe_control(prepared, *, alpha=0.0, m=5.0):
    """Run the same frozen builder on base and probe-expanded radii."""
    results = {}
    for label, radii in (("base_radii", prepared.base_radii),
                         ("base_plus_1p4A", prepared.expanded_radii)):
        results[label] = analyze_power_interface_curvature(
            prepared.points, radii, prepared.group_a, prepared.group_b,
            alpha=alpha, condition_beta_m=m,
        )
    rows = []
    edge_sets = {key: value.full_simplices[1] for key, value in results.items()}
    union = set.union(*map(set, edge_sets.values()))
    for pair in sorted(union):
        rows.append({
            "generator_i": pair[0], "generator_j": pair[1],
            "center_distance_A": float(np.linalg.norm(prepared.points[pair[0]] - prepared.points[pair[1]])),
            "in_base_full_regular": pair in edge_sets["base_radii"],
            "in_probe_full_regular": pair in edge_sets["base_plus_1p4A"],
            "base_alpha0": pair in results["base_radii"].alpha_simplices[1],
            "probe_alpha0": pair in results["base_plus_1p4A"].alpha_simplices[1],
            "base_bicolor": _is_bicolor(pair, results["base_radii"]),
            "probe_bicolor": _is_bicolor(pair, results["base_plus_1p4A"]),
            "base_final_M5": any(f.generator_ids == pair and f.final_selected for f in results["base_radii"].facets),
            "probe_final_M5": any(f.generator_ids == pair and f.final_selected for f in results["base_plus_1p4A"].facets),
        })
    return rows, {key: {
        "regular_edges": len(value.full_simplices[1]),
        "alpha0_edges": len(value.alpha_simplices[1]),
        "bicolor_alpha0_facets": len(value.alpha_zero_facets),
        "final_M5_facets": len(value.final_facets),
    } for key, value in results.items()}


def _is_bicolor(pair, result):
    a, b = result.group_a, result.group_b
    return (pair[0] in a and pair[1] in b) or (pair[1] in a and pair[0] in b)


def export_power_visualization(rows, atoms, output_dir, structure_path):
    """Write stage-separated PyMOL CGO objects using VorPy's draw_line mesh builder."""
    from vorpy.src.geometry.visualization.export import _edge_mesh
    from vorpy.src.output.mesh import write_mesh

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    xyz = {atom.index: np.asarray((atom.x, atom.y, atom.z), dtype=float) for atom in atoms}
    categories = {
        "2KAI_full_regular": rows,
        "2KAI_alpha0": [row for row in rows if row["alpha0_membership"]],
        "2KAI_cazals_final": [row for row in rows if row["final_cazals_membership"]],
        "2KAI_rejected_alpha0": [row for row in rows if row["stage"] == "rejected_alpha0"],
        "2KAI_rejected_M5": [row for row in rows if row["stage"] == "rejected_M5"],
    }
    colors = {"2KAI_full_regular": (0.50, 0.55, 0.62), "2KAI_alpha0": (0.95, 0.65, 0.1),
              "2KAI_cazals_final": (0.0, 0.85, 0.35), "2KAI_rejected_alpha0": (0.95, 0.2, 0.2),
              "2KAI_rejected_M5": (0.75, 0.1, 0.85)}
    files, counts = {}, {}
    for name, selected in categories.items():
        mesh = _edge_mesh([(row["generator_i"], row["generator_j"]) for row in selected], xyz, 0.045)
        if mesh is not None:
            filename = f"{name}.off"
            write_mesh(mesh, filename, directory=destination)
            files[name] = filename
        counts[name] = len(selected)
    # Use the same OFF-to-CGO loading approach as the shared visualization exporter.
    lines = ["from pymol import cmd", "from pymol.cgo import BEGIN, END, TRIANGLES, COLOR, VERTEX",
             "import os", f"cmd.load({str(Path(structure_path).resolve())!r}, '2KAI_atoms')",
             "cmd.hide('everything', '2KAI_atoms')", "cmd.show('spheres', '2KAI_atoms')",
             "cmd.set('sphere_scale', 0.22, '2KAI_atoms')", "cmd.color('gray70', '2KAI_atoms')",
             "def _load_off(path, name, color):", "    if not os.path.exists(path): return",
             "    with open(path, encoding='utf-8') as f: data=f.readlines()",
             "    nv,nf,_=map(int,data[1].split()); pts=[tuple(map(float,x.split())) for x in data[2:2+nv]]",
             "    obj=[BEGIN,TRIANGLES,COLOR,*color]",
             "    for line in data[2+nv:2+nv+nf]:",
             "        fields=line.split(); ids=list(map(int,fields[1:1+int(fields[0])]))",
             "        for k in range(1,len(ids)-1): obj.extend([VERTEX,*pts[ids[0]],VERTEX,*pts[ids[k]],VERTEX,*pts[ids[k+1]]])",
             "    obj.append(END); cmd.load_cgo(obj,name)"]
    for name, filename in files.items():
        lines.append(f"_load_off({str((destination / filename).resolve())!r}, {name!r}, {list(colors[name])!r})")
    pml = destination / "2KAI_dual_restriction_stages.pml"
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    metadata = {"structure": str(Path(structure_path).resolve()), "counts": counts,
                "objects": files, "pymol_script": str(pml),
                "distance_cutoff_used": False,
                "line_geometry": "VorPy draw_line-derived OFF meshes loaded as named PyMOL CGO objects"}
    (destination / "visualization_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def run_2kai(pdb, output_dir, *, aw_network_path=None):
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Production-compatible frozen result. Exports and summaries only read it.
    run = run_power_interface_pdb(pdb, {"A", "B"}, {"I"}, output / "production_source",
                                  probe_radius=1.4, alpha=0.0, condition_beta_m=5.0)
    result, prepared = run["result"], run["prepared"]
    rows = classify_power_edges(result, prepared.atoms)
    _write_csv(output / "2KAI_power_connections.csv", rows)
    distance_rows = summarize_distances(rows)
    _write_csv(output / "2KAI_power_distance_summary.csv", distance_rows)
    longest = sorted(rows, key=lambda row: row["center_distance_A"], reverse=True)[:25]
    _write_csv(output / "2KAI_longest_full_regular_connections.csv", longest)
    longest_stage_rows = []
    stage_masks = {
        "full_regular": lambda row: True,
        "rejected_alpha0": lambda row: row["stage"] == "rejected_alpha0",
        "retained_alpha0": lambda row: row["alpha0_membership"],
        "rejected_M5": lambda row: row["stage"] == "rejected_M5",
        "final_cazals": lambda row: row["final_cazals_membership"],
    }
    for stage_name, predicate in stage_masks.items():
        selected = sorted((row for row in rows if predicate(row)),
                          key=lambda row: row["center_distance_A"], reverse=True)[:25]
        longest_stage_rows.extend({**row, "rank_in_stage": rank, "distribution_stage": stage_name}
                                  for rank, row in enumerate(selected, 1))
    _write_csv(output / "2KAI_longest_connections_by_stage.csv", longest_stage_rows)
    control_rows, control_summary = build_probe_control(prepared)
    _write_csv(output / "2KAI_probe_radius_control.csv", control_rows)
    (output / "2KAI_probe_radius_control_summary.json").write_text(json.dumps(control_summary, indent=2) + "\n", encoding="utf-8")
    # Frozen-selection/curvature immutability check against exported rows.
    final_pairs = {tuple(sorted((row["generator_i"], row["generator_j"]))) for row in rows if row["final_cazals_membership"]}
    production_pairs = {tuple(f.generator_ids) for f in result.final_facets}
    if final_pairs != production_pairs or len(result.included_edges) != 864:
        raise AssertionError("Diagnostic export does not preserve the frozen 2KAI selection/edge count")
    expected = {"length": 953.5241609499, "signed": -127.0500153782, "unsigned": 607.3649734659}
    totals = result.edge_totals
    for key, value in expected.items():
        if not math.isclose(totals[key], value, abs_tol=2e-6):
            raise AssertionError(f"Frozen Cazals {key} changed: {totals[key]} != {value}")
    visualization = export_power_visualization(rows, prepared.atoms, output / "visualization", pdb)
    aw_summary = None
    if aw_network_path:
        aw_summary = export_aw_connections(aw_network_path, output, structure_path=pdb)
    summary = {
        "full_regular_edges": len(result.full_simplices[1]),
        "alpha0_edges": len(result.alpha_simplices[1]),
        "bicolor_regular_edges": sum(row["bicolor"] for row in rows),
        "bicolor_alpha0": len(result.alpha_zero_facets),
        "M5_rejected": sum(f.alpha_zero_selected and f.condition_beta_accepted is False for f in result.facets),
        "final_cazals": len(result.final_facets), "interface_edges": len(result.included_edges),
        "curvature": totals, "distance_summary": distance_rows,
        "longest_25": "2KAI_longest_full_regular_connections.csv",
        "probe_control": control_summary, "visualization": visualization,
        "aw_connections": aw_summary,
    }
    (output / "2KAI_dual_restriction_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")
    (output / "2KAI_dual_restriction_summary.txt").write_text(_format_summary(summary), encoding="utf-8")
    return summary


def export_aw_connections(network_path, output_dir, *, alpha=1.4, structure_path=None):
    """Export full and exact supported AW pair restriction from AW births.

    Pair births are taken from the experimental AW restricted-cell filtration,
    which tests/minimizes on the actual clipped network surface. Unsupported
    births remain unresolved; this function never substitutes Power alpha.
    """
    from vorpy.src.geometry.aw_alpha.births import surface_birth
    from vorpy.src.io import load_network
    from vorpy.src.geometry.visualization.export import _edge_mesh, _surface_mesh
    from vorpy.src.output.mesh import write_mesh
    from vorpy.src.analyze.apollonius import PrimalFeatureRef

    network = load_network(network_path)
    stable_atoms = _aw_stable_atom_map(Path(network_path).with_name("radii_audit.csv"))
    # Restriction is an interface question here, so only bicolor pair births
    # are evaluated. Every full Apollonius pair is still exported below.
    group_a_chains, group_b_chains = {"A", "B"}, {"I"}
    complete_cells = {int(i): bool(row.get("complete", True))
                      for i, row in network.balls.iterrows()}
    features_by_pair = {}
    for feature_id, row in network.surfs.iterrows():
        try:
            ids = tuple(sorted(int(value) for value in row.get("balls", ())))
        except (TypeError, ValueError):
            continue
        if len(ids) != 2 or ids[0] == ids[1] or any(i not in network.balls.index for i in ids):
            continue
        bounded = bool(row.get("edges", ())) and bool(row.get("verts", ()))
        cells_complete = all(complete_cells[i] for i in ids)
        feature = PrimalFeatureRef("surf", int(feature_id), ids, bounded, bounded,
                                   "supported" if cells_complete else "incomplete_network_geometry",
                                   generator_cells_complete=cells_complete)
        features_by_pair.setdefault(ids, []).append(feature)
    rows = []
    full_pairs, selected_pairs, rejected_pairs = [], [], []
    selected_surface_ids, unresolved_pairs, non_bicolor_pairs = [], [], []
    for pair, features in sorted(features_by_pair.items()):
        i, j = pair
        a, b = network.balls.loc[i], network.balls.loc[j]
        feature_supported = any(f.bounded and f.complete for f in features)
        chain_i = stable_atoms.get(i, "").split(":", 1)[0]
        chain_j = stable_atoms.get(j, "").split(":", 1)[0]
        bicolor = ((chain_i in group_a_chains and chain_j in group_b_chains)
                   or (chain_j in group_a_chains and chain_i in group_b_chains))
        birth, geometric_birth = None, None
        birth_diagnostic = {"birth_source": "not_evaluated_not_bicolor"}
        if bicolor:
            try:
                geometric_birth, birth_diagnostic = surface_birth(
                    network, SimpleNamespace(generator_ids=pair, primal_features=features), 1e-6
                )
                if geometric_birth is not None:
                    # This is the existing nerve face-maximum closure for a
                    # pair: max(pair birth, both formal generator births).
                    birth = max(float(geometric_birth), -float(a["rad"]), -float(b["rad"]))
            except (TypeError, ValueError, KeyError, IndexError, np.linalg.LinAlgError) as error:
                birth_diagnostic = {"birth_source": "unresolved", "reason": str(error)}
        supported_birth = bool(bicolor and birth is not None)
        selected = bool(supported_birth and birth <= float(alpha) + 1e-6)
        surface_area = _network_surface_area(network, features)
        if selected:
            selected_pairs.append(pair)
            selected_surface_ids.extend(f.feature_id for f in features)
        elif not bicolor:
            non_bicolor_pairs.append(pair)
        elif not supported_birth:
            unresolved_pairs.append(pair)
        else:
            rejected_pairs.append(pair)
        full_pairs.append(pair)
        rows.append({
            "generator_i": i, "generator_j": j,
            "atom_i": stable_atoms.get(i, _network_atom_id(a, i)),
            "atom_j": stable_atoms.get(j, _network_atom_id(b, j)),
            "center_distance_A": float(np.linalg.norm(np.asarray(a["loc"]) - np.asarray(b["loc"]))),
            "radius_i_A": float(a["rad"]), "radius_j_A": float(b["rad"]),
            "surface_ids": ";".join(str(f.feature_id) for f in features),
            "surface_count": len(features),
            "surface_finite_complete": any(f.bounded and f.complete for f in features),
            "surface_supported": feature_supported,
            "bicolor_interface_pair": bicolor,
            "geometric_birth_A": geometric_birth,
            "filtration_birth_A": birth,
            "birth_source": birth_diagnostic.get("birth_source", "unresolved"),
            "birth_reason": birth_diagnostic.get("reason", ""),
            "birth_supported": supported_birth,
            "alpha_native_A": float(alpha),
            "alpha_membership": selected,
            "restriction_stage": ("aw_alpha_selected" if selected else
                                  "not_evaluated_not_bicolor" if not bicolor else
                                  "unresolved_aw_birth" if not supported_birth else
                                  "rejected_aw_alpha"),
            "boundary_edge_count": sum(len(network.surfs.loc[f.feature_id].get("edges", ())) for f in features),
            "cell_i_complete": bool(a.get("complete", True)), "cell_j_complete": bool(b.get("complete", True)),
            "area_A2": surface_area,
            "area_status": "stored_scalar" if surface_area is not None
                          else "no_surface_area_scalar_in_network_record",
            "restriction_status": "full_AW_connection_unfiltered",
        })
    output_dir = Path(output_dir)
    _write_csv(output_dir / "2KAI_aw_full_connections.csv", rows)
    viz = output_dir / "visualization_aw"
    viz.mkdir(parents=True, exist_ok=True)
    coordinates = {int(i): np.asarray(row["loc"], dtype=float) for i, row in network.balls.iterrows()}
    layer_pairs = {"aw_full_dual": full_pairs, "aw_alpha_1p4": selected_pairs,
                   "aw_rejected_alpha": rejected_pairs, "aw_unresolved_birth": unresolved_pairs}
    layer_colors = {"aw_full_dual": (0.45, 0.5, 0.58), "aw_alpha_1p4": (0.0, 0.82, 0.4),
                    "aw_rejected_alpha": (0.95, 0.25, 0.15), "aw_unresolved_birth": (0.85, 0.55, 0.1)}
    layer_files = {}
    for name, pairs in layer_pairs.items():
        mesh = _edge_mesh(pairs, coordinates, 0.045)
        if mesh is not None:
            filename = f"{name}.off"
            write_mesh(mesh, filename, directory=viz)
            layer_files[name] = filename
    positions = [network.surfs.index.get_loc(index) for index in sorted(set(selected_surface_ids))
                 if index in network.surfs.index]
    selected_surface_mesh = _surface_mesh(network, positions)
    if selected_surface_mesh is not None:
        write_mesh(selected_surface_mesh, "aw_alpha_1p4_physical_surfaces.off", directory=viz)
    pml = viz / "2KAI_aw_restriction.pml"
    _write_pymol_stage_script(pml, layer_files, layer_colors, structure_path=structure_path)
    metadata = {"pair_count": len(rows),
                "bounded_complete_supported": sum(r["surface_finite_complete"] and r["surface_supported"] for r in rows),
                "aw_alpha_A": float(alpha), "alpha_selected_pairs": len(selected_pairs),
                "rejected_alpha_pairs": len(rejected_pairs), "unresolved_birth_pairs": len(unresolved_pairs),
                "not_evaluated_non_bicolor_pairs": len(non_bicolor_pairs),
                "selected_physical_surface_count": len(selected_surface_ids),
                "radius_model": "exact Network.balls.rad values in supplied AW network",
                "growth": "B_i(alpha)=B(p_i,r_i+alpha); pair minimum on actual network-clipped AW surface",
                "geometry_source": "vorpy.src.geometry.aw_alpha.births.surface_birth",
                "visualization": {"directory": str(viz), "layers": layer_files,
                                  "selected_surfaces_file": "aw_alpha_1p4_physical_surfaces.off" if selected_surface_mesh is not None else None,
                                  "pymol_script": str(pml)},
                "restriction": "experimental AW restricted-cell pair births on actual bounded network surface patches for bicolor molecular pairs; within-partner pairs are not evaluated and unresolved bicolor births remain unresolved"}
    (output_dir / "2KAI_aw_restriction_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (viz / "visualization_metadata.json").write_text(json.dumps({
        "scheme": "aw", "selection": "full Apollonius dual plus experimental AW pair restriction",
        "alpha_native": float(alpha), "alpha_units": "A",
        "growth": "r_i + alpha on actual AW pair surface",
        "full_pairs_rendered": len(full_pairs), "alpha_pairs_rendered": len(selected_pairs),
        "rejected_bicolor_pairs": len(rejected_pairs), "unresolved_bicolor_pairs": len(unresolved_pairs),
        "within_partner_pairs_not_evaluated": len(non_bicolor_pairs),
        "physical_selected_surfaces": len(selected_surface_ids),
        "unresolved_higher_simplices_rendered": 0,
        "pymol_script": str(pml),
    }, indent=2) + "\n", encoding="utf-8")
    return metadata


def _write_pymol_stage_script(path, files, colors, *, structure_path=None):
    lines = ["from pymol import cmd", "from pymol.cgo import BEGIN, END, TRIANGLES, COLOR, VERTEX", "import os"]
    if structure_path:
        lines.extend([f"cmd.load({str(Path(structure_path).resolve())!r}, '2KAI_atoms')",
                      "cmd.hide('everything','2KAI_atoms')", "cmd.show('spheres','2KAI_atoms')",
                      "cmd.set('sphere_scale',0.22,'2KAI_atoms')", "cmd.color('gray70','2KAI_atoms')"])
    lines += [
             "def _load_off(path,name,color):", "    if not os.path.exists(path): return",
             "    with open(path,encoding='utf-8') as f: data=f.readlines()",
             "    nv,nf,_=map(int,data[1].split()); pts=[tuple(map(float,x.split())) for x in data[2:2+nv]]",
             "    obj=[BEGIN,TRIANGLES,COLOR,*color]",
             "    for line in data[2+nv:2+nv+nf]:",
             "        fields=line.split(); ids=list(map(int,fields[1:1+int(fields[0])]))",
             "        for k in range(1,len(ids)-1): obj.extend([VERTEX,*pts[ids[0]],VERTEX,*pts[ids[k]],VERTEX,*pts[ids[k+1]]])",
             "    obj.append(END); cmd.load_cgo(obj,name)"]
    for name, filename in files.items():
        lines.append(f"_load_off({str((Path(path).parent / filename).resolve())!r}, {name!r}, {list(colors[name])!r})")
    physical = Path(path).parent / "aw_alpha_1p4_physical_surfaces.off"
    if physical.exists():
        lines.append(f"_load_off({str(physical.resolve())!r}, 'aw_alpha_1p4_physical_surfaces', [0.95,0.1,0.25])")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _network_atom_id(row, index):
    return f"{row.get('chn','')}:{row.get('num',index)}:{row.get('res','')}:{row.get('name','')}"


def _aw_stable_atom_map(path):
    if not path.is_file():
        return {}
    mapping = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                chain, resnum, icode, resname, atom = ast.literal_eval(row["atom_id"])
                mapping[int(row["aw_generator_id"])] = f"{chain}:{resnum}{icode}:{resname}:{atom}"
            except (KeyError, TypeError, ValueError, SyntaxError):
                continue
    return mapping


def _network_surface_area(network, features):
    values = []
    for feature in features:
        row = network.surfs.loc[feature.feature_id]
        for key in ("area", "area_A2", "surface_area"):
            try:
                value = float(row[key])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(value):
                values.append(value)
                break
    return sum(values) if values else None


def _format_summary(summary):
    lines = ["2KAI molecular dual restriction diagnostic", "", "Power stages:"]
    for key in ("full_regular_edges", "alpha0_edges", "bicolor_regular_edges", "bicolor_alpha0", "M5_rejected", "final_cazals", "interface_edges"):
        lines.append(f"  {key}: {summary[key]}")
    lines += ["", "Center-distance distributions are diagnostic only; no distance cutoff was used."]
    for row in summary["distance_summary"]:
        lines.append(f"  {row['stage']}: n={row['count']}, min={row['min_A']}, median={row['median_A']}, mean={row['mean_A']}, p95={row['p95_A']}, max={row['max_A']}")
    lines += ["", "Longest full-regular edges and exact stage decisions are in 2KAI_longest_full_regular_connections.csv.",
              "Power probe control (base vs base + 1.4 A) is in 2KAI_probe_radius_control_summary.json.",
              "AW full connections are retained, with exact experimental AW pair restriction reported separately for bicolor pairs."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, default=Path("vorpy/data/2KAI.pdb"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/dual_restriction_audit"))
    parser.add_argument("--aw-network", type=Path, default=None)
    args = parser.parse_args(argv)
    summary = run_2kai(args.pdb, args.output_dir, aw_network_path=args.aw_network)
    print(json.dumps({key: summary[key] for key in ("full_regular_edges", "alpha0_edges", "bicolor_alpha0", "M5_rejected", "final_cazals", "interface_edges")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
