"""User-facing Cazals-compatible Power/Laguerre interface runner."""
from __future__ import annotations

import csv
import math
import tempfile
from pathlib import Path

from vorpy.src.analyze.cazals_validation import (
    facet_area_statistics,
    interface_atom_asymmetry,
    interface_atom_ids,
    interface_components,
    pdb_inventory,
    power_facet_polygon_area,
)
from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import prepare_cazals_pdb


WATER_RESIDUE_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})


def aw_comparison_structure_without_waters(structure):
    """Make a temporary parser-safe PDB copy with crystallographic waters removed.

    The existing AW interface workflow builds from selected protein groups;
    water records do not belong to either requested group. The current legacy
    PDB loader cannot parse water records when its optional solvent container
    is absent, so the comparison command removes only water HETATM lines from
    this temporary input. The source PDB is never rewritten.
    """
    source = Path(structure)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    retained = []
    removed = 0
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip().upper() in WATER_RESIDUE_NAMES:
            removed += 1
        else:
            retained.append(line)
    handle = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", prefix="vorpy_aw_compare_",
                                         encoding="utf-8", delete=False)
    try:
        handle.writelines(retained)
    finally:
        handle.close()
    return Path(handle.name), removed


def extract_power_interface_options(args):
    """Remove Power-specific options while preserving the legacy argument list."""
    remaining=[]; settings={}; compare_aw=False; args=iter(args)
    names={"--power-interface":("preset",str), "--power-probe-radius":("probe_radius",float),
           "--power-alpha":("alpha",float), "--power-M":("condition_beta_m",float)}
    for arg in args:
        if arg == "--power-vs-aw":
            if compare_aw: raise SystemExit("--power-vs-aw may only be specified once")
            compare_aw=True; continue
        if arg not in names:
            remaining.append(arg); continue
        key,convert=names[arg]
        if key in settings:
            raise SystemExit(f"{arg} may only be specified once")
        try: settings[key]=convert(next(args))
        except StopIteration: raise SystemExit(f"{arg} requires a value") from None
        except ValueError: raise SystemExit(f"{arg} requires a valid value") from None
    if "preset" in settings and settings["preset"].lower() != "cazals":
        raise SystemExit("Supported Power-interface preset: cazals")
    if settings and "preset" not in settings and not compare_aw:
        raise SystemExit("Power parameters require --power-interface cazals")
    settings.setdefault("preset", "cazals" if compare_aw else None)
    settings["compare_aw"] = compare_aw
    settings.setdefault("probe_radius", 1.4)
    settings.setdefault("alpha", 0.0)
    settings.setdefault("condition_beta_m", 5.0)
    if (not math.isfinite(settings["probe_radius"]) or settings["probe_radius"] <= 0
            or not math.isfinite(settings["condition_beta_m"]) or settings["condition_beta_m"] <= 0
            or not math.isfinite(settings["alpha"])):
        raise SystemExit("Power probe radius and M must be positive; alpha must be finite")
    return remaining, settings


def chain_groups_from_legacy_args(args):
    """Read existing ``-g c CHAIN and c CHAIN`` group syntax.

    The first Power CLI release accepts chain-based selections, matching the
    existing group syntax. It rejects other legacy selectors instead of
    silently widening or changing their meaning.
    """
    groups=[]; current=None; tokens=[]; index=0
    while index < len(args):
        value=args[index]
        if value == "-g":
            if current is not None:
                groups.append(_chain_group(current))
            current=[]; index+=1; continue
        if value.startswith("-"):
            if current is not None:
                groups.append(_chain_group(current)); current=None
            index+=1
            # consume value of an export directory command
            if value in {"-e","--export"} and index<len(args):
                if args[index] in {"dir","directory"}:
                    index+=2
                else:
                    index+=1
            continue
        if current is not None:
            current.append(value)
        index+=1
    if current is not None: groups.append(_chain_group(current))
    if len(groups)!=2:
        raise SystemExit("Power-interface mode requires exactly two -g groups, for example: -g c A and c B -g c I")
    if groups[0] & groups[1]: raise SystemExit("Power-interface chain groups must be disjoint")
    return groups[0],groups[1]


def _chain_group(tokens):
    chains=set(); part=[]
    for token in [*tokens,"and"]:
        if token.lower()=="and":
            if not part: continue
            if len(part)<2 or part[0].lower() not in {"c","cs","chain"}:
                raise SystemExit("Cazals Power preset currently accepts existing chain selectors only: -g c CHAIN")
            chains.add(part[1]); part=[]
        else: part.append(token)
    return chains


def power_output_directory(args, structure):
    """Honor the standard ``-e dir BASE`` output-directory convention."""
    for i,token in enumerate(args[:-2]):
        if token in {"-e","--export"} and args[i+1] in {"dir","directory"}:
            return Path(args[i+2]).expanduser().resolve() / Path(structure).stem
    return Path("output") / Path(structure).stem


def run_power_interface_pdb(structure, group_a_chains, group_b_chains, output_dir, *,
                            probe_radius=1.4, alpha=0.0, condition_beta_m=5.0):
    """Construct, summarize, and export a Cazals-compatible Power interface."""
    structure=Path(structure); output_dir=Path(output_dir)
    if structure.suffix.lower()!=".pdb":
        raise ValueError("The Cazals Power preset currently requires a PDB input")
    prepared=prepare_cazals_pdb(structure,group_a_chains,group_b_chains,
                                probe_radius=probe_radius,include_hydrogens=False,
                                include_hetatm=False,strict_heavy_fallbacks=False,verbose=False)
    result=analyze_power_interface_curvature(prepared.points,prepared.expanded_radii,
             prepared.group_a,prepared.group_b,alpha=alpha,condition_beta_m=condition_beta_m,
             fallback_radius_assignments=len(prepared.fallback_atoms))
    output_dir.mkdir(parents=True,exist_ok=True)
    atoms={a.index:a for a in prepared.atoms}
    final={f.generator_ids for f in result.final_facets}
    areas={}; facet_rows=[]
    for facet in result.final_facets:
        poly=power_facet_polygon_area(result,facet); pair=facet.generator_ids
        areas[pair]=float(poly["area_A2"] or 0.)
        facet_rows.append({"facet_id":f"{pair[0]}-{pair[1]}","atom_i":pair[0],"atom_j":pair[1],
            "group_i":atoms[pair[0]].group,"group_j":atoms[pair[1]].group,
            "residue_i":_residue_label(atoms[pair[0]]),"residue_j":_residue_label(atoms[pair[1]]),
            "alpha_filtration":result.filtration.get(pair),"condition_b_m":facet.largest_orthogonal_ball_radius,
            "condition_b_status":facet.condition_beta_status,
            "condition_b_r":facet.smaller_expanded_radius,"condition_b_m_over_r":facet.m_over_r,
            "condition_b_accepted":facet.condition_beta_accepted,"area_A2":poly["area_A2"],
            "polygon_vertex_count":poly["polygon_vertex_count"],"polygon_vertices":repr(poly.get("polygon_vertices",())),
            "component_id":"","status":poly["status"]})
    components=interface_components(result,areas)
    component_by_facet={facet:component["component_id"] for component in components for facet in component["facets"]}
    for row in facet_rows: row["component_id"]=component_by_facet.get(tuple(sorted((row["atom_i"],row["atom_j"]))),"")
    interface_a,interface_b=interface_atom_ids(result.final_facets,result.group_a)
    atom_rows=[]
    for group_name,indices in (("A",sorted(interface_a)),("B",sorted(interface_b))):
        for atom_id in indices:
            a=atoms[atom_id]
            atom_rows.append({"atom_id":atom_id,"pdb_serial":a.pdb_serial,"group":group_name,
                "chain":a.chain,"residue":_residue_label(a),"atom_name":a.atom_name,
                "element":a.element,"x":a.x,"y":a.y,"z":a.z,"base_radius_A":a.base_radius,
                "expanded_radius_A":a.expanded_radius,"radius_class":a.radius_class})
    edge_rows=[]; boundary_edges=0
    for edge_id,edge in enumerate(result.edges,1):
        selected=[pair for pair in edge.bicolor_facets if pair in final]
        if edge.geometry_status!="finite": status="unbounded_or_unresolved"
        elif len(selected)==2: status="interior"
        elif len(selected)==1: status="boundary"; boundary_edges+=1
        else: status="outside_selected_interface"
        edge_rows.append({"edge_id":edge_id,"triangle_generators":";".join(map(str,edge.regular_triangle)),
            "generator_residues":";".join(_residue_label(atoms[i])+":"+atoms[i].atom_name for i in edge.regular_triangle),
            "group_pattern":edge.group_composition,"singleton_atom":edge.singleton_generator,
            "incident_facets":";".join(f"{i}-{j}" for i,j in edge.bicolor_facets),
            "selected_incident_facets":";".join(f"{i}-{j}" for i,j in selected),
            "selected_facet_count":len(selected),"length_A":edge.length,"beta_rad":edge.beta_radians,
            "beta_deg":edge.beta_degrees,"l_beta_A_rad":edge.length_beta,
            "l_abs_beta_A_rad":edge.length_abs_beta,"geometry_status":edge.geometry_status,"status":status})
    _write_csv(output_dir/"power_interface_atoms.csv",atom_rows)
    _write_csv(output_dir/"power_interface_facets.csv",facet_rows)
    _write_csv(output_dir/"power_interface_edges.csv",edge_rows)
    _write_csv(output_dir/"power_interface_components.csv",[
        {"component_id":c["component_id"],"facet_count":c["facet_count"],"VIA_A2":c["area_A2"],
         "VIA_fraction":c["via_fraction"],"significant":c["significant"],
         "facet_ids":";".join(f"{i}-{j}" for i,j in c["facets"])} for c in components])
    summary=_user_summary(structure,prepared,result,interface_a,interface_b,facet_rows,components,boundary_edges,
                          group_a_chains,group_b_chains)
    (output_dir/"power_interface_summary.txt").write_text(summary,encoding="utf-8")
    return {"result":result,"prepared":prepared,"summary":summary,"output_dir":output_dir,
            "facet_records":facet_rows,"facet_area_by_facet":areas}


def _residue_label(atom):
    return f"{atom.chain}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}"


def _write_csv(path,rows):
    rows=list(rows); path=Path(path)
    fields=list(rows[0]) if rows else []
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def _user_summary(structure,prepared,result,interface_a,interface_b,facets,components,boundary_edges,
                  group_a_chains,group_b_chains):
    inventory=pdb_inventory(structure,set(group_a_chains)|set(group_b_chains))
    stats=facet_area_statistics(float(row["area_A2"]) for row in facets if row["area_A2"] is not None)
    totals=result.edge_totals
    alpha_selected=[f for f in result.all_bicolor_facets
                    if f.regular_edge_filtration <= result.alpha + result.alpha_tolerance]
    beta_rejected=sum(f.condition_beta_status == "rejected_m_over_r" for f in alpha_selected)
    return ("Power/Laguerre Interface Analysis\n"
        f"Input: {structure}\nGeometry: weighted Power/Laguerre diagram (regular complex)\n"
        f"Chemistry/radius model: Cazals/Chothia; R=r+probe; w=R^2\nProbe: {prepared.probe_radius:.3f} A; alpha: {result.alpha:g}; M: {result.condition_beta_m:g}\n"
        f"Group A chains: {sorted(group_a_chains)}\nGroup B chains: {sorted(group_b_chains)}\n\n"
        f"INPUT\nAtoms retained: {len(prepared.atoms)}; group A/B atoms: {len(prepared.group_a)}/{len(prepared.group_b)}\n"
        f"Excluded water records: {inventory['excluded_water_atoms']}; excluded other HETATM: {inventory['excluded_other_hetatm_atoms']}\n"
        f"Radius fallbacks: {len(prepared.fallback_atoms)}\n\nINTERFACE\n"
        f"Interface atoms A/B/total: {len(interface_a)}/{len(interface_b)}/{len(interface_a)+len(interface_b)}; r_AB={interface_atom_asymmetry(interface_a,interface_b):.8g}\n"
        f"Bicolor regular facets: {len(result.all_bicolor_facets)}; alpha-selected: {len(alpha_selected)}; condition-b rejected: {beta_rejected}; final selected: {len(facets)}\n"
        f"Finite polygon facets: {sum(row['area_A2'] is not None for row in facets)}; unresolved polygons: {sum(row['area_A2'] is None for row in facets)}\n"
        f"VIA={stats['VIA_A2']:.9g} A^2; mean/median area={_format_stat(stats['mean_A2'])}/{_format_stat(stats['median_A2'])} A^2\n\n"
        f"TOPOLOGY\nConnected components: {len(components)}; significant components (>=7.5% VIA): {sum(c['significant'] for c in components)}\n\n"
        f"CURVATURE\nInterior edges: {len(result.included_edges)}; boundary edges: {boundary_edges}; total interior edge length L={totals['length']:.9g} A\n"
        f"Positive/negative edges: {totals['positive_edges']}/{totals['negative_edges']}\n"
        f"C_signed=sum(l beta)={totals['signed']:.9g} A rad; C_unsigned=sum(l |beta|)={totals['unsigned']:.9g} A rad\n"
        f"s_H signed={totals['signed_mean_rad']:.9g} rad / {totals['signed_mean_deg']:.9g} deg; unsigned={totals['unsigned_mean_rad']:.9g} rad / {totals['unsigned_mean_deg']:.9g} deg\n"
        f"Conventional H=(k1+k2)/2 edge measure C_H=C_signed/2={totals['signed']/2:.9g} A rad. C_signed is the raw Cazals turning-angle sum.\n")


def _format_stat(value):
    return "n/a" if value is None or not math.isfinite(float(value)) else f"{float(value):.9g}"
