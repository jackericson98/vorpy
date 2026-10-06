"""Generate an evidence-led 2KAI interface curvature/visualization checkpoint.

This exporter only reuses solved 2KAI geometry and saved matched-probe results.
It deliberately leaves union-of-balls contact patches, Gaussian interface
curvature, and finite-solvent tessellations unresolved where VorPy has no
validated construction for them.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import json
import math
import pickle
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

ROOT = Path("output/2KAI_curvature_showcase")
MATCHED = Path("output/comparative_study/matched_probe_control/2KAI")
NETWORK_PATH = Path("output/comparative_study/cases/2KAI/H/full_aw.vpy")
STATIC_PDB = Path("vorpy/data/2KAI.pdb")
SOLVATED_PDB = Path("vorpy/data/2KAI_MD_centered.pdb")
PROBE = 1.4
AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    names = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=names, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def _pdb_rows(path):
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[12:16].strip(); res = line[17:20].strip().upper()
        rows.append({"record": line[:6].strip(), "serial": line[6:11].strip(),
            "atom": name, "resname": res, "chain": line[21:22].strip(),
            "resnum": line[22:26].strip(), "icode": line[26:27].strip(),
            "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            "element": (line[76:78].strip().upper() if len(line) >= 78 else "") or
                       ("H" if name.upper().startswith(("H", "1H", "2H", "3H")) else name[0].upper()),
            "line": line})
    return rows


def audit_solvated(path=SOLVATED_PDB):
    rows = _pdb_rows(path)
    water_res = {"SOL", "HOH", "WAT", "TIP3", "TIP3P"}
    protein = [x for x in rows if x["resname"] in AA]
    waters = [x for x in rows if x["resname"] in water_res]
    ions = [x for x in rows if x["resname"] not in AA and x["resname"] not in water_res]
    # PDB chain fields are blank in this simulation export. Reconstruct chain
    # segmentation only when its protein residue sequence exactly matches the
    # ordered A/B/I sequence in the crystallographic PDB.
    static = _pdb_rows(STATIC_PDB)
    seq_by_chain = {}
    for chain in ("A", "B", "I"):
        residues = {}
        for x in static:
            if x["chain"] == chain and x["resname"] in AA:
                residues.setdefault((int(x["resnum"]), x["icode"], x["resname"]), None)
        seq_by_chain[chain] = [(r[2]) for r in sorted(residues)]
    md_res = []
    seen = None
    for x in protein:
        key = (x["resnum"], x["resname"])
        if key != seen:
            md_res.append(key); seen = key
    expected = [(res, "A") for res in seq_by_chain["A"]] + [(res, "B") for res in seq_by_chain["B"]] + [(res, "I") for res in seq_by_chain["I"]]
    matcher = difflib.SequenceMatcher(None, [r for r,_ in expected], [r for _,r in md_res], autojunk=False)
    chain_for_md_residue = {}
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            for si, mi in zip(range(i1,i2), range(j1,j2)):
                chain_for_md_residue[mi] = expected[si][1]
    # A chain map is accepted only if every MD protein residue is matched and
    # the only discrepancy is a static residue absent in the simulation model.
    sequence_match = (len(chain_for_md_residue) == len(md_res) and
        all(tag in {"equal","delete"} for tag, *_ in matcher.get_opcodes()))
    heavy = [x for x in protein if x["element"] != "H"]
    hydrogens = [x for x in protein if x["element"] == "H"]
    groups = {"A": 0, "B": 0, "I": 0}
    if sequence_match:
        residue_index = -1; last_key = None
        for atom in protein:
            key = (atom["resnum"], atom["resname"])
            if key != last_key:
                residue_index += 1; last_key = key
            chain = chain_for_md_residue[residue_index]
            atom["reconstructed_chain"] = chain
            groups[chain] += 1
    return {"path": str(Path(path)), "all_atom_records": len(rows),
        "ATOM": sum(x["record"] == "ATOM" for x in rows),
        "HETATM": sum(x["record"] == "HETATM" for x in rows),
        "protein_atoms": len(protein), "protein_heavy_atoms": len(heavy),
        "protein_hydrogens": len(hydrogens), "water_atoms": len(waters),
        "water_residues": len({(x["resnum"], x["resname"]) for x in waters}),
        "ion_or_other_atoms": len(ions), "other_residue_names": sorted({x["resname"] for x in ions}),
        "chains_in_file": sorted({x["chain"] for x in rows}),
        "protein_residues": len(md_res), "protein_chain_sequence_matches_static_A_B_I": sequence_match,
        "protein_residue_sequence_alignment": "exact, with at most a static-only missing residue" if sequence_match else "unresolved",
        "protein_atom_counts_by_reconstructed_chain": groups,
        "available_solvated_geometry_logs_atom_count": 4340,
        "available_logs_include_water_or_ions": False,
        "solvated_network_status": "NOT_SOLVED: existing network/log cache omits water and chloride generators"}


def _parse_off(path):
    lines = [line for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    nv, nf, _ = map(int, lines[1].split())
    points = np.array([[float(v) for v in line.split()[:3]] for line in lines[2:2+nv]])
    faces = []
    for line in lines[2+nv:2+nv+nf]:
        xs = line.split(); faces.append(list(map(int, xs[1:1+int(xs[0])])))
    return points, faces


def _write_off(path, points, faces):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write("OFF\n%d %d 0\n" % (len(points), len(faces)))
        for p in points: f.write("%.9f %.9f %.9f\n" % tuple(p))
        for face in faces: f.write(str(len(face)) + " " + " ".join(map(str, face)) + "\n")


def _power_faces(result, areas, pairs):
    points, faces = [], []
    for pair in sorted(pairs):
        entry = areas.get(f"{pair[0]}:{pair[1]}", {})
        poly = entry.get("polygon_vertices") or []
        if len(poly) < 3: continue
        start = len(points); points.extend(poly)
        faces.append(list(range(start, start + len(poly))))
    return points, faces


def _components(pairs, edge_pair_map):
    adj = {p: set() for p in pairs}
    for ps in edge_pair_map:
        selected = sorted(set(ps) & pairs)
        for p in selected:
            adj[p].update(x for x in selected if x != p)
    count = 0; unseen = set(adj)
    while unseen:
        count += 1; stack = [unseen.pop()]
        while stack:
            for n in adj[stack.pop()] & unseen:
                unseen.remove(n); stack.append(n)
    return count


def _edge_angle(center, prev, next_):
    a=np.asarray(prev)-np.asarray(center); b=np.asarray(next_)-np.asarray(center)
    den=float(np.linalg.norm(a)*np.linalg.norm(b))
    return math.acos(float(np.clip(a@b/den,-1.0,1.0))) if den else float("nan")


def _power_gaussian(power, area_cache, chosen_pairs):
    """Polyhedral Gauss-Bonnet terms on the finite selected Power facets."""
    chosen=set(chosen_pairs); by_vertex=defaultdict(float); face_count=0; bad_faces=0
    selected_edge_counts={}; edge_ends={}; edge_lookup={}
    for edge in power.edges:
        incident=set(tuple(sorted(p)) for p in edge.bicolor_facets)&chosen
        if not incident: continue
        key=tuple(sorted(edge.regular_triangle)); selected_edge_counts[key]=len(incident)
        tets=[tuple(sorted(t)) for t in edge.incident_tetrahedra if tuple(sorted(t)) in power.power_vertices]
        if len(tets)==2: edge_ends[key]=tuple(tets)
        else: bad_faces+=1
        edge_lookup[key]=edge
    for pair in chosen:
        item=area_cache.get(f"{pair[0]}:{pair[1]}",{})
        poly=np.asarray(item.get("polygon_vertices") or [],dtype=float)
        if item.get("status")!="finite_bounded" or len(poly)<3:
            bad_faces+=1; continue
        face_count+=1
        mapped=[]
        for xyz in poly:
            dists=[(float(np.linalg.norm(np.asarray(v.position)-xyz)),tet) for tet,v in power.power_vertices.items()]
            dist,tet=min(dists,key=lambda x:x[0]) if dists else (float("inf"),None)
            if dist>2e-5 or tet is None: mapped=[]; break
            mapped.append(tuple(tet))
        if len(mapped)!=len(poly) or len(set(mapped))<3:
            bad_faces+=1; continue
        for k,tet in enumerate(mapped):
            by_vertex[tet]+=_edge_angle(poly[k],poly[k-1],poly[(k+1)%len(poly)])
    boundary_vertices=set()
    for edge,count in selected_edge_counts.items():
        if count==1 and edge in edge_ends: boundary_vertices.update(edge_ends[edge])
    all_vertices=set(by_vertex)
    vertex_terms={v:((math.pi-a) if v in boundary_vertices else (2*math.pi-a)) for v,a in by_vertex.items()}
    boundary=sum(value for v,value in vertex_terms.items() if v in boundary_vertices)
    interior=sum(value for v,value in vertex_terms.items() if v not in boundary_vertices)
    edge_count=len(selected_edge_counts); vertex_count=len(all_vertices)
    chi=vertex_count-edge_count-face_count
    total=interior+boundary
    expected=2*math.pi*chi
    complete=(bad_faces==0 and all(k in edge_ends for k in selected_edge_counts) and
              all(c in (1,2) for c in selected_edge_counts.values()))
    bgraph=defaultdict(set)
    for edge,count in selected_edge_counts.items():
        if count==1 and edge in edge_ends:
            a,b=edge_ends[edge]; bgraph[a].add(b); bgraph[b].add(a)
    loops=None
    if bgraph and all(len(nbrs)==2 for nbrs in bgraph.values()):
        unseen=set(bgraph); loops=0
        while unseen:
            loops+=1; todo=[unseen.pop()]
            while todo:
                for v in bgraph[todo.pop()]&unseen: unseen.remove(v); todo.append(v)
    residual=total-expected if complete else None
    gb_valid=bool(complete and residual is not None and abs(residual)<=1e-6*max(1.0,abs(expected)))
    status=("available" if gb_valid else
            "UNRESOLVED: Gauss-Bonnet residual exceeds tolerance" if complete else
            "UNRESOLVED: incomplete/nonmanifold mapped Power incidence")
    return {"smooth":0.0,"edge":0.0,"vertex":interior,"boundary":boundary,
        "total":total if gb_valid else None,"calculated_total":total,"chi":chi,"expected":expected,
        "residual":residual,"vertices":vertex_count,"edges":edge_count,
        "faces":face_count,"bad_features":bad_faces,"unresolved_vertices":bad_faces,
        "components":_components(chosen,[e.bicolor_facets for e in power.edges]),
        "boundary_edges":sum(v==1 for v in selected_edge_counts.values()),
        "boundary_loops":loops,"vertex_terms":vertex_terms,"boundary_vertices":boundary_vertices,
        "status":status}


def _aw_gaussian(network, chosen_pairs, pair_supported, group_a, *, order=32):
    """AW face/edge/vertex Gaussian terms using existing VorPy integrators."""
    from vorpy.src.calculations.edge_geometry import oriented_aw_face_edge_geodesic_curvature
    from vorpy.src.calculations.vertex_geometry import aw_vertex_face_angle
    from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge

    surface_by_pair=defaultdict(list)
    for sid,row in network.surfs.iterrows():
        pair=tuple(sorted(map(int,row.get("balls",()))))
        if len(pair)==2 and pair in chosen_pairs: surface_by_pair[pair].append(int(sid))
    face_ids=[]; unsupported=[]
    for pair in chosen_pairs:
        for sid in surface_by_pair.get(pair,[]):
            if pair_supported.get(pair,False): face_ids.append((pair,sid))
            else: unsupported.append((pair,sid))
    # Resolve each boundary edge once and reuse it for all selected face terms.
    face_edge_ids={sid:{int(e) for e in network.surfs.loc[sid,"edges"]} for _,sid in face_ids}
    all_edge_ids=set().union(*face_edge_ids.values()) if face_edge_ids else set()
    resolved={}; unresolved_edges=set()
    for eid in all_edge_ids:
        try: resolved[eid]=resolve_aw_network_edge(network,eid)
        except Exception: unresolved_edges.add(eid)
    smooth=0.0; face_kg=defaultdict(float); angle_sums=defaultdict(float); face_count=0
    surface_rows=[]; face_supported_ids={sid for _,sid in face_ids}
    for pair,sid in face_ids:
        row=network.surfs.loc[sid]
        a=next((i for i in pair if i in group_a),None)
        b=next((i for i in pair if i!=a),None)
        if a is None or b is None:
            unsupported.append((pair,sid)); face_supported_ids.discard(sid); continue
        k=float(row.get("int_gauss_curv",float("nan")))
        if not math.isfinite(k): unsupported.append((pair,sid)); face_supported_ids.discard(sid); continue
        smooth+=k; face_count+=1
        kg=0.0; bad=False
        for eid in face_edge_ids[sid]:
            if eid not in resolved: bad=True; continue
            er=resolved[eid]
            try:
                kg += oriented_aw_face_edge_geodesic_curvature(er.geometry,er.ball_indices,er.locations,a,b,order=order)
            except Exception: bad=True
        face_kg[sid]=kg
        face_angles={}
        for vid in map(int,row.get("verts",())):
            try: face_angles[vid]=float(aw_vertex_face_angle(network,vid,a,sid,resolved_edges=resolved))
            except Exception: bad=True
        if not bad:
            for vid,angle in face_angles.items(): angle_sums[vid]+=angle
        if bad: unsupported.append((pair,sid)); face_supported_ids.discard(sid)
        surface_rows.append({"surface_id":sid,"generator_tuple":pair,"smooth_gaussian":k,
            "face_boundary_geodesic_gaussian":kg,"status":"supported" if not bad else "unresolved_geometry"})
    # Build selected face incidence from network edge/surface tables.
    edge_selected=defaultdict(set); vertex_selected=defaultdict(set)
    for pair,sid in face_ids:
        if sid not in face_supported_ids: continue
        for eid in map(int,network.surfs.loc[sid,"edges"]): edge_selected[eid].add(sid)
        for vid in map(int,network.surfs.loc[sid,"verts"]): vertex_selected[vid].add(sid)
    boundary_edges={eid for eid,sids in edge_selected.items() if len(sids)==1}
    nonmanifold={eid for eid,sids in edge_selected.items() if len(sids)>2}
    edge_interior=sum(face_kg.get(sid,0.0) for eid,sids in edge_selected.items() if len(sids)==2 for sid in sids)
    edge_boundary=sum(face_kg.get(sid,0.0) for eid,sids in edge_selected.items() if len(sids)==1 for sid in sids)
    edge_contributions={eid:sum(face_kg.get(sid,0.0) for sid in sids) for eid,sids in edge_selected.items()}
    interior_vertex=boundary_vertex=0.0
    boundary_vertex_ids=set()
    for eid in boundary_edges:
        boundary_vertex_ids.update(map(int,network.edges.loc[eid,"verts"]))
    vertex_terms={}
    for vid,sids in vertex_selected.items():
        angle=float(angle_sums.get(vid,float("nan")))
        if not math.isfinite(angle): continue
        term=(math.pi-angle) if vid in boundary_vertex_ids else (2*math.pi-angle)
        vertex_terms[vid]=term
        if vid in boundary_vertex_ids: boundary_vertex+=term
        else: interior_vertex+=term
    # Cell-complex Euler characteristic for the supported selected surface subcomplex.
    vertices=set(vertex_selected); edges=set(edge_selected); faces=len(face_supported_ids)
    chi=len(vertices)-len(edges)+faces
    total=smooth+edge_interior+edge_boundary+interior_vertex+boundary_vertex
    expected=2*math.pi*chi
    complete=(not unsupported and not unresolved_edges and not nonmanifold and
              all(len(sids) in (1,2) for sids in edge_selected.values()) and len(vertex_terms)==len(vertices))
    residual=total-expected if complete else None
    gb_valid=bool(complete and residual is not None and abs(residual)<=1e-6*max(1.0,abs(expected)))
    # A valid boundary is a union of loops: each boundary vertex has degree two.
    boundary_degree=defaultdict(int)
    for eid in boundary_edges:
        for vid in map(int,network.edges.loc[eid,"verts"]): boundary_degree[vid]+=1
    loops=None
    if boundary_edges and all(deg==2 for deg in boundary_degree.values()):
        graph=defaultdict(set)
        for eid in boundary_edges:
            v=tuple(map(int,network.edges.loc[eid,"verts"])); graph[v[0]].add(v[1]); graph[v[1]].add(v[0])
        unseen=set(graph); loops=0
        while unseen:
            loops+=1; todo=[unseen.pop()]
            while todo:
                for v in graph[todo.pop()]&unseen: unseen.remove(v); todo.append(v)
    return {"smooth":smooth,"edge":edge_interior+edge_boundary,"edge_interior":edge_interior,
        "edge_boundary":edge_boundary,"vertex":interior_vertex,"boundary":boundary_vertex,
        "total":total if gb_valid else None,"partial_total":total,"calculated_total":total,"chi":chi,"expected":expected,
        "residual":residual,"vertices":len(vertices),"edges":len(edges),
        "faces":faces,"boundary_edges":len(boundary_edges),"boundary_loops":loops,
        "unresolved_surfaces":len(set(sid for _,sid in unsupported)),"unresolved_edges":len(unresolved_edges),
        "unresolved_vertices":len(vertices)-len(vertex_terms),"nonmanifold_edges":len(nonmanifold),
        "components":None,"surface_rows":surface_rows,"vertex_rows":[{"vertex_id":v,"angle_defect_rad":value,"status":"supported"} for v,value in vertex_terms.items()],
        "edge_rows":[{"edge_id":eid,"geodesic_gaussian":value,"selected_surface_count":len(edge_selected[eid]),"status":"supported"} for eid,value in edge_contributions.items()],
        "status":"available" if gb_valid else
            "UNRESOLVED: Gauss-Bonnet residual exceeds tolerance" if complete else
            "PARTIAL/UNRESOLVED: unsupported geometry or incidence"}


def _project_off(points, faces, basis, center):
    xy = (np.asarray(points) - center) @ basis[:2].T
    projected_faces = [xy[f] for f in faces if len(f) >= 3]
    return xy, projected_faces


def _export_pymol(root, pdb, colors, layers, title):
    from vorpy.src.geometry.visualization.export import _edge_mesh
    del _edge_mesh  # assert/reuse exporter ecosystem; geometry already uses its OFF output
    py = ["from pymol import cmd", "from pymol.cgo import BEGIN, TRIANGLES, COLOR, VERTEX, END", "import os"]
    py.append(f"cmd.load({str(Path(pdb).resolve())!r}, 'molecule')")
    for sel, color in (("chain A or chain B", "blue"), ("chain I", "red")):
        py.append(f"cmd.color({color!r}, 'molecule and ({sel})')")
    py += ["cmd.hide('everything', 'molecule')", "cmd.show('cartoon', 'molecule and polymer.protein')",
           "cmd.set('cartoon_transparency', 0.35)", "cmd.show('spheres', 'molecule and polymer.protein and name CA')",
           "cmd.set('sphere_scale', 0.20)"]
    py += ["def _off(path, name, color):", "    if not os.path.exists(path): return",
           "    rows=[line for line in open(path,encoding='utf8').read().splitlines() if line.strip()]", "    nv,nf,_=map(int,rows[1].split())",
           "    pts=[tuple(map(float,x.split())) for x in rows[2:2+nv]]", "    obj=[BEGIN,TRIANGLES,COLOR,*color]",
           "    for line in rows[2+nv:2+nv+nf]:", "        fields=line.split(); n=int(fields[0]); ids=list(map(int,fields[1:1+n]))",
           "        for k in range(1,len(ids)-1): obj += [VERTEX,*pts[ids[0]],VERTEX,*pts[ids[k]],VERTEX,*pts[ids[k+1]]]",
           "    obj.append(END); cmd.load_cgo(obj,name)"]
    for name, file, color in layers:
        py.append(f"_off({str((root/file).resolve())!r}, {name!r}, {color!r})")
    py += ["cmd.disable('all')", "cmd.enable('molecule')"]
    for name, _, _ in layers: py.append(f"cmd.enable({name!r})")
    py += ["cmd.orient('molecule')", "cmd.bg_color('white')", "cmd.set('ray_opaque_background', 0)",
           "cmd.set('antialias', 2)", "cmd.set('orthoscopic', 1)",
           f"cmd.set_title('molecule', {title!r})"]
    (root / "showcase_layers.py").write_text("\n".join(py) + "\n", encoding="utf-8")
    (root / "showcase_layers.pml").write_text(f'run "{(root / "showcase_layers.py").resolve().as_posix()}"\n', encoding="utf-8")


def run(output=ROOT):
    from vorpy.src.io import load_network
    from vorpy.src.geometry.visualization.export import _surface_mesh
    from vorpy.src.output.mesh import write_mesh

    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=True)
    viz = MATCHED / "visualization"
    summary = json.loads((MATCHED / "comparison_summary.json").read_text(encoding="utf-8"))
    comparison_rows = read_csv(MATCHED / "shared_pair_surfaces.csv")
    shared = {tuple(sorted((int(r["generator_i"]), int(r["generator_j"])))) for r in comparison_rows}
    aw_selected = {tuple(sorted((int(r["generator_i"]), int(r["generator_j"])))) for r in read_csv(MATCHED / "bicolor_pair_comparison.csv") if r["in_AW"] == "True"}
    power_selected = {tuple(sorted((int(r["generator_i"]), int(r["generator_j"])))) for r in read_csv(MATCHED / "bicolor_pair_comparison.csv") if r["in_Power"] == "True"}
    net = load_network(NETWORK_PATH)
    # The selected geometry uses the exact same saved network and pair IDs as
    # the matched-probe numerical checkpoint.
    pos_by_pair = {}
    for sid, row in net.surfs.iterrows():
        pair = tuple(sorted(map(int, row.get("balls", ()))))
        if len(pair) == 2: pos_by_pair.setdefault(pair, []).append(net.surfs.index.get_loc(sid))
    surface_positions = [pos for pair in sorted(aw_selected) for pos in pos_by_pair.get(pair, [])]
    aw_mesh = _surface_mesh(net, surface_positions)
    if aw_mesh is not None: write_mesh(aw_mesh, "static_aw_alpha_physical.off", directory=output)
    area_cache = json.loads((MATCHED / "power_facet_area_cache.json").read_text(encoding="utf-8"))["areas"]
    with (MATCHED / "power_geometry_cache.pkl").open("rb") as f: _, power = pickle.load(f)
    atom_audit=read_csv(MATCHED/"input_identity_audit.csv")
    group_a={int(r["generator_id"]) for r in atom_audit if r["chain"] in {"A","B"}}
    pair_rows=read_csv(MATCHED/"bicolor_pair_comparison.csv")
    aw_support={tuple(sorted((int(r["generator_i"]),int(r["generator_j"])))):r["AW_surface_status"]=="supported"
                for r in pair_rows if r["in_AW"]=="True"}
    gaussian_cache=output/"gaussian_components_cache.pkl"
    gaussian_key=("curvature-showcase-gaussian-v2-gb-validated",tuple(sorted(power_selected)),tuple(sorted(aw_selected)))
    cached_gaussian=None
    if gaussian_cache.exists():
        try:
            with gaussian_cache.open("rb") as f: cached_key,cached_gaussian=pickle.load(f)
            if cached_key!=gaussian_key: cached_gaussian=None
        except (OSError,EOFError,pickle.PickleError,ValueError): cached_gaussian=None
    if cached_gaussian is None:
        power_gaussian=_power_gaussian(power,area_cache,power_selected)
        aw_gaussian=_aw_gaussian(net,aw_selected,aw_support,group_a)
        cached_gaussian=(power_gaussian,aw_gaussian)
        with gaussian_cache.open("wb") as f: pickle.dump((gaussian_key,cached_gaussian),f,pickle.HIGHEST_PROTOCOL)
    else:
        power_gaussian,aw_gaussian=cached_gaussian
    p_pts, p_faces = _power_faces(None, area_cache, power_selected)
    _write_off(output / "static_power_alpha_physical.off", p_pts, p_faces)
    s_pts, s_faces = _power_faces(None, area_cache, shared)
    _write_off(output / "shared_power_physical.off", s_pts, s_faces)
    # Reuse exact AW tessellated surface representation already created by the
    # Phase-5 exporter for the shared pair set.
    (output / "shared_aw_physical.off").write_bytes((viz / "shared_AW_physical_surfaces.off").read_bytes())

    # Feature-level mean-curvature records and explicit null Gaussian terms.
    aw_curv_surfs = {int(r["surface_id"]): r for r in read_csv(Path("output/comparative_study/cases/2KAI/H/aw/aw_interface_surfaces.csv"))}
    shared_rows = {tuple(sorted((int(r["generator_i"]), int(r["generator_j"])))): r for r in comparison_rows}
    aw_surface_features = []
    for pair in sorted(aw_selected):
        cr = shared_rows.get(pair)
        if cr:
            # Shared-pair rows have a direct AW ID; nonshared selected rows map
            # through the solved network's surface incidence.
            ids = [int(x) for x in cr["AW_surface_id"].split(";") if x.isdigit()]
        else:
            ids = [int(sid) for sid, row in net.surfs.iterrows() if tuple(sorted(map(int,row.get("balls",())))) == pair]
        for sid in ids:
            src = aw_curv_surfs.get(sid, {})
            aw_surface_features.append({"surface_id": sid, "generator_i": pair[0], "generator_j": pair[1],
                "area_A2": src.get("area", cr.get("AW_area_A2") if cr else None),
                "integral_H_dA_A": src.get("integrated_mean_curvature"),
                "integral_K_dA": float(net.surfs.loc[sid,"int_gauss_curv"]), "support_status": src.get("status", "unavailable"),
                "gaussian_status": "existing AW analytic-surface quadrature"})
    write_csv(output / "aw_surface_curvature_components.csv", aw_surface_features)
    power_surface_features = []
    for pair in sorted(power_selected):
        entry = area_cache.get(f"{pair[0]}:{pair[1]}", {})
        power_surface_features.append({"surface_id": f"power_pair:{pair[0]}:{pair[1]}", "generator_i": pair[0], "generator_j": pair[1],
            "area_A2": entry.get("area_A2"), "integral_H_dA_A": 0.0 if entry.get("status") == "finite_bounded" else None,
            "integral_K_dA": 0.0 if entry.get("status") == "finite_bounded" else None, "support_status": entry.get("status"),
            "gaussian_status": "planar Power facet interior; concentrated vertex/boundary terms audited separately"})
    write_csv(output / "power_surface_curvature_components.csv", power_surface_features)

    # Per-edge curvature data are sourced from the frozen analysis cache/logs.
    power_edges = []
    for e in power.edges:
        incident = [tuple(sorted(p)) for p in e.bicolor_facets if tuple(sorted(p)) in power_selected]
        if not incident: continue
        value = e.length_beta
        power_edges.append({"edge_id": ":".join(map(str,e.regular_triangle)), "generator_tuple": ";".join(map(str,e.regular_triangle)),
            "length_A": e.length, "signed_mean_turn_rad": e.beta_radians if value is not None else None,
            "raw_signed_contribution_A_rad": value, "conventional_mean_contribution_A": None if value is None else .5*value,
            "gaussian_geodesic_contribution": 0.0, "support_status": e.geometry_status,
            "selected_surface_count": len(incident), "classification": "interior_selected" if len(incident)==2 else "boundary_selected"})
    write_csv(output / "power_edge_curvature_components.csv", power_edges)
    aw_edge_source = read_csv(Path("output/comparative_study/cases/2KAI/H/aw/aw_interface_edges.csv"))
    aw_pair_to_sid = defaultdict(set)
    for r in aw_surface_features: aw_pair_to_sid[tuple(sorted((int(r["generator_i"]),int(r["generator_j"]))))].add(int(r["surface_id"]))
    aw_edge_features = []
    aw_gaussian_by_edge={int(r["edge_id"]):r["geodesic_gaussian"] for r in aw_gaussian["edge_rows"]}
    aw_selected_surface_ids=set().union(*(aw_pair_to_sid[p] for p in aw_selected if p in aw_pair_to_sid))
    for e in aw_edge_source:
        surf_ids = {int(x) for x in e["surface_ids"].split(";") if x.isdigit()}
        incident = surf_ids & aw_selected_surface_ids
        if not incident: continue
        aw_edge_features.append({"edge_id": e["edge_id"], "generator_tuple": e["generator_ids"], "length_A": e["length"],
            "raw_signed_contribution_A_rad": e["signed_contribution"],
            "conventional_mean_contribution_A": None if not e["signed_contribution"] else .5*float(e["signed_contribution"]),
            "gaussian_geodesic_contribution": aw_gaussian_by_edge.get(int(e["edge_id"])), "support_status": e["status"],
            "selected_surface_count": len(incident), "classification": "interior_selected" if len(incident)==2 else "boundary_selected" if len(incident)==1 else "nonmanifold"})
    write_csv(output / "aw_edge_curvature_components.csv", aw_edge_features)
    write_csv(output / "vertex_gaussian_components.csv",
        [{"vertex_id":r["vertex_id"],"generator_tuple":";".join(map(str,net.verts.loc[r["vertex_id"],"balls"])),"angle_defect_rad":r["angle_defect_rad"],"support_status":r["status"]} for r in aw_gaussian["vertex_rows"]] +
        [{"vertex_id":t,"generator_tuple":";".join(map(str,t)),"angle_defect_rad":g,"support_status":"supported"} for t,g in power_gaussian.get("vertex_terms",{}).items()],
        fields=["vertex_id","generator_tuple","angle_defect_rad","support_status"])
    write_csv(output / "boundary_curvature_components.csv",
        [{"boundary_id":r["edge_id"],"length_A":None,"geodesic_curvature":r["geodesic_gaussian"],"support_status":r["status"]}
         for r in aw_gaussian["edge_rows"] if r["selected_surface_count"]==1] +
        [{"boundary_id":"power_vertex:"+";".join(map(str,v)),"length_A":None,"geodesic_curvature":power_gaussian["vertex_terms"][v],"support_status":"polyhedral boundary vertex turn"}
         for v in power_gaussian.get("boundary_vertices",set()) if v in power_gaussian.get("vertex_terms",{})],
        fields=["boundary_id","length_A","geodesic_curvature","support_status"])

    # Resolve a local pair by quality criteria: shared, finite, supported, and
    # maximum of its smaller Power/AW patch areas.
    candidates = [r for r in comparison_rows if r["AW_status"] == "supported" and r["Power_status"] == "finite_bounded"]
    example = max(candidates, key=lambda r: min(float(r["Power_area_A2"]), float(r["AW_area_A2"])))
    ex_pair = tuple(sorted((int(example["generator_i"]), int(example["generator_j"]))))
    ex_sids = pos_by_pair.get(ex_pair, [])
    local_aw = _surface_mesh(net, ex_sids)
    if local_aw is not None: write_mesh(local_aw, "local_example_aw_surface.off", directory=output)
    ex_entry = area_cache.get(f"{ex_pair[0]}:{ex_pair[1]}", {})
    ep, ef = _power_faces(None, {f"{ex_pair[0]}:{ex_pair[1]}": ex_entry}, {ex_pair})
    _write_off(output / "local_example_power_facet.off", ep, ef)
    (output / "local_example.json").write_text(json.dumps({"stable_atom_i":example["stable_atom_i"],
        "stable_atom_j":example["stable_atom_j"], "generator_i":ex_pair[0], "generator_j":ex_pair[1],
        "Power_area_A2":float(example["Power_area_A2"]), "AW_area_A2":float(example["AW_area_A2"]),
        "AW_status":example["AW_status"], "selection":"shared Power/AW alpha-selected bicolor contact",
        "selection_rule":"largest min(Power area, AW area) among supported shared pairs"}, indent=2)+"\n", encoding="utf-8")

    # Solvent structure audit and exact protein-to-static sequence crosswalk.
    solvent = audit_solvated()
    (output / "solvated_input_audit.json").write_text(json.dumps(solvent, indent=2)+"\n", encoding="utf-8")
    mdrows = _pdb_rows(SOLVATED_PDB)
    write_csv(output / "solvated_atom_audit.csv", [{k: v for k,v in x.items() if k not in {"line","xyz"}} |
        {"x_A":x["xyz"][0],"y_A":x["xyz"][1],"z_A":x["xyz"][2]} for x in mdrows])

    # Comparison rows contain only certified partials. Missing geometric terms
    # stay blank and status tells downstream plotting not to replace them.
    metrics = { (r["comparison"],r["scheme"]):r for r in read_csv(MATCHED/"total_interface_metrics.csv") }
    table = []
    for representation in ("alpha_contact_A", "alpha_contact_B", "alpha_contact_combined"):
        table.append({"representation":representation,"certified_mean_total":False,"certified_gaussian_total":False,
            "status":"UNRESOLVED: no exact contact-patch decomposition of the union-of-balls boundary is implemented"})
    for label, key in (("static_power_alpha_interface",("total_restricted_interface","power")),
                       ("static_aw_alpha_interface",("total_restricted_interface","aw"))):
        r=metrics[key]
        chosen_pairs = power_selected if "power" in label else aw_selected
        source_surfs = power_surface_features if "power" in label else aw_surface_features
        chosen_atoms = {int(i) for pair in chosen_pairs for i in pair}
        components = None
        if "power" in label:
            components = _components(chosen_pairs, [e.bicolor_facets for e in power.edges])
        else:
            surface_pair_by_id = {int(sid): tuple(sorted(map(int, row.get("balls", ()))))
                                  for sid, row in net.surfs.iterrows() if len(row.get("balls", ())) == 2}
            chosen_surface_ids = {sid for sid, pair in surface_pair_by_id.items() if pair in chosen_pairs}
            edge_pair_map = []
            for _, edge_row in net.edges.iterrows():
                incident = [surface_pair_by_id[int(sid)] for sid in edge_row.get("surfs", ())
                            if int(sid) in chosen_surface_ids]
                if incident: edge_pair_map.append(incident)
            components = _components(chosen_pairs, edge_pair_map)
        edge_rows_for_scheme = power_edges if "power" in label else aw_edge_features
        edge_known_len = sum(float(e["length_A"]) for e in edge_rows_for_scheme
                             if e.get("length_A") not in (None, "") and e.get("classification") == "interior_selected"
                             and e.get("raw_signed_contribution_A_rad") not in (None, ""))
        edge_interior_len = sum(float(e["length_A"]) for e in edge_rows_for_scheme
                                if e.get("length_A") not in (None, "") and e.get("classification") == "interior_selected")
        edge_boundary_len = float(r["boundary_length_A"])
        edge_coverage = edge_known_len/(edge_interior_len+edge_boundary_len) if edge_interior_len+edge_boundary_len else None
        smooth_partial = (sum(float(x["integral_H_dA_A"]) for x in source_surfs
                               if x.get("integral_H_dA_A") not in (None,"") and x.get("support_status") == "included") if "aw" in label else 0.0)
        smooth_known_area = (sum(float(x["area_A2"]) for x in source_surfs
                                 if x.get("integral_H_dA_A") not in (None,"") and x.get("area_A2") not in (None,"") and x.get("support_status") == "included") if "aw" in label else float(r["area_A2"]))
        smooth_coverage = smooth_known_area/float(r["area_A2"]) if r["area_A2"] else None
        partial_combined = (smooth_partial+float(r["conventional_edge_H_A"]) if "aw" in label and r["conventional_edge_H_A"] not in (None,"")
                            else r["combined_H_A"])
        g=power_gaussian if "power" in label else aw_gaussian
        g_smooth=float(g["smooth"])
        g_edge=float(g.get("edge_interior",0.0))
        g_vertex=float(g.get("vertex",0.0))
        g_boundary=float(g.get("boundary",0.0))+(float(g.get("edge_boundary",0.0)) if "aw" in label else 0.0)
        g_smooth_coverage=(sum(float(x["area_A2"]) for x in source_surfs if x.get("integral_K_dA") not in (None,"") and x.get("support_status")=="included")/float(r["area_A2"]) if "aw" in label and r["area_A2"] else (1.0 if "power" in label else None))
        gaussian_coverage=(1.0 if g.get("total") is not None else g_smooth_coverage)
        table.append({"representation":label,"area_A2":float(r["area_A2"]),"smooth_mean_A":smooth_partial,
            "edge_mean_A":r["conventional_edge_H_A"],"integrated_mean_A":partial_combined,
            "smooth_gaussian":g_smooth,"edge_gaussian":g_edge,"vertex_gaussian":g_vertex,"boundary_gaussian":g_boundary,
            "integrated_gaussian":g.get("total"),"euler_characteristic":g.get("chi"),"components":components,"boundary_loops":g.get("boundary_loops"),
            "boundary_length_A":float(r["boundary_length_A"]),"pair_count":len(chosen_pairs),
            "interface_atom_count":len(chosen_atoms),"components":components,
            "smooth_mean_coverage_area":smooth_coverage,"edge_mean_coverage_length":edge_coverage,
            "mean_coverage":min(x for x in (smooth_coverage,edge_coverage) if x is not None) if any(x is not None for x in (smooth_coverage,edge_coverage)) else None,
            "gaussian_smooth_coverage_area":g_smooth_coverage,"gaussian_coverage":gaussian_coverage,
            "unresolved_surfaces":int(r["unresolved_smooth_surfaces"]),"unresolved_edges":int(r["unresolved_edges"]),
            "unresolved_vertices":g.get("unresolved_vertices"),"certified_mean_total":False,
            "certified_gaussian_total":g.get("total") is not None,
            "gauss_bonnet_expected":g.get("expected"),"gauss_bonnet_calculated":g.get("calculated_total",g.get("partial_total",g.get("total"))),"gauss_bonnet_residual":g.get("residual"),
            "status":"PARTIAL mean; Gaussian "+str(g.get("status"))})
    for label in ("solvated_power_alpha_interface","solvated_aw_alpha_interface"):
        table.append({"representation":label,"certified_mean_total":False,"certified_gaussian_total":False,
            "status":"UNRESOLVED: available cached Power/AW logs contain protein generators only; water/ion-inclusive network not solved"})
    write_csv(output/"curvature_comparison.csv", table)
    write_csv(output/"gauss_bonnet_audit.csv", [
        {"representation":r["representation"],"euler_characteristic":r.get("euler_characteristic"),
         "expected_topological_total":r.get("gauss_bonnet_expected"),"calculated_total":r.get("gauss_bonnet_calculated"),"residual":r.get("gauss_bonnet_residual"),
         "status":r.get("status","UNRESOLVED")}
        for r in table])
    (output/"alpha_contact_patch_status.txt").write_text(
        "Construction audit: unresolved. Existing Power/AW alpha data select generator pairs and Voronoi surfaces, but do not supply the exact spherical arrangement required to partition boundary(A isolated), boundary(B isolated), and boundary(A+B). No entire atom spheres are exported as contact patches.\n"
        "Required geometry: per-generator spherical patches clipped by same-group and opposite-group union occlusion; resolve circle-arc arrangements and patch adjacency, then integrate smooth sphere curvature and boundary terms with a validated Gauss-Bonnet accounting.\n",encoding="utf-8")

    # Actual geometry projection figure (same PCA frame for both schemes).
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    aw_pts, aw_faces = _parse_off(output/"static_aw_alpha_physical.off")
    p_points=np.asarray(p_pts,float); aw_points=np.asarray(aw_pts,float)
    allpts=np.vstack([p_points,aw_points]); center=allpts.mean(axis=0)
    _u,_s,vt=np.linalg.svd(allpts-center,full_matrices=False); basis=vt[:2]
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for ax,pts,faces,title,color in ((axes[0],p_points,p_faces,"Power alpha-selected physical facets","#2962a3"),
                                     (axes[1],aw_points,aw_faces,"AW alpha-selected physical surfaces","#c24949")):
        xy,polys=_project_off(pts,faces,basis,center)
        ax.add_collection(PolyCollection(polys,facecolors=color,edgecolors="none",alpha=.88))
        ax.autoscale(); ax.set_aspect("equal"); ax.set_title(title); ax.set_xlabel("PCA 1 (Å)"); ax.set_ylabel("PCA 2 (Å)")
    fig.suptitle("2KAI matched-probe alpha-selected interfaces — projection of exported physical geometry")
    fig.savefig(output/"static_interface_projection.png",dpi=300); fig.savefig(output/"static_interface_projection.svg"); plt.close(fig)
    camera={"projection":"PCA of combined actual Power/AW physical surface mesh vertices",
            "basis_vectors_xyz":basis.tolist(),"center_xyz_A":center.tolist(),
            "same_projection_for_both_panels":True}
    (output/"camera_metadata.json").write_text(json.dumps(camera,indent=2)+"\n",encoding="utf-8")
    layers=[("static_power_interface","static_power_alpha_physical.off",[0.12,0.35,0.85]),
            ("static_aw_interface","static_aw_alpha_physical.off",[0.85,0.20,0.20]),
            ("shared_power_interface","shared_power_physical.off",[0.10,0.40,0.95]),
            ("shared_aw_interface","shared_aw_physical.off",[0.90,0.25,0.20]),
            ("local_power_facet","local_example_power_facet.off",[0.05,0.25,1.0]),
            ("local_aw_surface","local_example_aw_surface.off",[1.0,0.15,0.12])]
    _export_pymol(output,STATIC_PDB,"",layers,"2KAI real solved interface geometry")
    local = json.loads((output / "local_example.json").read_text(encoding="utf-8"))
    def _stable_parts(stable):
        chain, resnum, icode, resname, atom = stable.split("|")
        return chain, resnum, atom
    ci,ri,ai=_stable_parts(local["stable_atom_i"]); cj,rj,aj=_stable_parts(local["stable_atom_j"])
    local_pml=[f'run "{(output / "showcase_layers.py").resolve().as_posix()}"',
        "cmd.disable('all')", "cmd.enable('molecule')", "cmd.enable('local_power_facet')", "cmd.enable('local_aw_surface')",
        f"cmd.select('example_i', 'molecule and chain {ci} and resi {ri} and name {ai}')",
        f"cmd.select('example_j', 'molecule and chain {cj} and resi {rj} and name {aj}')",
        "cmd.select('example_neighbors', 'byres (molecule within 8 of (example_i or example_j))')",
        "cmd.show('sticks','example_neighbors')",
        "cmd.show('spheres','example_i or example_j')", "cmd.set('sphere_scale',0.45,'example_i or example_j')",
        "cmd.zoom('example_i or example_j', buffer=10)"]
    (output/"local_example.pml").write_text("\n".join(local_pml)+"\n",encoding="utf-8")
    for filename, active in (("static_power_alpha.pml","static_power_interface"),
                             ("static_aw_alpha.pml","static_aw_interface"),
                             ("static_overlay.pml",None)):
        commands=[f'run "{(output / "showcase_layers.py").resolve().as_posix()}"',"cmd.disable('all')","cmd.enable('molecule')"]
        if active: commands.append(f"cmd.enable('{active}')")
        else: commands += ["cmd.enable('static_power_interface')","cmd.enable('static_aw_interface')"]
        commands.append("cmd.orient('molecule')")
        (output/filename).write_text("\n".join(commands)+"\n",encoding="utf-8")
    md_pml=[f"cmd.load({str(SOLVATED_PDB.resolve())!r},'finite_solvent_context')",
        "cmd.hide('everything','finite_solvent_context')",
        "cmd.show('cartoon','finite_solvent_context and not resn SOL and not resn NA+CL')",
        "cmd.show('spheres','finite_solvent_context and resn SOL and name OW')",
        "cmd.set('sphere_scale',0.18,'finite_solvent_context and resn SOL and name OW')",
        "cmd.color('gray70','finite_solvent_context and resn SOL')",
        "cmd.orient('finite_solvent_context')"]
    (output/"finite_solvent_context.pml").write_text("\n".join(md_pml)+"\n",encoding="utf-8")

    # Real-data contact graph and solvent-coordinate projection, not fabricated
    # boundary patches or an alleged solvated Voronoi interface.
    gen = {int(r["generator_id"]):r for r in read_csv(MATCHED/"input_identity_audit.csv")}
    graph_rows=[]
    comp=read_csv(MATCHED/"bicolor_pair_comparison.csv")
    for r in comp:
        i,j=int(r["generator_i"]),int(r["generator_j"]); ai,aj=gen[i],gen[j]
        graph_rows.append({"pair_id":f"{i}:{j}","stable_atom_i":ai["stable_atom_id"],"stable_atom_j":aj["stable_atom_id"],
            "x_i_A":ai["x_A"],"y_i_A":ai["y_A"],"z_i_A":ai["z_A"],"x_j_A":aj["x_A"],"y_j_A":aj["y_A"],"z_j_A":aj["z_A"],
            "Power_selected":r["in_Power"],"AW_selected":r["in_AW"],"classification":r["classification"]})
    write_csv(output/"alpha_bicolor_contact_dual_edges.csv",graph_rows)
    fig,ax=plt.subplots(figsize=(7,6),constrained_layout=True)
    edge_colors={"shared":"#667788","Power_only":"#e68632","AW_only":"#9b59b6","AW_unresolved_for_Power_pair":"#d2ad32"}
    for row in graph_rows:
        ax.plot([float(row["x_i_A"]),float(row["x_j_A"])],[float(row["y_i_A"]),float(row["y_j_A"])],
                color=edge_colors.get(row["classification"],"#667788"),alpha=.5,lw=.55)
    atom_rows=list(gen.values())
    for chain,color in (("A","#3979b9"),("B","#3979b9"),("I","#c34e4e")):
        selected=[x for x in atom_rows if x["chain"]==chain]
        ax.scatter([float(x["x_A"]) for x in selected],[float(x["y_A"]) for x in selected],s=3,c=color,alpha=.35)
    ax.set_aspect("equal"); ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)")
    ax.set_title("2KAI alpha-selected bicolor dual contacts (real centers and connections; no ball patches)")
    fig.savefig(output/"alpha_dual_contacts_2d.png",dpi=300); fig.savefig(output/"alpha_dual_contacts_2d.svg"); plt.close(fig)
    mdxyz=np.asarray([x["xyz"] for x in mdrows]); mdcenter=mdxyz.mean(axis=0)
    _u,_s,mdvt=np.linalg.svd(mdxyz-mdcenter,full_matrices=False); mdxy=(mdxyz-mdcenter)@mdvt[:2].T
    protein_ix=[i for i,x in enumerate(mdrows) if x["resname"] in AA and x["element"]!="H"]
    water_ix=[i for i,x in enumerate(mdrows) if x["resname"] in {"SOL","HOH","WAT","TIP3","TIP3P"} and x["atom"].upper() in {"OW","O"}]
    ion_ix=[i for i,x in enumerate(mdrows) if x["resname"] in {"NA","CL"}]
    fig,ax=plt.subplots(figsize=(7,6),constrained_layout=True)
    ax.scatter(mdxy[water_ix,0],mdxy[water_ix,1],s=1,c="#8797a6",alpha=.24,label="water oxygen")
    ax.scatter(mdxy[protein_ix,0],mdxy[protein_ix,1],s=3,c="#444444",alpha=.55,label="protein heavy atoms")
    if ion_ix: ax.scatter(mdxy[ion_ix,0],mdxy[ion_ix,1],s=9,c="#e0a030",alpha=.8,label="ions")
    ax.set_aspect("equal"); ax.set_xlabel("MD PCA 1 (Å)"); ax.set_ylabel("MD PCA 2 (Å)")
    ax.set_title("2KAI finite explicit-solvent coordinates — context only"); ax.legend(markerscale=2)
    fig.savefig(output/"finite_solvent_context_2d.png",dpi=300); fig.savefig(output/"finite_solvent_context_2d.svg"); plt.close(fig)
    lookup={r["representation"]:r for r in table}
    display=[("Alpha A","alpha_contact_A"),("Alpha B","alpha_contact_B"),("Power","static_power_alpha_interface"),
             ("AW","static_aw_alpha_interface"),("Solvated Power","solvated_power_alpha_interface"),("Solvated AW","solvated_aw_alpha_interface")]
    table_text="Representation         Area A^2      int H (A)    int K       Status\n"+"-"*88+"\n"
    for name,key in display:
        r=lookup[key]
        def fmt(v):
            if v in (None,""): return "UNRESOLVED"
            try: return f"{float(v):.4f}"
            except (TypeError,ValueError): return str(v)
        table_text+=f"{name:<20} {fmt(r.get('area_A2')):>12} {fmt(r.get('integrated_mean_A')):>12} {fmt(r.get('integrated_gaussian')):>12} {r.get('status','UNRESOLVED')}\n"
    report=("2KAI curvature showcase checkpoint\n\n"+table_text+"\n"
        "The static Power/AW surfaces in OFF, PyMOL, and the 2D projection are the same physical geometries used by the matched-probe analysis. "
        "Power uses expanded radii r+1.4 and alpha=0; AW uses base radii and actual clipped pair births <=1.4 A.\n\n"
        f"Static area: Power {summary['total_interface_area_A2']['Power']:.6f} A^2; AW {summary['total_interface_area_A2']['AW']:.6f} A^2.\n"
        f"Static pair counts: Power {len(power_selected)}, AW {len(aw_selected)}.\n"
        f"Representative resolved shared contact: {example['stable_atom_i']} -- {example['stable_atom_j']} "
        f"(Power {float(example['Power_area_A2']):.4f} A^2; AW {float(example['AW_area_A2']):.4f} A^2).\n\n"
        "Alpha contact patches: unresolved. VorPy has alpha-selected dual pairs, but no exact exposed/ buried spherical-patch arrangement that partitions boundary(A), boundary(B), and boundary(A+B). Whole spheres were not substituted.\n"
        "Gaussian curvature: Power smooth/edge terms and polyhedral vertex/boundary angle terms are tabulated; AW smooth-face, oriented face-edge geodesic, and vertex-angle terms use existing VorPy geometry routines. Totals are certified only when incidence is complete and the Gauss-Bonnet residual is available; otherwise the total remains UNRESOLVED. See gauss_bonnet_audit.csv.\n"
        "Mean curvature: cached signed edge integrals are interior-seam partials. Selected boundary edges have no invented curvature term. AW also has unsupported smooth patches and unresolved edge quadratures; its integrated total remains unresolved.\n\n"
        f"Solvated input audit: {solvent['all_atom_records']} coordinate records: {solvent['protein_atoms']} protein atoms "
        f"({solvent['protein_heavy_atoms']} heavy, {solvent['protein_hydrogens']} H), {solvent['water_atoms']} water atoms, "
        f"{solvent['ion_or_other_atoms']} ion/other atoms. Chain field is blank; ordered protein sequence maps to static A/B/I={solvent['protein_chain_sequence_matches_static_A_B_I']}.\n"
        "The existing cached Power/AW logs have 4,340 protein generators and omit all 40,086 water atoms and the ion/other atoms. Therefore the requested finite explicit-solvent protein-protein geometry, its curvature, and static-vs-solvated pair comparison are unresolved; the cached protein-only network is not relabeled solvated.\n\n"
        "PyMOL is unavailable in this environment; PML launchers and physical OFF meshes were emitted for local rendering. The supplied 2KAI_MD_centered coordinates represent a distinct MD conformation, so they were not overlaid with static-crystal interface patches.\n")
    (output/"curvature_comparison_report.txt").write_text(report,encoding="utf-8")
    (output/"showcase_metadata.json").write_text(json.dumps({"system":"2KAI","input_geometry":"static solved H-tier plus separate centered MD PDB audit",
        "probe_A":PROBE,"group_A":"chains A+B","group_B":"chain I","static_geometry_source":str(NETWORK_PATH),
        "input_audit":solvent,"example_pair":json.loads((output/"local_example.json").read_text()),
        "gaussian_curvature_status":{r["representation"]:r.get("status") for r in table},
        "alpha_contact_patch_status":"unresolved; no exact union-of-balls patch arrangement exists in current pipeline",
        "solvated_interface_status":"unresolved; full finite solvent network not present in cached geometry"},indent=2)+"\n",encoding="utf-8")
    return output


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",default=str(ROOT))
    args=p.parse_args(argv)
    out=run(args.output_dir); print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
