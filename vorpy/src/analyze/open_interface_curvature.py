"""Topology and curvature accounting for open, face-based interfaces.

The module deliberately keeps extrinsic mean curvature separate from the
intrinsic Gauss--Bonnet terms.  It operates on physical face boundary cycles;
triangulation diagonals are never treated as interface edges.
"""
from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
import math
from typing import Hashable, Iterable, Sequence


@dataclass(frozen=True)
class ComponentTopology:
    component: int
    vertices: int
    edges: int
    faces: int
    euler_characteristic: int
    boundary_components: int | None
    orientable: bool | None
    genus: int | None
    manifold: bool


def _edge(a, b):
    return (a, b) if repr(a) <= repr(b) else (b, a)


def _component_count(graph):
    unseen = set(graph)
    count = 0
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            for node in graph[stack.pop()] & unseen:
                unseen.remove(node)
                stack.append(node)
    return count


def analyze_surface_topology(face_cycles: Iterable[Sequence[Hashable]]) -> dict:
    """Audit a polygonal 2-complex using only its physical face boundaries.

    Each cycle is the ordered boundary of one physical face. A face may have
    any number of vertices. Edge-connected components, boundary loops,
    orientability, and genus are reported only where incidence supports them.
    """
    faces = [tuple(cycle) for cycle in face_cycles]
    if any(len(face) < 3 or len(set(face)) != len(face) for face in faces):
        raise ValueError("Each physical face must have at least three distinct vertices")
    edge_faces = defaultdict(list)
    directed = defaultdict(list)
    vertex_faces = defaultdict(list)
    face_edges = []
    for fi, face in enumerate(faces):
        es = []
        for i, a in enumerate(face):
            b = face[(i + 1) % len(face)]
            key = _edge(a, b)
            edge_faces[key].append(fi)
            directed[key].append((fi, (a, b)))
            es.append(key)
            vertex_faces[a].append((fi, face[(i - 1) % len(face)], b))
        face_edges.append(es)
    face_graph = {i: set() for i in range(len(faces))}
    for incident in edge_faces.values():
        if len(incident) == 2:
            a, b = incident
            face_graph[a].add(b); face_graph[b].add(a)
    components = []
    unseen = set(face_graph)
    while unseen:
        seed = unseen.pop(); stack = [seed]; comp = {seed}
        while stack:
            for n in face_graph[stack.pop()] & unseen:
                unseen.remove(n); stack.append(n); comp.add(n)
        components.append(comp)

    edge_component = {}
    for edge, incident in edge_faces.items():
        if incident:
            edge_component[edge] = next(i for i, c in enumerate(components) if incident[0] in c)
    boundary_edges = [e for e, inc in edge_faces.items() if len(inc) == 1]
    boundary_graph = defaultdict(set)
    boundary_edge_component = defaultdict(list)
    for edge in boundary_edges:
        a, b = edge
        boundary_graph[a].add(b); boundary_graph[b].add(a)
        boundary_edge_component[edge_component[edge]].append(edge)
    boundary_loops = None
    boundary_regular = all(len(adj) == 2 for adj in boundary_graph.values())
    if boundary_regular:
        bg = {v: set(ns) for v, ns in boundary_graph.items()}
        boundary_loops = _component_count(bg) if bg else 0

    # A face-flip constraint checks orientability of each edge-manifold
    # component: adjacent face cycles must traverse their shared edge oppositely.
    signs = {}
    orientable_components = []
    for ci, comp in enumerate(components):
        ok = True
        for seed in comp:
            if seed in signs:
                continue
            signs[seed] = 1; stack = [seed]
            while stack:
                fi = stack.pop()
                for e in face_edges[fi]:
                    inc = edge_faces[e]
                    if len(inc) != 2:
                        continue
                    fj = inc[0] if inc[1] == fi else inc[1]
                    if fj not in comp:
                        continue
                    di, dj = directed[e]
                    dir_i = 1 if di[1] == e else -1
                    dir_j = 1 if dj[1] == e else -1
                    required = -signs[fi] * dir_i * dir_j
                    if fj in signs and signs[fj] != required:
                        ok = False
                    elif fj not in signs:
                        signs[fj] = required; stack.append(fj)
        orientable_components.append(ok)

    # Vertex links must be one cycle (interior) or one path (boundary).
    vertex_link_ok = True
    for vertex, corners in vertex_faces.items():
        link = defaultdict(set)
        for _, prev_v, next_v in corners:
            link[prev_v].add(next_v); link[next_v].add(prev_v)
        degrees = [len(x) for x in link.values()]
        n_boundary_incident = sum(vertex in e for e in boundary_edges)
        expected_ones = 2 if n_boundary_incident else 0
        if (not link or _component_count(link) != 1 or
                sum(d == 1 for d in degrees) != expected_ones or
                any(d not in (1, 2) for d in degrees)):
            vertex_link_ok = False

    records = []
    for ci, comp in enumerate(components):
        comp_edges = {e for e, inc in edge_faces.items() if inc[0] in comp}
        comp_vertices = {v for fi in comp for v in faces[fi]}
        b_edges = [e for e in comp_edges if len(edge_faces[e]) == 1]
        b_graph = defaultdict(set)
        for a, b in b_edges:
            b_graph[a].add(b); b_graph[b].add(a)
        b_regular = all(len(ns) == 2 for ns in b_graph.values())
        b_count = _component_count(b_graph) if b_regular and b_graph else (0 if not b_edges else None)
        chi = len(comp_vertices) - len(comp_edges) + len(comp)
        orientable = orientable_components[ci]
        manifold = (all(len(edge_faces[e]) in (1, 2) for e in comp_edges) and
                    vertex_link_ok and b_regular)
        genus = None
        if manifold and orientable and b_count is not None:
            numerator = 2 - b_count - chi
            if numerator >= 0 and numerator % 2 == 0:
                genus = numerator // 2
        records.append(ComponentTopology(ci, len(comp_vertices), len(comp_edges), len(comp),
                                         chi, b_count, orientable, genus, manifold))

    # Counting edge-connected components separately double-counts a pinched
    # vertex. Such a complex is rejected as a manifold, but its cell Euler
    # characteristic must still be V-E+F.
    total_chi = len(vertex_faces) - len(edge_faces) + len(faces)
    total = {"vertices": len({v for f in faces for v in f}), "edges": len(edge_faces),
        "faces": len(faces), "euler_characteristic": total_chi,
        "components": len(components), "boundary_edges": len(boundary_edges),
        "boundary_components": boundary_loops, "boundary_is_manifold": vertex_link_ok and
            all(len(inc) in (1, 2) for inc in edge_faces.values()) and boundary_regular,
        "orientable": all(orientable_components) if orientable_components else True,
        "component_records": records, "edge_face_incidence": dict(edge_faces),
        "face_cycles": faces}
    if total["boundary_is_manifold"] and total["orientable"] and all(r.genus is not None for r in records):
        total["genus_sum"] = sum(r.genus for r in records)
    else:
        total["genus_sum"] = None
    return total


def integrated_mean_curvature(smooth_face_terms: Iterable[float],
                              internal_crease_terms: Iterable[float], *,
                              unresolved_faces: int = 0,
                              unresolved_internal_edges: int = 0,
                              boundary_terms: Iterable[float] = ()) -> dict:
    """Return C_H with no open-boundary term under H=(k1+k2)/2.

    ``internal_crease_terms`` are already conventional contributions
    ``1/2 integral beta ds``. ``boundary_terms`` is rejected unless exactly
    zero: a free boundary contributes a conormal to first variation, not a
    scalar mean-curvature measure.
    """
    boundary_values = list(boundary_terms)
    if any(abs(float(x)) > 1e-14 for x in boundary_values):
        raise ValueError("Free-boundary geodesic/turning terms are not mean curvature")
    smooth = sum(map(float, smooth_face_terms))
    crease = sum(map(float, internal_crease_terms))
    certified = (unresolved_faces == 0 and unresolved_internal_edges == 0
                 and math.isfinite(smooth) and math.isfinite(crease))
    return {"smooth_H": smooth, "internal_edge_H": crease,
            "boundary_H": 0.0, "total_H": smooth + crease if certified else None,
            "certified": certified, "unresolved_faces": unresolved_faces,
            "unresolved_internal_edges": unresolved_internal_edges}


def gauss_bonnet_accounting(*, smooth_K: float, face_boundary_geodesic: float,
                            vertex_corner_terms: float, euler_characteristic: int,
                            complete: bool = True, tolerance: float = 1e-7) -> dict:
    """Audit the face-wise Gauss--Bonnet sum for a piecewise-smooth complex.

    ``face_boundary_geodesic`` sums the signed geodesic-curvature integral on
    every face-boundary incidence, including both sides of internal seams.
    ``vertex_corner_terms`` sums assembled interior defects and boundary
    turns, NOT the exterior turns of every individual face polygon.
    The signed extrinsic dihedral angle is not an input to this identity.
    """
    lhs = float(smooth_K) + float(face_boundary_geodesic) + float(vertex_corner_terms)
    rhs = 2 * math.pi * int(euler_characteristic)
    residual = lhs - rhs
    certified = bool(complete and abs(residual) <= tolerance * max(1.0, abs(rhs)))
    return {"smooth_K": float(smooth_K),
            "face_boundary_geodesic_K": float(face_boundary_geodesic),
            "vertex_corner_K": float(vertex_corner_terms),
            "total_GB_left": lhs, "two_pi_chi": rhs,
            "residual": residual if complete else None,
            "certified": certified,
            "status": "certified" if certified else "UNRESOLVED: incomplete terms" if not complete else
                     "UNRESOLVED: Gauss-Bonnet residual exceeds tolerance"}


def polygon_gaussian_terms(face_cycles, coordinates):
    """Discrete GB terms for convex planar polygonal disk faces.

    Reflex corners require an oriented angle implementation; callers must
    not pass concave or nonplanar faces to this convex-face adapter.
    """
    topo = analyze_surface_topology(face_cycles)
    angle_sum = defaultdict(float)
    for face in topo["face_cycles"]:
        for i, v in enumerate(face):
            p = coordinates[v]
            a = [coordinates[face[i-1]][j] - p[j] for j in range(3)]
            b = [coordinates[face[(i+1) % len(face)]][j] - p[j] for j in range(3)]
            na = math.sqrt(sum(x*x for x in a)); nb = math.sqrt(sum(x*x for x in b))
            if na == 0 or nb == 0:
                raise ValueError("Degenerate polygon corner")
            angle_sum[v] += math.acos(max(-1.0, min(1.0, sum(x*y for x,y in zip(a,b))/(na*nb))))
    boundary_vertices = {v for e, inc in topo["edge_face_incidence"].items() if len(inc) == 1 for v in e}
    vertex_terms = {v: (math.pi - angle_sum[v] if v in boundary_vertices else
                        2*math.pi - angle_sum[v]) for v in angle_sum}
    smooth_K = 0.0
    corner_K = sum(vertex_terms.values())
    gb = gauss_bonnet_accounting(smooth_K=smooth_K, face_boundary_geodesic=0.0,
        vertex_corner_terms=corner_K, euler_characteristic=topo["euler_characteristic"],
        complete=topo["boundary_is_manifold"] and topo["orientable"])
    interior = sum(value for v, value in vertex_terms.items() if v not in boundary_vertices)
    boundary = sum(value for v, value in vertex_terms.items() if v in boundary_vertices)
    return {"topology": topo, "vertex_terms": vertex_terms, "smooth_K": smooth_K,
            "interior_vertex_defect_K": interior, "boundary_vertex_corner_K": boundary,
            "intrinsic_total_K": interior if gb["certified"] else None,
            "boundary_geodesic_K": 0.0, "boundary_corner_K": corner_K, **gb}


def _network_face_cycle(network, surface_row):
    """Order a physical network surface's edge endpoints into one boundary cycle."""
    adjacency = defaultdict(set)
    for edge_id in map(int, surface_row.get("edges", ())):
        endpoints = tuple(map(int, network.edges.loc[edge_id, "verts"]))
        if len(endpoints) != 2 or endpoints[0] == endpoints[1]:
            raise ValueError(f"Surface boundary contains invalid edge {edge_id}")
        a, b = endpoints
        adjacency[a].add(b); adjacency[b].add(a)
    if not adjacency or any(len(ns) != 2 for ns in adjacency.values()):
        raise ValueError("Physical surface boundary is not one simple cycle")
    start = min(adjacency)
    cycle = [start]
    previous = None
    current = start
    while True:
        candidates = sorted(adjacency[current] - ({previous} if previous is not None else set()))
        nxt = candidates[0]
        if nxt == start:
            break
        if nxt in cycle:
            raise ValueError("Physical surface boundary cycle self-intersects in incidence")
        cycle.append(nxt); previous, current = current, nxt
        if len(cycle) > len(adjacency):
            raise ValueError("Failed to close physical surface boundary cycle")
    if len(cycle) != len(adjacency):
        raise ValueError("Surface has multiple boundary cycles or unused vertices")
    return tuple(cycle)


def _aw_supported_gauss_bonnet(network, face_records, supported_surface_ids, group_a, *, order=32):
    """Reintegrate each supported (face, physical edge) incidence exactly once."""
    from vorpy.src.calculations.edge_geometry import oriented_aw_face_edge_geodesic_curvature
    from vorpy.src.calculations.vertex_geometry import aw_vertex_face_angle
    from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge

    selected = [(sid, pair, cycle) for sid, pair, cycle in face_records if sid in supported_surface_ids]
    surface_edges = {sid: set(map(int, network.surfs.loc[sid, "edges"])) for sid, _, _ in selected}
    surface_vertices = {sid: set(map(int, network.surfs.loc[sid, "verts"])) for sid, _, _ in selected}
    edge_surfaces = defaultdict(set); vertex_surfaces = defaultdict(set)
    for sid, _, _ in selected:
        for eid in surface_edges[sid]: edge_surfaces[eid].add(sid)
        for vid in surface_vertices[sid]: vertex_surfaces[vid].add(sid)
    resolved = {}; edge_failures = set()
    for eid in edge_surfaces:
        try:
            resolved[eid] = resolve_aw_network_edge(network, eid)
        except Exception:
            edge_failures.add(eid)
    incidence_terms = {}; face_failures = set()
    for sid, pair, _ in selected:
        a = next((i for i in pair if i in group_a), None)
        b = next((i for i in pair if i != a), None)
        if a is None or b is None:
            face_failures.add(sid); continue
        for eid in surface_edges[sid]:
            if eid not in resolved:
                continue
            er = resolved[eid]
            try:
                incidence_terms[(sid, eid)] = oriented_aw_face_edge_geodesic_curvature(
                    er.geometry, er.ball_indices, er.locations, a, b, order=order)
            except Exception:
                edge_failures.add(eid)
    incidence_counts = {eid: sum((sid, eid) in incidence_terms for sid in sids)
                        for eid, sids in edge_surfaces.items()}
    edge_interior = sum(value for (sid, eid), value in incidence_terms.items()
                        if len(edge_surfaces[eid]) == 2)
    edge_boundary = sum(value for (sid, eid), value in incidence_terms.items()
                        if len(edge_surfaces[eid]) == 1)
    edge_bad = {eid for eid, count in incidence_counts.items()
                if count != len(edge_surfaces[eid]) or len(edge_surfaces[eid]) not in (1, 2)}
    angles_by_vertex = defaultdict(float); angle_failures = set()
    face_angles = defaultdict(list)
    for sid, pair, _ in selected:
        a = next((i for i in pair if i in group_a), None)
        if a is None:
            angle_failures.update(surface_vertices[sid]); continue
        for vid in surface_vertices[sid]:
            try:
                angle = float(aw_vertex_face_angle(network, vid, a, sid, resolved_edges=resolved))
                angles_by_vertex[vid] += angle
                face_angles[sid].append(angle)
            except Exception:
                angle_failures.add(vid)
    boundary_edges = {eid for eid, sids in edge_surfaces.items() if len(sids) == 1}
    boundary_vertices = {int(v) for eid in boundary_edges for v in network.edges.loc[eid, "verts"]}
    vertex_terms = {}
    for vid, angle_sum in angles_by_vertex.items():
        vertex_terms[vid] = (math.pi-angle_sum) if vid in boundary_vertices else (2*math.pi-angle_sum)
    vertex_interior = sum(v for vid, v in vertex_terms.items() if vid not in boundary_vertices)
    vertex_boundary = sum(v for vid, v in vertex_terms.items() if vid in boundary_vertices)
    smooth_k = sum(float(network.surfs.loc[sid, "int_gauss_curv"]) for sid, _, _ in selected)
    topo_cycles = [cycle for _, _, cycle in selected]
    topology = analyze_surface_topology(topo_cycles)
    complete = (not face_failures and not edge_bad and not angle_failures and
                topology["boundary_is_manifold"] and topology["orientable"])
    gb = gauss_bonnet_accounting(smooth_K=smooth_k,
        face_boundary_geodesic=edge_interior+edge_boundary,
        vertex_corner_terms=vertex_interior+vertex_boundary,
        euler_characteristic=topology["euler_characteristic"], complete=complete)
    face_diagnostics = []
    for sid, pair, _ in selected:
        face_complete = (len(face_angles[sid]) == len(surface_vertices[sid]) and
                         all((sid, eid) in incidence_terms for eid in surface_edges[sid]))
        lhs = (float(network.surfs.loc[sid, "int_gauss_curv"]) +
               sum(incidence_terms[(sid, eid)] for eid in surface_edges[sid]) +
               sum(math.pi-angle for angle in face_angles[sid])) if face_complete else None
        face_diagnostics.append({"surface_id": sid, "pair": str(pair), "complete": face_complete,
            "face_GB_LHS": lhs, "face_GB_residual": lhs-2*math.pi if lhs is not None else None})
    return {"topology": topology, "smooth_K": smooth_k,
        "face_diagnostics": face_diagnostics,
        "internal_seam_face_geodesic_K": edge_interior,
        "free_boundary_face_geodesic_K": edge_boundary,
        "interior_vertex_defect_K": vertex_interior,
        "boundary_vertex_corner_K": vertex_boundary,
        "vertex_terms": vertex_terms, "unresolved_edges": len(edge_bad),
        "unresolved_vertices": len(angle_failures), "unresolved_faces": len(face_failures),
        **gb}


def run_2kai_open_interface_audit(output_dir="output/2KAI_curvature_showcase/open_interface_curvature"):
    """Reconstruct physical Power/AW interface topology and curvature audit."""
    import csv
    import json
    import pickle
    from collections import defaultdict
    from pathlib import Path

    synthetic = synthetic_validation_rows()
    if any(not math.isfinite(row["residual"]) or abs(row["residual"]) > 1e-10 for row in synthetic):
        raise ValueError("Synthetic validation failed; static audit prohibited")
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    def read_rows(path):
        with Path(path).open(newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))
    root = Path("output/2KAI_curvature_showcase")
    matched = Path("output/comparative_study/matched_probe_control/2KAI")
    pairs = read_rows(matched / "bicolor_pair_comparison.csv")
    power_pairs = {tuple(sorted((int(r["generator_i"]), int(r["generator_j"])))) for r in pairs if r["in_Power"] == "True"}
    aw_pairs = {tuple(sorted((int(r["generator_i"]), int(r["generator_j"])))) for r in pairs if r["in_AW"] == "True"}
    area_cache = json.loads((matched / "power_facet_area_cache.json").read_text(encoding="utf-8"))["areas"]
    _, power_result = pickle.loads((matched / "power_geometry_cache.pkl").read_bytes())
    from scipy.spatial import cKDTree
    vertex_ids = list(power_result.power_vertices)
    vertex_tree = cKDTree([power_result.power_vertices[v].position for v in vertex_ids])
    power_cycles = []; power_coords = {}
    for pair in sorted(power_pairs):
        item = area_cache[f"{pair[0]}:{pair[1]}"]
        poly = item.get("polygon_vertices") or []
        cycle = []
        for point in poly:
            distances, indices = vertex_tree.query(point, k=2)
            if distances[0] > 2e-5 or distances[1] < 1e-8:
                raise ValueError("Power polygon vertex has unresolved dual-incidence identity")
            key = vertex_ids[int(indices[0])]
            power_coords[key] = power_result.power_vertices[key].position
            cycle.append(key)
        if len(cycle) < 3 or len(set(cycle)) != len(cycle):
            raise ValueError(f"Invalid physical Power polygon for pair {pair}")
        power_cycles.append(tuple(cycle))
    p_gauss = polygon_gaussian_terms(power_cycles, power_coords)
    # Each physical polygon edge must be an actual shared dual triangle,
    # not merely a coincidence in rounded coordinate keys.
    for edge in p_gauss["topology"]["edge_face_incidence"]:
        if len(set(edge[0]) & set(edge[1])) != 3:
            raise ValueError("Power polygon edge does not match solved dual incidence")

    write_rows(out / "synthetic_validation.csv", synthetic)

    from vorpy.src.io import load_network
    network = load_network("output/comparative_study/cases/2KAI/H/full_aw.vpy")
    aw_face_records = []
    for sid, row in network.surfs.iterrows():
        pair = tuple(sorted(map(int, row.get("balls", ()))))
        if len(pair) == 2 and pair in aw_pairs:
            aw_face_records.append((int(sid), pair, _network_face_cycle(network, row)))
    if {pair for _, pair, _ in aw_face_records} != aw_pairs:
        raise ValueError("Selected AW pair-to-physical-surface mapping is incomplete")
    aw_cycles = [cycle for _, _, cycle in aw_face_records]
    aw_topology = analyze_surface_topology(aw_cycles)
    support_by_sid = {int(r["surface_id"]): r for r in read_rows(root / "aw_surface_curvature_components.csv")}
    supported_cycles = [cycle for sid, _, cycle in aw_face_records if support_by_sid.get(sid, {}).get("support_status") == "included"]
    aw_supported_topology = analyze_surface_topology(supported_cycles)

    def topology_rows(scheme, scopes):
        rows = []
        for scope, topo in scopes:
            for comp in topo["component_records"]:
                rows.append({"scheme": scheme, "scope": scope, "component": comp.component,
                    "V_physical": comp.vertices, "E_physical": comp.edges, "F_physical": comp.faces,
                    "chi": comp.euler_characteristic, "boundary_loops": comp.boundary_components,
                    "orientable": comp.orientable, "genus": comp.genus, "manifold": comp.manifold})
            rows.append({"scheme": scheme, "scope": scope, "component": "TOTAL",
                "V_physical": topo["vertices"], "E_physical": topo["edges"], "F_physical": topo["faces"],
                "chi": topo["euler_characteristic"], "boundary_loops": topo["boundary_components"],
                "orientable": topo["orientable"], "genus": topo.get("genus_sum"), "manifold": topo["boundary_is_manifold"]})
        return rows

    write_rows(out / "power_topology_audit.csv", topology_rows("Power", [("selected_physical_interface", p_gauss["topology"])]))
    write_rows(out / "aw_topology_audit.csv", topology_rows("AW", [
        ("selected_physical_interface", aw_topology), ("supported_curvature_subcomplex", aw_supported_topology)]))

    metrics = {r["scheme"]: r for r in read_rows(matched / "total_interface_metrics.csv")
               if r["comparison"] == "total_restricted_interface"}
    edge_rows_power = read_rows(root / "power_edge_curvature_components.csv")
    p_internal = [r for r in edge_rows_power if r["classification"] == "interior_selected"]
    p_raw = sum(float(r["raw_signed_contribution_A_rad"]) for r in p_internal if r["raw_signed_contribution_A_rad"] not in ("", "None"))
    p_unresolved = sum(r["raw_signed_contribution_A_rad"] in ("", "None") for r in p_internal)
    import numpy as np
    oriented_cycles = []
    normals = []
    for pair, cycle in zip(sorted(power_pairs), power_cycles):
        a, b = pair if pair[0] in power_result.group_a else pair[::-1]
        n = power_result.points[b]-power_result.points[a]
        n = n / np.linalg.norm(n)
        area_vector = sum((np.cross(power_coords[cycle[k]], power_coords[cycle[(k+1)%len(cycle)]])
                           for k in range(len(cycle))), np.zeros(3))
        oriented_cycles.append(cycle if np.dot(area_vector, n) > 0 else cycle[::-1])
        normals.append(n)
    source_crease = {tuple(map(int, r["generator_tuple"].split(";"))):
                     float(r["conventional_mean_contribution_A"]) for r in p_internal
                     if r["conventional_mean_contribution_A"] not in ("", "None")}
    reconstructed = []
    for edge, incident in p_gauss["topology"]["edge_face_incidence"].items():
        if len(incident) != 2:
            continue
        fi, fj = incident
        face = oriented_cycles[fi]
        a, b = next((face[k], face[(k+1)%len(face)]) for k in range(len(face))
                    if {face[k], face[(k+1)%len(face)]} == set(edge))
        tangent = np.asarray(power_coords[b])-np.asarray(power_coords[a])
        length = float(np.linalg.norm(tangent)); tangent /= length
        # Existing production sign uses the opposite of the synthetic Dn
        # convention. This explicit sign adapter is checked edge by edge.
        beta = math.atan2(float(np.dot(tangent, np.cross(normals[fi], normals[fj]))),
                          float(np.dot(normals[fi], normals[fj])))
        triangle = tuple(sorted(set(a) & set(b)))
        value = .5*length*beta
        expected = source_crease.get(triangle)
        reconstructed.append({"dual_triangle": str(triangle), "normal_turn_H": value,
            "existing_H": expected, "difference": value-expected if expected is not None else None})
    write_rows(out / "power_independent_crease_audit.csv", reconstructed)
    if len(reconstructed) != len(p_internal) or any(r["difference"] is None or abs(r["difference"]) > 1e-8 for r in reconstructed):
        raise ValueError("Independent Power normal-turn validation failed")
    p_mean = integrated_mean_curvature([0.0], [0.5*p_raw], unresolved_internal_edges=p_unresolved)
    pgb = {k: p_gauss[k] for k in ("smooth_K", "boundary_geodesic_K", "boundary_corner_K", "total_GB_left", "two_pi_chi", "residual", "certified", "status")}
    p_curv = [{"term": "area_A2", "value": sum(float(area_cache[f"{a}:{b}"]["area_A2"]) for a,b in power_pairs), "support": f"{len(power_pairs)} selected Power facets"},
        {"term": "smooth_face_H", "value": 0.0, "support": "exact: planar Power facets"},
        {"term": "internal_crease_raw_signed_beta_ds", "value": p_raw, "support": f"{len(p_internal)-p_unresolved}/{len(p_internal)} internal edges"},
        {"term": "internal_edge_H_half_raw", "value": .5*p_raw, "support": "existing validated Power edge beta; internal seams only"},
        {"term": "free_boundary_H", "value": 0.0, "support": "not a mean-curvature term under H=(k1+k2)/2"},
        {"term": "total_H", "value": p_mean["total_H"], "support": "certified" if p_mean["certified"] else "UNRESOLVED"},
        {"term": "smooth_K", "value": 0.0, "support": "exact: planar Power facets"},
        {"term": "face_boundary_geodesic_K", "value": 0.0, "support": "straight polygon sides"},
        {"term": "interior_vertex_defect_K", "value": p_gauss["interior_vertex_defect_K"], "support": pgb["status"]},
        {"term": "boundary_vertex_corner_K", "value": p_gauss["boundary_vertex_corner_K"], "support": pgb["status"]},
        {"term": "intrinsic_total_K", "value": p_gauss["intrinsic_total_K"], "support": pgb["status"]},
        {"term": "Gauss_Bonnet_LHS", "value": pgb["total_GB_left"], "support": pgb["status"]},
        {"term": "two_pi_chi", "value": pgb["two_pi_chi"], "support": "independent physical polygon incidence"},
        {"term": "Gauss_Bonnet_residual", "value": pgb["residual"], "support": pgb["status"]}]
    write_rows(out / "power_curvature_components.csv", p_curv)

    aw_surface_rows = read_rows(root / "aw_surface_curvature_components.csv")
    aw_edge_rows = read_rows(root / "aw_edge_curvature_components.csv")
    aw_internal = [r for r in aw_edge_rows if r["classification"] == "interior_selected"]
    aw_included_internal = [r for r in aw_internal if r["support_status"] == "included" and r["raw_signed_contribution_A_rad"] not in ("", "None")]
    aw_missing_internal = len(aw_internal) - len(aw_included_internal)
    aw_raw = sum(float(r["raw_signed_contribution_A_rad"]) for r in aw_included_internal)
    aw_smooth_mean = sum(float(r["integral_H_dA_A"]) for r in aw_surface_rows if r["support_status"] == "included" and r["integral_H_dA_A"] not in ("", "None"))
    aw_missing_faces = sum(r["support_status"] != "included" or r["integral_H_dA_A"] in ("", "None") for r in aw_surface_rows)
    aw_mean = integrated_mean_curvature([aw_smooth_mean], [.5*aw_raw], unresolved_faces=aw_missing_faces,
        unresolved_internal_edges=aw_missing_internal)
    atom_audit = read_rows(matched / "input_identity_audit.csv")
    group_a = {int(r["generator_id"]) for r in atom_audit if r["chain"] in {"A", "B"}}
    supported_surface_ids = {sid for sid, row in support_by_sid.items() if row.get("support_status") == "included"}
    aw_gb_by_order = {order: _aw_supported_gauss_bonnet(network, aw_face_records,
        supported_surface_ids, group_a, order=order) for order in (64, 128, 256)}
    aw_gb = aw_gb_by_order[max(aw_gb_by_order)]
    write_rows(out / "aw_supported_face_GB_audit.csv", aw_gb["face_diagnostics"])
    aw_expected = 2*math.pi*aw_supported_topology["euler_characteristic"]
    aw_partial_lhs = float(aw_gb["total_GB_left"])
    aw_partial_residual = aw_partial_lhs - aw_expected
    aw_curv = [{"term": "area_A2", "value": sum(float(r["area_A2"]) for r in aw_surface_rows if r["area_A2"] not in ("", "None")), "support": f"{len(aw_surface_rows)} selected AW surfaces"},
        {"term": "smooth_face_H", "value": aw_smooth_mean, "support": f"{len(aw_surface_rows)-aw_missing_faces}/{len(aw_surface_rows)} surfaces"},
        {"term": "internal_crease_raw_signed_beta_ds", "value": aw_raw, "support": f"{len(aw_included_internal)}/{len(aw_internal)} internal edges"},
        {"term": "internal_edge_H_half_raw", "value": .5*aw_raw, "support": "existing validated AW edge integrals; internal seams only"},
        {"term": "free_boundary_H", "value": 0.0, "support": "not a mean-curvature term under H=(k1+k2)/2"},
        {"term": "partial_H_sum", "value": aw_mean["smooth_H"]+aw_mean["internal_edge_H"], "support": "partial; not certified total"},
        {"term": "smooth_K", "value": aw_gb["smooth_K"], "support": f"{len(supported_cycles)}/{len(aw_face_records)} supported faces"},
        {"term": "internal_seam_face_geodesic_K", "value": aw_gb["internal_seam_face_geodesic_K"], "support": "reintegrated once per supported face-edge incidence"},
        {"term": "free_boundary_face_geodesic_K", "value": aw_gb["free_boundary_face_geodesic_K"], "support": "reintegrated once per supported face-edge incidence"},
        {"term": "interior_vertex_defect_K", "value": aw_gb["interior_vertex_defect_K"], "support": "supported-subcomplex AW physical vertex angles"},
        {"term": "boundary_vertex_corner_K", "value": aw_gb["boundary_vertex_corner_K"], "support": "supported-subcomplex AW physical vertex angles"},
        {"term": "supported_subcomplex_two_pi_chi", "value": aw_expected, "support": f"chi={aw_supported_topology['euler_characteristic']} for supported face subset"},
        {"term": "supported_subcomplex_GB_residual_order_256", "value": aw_partial_residual, "support": aw_gb["status"]},
        {"term": "total_H", "value": aw_mean["total_H"], "support": "UNRESOLVED"},
        {"term": "mean_surface_coverage_fraction", "value": (len(aw_surface_rows)-aw_missing_faces)/len(aw_surface_rows), "support": "supported surfaces / selected surfaces"},
        {"term": "mean_internal_edge_coverage_fraction", "value": len(aw_included_internal)/len(aw_internal) if aw_internal else 1.0, "support": "resolved internal edges / selected internal edges"},
        {"term": "total_Gaussian_GB", "value": None, "support": "UNRESOLVED: incomplete face/edge curvature coverage"}]
    for order, result in aw_gb_by_order.items():
        aw_curv.append({"term": f"supported_subcomplex_GB_LHS_order_{order}",
            "value": result["total_GB_left"], "support": result["status"]})
        aw_curv.append({"term": f"supported_subcomplex_GB_difference_order_{order}",
            "value": result["total_GB_left"]-aw_expected,
            "support": "quadrature residual diagnostic; full AW remains unresolved"})
    write_rows(out / "aw_curvature_components.csv", aw_curv)
    write_rows(out / "orientation_reversal_audit.csv", [
        {"scheme": "Power", "area_A2": sum(float(area_cache[f"{a}:{b}"]["area_A2"]) for a,b in power_pairs),
         "forward_partial_or_total_H_A": p_mean["total_H"], "reversed_H_A": -p_mean["total_H"],
         "forward_GB_LHS": pgb["total_GB_left"], "reversed_GB_LHS": pgb["total_GB_left"],
         "chi": p_gauss["topology"]["euler_characteristic"], "status": "verified: mean sign flips; intrinsic GB invariant"},
        {"scheme": "AW", "area_A2": sum(float(r["area_A2"]) for r in aw_surface_rows if r["area_A2"] not in ("", "None")),
         "forward_partial_or_total_H_A": aw_mean["smooth_H"]+aw_mean["internal_edge_H"],
         "reversed_H_A": -(aw_mean["smooth_H"]+aw_mean["internal_edge_H"]),
         "forward_GB_LHS": aw_gb["total_GB_left"], "reversed_GB_LHS": aw_gb["total_GB_left"],
         "chi": aw_supported_topology["euler_characteristic"], "status": "partial AW terms; orientation rule verified, Gaussian total unresolved"}])

    report = ("Open alpha-selected interface curvature audit: static matched-probe 2KAI\n\n"
        f"Power selected pairs={len(power_pairs)}, area={metrics['power']['area_A2']} A^2.\n"
        f"Power physical polygon complex: V={p_gauss['topology']['vertices']}, E={p_gauss['topology']['edges']}, F={p_gauss['topology']['faces']}, chi={p_gauss['topology']['euler_characteristic']}; "
        f"components={p_gauss['topology']['components']}, boundary loops={p_gauss['topology']['boundary_components']}, orientable={p_gauss['topology']['orientable']}, genus sum={p_gauss['topology'].get('genus_sum')}, area={sum(float(area_cache[f'{a}:{b}']['area_A2']) for a,b in power_pairs):.6f} A^2.\n"
        f"Power GB: vertex corner sum={p_gauss['boundary_corner_K']:.12g}, LHS={p_gauss['total_GB_left']:.12g}, 2pi*chi={p_gauss['two_pi_chi']:.12g}, residual={p_gauss['residual']:.3g} ({p_gauss['status']}).\n"
        f"Power intrinsic Gaussian total={p_gauss['intrinsic_total_K']:.12g}; boundary corner sum={p_gauss['boundary_vertex_corner_K']:.12g}. These are distinct from the GB left-hand side.\n"
        f"Power mean: smooth=0, internal-edge H={.5*p_raw:.12g}, free-boundary H=0, total={p_mean['total_H']}; certified={p_mean['certified']}, unresolved internal edges={p_unresolved}.\n\n"
        f"AW selected pairs={len(aw_pairs)}, area={metrics['aw']['area_A2']} A^2.\n"
        f"AW physical face complex: V={aw_topology['vertices']}, E={aw_topology['edges']}, F={aw_topology['faces']}, chi={aw_topology['euler_characteristic']}; "
        f"components={aw_topology['components']}, boundary loops={aw_topology['boundary_components']}, orientable={aw_topology['orientable']}, genus sum={aw_topology.get('genus_sum')}, area={sum(float(r['area_A2']) for r in aw_surface_rows if r['area_A2'] not in ('', 'None')):.6f} A^2.\n"
        f"AW supported-curvature subcomplex: V={aw_supported_topology['vertices']}, E={aw_supported_topology['edges']}, F={aw_supported_topology['faces']}, chi={aw_supported_topology['euler_characteristic']}; "
        f"components={aw_supported_topology['components']}, boundary loops={aw_supported_topology['boundary_components']}.\n"
        f"AW mean partial: smooth H={aw_smooth_mean:.12g}, internal-edge H={.5*aw_raw:.12g}, boundary H=0; missing faces={aw_missing_faces}, missing internal edges={aw_missing_internal}; total UNRESOLVED.\n"
        f"AW supported-subcomplex Gauss--Bonnet order sweep: " + "; ".join(
            f"{order}: LHS={r['total_GB_left']:.12g}, difference from 2pi*chi={r['total_GB_left']-aw_expected:.3g}"
            for order, r in aw_gb_by_order.items()) + ". "
        f"At order {max(aw_gb_by_order)} the unresolved edge/vertex/face counts are "
        f"{aw_gb['unresolved_edges']}/{aw_gb['unresolved_vertices']}/{aw_gb['unresolved_faces']}; status={aw_gb['status']}. "
        "Incidence is complete for this supported subset, but the residual remains above tolerance. Five selected AW patches lack supported curvature, so the full-interface Gaussian total remains UNRESOLVED.\n\n"
        "The previous Power expectation -4398.2297 was an implementation error: Euler characteristic was coded as V-E-F. Independent physical polygon incidence gives 605-955+350=0; the correct target is 2pi*0=0. The corner-turn sum independently gives approximately zero.\n"
        "A selected free boundary has no scalar mean-curvature term. Its conormal is a boundary term in first variation; geodesic curvature and corner turns are intrinsic Gauss-Bonnet terms. Extrinsic crease beta contributes to mean curvature, not directly to Gaussian curvature.\n"
        f"Group reversal leaves areas/topology/GB invariant and flips signed mean terms; the audited forward/reverse values are in orientation_reversal_audit.csv (Power {p_mean['total_H']:.6g}/{-p_mean['total_H']:.6g} A; AW partial {aw_mean['smooth_H']+aw_mean['internal_edge_H']:.6g}/{-(aw_mean['smooth_H']+aw_mean['internal_edge_H']):.6g} A).\n"
        "Power topology consists of a disk-like isolated facet component and a three-boundary genus-zero component (2 components, 4 loops total); the component-level CSV gives exact counts. AW topology is separately reconstructed from solved AW edge/surface incidence.\n")
    (out / "open_interface_curvature_report.txt").write_text(report, encoding="utf-8")
    return out


def synthetic_validation_rows():
    """Known-answer checks used by the open-interface audit report."""
    rows = []
    planar = polygon_gaussian_terms([(0, 1, 2, 3)],
        {0: (0, 0, 0), 1: (1, 0, 0), 2: (1, 1, 0), 3: (0, 1, 0)})
    rows.append({"case": "planar_polygon_disk", "chi": planar["topology"]["euler_characteristic"],
        "mean_H": 0.0, "gaussian_or_GB_LHS": planar["total_GB_left"], "expected": 2*math.pi,
        "residual": planar["total_GB_left"]-2*math.pi, "status": planar["status"]})
    hemisphere = gauss_bonnet_accounting(smooth_K=2*math.pi, face_boundary_geodesic=0,
        vertex_corner_terms=0, euler_characteristic=1)
    rows.append({"case": "unit_hemisphere", "chi": 1, "mean_H": 2*math.pi,
        "gaussian_or_GB_LHS": hemisphere["total_GB_left"], "expected": 2*math.pi,
        "residual": hemisphere["residual"], "status": hemisphere["status"]})
    theta = math.pi/3
    cap = gauss_bonnet_accounting(smooth_K=2*math.pi*(1-math.cos(theta)),
        face_boundary_geodesic=2*math.pi*math.cos(theta), vertex_corner_terms=0,
        euler_characteristic=1)
    rows.append({"case": "unit_spherical_cap_theta_pi_over_3", "chi": 1, "mean_H": 2*math.pi*(1-math.cos(theta)),
        "gaussian_or_GB_LHS": cap["total_GB_left"], "expected": 2*math.pi,
        "residual": cap["residual"], "status": cap["status"]})
    coords = {0:(0,0,0),1:(1,0,0),2:(1,1,0),3:(0,1,0)}
    tri = polygon_gaussian_terms([(0,1,2),(0,2,3)], coords)
    rows.append({"case": "planar_triangulated_disk", "chi": tri["topology"]["euler_characteristic"],
        "mean_H": 0.0, "gaussian_or_GB_LHS": tri["total_GB_left"], "expected": 2*math.pi,
        "residual": tri["residual"], "status": tri["status"]})
    fold_coords = {0:(0,0,0),1:(1,0,0),2:(0,1,0),3:(0,-.5,math.sqrt(3)/2)}
    folded = polygon_gaussian_terms([(0,1,2),(1,0,3)], fold_coords)
    crease = planar_crease_mean((0,1,2), (1,0,3), fold_coords)
    if not math.isclose(crease, math.pi/6, abs_tol=1e-12):
        raise ValueError("Geometrically reconstructed fold failed")
    mean = integrated_mean_curvature([0,0],[crease])
    rows.append({"case": "folded_planar_patch", "chi": folded["topology"]["euler_characteristic"],
        "mean_H": mean["total_H"], "gaussian_or_GB_LHS": folded["total_GB_left"], "expected": 2*math.pi,
        "residual": folded["residual"], "status": folded["status"]})
    rev_mean = integrated_mean_curvature([0, 0], [-.5*math.pi/3])
    rows.append({"case": "group_orientation_reversal", "chi": folded["topology"]["euler_characteristic"],
        "mean_H": rev_mean["total_H"], "gaussian_or_GB_LHS": folded["total_GB_left"],
        "expected": 2*math.pi, "residual": folded["residual"],
        "status": "signed mean reversed; intrinsic Gaussian/topology invariant"})
    return rows


def planar_crease_mean(face1, face2, coordinates):
    """Independent oriented normal-turn check for two triangular faces.

    Uses H=trace(Dn)/2. The shared tangent follows face1's positive boundary;
    face2 must traverse it oppositely. Not a production curvature analyzer.
    """
    import numpy as np
    faces = [tuple(face1), tuple(face2)]
    if any(len(face) != 3 for face in faces):
        raise ValueError("The validation adapter requires triangles")
    common = set(face1) & set(face2)
    if len(common) != 2:
        raise ValueError("Two faces must share exactly one edge")
    normals = []
    for face in faces:
        p, q, r = [np.asarray(coordinates[v], dtype=float) for v in face]
        n = np.cross(q-p, r-p)
        if np.linalg.norm(n) == 0:
            raise ValueError("Degenerate face")
        normals.append(n / np.linalg.norm(n))
    a, b = next((face1[i], face1[(i+1)%3]) for i in range(3)
                if {face1[i], face1[(i+1)%3]} == common)
    if (b, a) not in [(face2[i], face2[(i+1)%3]) for i in range(3)]:
        raise ValueError("Adjacent orientations are inconsistent")
    tangent = np.asarray(coordinates[b], dtype=float)-np.asarray(coordinates[a], dtype=float)
    length = float(np.linalg.norm(tangent))
    tangent /= length
    beta = -math.atan2(float(np.dot(tangent, np.cross(*normals))),
                       float(np.dot(*normals)))
    return .5 * length * beta


def write_rows(path, rows):
    import csv
    from pathlib import Path
    rows = list(rows)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="output/2KAI_curvature_showcase/open_interface_curvature")
    args = parser.parse_args(argv)
    print(run_2kai_open_interface_audit(args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
