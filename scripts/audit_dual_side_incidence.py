"""Read-only Power dual-side incidence audit for dry 2KAI.

This writes audit artifacts only. It does not build or export a Dual A/B
surface.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "dual_side_incidence_audit" / "2kai_power"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def ints(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.replace(",", " ").split())


def simplex_key(row: dict[str, str]) -> tuple[int, tuple[int, ...]]:
    return int(row["dimension"]), tuple(sorted(ints(row["simplex"])))


def composition(vertices: tuple[int, ...], labels: dict[int, str]) -> str:
    return "".join(sorted(labels[index] for index in vertices))


def state(row: dict[str, str] | None) -> str:
    if row is None:
        return "unresolved_missing"
    return "alpha_selected" if row["in_alpha"] == "1" else "not_alpha_selected"


def atom_id(index: int, atoms: dict[int, dict[str, str]]) -> str:
    row = atoms[index]
    return f"{index}|pdb={row['pdb_serial']}|{row['chain']}:{row['atom_name']}"


def provenance(dimension: int) -> dict[str, str]:
    # The Cazals cache has generator tuples and simplex dimensions, but no
    # per-simplex VorPy primal-feature IDs or boundedness flags.
    return {
        "dual_dimension": str(dimension),
        "primal_feature_kind": {1: "surface", 2: "edge", 3: "vertex"}.get(
            dimension, "unknown"
        ),
        "primal_feature_ref": "not_recorded_in_cazals_cache",
        "bounded_complete_supported": "not_recorded_in_cazals_cache",
    }


def vector_sub(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(x - y for x, y in zip(a, b))


def dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sum(x * y for x, y in zip(a, b))


def cross(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def length(a: tuple[float, ...]) -> float:
    return math.sqrt(dot(a, a))


def tetra_geometry(vertices: tuple[int, ...], coordinates: dict[int, tuple[float, ...]]) -> dict:
    points = [coordinates[index] for index in vertices]
    edges = {
        f"{vertices[i]}-{vertices[j]}": length(vector_sub(points[i], points[j]))
        for i, j in combinations(range(4), 2)
    }
    a, b, c, d = points
    volume = abs(dot(vector_sub(b, a), cross(vector_sub(c, a), vector_sub(d, a)))) / 6.0
    faces = {}
    for face in combinations(range(4), 3):
        p, q, r = (points[index] for index in face)
        faces[" ".join(str(vertices[index]) for index in face)] = 0.5 * length(
            cross(vector_sub(q, p), vector_sub(r, p))
        )
    return {"vertices": list(vertices), "coordinates": points, "edge_lengths_A": edges,
            "face_areas_A2": faces, "volume_A3": volume}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    atoms_rows = read_csv(ROOT / "2kai_cazals_atoms.csv")
    atoms = {int(row["index"]): row for row in atoms_rows}
    labels = {index: row["group"] for index, row in atoms.items()}
    coordinates = {
        index: (float(row["x"]), float(row["y"]), float(row["z"]))
        for index, row in atoms.items()
    }

    simplex_rows = read_csv(ROOT / "2kai_cazals_alpha0_audit" / "simplex_filtration.csv")
    simplices = {simplex_key(row): row for row in simplex_rows}
    by_dimension: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    for dimension, simplex in simplices:
        by_dimension[dimension].append(simplex)

    selected_rows = read_csv(ROOT / "2kai_alpha0_edges.csv")
    selected_triangles = {
        tuple(sorted((int(row["triangle_i"]), int(row["triangle_j"]), int(row["triangle_k"]))))
        for row in selected_rows
    }
    selected_ab = set()
    for row in selected_rows:
        selected_ab.add(tuple(sorted((int(row["facet_1_i"]), int(row["facet_1_j"])))) )
        selected_ab.add(tuple(sorted((int(row["facet_2_i"]), int(row["facet_2_j"])))) )

    incident_triangles = {
        edge: [simplex for simplex in by_dimension[2] if set(edge) <= set(simplex)]
        for edge in selected_ab
    }
    incident_tetrahedra = {
        edge: [simplex for simplex in by_dimension[3] if set(edge) <= set(simplex)]
        for edge in selected_ab
    }
    interface_tetrahedra = sorted({tet for tets in incident_tetrahedra.values() for tet in tets})

    edge_rows = []
    for edge in sorted(selected_ab):
        triangles = incident_triangles[edge]
        tetrahedra = incident_tetrahedra[edge]
        tri_comps = Counter(composition(tri, labels) for tri in triangles)
        tri_states = Counter(state(simplices.get((2, tri))) for tri in triangles)
        tri_alpha_by_comp = Counter(
            (composition(tri, labels), state(simplices.get((2, tri)))) for tri in triangles
        )
        selected_tri_count = sum(tri in selected_triangles for tri in triangles)
        selected_by_comp = Counter(composition(tri, labels) for tri in triangles if tri in selected_triangles)
        tet_comps = Counter(composition(tet, labels) for tet in tetrahedra)
        edge_rows.append({
            "a": atom_id(edge[0], atoms),
            "b": atom_id(edge[1], atoms),
            "generator_indices": f"{edge[0]} {edge[1]}",
            "endpoint_groups": f"{labels[edge[0]]}{labels[edge[1]]}",
            "incident_triangles": len(triangles),
            "incident_AAB": tri_comps["AAB"],
            "incident_ABB": tri_comps["ABB"],
            "selected_contact_triangles": selected_tri_count,
            "selected_AAB_triangles": selected_by_comp["AAB"],
            "selected_ABB_triangles": selected_by_comp["ABB"],
            "alpha_selected_triangles": tri_states["alpha_selected"],
            "alpha_selected_AAB": tri_alpha_by_comp["AAB", "alpha_selected"],
            "alpha_selected_ABB": tri_alpha_by_comp["ABB", "alpha_selected"],
            "not_alpha_triangles": tri_states["not_alpha_selected"],
            "not_alpha_AAB": tri_alpha_by_comp["AAB", "not_alpha_selected"],
            "not_alpha_ABB": tri_alpha_by_comp["ABB", "not_alpha_selected"],
            "unresolved_triangles": tri_states["unresolved_missing"],
            "incident_tetrahedra": len(tetrahedra),
            "tetra_compositions": json.dumps(dict(sorted(tet_comps.items())), sort_keys=True),
            **provenance(1),
        })

    triangle_rows = []
    neighborhood_triangles = sorted({tri for tris in incident_triangles.values() for tri in tris})
    for tri in neighborhood_triangles:
        row = simplices.get((2, tri))
        triangle_rows.append({
            "generator_indices": " ".join(map(str, tri)),
            "generator_ids": " | ".join(atom_id(index, atoms) for index in tri),
            "composition": composition(tri, labels),
            "selected_contact": tri in selected_triangles,
            "alpha_status": state(row),
            "filtration": row["filtration"] if row else "",
            **provenance(2),
        })

    tetra_rows = []
    for tet in interface_tetrahedra:
        row = simplices.get((3, tet))
        tet_edges = [tuple(sorted(edge)) for edge in combinations(tet, 2)]
        tet_faces = [tuple(sorted(face)) for face in combinations(tet, 3)]
        selected_edges = [edge for edge in tet_edges if edge in selected_ab]
        mixed_faces = [face for face in tet_faces if composition(face, labels) in ("AAB", "ABB")]
        mono_faces = [face for face in tet_faces if composition(face, labels) in ("AAA", "BBB")]
        tetra_rows.append({
            "generator_indices": " ".join(map(str, tet)),
            "generator_ids": " | ".join(atom_id(index, atoms) for index in tet),
            "composition": composition(tet, labels),
            "all_edges": " | ".join(" ".join(map(str, edge)) + ":" + composition(edge, labels) for edge in tet_edges),
            "all_faces": " | ".join(" ".join(map(str, face)) + ":" + composition(face, labels) for face in tet_faces),
            "selected_AB_edges": " | ".join(" ".join(map(str, edge)) for edge in selected_edges),
            "mixed_faces": " | ".join(" ".join(map(str, face)) + ":" + composition(face, labels) for face in mixed_faces),
            "monochromatic_faces": " | ".join(" ".join(map(str, face)) + ":" + composition(face, labels) for face in mono_faces),
            "alpha_status": state(row),
            "filtration": row["filtration"] if row else "",
            **provenance(3),
        })

    aabb_rows = [row for row in tetra_rows if row["composition"] == "AABB"]
    representatives = {}
    for target in ("AAAB", "AABB", "ABBB"):
        candidates = [tet for tet in interface_tetrahedra if composition(tet, labels) == target]
        candidate = next(
            (tet for tet in candidates if state(simplices.get((3, tet))) == "alpha_selected"),
            candidates[0] if candidates else None,
        )
        if candidate is not None:
            representatives[target] = {
                "generator_ids": [atom_id(index, atoms) for index in candidate],
                "composition": target,
                "selected_AB_edges": [list(edge) for edge in combinations(candidate, 2)
                                       if tuple(sorted(edge)) in selected_ab],
                "geometry": tetra_geometry(candidate, coordinates),
                "alpha_status": state(simplices.get((3, candidate))),
            }

    synthetic = {}
    for name, groups in {
        "AAAB": "AAAB", "AABB": "AABB", "ABBB": "ABBB"
    }.items():
        vertices = tuple(f"{group}{index + 1}" for index, group in enumerate(groups))
        labels_s = dict(zip(vertices, groups))
        edges = [
            {"vertices": list(edge), "composition": "".join(sorted(labels_s[v] for v in edge))}
            for edge in combinations(vertices, 2)
        ]
        faces = [
            {"vertices": list(face), "composition": "".join(sorted(labels_s[v] for v in face))}
            for face in combinations(vertices, 3)
        ]
        reverse = {vertex: vertex.replace("A", "X").replace("B", "A").replace("X", "B") for vertex in vertices}
        synthetic[name] = {"vertices": list(vertices), "edges": edges, "faces": faces,
                           "reversal_vertices": [reverse[v] for v in vertices]}

    global_counts = Counter((int(row["dimension"]), row["composition"]) for row in simplex_rows)
    neighborhood_tri_counts = Counter(row["composition"] for row in triangle_rows)
    neighborhood_tet_counts = Counter(row["composition"] for row in tetra_rows)
    summary = {
        "selection_source": "2kai_alpha0_edges.csv",
        "selected_mixed_triangles": len(selected_triangles),
        "selected_unique_AB_contacts": len(selected_ab),
        "neighborhood_unique_triangles": len(triangle_rows),
        "neighborhood_unique_tetrahedra": len(tetra_rows),
        "global_simplex_counts": {f"d{d}:{c}": count for (d, c), count in sorted(global_counts.items())},
        "neighborhood_triangle_counts": dict(sorted(neighborhood_tri_counts.items())),
        "neighborhood_tetrahedron_counts": dict(sorted(neighborhood_tet_counts.items())),
        "aabb_tetrahedra": len(aabb_rows),
        "alpha_membership_coverage": {
            "d0_d3_rows": len(simplex_rows),
            "missing_neighborhood_triangles": sum(row["alpha_status"] == "unresolved_missing" for row in triangle_rows),
            "missing_interface_tetrahedra": sum(row["alpha_status"] == "unresolved_missing" for row in tetra_rows),
        },
        "provenance_limit": "Cazals CSV has generator tuples and alpha state, but not VorPy primal IDs or bounded/complete/support flags.",
    }

    write_csv(OUT / "selected_ab_edges.csv", edge_rows)
    write_csv(OUT / "interface_triangles.csv", triangle_rows)
    write_csv(OUT / "interface_tetrahedra.csv", tetra_rows)
    write_csv(OUT / "aabb_tetrahedra.csv", aabb_rows)
    for name, value in (("summary.json", summary), ("representatives.json", representatives), ("synthetic_cases.json", synthetic)):
        (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
