"""Read-only audit of Power alpha-boundary contact restrictions for dry 2KAI.

This is deliberately not a production molecular-surface exporter.  It audits
the existing weighted alpha data and reports the PL alpha-shape surrogate.
The exact union-of-balls contact surface is curved and is not constructed here.
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


def simplex_key(row: dict[str, str]) -> tuple[int, tuple[int, ...]]:
    simplex = tuple(sorted(int(part) for part in row["simplex"].replace(",", " ").split()))
    return int(row["dimension"]), simplex


def triangle_area(points: list[tuple[float, float, float]]) -> float:
    a, b, c = points
    u = tuple(b[i] - a[i] for i in range(3))
    v = tuple(c[i] - a[i] for i in range(3))
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return 0.5 * math.sqrt(sum(value * value for value in cross))


def make_contact_complex(
    simplices: dict[tuple[int, tuple[int, ...]], dict[str, str]],
    selected_ab: set[tuple[int, int]],
) -> dict[int, set[tuple[int, ...]]]:
    full = {
        dimension: {simplex for dim, simplex in simplices if dim == dimension}
        for dimension in range(4)
    }
    retained = {dimension: set() for dimension in range(4)}
    for dimension in range(1, 4):
        for simplex in full[dimension]:
            if any(set(edge) <= set(simplex) for edge in selected_ab):
                retained[dimension].add(simplex)
                for face_dimension in range(dimension):
                    retained[face_dimension].update(
                        tuple(sorted(face))
                        for face in combinations(simplex, face_dimension + 1)
                    )
    retained[1].update(selected_ab)
    retained[0].update((vertex,) for edge in selected_ab for vertex in edge)
    return retained


def induced(complex_: dict[int, set[tuple[int, ...]]], labels: dict[int, str], group: str) -> dict[int, set[tuple[int, ...]]]:
    return {
        dimension: {
            simplex for simplex in values
            if all(labels[vertex] == group for vertex in simplex)
        }
        for dimension, values in complex_.items()
    }


def alpha_boundary_triangles(
    simplices: dict[tuple[int, tuple[int, ...]], dict[str, str]],
    labels: dict[int, str],
    group: str,
) -> tuple[set[tuple[int, int, int]], dict[tuple[int, int, int], int]]:
    alpha_triangles = {
        simplex for (dimension, simplex), row in simplices.items()
        if dimension == 2 and row["in_alpha"] == "1"
        and all(labels[vertex] == group for vertex in simplex)
    }
    alpha_tetrahedra = {
        simplex for (dimension, simplex), row in simplices.items()
        if dimension == 3 and row["in_alpha"] == "1"
        and all(labels[vertex] == group for vertex in simplex)
    }
    cofaces = Counter()
    for tetrahedron in alpha_tetrahedra:
        for face in combinations(tetrahedron, 3):
            cofaces[tuple(sorted(face))] += 1
    # Maximal 2-simplices of the induced alpha complex: true boundary facets
    # have one tetrahedral coface; standalone triangles have none.
    boundary = {triangle for triangle in alpha_triangles if cofaces[triangle] <= 1}
    return boundary, dict(cofaces)


def graph_components(vertices: set[int], edges: set[tuple[int, int]]) -> int:
    graph = {vertex: set() for vertex in vertices}
    for a, b in edges:
        graph[a].add(b)
        graph[b].add(a)
    unseen = set(vertices)
    count = 0
    while unseen:
        count += 1
        root = unseen.pop()
        stack = [root]
        while stack:
            vertex = stack.pop()
            for neighbor in graph[vertex] & unseen:
                unseen.remove(neighbor)
                stack.append(neighbor)
    return count


def surface_metrics(
    triangles: set[tuple[int, int, int]],
    points: dict[int, tuple[float, float, float]],
) -> dict[str, object]:
    vertices = {vertex for triangle in triangles for vertex in triangle}
    edges = {tuple(sorted(edge)) for triangle in triangles for edge in combinations(triangle, 2)}
    edge_triangles: dict[tuple[int, int], int] = Counter(
        tuple(sorted(edge)) for triangle in triangles for edge in combinations(triangle, 2)
    )
    return {
        "triangles": len(triangles),
        "vertices": len(vertices),
        "edges": len(edges),
        "components": graph_components(vertices, edges) if vertices else 0,
        "boundary_edges": sum(value == 1 for value in edge_triangles.values()),
        "nonmanifold_edges": sum(value > 2 for value in edge_triangles.values()),
        "area_A2_center_triangle_surrogate": sum(
            triangle_area([points[vertex] for vertex in triangle])
            for triangle in triangles
        ),
        "coface_zero_triangles": 0,
    }


def main() -> None:
    atom_rows = read_csv(ROOT / "2kai_cazals_atoms.csv")
    labels = {int(row["index"]): row["group"] for row in atom_rows}
    points = {
        int(row["index"]): (float(row["x"]), float(row["y"]), float(row["z"]))
        for row in atom_rows
    }
    simplex_rows = read_csv(ROOT / "2kai_cazals_alpha0_audit" / "simplex_filtration.csv")
    simplices = {simplex_key(row): row for row in simplex_rows}
    selected_rows = read_csv(ROOT / "2kai_alpha0_edges.csv")
    selected_ab = {
        tuple(sorted((int(row["facet_1_i"]), int(row["facet_1_j"]))))
        for row in selected_rows
    } | {
        tuple(sorted((int(row["facet_2_i"]), int(row["facet_2_j"]))))
        for row in selected_rows
    }
    contact = make_contact_complex(simplices, selected_ab)
    contact_sides = {
        group: induced(contact, labels, group)
        for group in ("A", "B")
    }
    direct_atoms = {
        group: {vertex for edge in selected_ab for vertex in edge if labels[vertex] == group}
        for group in ("A", "B")
    }

    candidates = {}
    for group in ("A", "B"):
        boundary, cofaces = alpha_boundary_triangles(simplices, labels, group)
        contact_triangles = contact_sides[group][2]
        contact_vertices = {simplex[0] for simplex in contact_sides[group][0]}
        generator_restricted = {
            triangle for triangle in boundary if set(triangle) <= contact_vertices
        }
        direct_restricted = {
            triangle for triangle in boundary if set(triangle) <= direct_atoms[group]
        }
        incidence_restricted = boundary & contact_triangles
        candidates[group] = {
            "all_group_alpha_boundary": surface_metrics(boundary, points),
            "generator_restriction_N_contact_vertices": surface_metrics(generator_restricted, points),
            "direct_contact_atom_restriction": surface_metrics(direct_restricted, points),
            "incidence_restriction_N_contact_triangles": surface_metrics(incidence_restricted, points),
            "alpha_boundary_total_group_triangles": len(boundary),
            "alpha_boundary_standalone_triangles": sum(cofaces.get(triangle, 0) == 0 for triangle in boundary),
            "generator_restriction_not_N_contact_incidence": len(generator_restricted - contact_triangles),
            "direct_restriction_not_N_contact_incidence": len(direct_restricted - contact_triangles),
            "N_contact_side_triangles": len(contact_triangles),
            "N_contact_side_vertices": len(contact_vertices),
            "direct_contact_atoms": len(direct_atoms[group]),
            "N_contact_extra_atoms": len(contact_vertices - direct_atoms[group]),
        }

    summary = {
        "status": "read-only PL alpha-boundary surrogate; not exact curved union-of-balls surface",
        "selection": {
            "selected_ab_contacts": len(selected_ab),
            "direct_contact_atoms": {group: len(values) for group, values in direct_atoms.items()},
            "N_contact_generator_atoms": {
                group: len({simplex[0] for simplex in contact_sides[group][0]})
                for group in ("A", "B")
            },
        },
        "groups": candidates,
        "interpretation": {
            "back_facing_proxy": "A boundary triangle admitted by atom restriction but absent from N_contact side triangles is not supported by selected contact coface incidence; this is an exact incidence proxy, not a claim that every such triangle is geometrically back-facing.",
            "not_implemented": "spherical clipping of exact union-of-balls boundary by selected partner balls",
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "restricted_molecular_surface_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
