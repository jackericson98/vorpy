"""Read-only induced-generator-side complex audit for dry 2KAI Power data."""

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


def tuple_key(text: str) -> tuple[int, ...]:
    return tuple(sorted(int(part) for part in text.replace(",", " ").split()))


def simplex_key(row: dict[str, str]) -> tuple[int, tuple[int, ...]]:
    return int(row["dimension"]), tuple_key(row["simplex"])


def coords(row: dict[str, str]) -> tuple[float, float, float]:
    return float(row["x"]), float(row["y"]), float(row["z"])


def stable_id(index: int, row: dict[str, str]) -> str:
    return f"{index}|pdb={row['pdb_serial']}|group={row['group']}|{row['chain']}:{row['residue_name']}{row['residue_number']}:{row['atom_name']}"


def point_distance(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def triangle_area(points: list[tuple[float, ...]]) -> float:
    a, b, c = points
    u = tuple(b[i] - a[i] for i in range(3))
    v = tuple(c[i] - a[i] for i in range(3))
    cross = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    return 0.5 * math.sqrt(sum(value * value for value in cross))


def rank_mod2(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = row
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def components(vertices: set[int], edges: set[tuple[int, int]]) -> list[set[int]]:
    graph = {vertex: set() for vertex in vertices}
    for a, b in edges:
        graph[a].add(b)
        graph[b].add(a)
    result = []
    unseen = set(vertices)
    while unseen:
        root = unseen.pop()
        component = {root}
        stack = [root]
        while stack:
            vertex = stack.pop()
            for neighbor in graph[vertex] & unseen:
                unseen.remove(neighbor)
                component.add(neighbor)
                stack.append(neighbor)
        result.append(component)
    return result


def triangle_components(triangles: set[tuple[int, int, int]]) -> list[set[tuple[int, int, int]]]:
    edge_to_triangles: dict[tuple[int, int], set[tuple[int, int, int]]] = defaultdict(set)
    for triangle in triangles:
        for edge in combinations(triangle, 2):
            edge_to_triangles[tuple(sorted(edge))].add(triangle)
    graph = {triangle: set() for triangle in triangles}
    for incident in edge_to_triangles.values():
        for left in incident:
            graph[left].update(incident - {left})
    result = []
    unseen = set(triangles)
    while unseen:
        root = unseen.pop()
        component = {root}
        stack = [root]
        while stack:
            triangle = stack.pop()
            for neighbor in graph[triangle] & unseen:
                unseen.remove(neighbor)
                component.add(neighbor)
                stack.append(neighbor)
        result.append(component)
    return result


def make_complex(
    simplices: dict[tuple[int, tuple[int, ...]], dict[str, str]],
    selected_ab: set[tuple[int, int]],
    alpha_only: bool,
) -> dict[int, set[tuple[int, ...]]]:
    full = {dimension: {simplex for dim, simplex in simplices if dim == dimension}
            for dimension in range(4)}
    cofaces = {
        dimension: {
            simplex for simplex in full[dimension]
            if any(set(edge) <= set(simplex) for edge in selected_ab)
        }
        for dimension in range(1, 4)
    }
    retained: dict[int, set[tuple[int, ...]]] = {3: set(), 2: set(), 1: set(), 0: set()}
    for dimension, values in cofaces.items():
        retained[dimension].update(values)
        for simplex in values:
            for face_dimension in range(dimension):
                retained[face_dimension].update(
                    tuple(sorted(face)) for face in combinations(simplex, face_dimension + 1)
                )
    retained[1].update(selected_ab)
    for edge in selected_ab:
        retained[0].update((vertex,) for vertex in edge)
    if alpha_only:
        retained = {
            dimension: {simplex for simplex in values
                        if simplices[(dimension, simplex)]["in_alpha"] == "1"}
            for dimension, values in retained.items()
        }
    return retained


def induced(complex_: dict[int, set[tuple[int, ...]]], labels: dict[int, str], group: str) -> dict[int, set[tuple[int, ...]]]:
    return {
        dimension: {simplex for simplex in values if all(labels[vertex] == group for vertex in simplex)}
        for dimension, values in complex_.items()
    }


def topology(
    complex_: dict[int, set[tuple[int, ...]]],
    atom_coordinates: dict[int, tuple[float, ...]],
) -> dict:
    vertices = {simplex[0] for simplex in complex_[0]}
    edges = complex_[1]
    triangles = complex_[2]
    tetrahedra = complex_[3]
    edge_triangle_counts = Counter(
        edge
        for triangle in triangles
        for edge in combinations(triangle, 2)
    )
    edge_triangle_counts = Counter(tuple(sorted(edge)) for edge in edge_triangle_counts.elements())
    boundary_edges = {edge for edge, count in edge_triangle_counts.items() if count == 1}
    nonmanifold_edges = {edge for edge, count in edge_triangle_counts.items() if count > 2}
    triangle_components_ = triangle_components(triangles)
    graph_components = components(vertices, edges)
    edge_index = {edge: index for index, edge in enumerate(sorted(edges))}
    boundary_rows = []
    for triangle in sorted(triangles):
        bits = 0
        for edge in combinations(triangle, 2):
            bits ^= 1 << edge_index[tuple(sorted(edge))]
        boundary_rows.append(bits)
    rank_boundary_2 = rank_mod2(boundary_rows)
    rank_boundary_1 = len(vertices) - len(graph_components)
    betti_0 = len(graph_components)
    betti_1 = len(edges) - rank_boundary_1 - rank_boundary_2
    betti_2 = len(triangles) - rank_boundary_2
    dangling_edges = {edge for edge in edges if edge_triangle_counts[edge] == 0}
    degree = Counter(vertex for edge in edges for vertex in edge)
    graph_leaf_edges = {edge for edge in edges if degree[edge[0]] == 1 or degree[edge[1]] == 1}
    return {
        "vertices": len(vertices),
        "edges": len(edges),
        "triangles": len(triangles),
        "tetrahedra": len(tetrahedra),
        "maximum_dimension": max((dimension for dimension, values in complex_.items() if values), default=-1),
        "euler_characteristic": len(vertices) - len(edges) + len(triangles) - len(tetrahedra),
        "connected_components": len(graph_components),
        "betti_numbers_mod2": [betti_0, betti_1, betti_2, 0],
        "isolated_vertices": len(vertices - {vertex for edge in edges for vertex in edge}),
        "dangling_edges_no_triangle": len(dangling_edges),
        "graph_leaf_edges": len(graph_leaf_edges),
        "triangle_components": len(triangle_components_),
        "triangle_component_sizes": sorted((len(component) for component in triangle_components_), reverse=True),
        "boundary_edges_of_2d_portion": len(boundary_edges),
        "nonmanifold_edges_of_2d_portion": len(nonmanifold_edges),
        "total_generator_center_edge_length_A": sum(
            point_distance(atom_coordinates[a], atom_coordinates[b]) for a, b in edges
        ),
        "total_generator_center_triangle_area_A2": sum(
            triangle_area([atom_coordinates[vertex] for vertex in triangle]) for triangle in triangles
        ),
    }


def write_complex_csv(path: Path, complexes: dict[str, dict[int, set[tuple[int, ...]]]], labels: dict[int, str]) -> None:
    rows = []
    for name, complex_ in complexes.items():
        for dimension, values in complex_.items():
            for simplex in sorted(values):
                rows.append({
                    "neighborhood": name,
                    "dimension": dimension,
                    "generator_indices": " ".join(map(str, simplex)),
                    "composition": "".join(sorted(labels[vertex] for vertex in simplex)),
                })
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    atoms_rows = read_csv(ROOT / "2kai_cazals_atoms.csv")
    atoms = {int(row["index"]): row for row in atoms_rows}
    labels = {index: row["group"] for index, row in atoms.items()}
    atom_coordinates = {index: coords(row) for index, row in atoms.items()}
    simplex_rows = read_csv(ROOT / "2kai_cazals_alpha0_audit" / "simplex_filtration.csv")
    simplices = {simplex_key(row): row for row in simplex_rows}
    selected_rows = read_csv(ROOT / "2kai_alpha0_edges.csv")
    selected_ab = set()
    selected_triangles = set()
    for row in selected_rows:
        selected_ab.add(tuple(sorted((int(row["facet_1_i"]), int(row["facet_1_j"])))) )
        selected_ab.add(tuple(sorted((int(row["facet_2_i"]), int(row["facet_2_j"])))) )
        selected_triangles.add(tuple(sorted((int(row["triangle_i"]), int(row["triangle_j"]), int(row["triangle_k"])))) )

    n_contact = make_complex(simplices, selected_ab, alpha_only=False)
    n_alpha = make_complex(simplices, selected_ab, alpha_only=True)
    complexes = {
        "N_contact": n_contact,
        "N_alpha": n_alpha,
        "N_contact_A": induced(n_contact, labels, "A"),
        "N_contact_B": induced(n_contact, labels, "B"),
        "N_alpha_A": induced(n_alpha, labels, "A"),
        "N_alpha_B": induced(n_alpha, labels, "B"),
    }
    side_complexes = {name: value for name, value in complexes.items() if name.endswith(("_A", "_B"))}
    summary = {
        "definitions": {
            "N_contact": "smallest closed dual subcomplex containing every simplex coface of a selected AB edge and all faces of those cofaces",
            "N_alpha": "N_contact intersected with the Power alpha-selected simplicial complex",
            "selected_ab_contacts": len(selected_ab),
            "selected_mixed_triangles": len(selected_triangles),
            "selected_mixed_triangles_alpha": sum(
                simplices[(2, triangle)]["in_alpha"] == "1" for triangle in selected_triangles
            ),
        },
        "neighborhood_counts": {
            name: {str(dimension): len(values) for dimension, values in complex_.items()}
            for name, complex_ in complexes.items() if name in ("N_contact", "N_alpha")
        },
        "induced_topology": {
            name: topology(complex_, atom_coordinates) for name, complex_ in side_complexes.items()
        },
    }

    contact_atom_set = {vertex for edge in selected_ab for vertex in edge}
    atom_relationship = {}
    for name, complex_ in side_complexes.items():
        side_vertices = {simplex[0] for simplex in complex_[0]}
        extra = sorted(side_vertices - contact_atom_set)
        atom_relationship[name] = {
            "vertices": len(side_vertices),
            "direct_selected_contact_atoms": len(side_vertices & contact_atom_set),
            "coface_only_atoms": len(extra),
            "coface_only_atom_indices": extra,
            "coface_only_atom_ids": [stable_id(index, atoms[index]) for index in extra],
        }
    summary["atom_relationship"] = atom_relationship
    summary["closure_checks"] = {
        name: all(
            tuple(sorted(face)) in complex_[dimension - 1]
            for dimension, values in complex_.items() if dimension > 0
            for simplex in values for face in combinations(simplex, dimension)
        )
        for name, complex_ in side_complexes.items()
    }
    OUT.mkdir(parents=True, exist_ok=True)
    write_complex_csv(OUT / "induced_generator_side_simplices.csv", side_complexes, labels)
    (OUT / "induced_generator_side_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
