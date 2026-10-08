"""Dependency-free exact combinatorial topology checks."""

from __future__ import annotations

from collections import Counter, defaultdict
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vorpy.src.geometry.molecular_contact import (
    AtomBall,
    CircularArc,
    GeometryTolerance,
    JunctionVertex,
    SphericalPatch,
    TopologicalCut,
    _canonical_shared_circle,
    _cell_edge_incidence_distribution,
    _refine_shared_circle,
    _split_common_arc,
)


def common_refinement_check():
    """Regression for unequal subdivision, wraparound, and pair reversal."""
    tolerance = GeometryTolerance()
    first = AtomBall(0, "A", "A", (0.0, 0.0, 0.0), 2.0)
    second = AtomBall(1, "B", "B", (2.0, 0.0, 0.0), 2.0)
    geometry = _canonical_shared_circle(first, second, tolerance)
    reverse_geometry = _canonical_shared_circle(second, first, tolerance)
    assert geometry is not None and reverse_geometry is not None
    assert geometry[0] == reverse_geometry[0]
    assert geometry[1:5] == reverse_geometry[1:5]
    _, center, radial, _, basis_1, basis_2 = geometry

    def point(theta):
        return tuple(
            center[index] + radial * (
                basis_1[index] * math.cos(theta)
                + basis_2[index] * math.sin(theta)
            )
            for index in range(3)
        )

    def arc(arc_id, carrier, source, start, end):
        span = (end - start) % (2.0 * math.pi)
        return CircularArc(
            arc_id, carrier, "local", arc_id + ":s", arc_id + ":e",
            point(start), point(end), (1.0, 0.0, 0.0), 0.0, span,
            (source,), (), 1,
        )

    # First carrier has a subdivision at 1.4; the second has one at 2.2.
    # Both copies wrap through zero using the same geometric circle.
    arcs = [
        arc("a0", "A", "B", 0.2, 1.4),
        arc("a1", "A", "B", 1.4, 3.8),
        arc("a2", "A", "B", 3.8, 0.2),
        arc("b0", "B", "A", 0.2, 2.2),
        arc("b1", "B", "A", 2.2, 3.8),
        arc("b2", "B", "A", 3.8, 0.2),
    ]
    patch_a = SphericalPatch(
        "patch-a", "A", first.center, first.radius, "A", ("B",), (), ("contact",),
        "network", "interface", "test", ("a0", "a1", "a2"), (), "patch-a", 1.0, 0.0,
        "arrangement_component",
    )
    patch_b = SphericalPatch(
        "patch-b", "B", second.center, second.radius, "B", ("A",), (), ("contact",),
        "network", "interface", "test", ("b0", "b1", "b2"), (), "patch-b", 1.0, 0.0,
        "arrangement_component",
    )
    pieces, diagnostics = _refine_shared_circle(
        first,
        second,
        arcs,
        {"a0": "patch-a", "a1": "patch-a", "a2": "patch-a", "b0": "patch-b", "b1": "patch-b", "b2": "patch-b"},
        {"patch-a": patch_a, "patch-b": patch_b},
        [],
        tolerance,
    )
    assert not diagnostics
    assert len(pieces) == 5
    assert all(piece.classification == "internal" for piece in pieces)
    assert all(len(piece.incident_patch_ids) == 2 for piece in pieces)
    reversed_pieces, reversed_diagnostics = _refine_shared_circle(
        second,
        first,
        list(reversed(arcs)),
        {"a0": "patch-a", "a1": "patch-a", "a2": "patch-a", "b0": "patch-b", "b1": "patch-b", "b2": "patch-b"},
        {"patch-a": patch_a, "patch-b": patch_b},
        [],
        tolerance,
    )
    assert not reversed_diagnostics
    assert [piece.seam_id for piece in pieces] == [piece.seam_id for piece in reversed_pieces]


def topology(vertices, edges, faces):
    edge_faces = Counter()
    face_degrees = Counter()
    for face in faces:
        for vertex in face:
            face_degrees[vertex] += 1
        for edge in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_faces[tuple(sorted(edge))] += 1
    boundary = [edge for edge, count in edge_faces.items() if count == 1]
    graph = defaultdict(set)
    surface_graph = defaultdict(set)
    for left, right in edge_faces:
        surface_graph[left].add(right)
        surface_graph[right].add(left)
    for left, right in boundary:
        graph[left].add(right)
        graph[right].add(left)
    components = 0
    seen = set()
    for vertex in vertices:
        if vertex in seen:
            continue
        components += 1
        todo = [vertex]
        while todo:
            current = todo.pop()
            if current in seen:
                continue
            seen.add(current)
            todo.extend(surface_graph[current] - seen)
    boundary_components = 0
    seen_boundary = set()
    for vertex in list(graph):
        if not graph[vertex]:
            continue
        if vertex in seen_boundary:
            continue
        boundary_components += 1
        todo = [vertex]
        while todo:
            current = todo.pop()
            if current in seen_boundary:
                continue
            seen_boundary.add(current)
            todo.extend(graph[current] - seen_boundary)
    return {
        "V": len(vertices),
        "E": len(edge_faces),
        "F": len(faces),
        "chi": len(vertices) - len(edge_faces) + len(faces),
        "components": components,
        "boundary_edges": len(boundary),
        "boundary_components": boundary_components,
        "boundary_degrees": {vertex: len(graph[vertex]) for vertex in graph},
        "face_degrees": dict(face_degrees),
    }


def rank_mod2(rows, width):
    rows = [sum(1 << index for index in row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((index for index in range(rank, len(rows)) if rows[index] & (1 << column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index] & (1 << column):
                rows[index] ^= rows[rank]
        rank += 1
    return rank


def chain_stats(vertices, edges, face_words):
    """Exact integer/GF(2) audit for small oriented cell complexes."""
    edge_index = {edge_id: index for index, (edge_id, _start, _end) in enumerate(edges)}
    vertex_index = {vertex: index for index, vertex in enumerate(sorted(vertices))}
    d1_rows = [[] for _ in vertices]
    for column, (_edge_id, start, end) in enumerate(edges):
        d1_rows[vertex_index[start]].append(column)
        d1_rows[vertex_index[end]].append(column)
    d2_rows = []
    defects = []
    incidence = Counter()
    for face_index, word in enumerate(face_words):
        coefficients = Counter()
        boundary = Counter()
        for edge_id, start, end in word:
            if edge_id not in edge_index:
                defects.append(f"missing_edge:{face_index}:{edge_id}")
                continue
            edge_start, edge_end = edges[edge_index[edge_id]][1:]
            sign = 1 if (start, end) == (edge_start, edge_end) else -1
            coefficients[edge_id] += sign
            boundary[edge_start] -= sign
            boundary[edge_end] += sign
            incidence[edge_id] += 1
        if any(value for value in boundary.values()):
            defects.append(f"boundary_of_boundary:{face_index}")
        d2_rows.append([edge_index[edge_id] for edge_id, value in coefficients.items() if value % 2])
    for edge_id in edge_index:
        count = incidence[edge_id]
        incidence[edge_id] = count
    rank_1 = rank_mod2(d1_rows, len(edges))
    rank_2 = rank_mod2(d2_rows, len(edges))
    betti = (len(vertices) - rank_1, len(edges) - rank_1 - rank_2, len(face_words) - rank_2)
    return {
        "V": len(vertices),
        "E": len(edges),
        "F": len(face_words),
        "chi": len(vertices) - len(edges) + len(face_words),
        "betti": betti,
        "incidence": dict(incidence),
        "defects": tuple(defects),
    }


def topology_audit_regression():
    # One disk, one annulus, and a pair of pants.  The auxiliary cut is
    # attached twice with opposite orientation to the same face word.
    disk_edges = [("ab", "a", "b"), ("bc", "b", "c"), ("ca", "c", "a")]
    disk = chain_stats({"a", "b", "c"}, disk_edges, [[
        ("ab", "a", "b"), ("bc", "b", "c"), ("ca", "c", "a")]])
    assert disk["defects"] == () and disk["betti"] == (1, 0, 0) and disk["chi"] == 1

    annulus_edges = [
        ("o0", "o0", "o1"), ("o1", "o1", "o2"), ("o2", "o2", "o3"), ("o3", "o3", "o0"),
        ("i0", "i0", "i1"), ("i1", "i1", "i2"), ("i2", "i2", "i3"), ("i3", "i3", "i0"),
        ("cut", "o0", "i0"),
    ]
    annulus = chain_stats(
        {"o0", "o1", "o2", "o3", "i0", "i1", "i2", "i3"},
        annulus_edges,
        [[
            ("o0", "o0", "o1"), ("o1", "o1", "o2"), ("o2", "o2", "o3"), ("o3", "o3", "o0"),
            ("cut", "o0", "i0"), ("i3", "i0", "i3"), ("i2", "i3", "i2"), ("i1", "i2", "i1"), ("i0", "i1", "i0"),
            ("cut", "i0", "o0"),
        ]],
    )
    assert annulus["defects"] == () and annulus["betti"] == (1, 1, 0)
    assert annulus["chi"] == 0 and annulus["incidence"]["cut"] == 2

    # Three boundary cycles require two independent cuts.
    pants = chain_stats(
        {"a0", "a1", "a2", "b0", "b1", "b2", "c0", "c1", "c2"},
        [("a0", "a0", "a1"), ("a1", "a1", "a2"), ("a2", "a2", "a0"),
         ("b0", "b0", "b1"), ("b1", "b1", "b2"), ("b2", "b2", "b0"),
         ("c0", "c0", "c1"), ("c1", "c1", "c2"), ("c2", "c2", "c0"),
         ("cut1", "a0", "b0"), ("cut2", "a0", "c0")],
        [[("a0", "a0", "a1"), ("a1", "a1", "a2"), ("a2", "a2", "a0"),
          ("cut1", "a0", "b0"), ("b2", "b0", "b2"), ("b1", "b2", "b1"), ("b0", "b1", "b0"),
          ("cut1", "b0", "a0"), ("cut2", "a0", "c0"), ("c2", "c0", "c2"), ("c1", "c2", "c1"), ("c0", "c1", "c0"),
          ("cut2", "c0", "a0")]],
    )
    assert pants["defects"] == () and pants["betti"] == (1, 2, 0) and pants["chi"] == -1

    disconnected = chain_stats(
        {"a", "b", "c", "d", "e", "f"},
        [("ab", "a", "b"), ("bc", "b", "c"), ("ca", "c", "a"),
         ("de", "d", "e"), ("ef", "e", "f"), ("fd", "f", "d")],
        [[("ab", "a", "b"), ("bc", "b", "c"), ("ca", "c", "a")],
         [("de", "d", "e"), ("ef", "e", "f"), ("fd", "f", "d")]],
    )
    assert disconnected["betti"] == (2, 0, 0) and disconnected["chi"] == 2


def orientation_parity_regression():
    def consistent(edges):
        labels = {}
        for left, right, required in edges:
            if left in labels and right in labels and labels[left] ^ labels[right] != required:
                return False
            labels.setdefault(left, 0)
            labels.setdefault(right, labels[left] ^ required)
        return True

    assert consistent([("f0", "f1", 0), ("f1", "f2", 0), ("f2", "f0", 0)])
    assert not consistent([("f0", "f1", 0), ("f1", "f2", 0), ("f2", "f0", 1)])


def same_face_cut_regression():
    edges = [
        ("e0", "a", "b", ("face",)),
        ("e1", "b", "c", ("face",)),
        ("e2", "c", "a", ("face",)),
        ("cut", "a", "b", ("face", "face")),
    ]
    cut = TopologicalCut("cut", "face", "carrier", (0, 1), ("a", "b"), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    distribution, issues = _cell_edge_incidence_distribution(
        edges,
        {"face": ((("e0", "e1", "e2"), ("a", "b", "c")),)},
        (cut,),
    )
    assert issues == ()
    assert dict(distribution) == {"one_face": 3, "same_face_double": 1}


def check():
    common_refinement_check()
    topology_audit_regression()
    orientation_parity_regression()
    same_face_cut_regression()
    disk = topology({0, 1, 2}, {(0, 1), (1, 2), (0, 2)}, [(0, 1, 2)])
    assert (disk["V"], disk["E"], disk["F"], disk["chi"], disk["components"], disk["boundary_components"]) == (3, 3, 1, 1, 1, 1)

    unioned_caps = topology({0, 1, 2, 3}, set(), [(0, 1, 2), (0, 1, 3)])
    assert (unioned_caps["chi"], unioned_caps["components"], unioned_caps["boundary_components"]) == (1, 1, 1)

    annulus_faces = []
    for index in range(4):
        outer = index
        next_outer = (index + 1) % 4
        inner = 4 + index
        next_inner = 4 + (index + 1) % 4
        annulus_faces.extend(((outer, next_outer, inner), (next_outer, next_inner, inner)))
    annulus = topology(set(range(8)), set(), annulus_faces)
    assert (annulus["V"], annulus["E"], annulus["F"], annulus["chi"], annulus["boundary_components"]) == (8, 16, 8, 0, 2)

    seam_joined = topology({0, 1, 2, 3}, set(), [(0, 1, 2), (0, 1, 3)])
    assert seam_joined["boundary_components"] == 1

    disconnected = topology({0, 1, 2, 3, 4, 5}, set(), [(0, 1, 2), (3, 4, 5)])
    assert (disconnected["chi"], disconnected["components"], disconnected["boundary_components"]) == (2, 2, 2)

    triple_junction = topology({0, 1, 2, 3}, set(), [(0, 1, 2), (0, 2, 3), (0, 3, 1)])
    assert triple_junction["chi"] == 1
    assert triple_junction["face_degrees"][0] == 3
    assert triple_junction["boundary_components"] == 1

    degenerate = topology({0, 1, 2, 3, 4}, set(), [(0, 1, 2), (0, 3, 4)])
    assert degenerate["boundary_degrees"][0] == 4
    assert degenerate["boundary_components"] == 1


if __name__ == "__main__":
    check()
    print("molecular contact topology synthetic checks: PASS")
