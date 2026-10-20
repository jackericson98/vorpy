"""Integrated Gaussian curvature for closed planar polygon boundaries.

The AW network has intrinsic curvature on its analytic faces and edges.  POW
and PRM faces are planar instead: their Gaussian curvature is concentrated at
polyhedral vertices.  This module deliberately works from the stored face /
edge / vertex topology rather than from the surface mesh triangulation.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field
import math

import numpy as np


POLYHEDRAL_STATUS_CALCULATED = "CALCULATED"
POLYHEDRAL_STATUS_UNAVAILABLE = "UNAVAILABLE"


@dataclass
class PolyhedralGaussianCurvatureResult:
    """Measured polyhedral curvature and independent topology evidence."""

    status: str
    total: float | None
    vertex_curvatures: dict[int, float | None] = field(default_factory=dict)
    face_loops: dict[int, tuple[int, ...]] = field(default_factory=dict)
    face_angles: dict[int, tuple[float, ...]] = field(default_factory=dict)
    euler_characteristic: int | None = None
    component_count: int = 0
    closed: bool = False
    manifold: bool = False
    orientable: bool = False
    gauss_bonnet_expected: float | None = None
    gauss_bonnet_residual: float | None = None
    reasons: tuple[str, ...] = ()

    @property
    def available(self) -> bool:
        return self.status == POLYHEDRAL_STATUS_CALCULATED


def _result_unavailable(reason, *, vertex_curvatures=None, **kwargs):
    reasons = tuple(reason) if isinstance(reason, (tuple, list)) else (str(reason),)
    return PolyhedralGaussianCurvatureResult(
        status=POLYHEDRAL_STATUS_UNAVAILABLE,
        total=None,
        vertex_curvatures=dict(vertex_curvatures or {}),
        reasons=reasons,
        **kwargs,
    )


def _row(table, index, label):
    try:
        return table.iloc[int(index)]
    except (AttributeError, IndexError, KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label} {index} is not present in the network table") from error


def _int_values(values, label):
    if values is None:
        raise ValueError(f"{label} is missing")
    try:
        result = tuple(int(value) for value in values)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} contains a non-integer topology ID") from error
    if not result:
        raise ValueError(f"{label} is empty")
    return result


def _edge_endpoints(net, edge_id):
    edge = _row(net.edges, edge_id, "Edge")
    endpoints = _int_values(edge.get("verts"), f"Edge {edge_id} vertices")
    if len(endpoints) != 2 or endpoints[0] == endpoints[1]:
        raise ValueError(
            f"Edge {edge_id} is not a simple two-endpoint polygon edge"
        )
    return endpoints


def _trace_face_loop(net, face_id, face):
    edge_ids = _int_values(face.get("edges"), f"Face {face_id} edges")
    if len(set(edge_ids)) != len(edge_ids):
        raise ValueError(f"Face {face_id} repeats a boundary edge")

    listed_vertices = _int_values(face.get("verts"), f"Face {face_id} vertices")
    if len(set(listed_vertices)) != len(listed_vertices):
        raise ValueError(f"Face {face_id} repeats a boundary vertex")

    edge_endpoints = {
        edge_id: _edge_endpoints(net, edge_id)
        for edge_id in edge_ids
    }
    expected_vertices = {
        vertex for endpoints in edge_endpoints.values() for vertex in endpoints
    }
    if set(listed_vertices) != expected_vertices:
        raise ValueError(
            f"Face {face_id} vertex/edge incidence is inconsistent"
        )

    incident_edges = defaultdict(list)
    for edge_id, endpoints in edge_endpoints.items():
        for vertex in endpoints:
            incident_edges[vertex].append(edge_id)
    if any(len(edges) != 2 for edges in incident_edges.values()):
        raise ValueError(f"Face {face_id} does not form a single polygon loop")

    start = min(expected_vertices)
    current = start
    edge_id = min(incident_edges[start])
    used_edges = set()
    ordered_edges = []
    loop = [start]

    while True:
        if edge_id in used_edges:
            raise ValueError(f"Face {face_id} repeats an edge before closing")
        used_edges.add(edge_id)
        ordered_edges.append(edge_id)
        first, second = edge_endpoints[edge_id]
        if current == first:
            next_vertex = second
        elif current == second:
            next_vertex = first
        else:
            raise ValueError(f"Face {face_id} edge walk lost incidence")

        if next_vertex == start:
            break
        loop.append(next_vertex)
        candidates = [
            candidate for candidate in incident_edges[next_vertex]
            if candidate not in used_edges
        ]
        if len(candidates) != 1:
            raise ValueError(f"Face {face_id} is disconnected or ambiguous")
        current = next_vertex
        edge_id = candidates[0]

    if used_edges != set(edge_ids) or set(loop) != expected_vertices:
        raise ValueError(f"Face {face_id} is not one simple polygon loop")
    if len(loop) < 3:
        raise ValueError(f"Face {face_id} has fewer than three polygon vertices")
    return tuple(loop), tuple(ordered_edges)


def _plane_and_scale(points, face_id, tolerance):
    origin = points[0]
    scale = max(float(np.max(np.linalg.norm(points - origin, axis=1))), 1.0)
    normal = None
    for index in range(1, len(points) - 1):
        candidate = np.cross(points[index] - origin, points[index + 1] - origin)
        norm = float(np.linalg.norm(candidate))
        if np.isfinite(norm) and norm > tolerance * scale * scale:
            normal = candidate / norm
            break
    if normal is None:
        raise ValueError(f"Face {face_id} is degenerate or collinear")

    polygon_normal = np.sum(np.cross(points, np.roll(points, -1, axis=0)), axis=0)
    polygon_norm = float(np.linalg.norm(polygon_normal))
    if not np.isfinite(polygon_norm) or polygon_norm <= tolerance * scale * scale:
        raise ValueError(f"Face {face_id} has zero projected area")
    if float(np.dot(polygon_normal, normal)) < 0.0:
        normal = -normal

    distances = np.abs((points - origin) @ normal)
    if not np.all(np.isfinite(distances)):
        raise ValueError(f"Face {face_id} has non-finite coordinates")
    if float(np.max(distances)) > tolerance * scale:
        raise ValueError(f"Face {face_id} is not planar")
    return normal, scale


def _project(points, normal):
    axis = int(np.argmax(np.abs(normal)))
    return np.delete(points, axis, axis=1)


def _orientation_2d(a, b, c):
    first = b - a
    second = c - a
    return float(first[0] * second[1] - first[1] * second[0])


def _segments_intersect(a, b, c, d, tolerance):
    scale = max(
        float(np.linalg.norm(b - a)),
        float(np.linalg.norm(d - c)),
        float(np.linalg.norm(c - a)),
        float(np.linalg.norm(d - a)),
        1.0,
    )
    epsilon = tolerance * scale * scale

    def on_segment(p, q, r):
        return (
            min(p[0], r[0]) - epsilon <= q[0] <= max(p[0], r[0]) + epsilon
            and min(p[1], r[1]) - epsilon <= q[1] <= max(p[1], r[1]) + epsilon
        )

    ab_c = _orientation_2d(a, b, c)
    ab_d = _orientation_2d(a, b, d)
    cd_a = _orientation_2d(c, d, a)
    cd_b = _orientation_2d(c, d, b)
    if (ab_c > epsilon and ab_d < -epsilon or ab_c < -epsilon and ab_d > epsilon) \
            and (cd_a > epsilon and cd_b < -epsilon or cd_a < -epsilon and cd_b > epsilon):
        return True
    return (
        abs(ab_c) <= epsilon and on_segment(a, c, b)
        or abs(ab_d) <= epsilon and on_segment(a, d, b)
        or abs(cd_a) <= epsilon and on_segment(c, a, d)
        or abs(cd_b) <= epsilon and on_segment(c, b, d)
    )


def _validate_simple_polygon(points, normal, face_id, tolerance):
    projected = _project(points, normal)
    for index in range(len(projected)):
        for other in range(index):
            if np.linalg.norm(projected[index] - projected[other]) <= tolerance:
                raise ValueError(f"Face {face_id} repeats or collapses a vertex")

    count = len(projected)
    for first in range(count):
        first_next = (first + 1) % count
        for second in range(first + 1, count):
            second_next = (second + 1) % count
            if first == second or first_next == second or second_next == first:
                continue
            if _segments_intersect(
                projected[first], projected[first_next],
                projected[second], projected[second_next], tolerance,
            ):
                raise ValueError(f"Face {face_id} polygon edges intersect")


def _face_angles(points, normal, face_id, tolerance):
    _validate_simple_polygon(points, normal, face_id, tolerance)
    angles = []
    normal = np.asarray(normal, dtype=float)
    for index, current in enumerate(points):
        previous = points[index - 1]
        following = points[(index + 1) % len(points)]
        incoming = current - previous
        outgoing = following - current
        left = previous - current
        right = following - current
        left_norm = float(np.linalg.norm(left))
        right_norm = float(np.linalg.norm(right))
        if (
            not np.isfinite(left_norm)
            or not np.isfinite(right_norm)
            or left_norm <= tolerance
            or right_norm <= tolerance
        ):
            raise ValueError(f"Face {face_id} has a zero-length corner")
        cosine = float(np.dot(left, right) / (left_norm * right_norm))
        if not np.isfinite(cosine):
            raise ValueError(f"Face {face_id} has a non-finite corner angle")
        cosine = min(1.0, max(-1.0, cosine))
        angle = math.acos(cosine)

        turn = float(np.dot(np.cross(incoming, outgoing), normal))
        turn_tolerance = tolerance * left_norm * right_norm
        if turn < -turn_tolerance:
            angle = 2.0 * math.pi - angle
        angles.append(float(angle))
    return tuple(angles)


def _orientability(face_loops, edge_faces, edge_directions):
    orientation = {}
    for start in face_loops:
        if start in orientation:
            continue
        orientation[start] = 1
        queue = deque([start])
        while queue:
            face_id = queue.popleft()
            for edge_id in face_loops[face_id][1]:
                neighbors = edge_faces[edge_id]
                if len(neighbors) != 2:
                    continue
                other = neighbors[0] if neighbors[1] == face_id else neighbors[1]
                first = edge_directions[(face_id, edge_id)]
                second = edge_directions[(other, edge_id)]
                same_direction = first == second
                required = -orientation[face_id] if same_direction else orientation[face_id]
                if other in orientation:
                    if orientation[other] != required:
                        return False
                else:
                    orientation[other] = required
                    queue.append(other)
    return True


def calculate_polyhedral_gaussian_curvature(
        net,
        surface_ids,
        *,
        tolerance=1e-8,
        gauss_bonnet_tolerance=1e-6):
    """Calculate angle-defect curvature for a closed polygonal boundary.

    ``surface_ids`` are network surface-table positions.  Polygon order is
    reconstructed from each surface's edge endpoint topology, so the stored
    ``verts`` column need not be ordered.

    The returned total is always the sum of measured vertex defects.  The
    Gauss--Bonnet value is retained only as an independent validation field;
    it is never substituted for the measured result.
    """
    settings = getattr(net, "settings", None) or {}
    net_type = str(settings.get("net_type", "")).lower()
    if net_type not in {"pow", "prm"}:
        raise ValueError(
            "Polyhedral Gaussian curvature requires a POW or PRM network."
        )
    if getattr(net, "verts", None) is None or getattr(net, "edges", None) is None \
            or getattr(net, "surfs", None) is None:
        return _result_unavailable("network vertex, edge, or surface table is missing")

    try:
        boundary_faces = tuple(dict.fromkeys(int(value) for value in surface_ids))
    except (TypeError, ValueError):
        return _result_unavailable("boundary surface IDs are invalid")
    if not boundary_faces:
        return _result_unavailable("boundary has no faces")

    face_records = {}
    boundary_vertex_ids = set()
    vertex_curvatures = {}
    reasons = []
    for face_id in boundary_faces:
        try:
            face = _row(net.surfs, face_id, "Surface")
            try:
                listed_vertices = face.get("verts")
                if listed_vertices is not None:
                    boundary_vertex_ids.update(int(value) for value in listed_vertices)
            except (TypeError, ValueError):
                pass
            try:
                listed_edges = face.get("edges")
                if listed_edges is not None:
                    for edge_id in listed_edges:
                        boundary_vertex_ids.update(_edge_endpoints(net, int(edge_id)))
            except (IndexError, KeyError, TypeError, ValueError):
                pass
            loop, ordered_edge_ids = _trace_face_loop(net, face_id, face)
            points = np.asarray([
                np.asarray(_row(net.verts, vertex_id, "Vertex").get("loc"), dtype=float)
                for vertex_id in loop
            ], dtype=float)
            if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
                raise ValueError(f"Face {face_id} has invalid vertex coordinates")
            if (
                "flat" in face
                and face.get("flat") is not None
                and not bool(face.get("flat"))
            ):
                raise ValueError(f"Face {face_id} is not marked planar")
            normal, _ = _plane_and_scale(points, face_id, tolerance)
            angles = _face_angles(points, normal, face_id, tolerance)
            face_records[face_id] = {
                "loop": loop,
                "angles": angles,
                "edge_ids": ordered_edge_ids,
            }
        except (IndexError, KeyError, TypeError, ValueError) as error:
            reasons.append(str(error))

    if reasons:
        return _result_unavailable(
            reasons,
            face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
            face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
            vertex_curvatures={vertex_id: None for vertex_id in boundary_vertex_ids},
        )

    edge_faces = defaultdict(list)
    edge_directions = {}
    vertex_faces = defaultdict(set)
    boundary_edges = set()
    for face_id, record in face_records.items():
        loop = record["loop"]
        for index, edge_id in enumerate(record["edge_ids"]):
            edge_faces[edge_id].append(face_id)
            boundary_edges.add(edge_id)
            start = loop[index]
            end = loop[(index + 1) % len(loop)]
            edge_directions[(face_id, edge_id)] = (start, end)
        for vertex_id in loop:
            vertex_faces[vertex_id].add(face_id)

    nonmanifold_edges = [edge_id for edge_id, faces in edge_faces.items() if len(faces) != 2]
    if nonmanifold_edges:
        return _result_unavailable(
            f"boundary edges do not have exactly two faces: {sorted(nonmanifold_edges)}",
            face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
            face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
            vertex_curvatures={vertex_id: None for vertex_id in boundary_vertex_ids},
            closed=False,
            manifold=False,
        )

    vertex_edges = defaultdict(set)
    for edge_id in boundary_edges:
        try:
            first, second = _edge_endpoints(net, edge_id)
        except ValueError as error:
            return _result_unavailable(
                str(error),
                face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
                face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
                vertex_curvatures={vertex_id: None for vertex_id in boundary_vertex_ids},
            )
        vertex_edges[first].add(edge_id)
        vertex_edges[second].add(edge_id)

    bad_vertices = []
    for vertex_id, faces in vertex_faces.items():
        incident_edges = vertex_edges.get(vertex_id, set())
        link = defaultdict(set)
        for edge_id in incident_edges:
            edge_incident_faces = edge_faces[edge_id]
            for face_id in edge_incident_faces:
                if face_id in faces:
                    link[face_id].update(
                        other for other in edge_incident_faces
                        if other != face_id and other in faces
                    )
        if len(incident_edges) < 3 or any(len(neighbors) != 2 for neighbors in link.values()):
            bad_vertices.append(vertex_id)
            continue
        unseen = set(faces)
        queue = [next(iter(unseen))]
        unseen.remove(queue[0])
        while queue:
            face_id = queue.pop()
            for neighbor in link[face_id] & unseen:
                unseen.remove(neighbor)
                queue.append(neighbor)
        if unseen:
            bad_vertices.append(vertex_id)
    if bad_vertices:
        return _result_unavailable(
            f"boundary vertex links are not single cycles: {sorted(bad_vertices)}",
            face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
            face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
            vertex_curvatures={vertex_id: None for vertex_id in boundary_vertex_ids},
            closed=True,
            manifold=False,
        )

    face_neighbors = defaultdict(set)
    for faces in edge_faces.values():
        first, second = faces
        face_neighbors[first].add(second)
        face_neighbors[second].add(first)
    components = 0
    unseen_faces = set(face_records)
    while unseen_faces:
        components += 1
        queue = [unseen_faces.pop()]
        while queue:
            face_id = queue.pop()
            for neighbor in face_neighbors[face_id] & unseen_faces:
                unseen_faces.remove(neighbor)
                queue.append(neighbor)

    orientable = _orientability(
        {face_id: (record["loop"], record["edge_ids"]) for face_id, record in face_records.items()},
        edge_faces,
        edge_directions,
    )
    if not orientable:
        return _result_unavailable(
            "boundary orientation is not globally consistent",
            face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
            face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
            vertex_curvatures={vertex_id: None for vertex_id in boundary_vertex_ids},
            euler_characteristic=len(vertex_faces) - len(boundary_edges) + len(face_records),
            component_count=components,
            closed=True,
            manifold=True,
            orientable=False,
        )

    for face_id, record in face_records.items():
        for vertex_id, angle in zip(record["loop"], record["angles"]):
            vertex_curvatures[vertex_id] = vertex_curvatures.get(vertex_id, 2.0 * math.pi) - angle

    total = float(sum(value for value in vertex_curvatures.values() if value is not None))
    euler_characteristic = len(vertex_faces) - len(boundary_edges) + len(face_records)
    expected = 2.0 * math.pi * euler_characteristic
    residual = total - expected
    if not np.isfinite(total) or not np.isfinite(residual):
        return _result_unavailable(
            "polyhedral curvature is non-finite",
            vertex_curvatures=vertex_curvatures,
            face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
            face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
            euler_characteristic=euler_characteristic,
            component_count=components,
            closed=True,
            manifold=True,
            orientable=True,
            gauss_bonnet_expected=expected,
            gauss_bonnet_residual=residual,
        )
    if abs(residual) >= gauss_bonnet_tolerance:
        return _result_unavailable(
            f"polyhedral Gauss-Bonnet residual {residual:.12g} exceeds tolerance",
            vertex_curvatures=vertex_curvatures,
            face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
            face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
            euler_characteristic=euler_characteristic,
            component_count=components,
            closed=True,
            manifold=True,
            orientable=True,
            gauss_bonnet_expected=expected,
            gauss_bonnet_residual=residual,
        )

    return PolyhedralGaussianCurvatureResult(
        status=POLYHEDRAL_STATUS_CALCULATED,
        total=total,
        vertex_curvatures=vertex_curvatures,
        face_loops={face_id: record["loop"] for face_id, record in face_records.items()},
        face_angles={face_id: record["angles"] for face_id, record in face_records.items()},
        euler_characteristic=euler_characteristic,
        component_count=components,
        closed=True,
        manifold=True,
        orientable=True,
        gauss_bonnet_expected=expected,
        gauss_bonnet_residual=residual,
    )


def calculate_polyhedral_network_vertex_gaussian_curvatures(
        net, surface_ids, *, tolerance=1e-8, gauss_bonnet_tolerance=1e-6):
    """Compatibility-named entry point for network/group callers."""
    return calculate_polyhedral_gaussian_curvature(
        net,
        surface_ids,
        tolerance=tolerance,
        gauss_bonnet_tolerance=gauss_bonnet_tolerance,
    )


__all__ = [
    "POLYHEDRAL_STATUS_CALCULATED",
    "POLYHEDRAL_STATUS_UNAVAILABLE",
    "PolyhedralGaussianCurvatureResult",
    "calculate_polyhedral_gaussian_curvature",
    "calculate_polyhedral_network_vertex_gaussian_curvatures",
]
