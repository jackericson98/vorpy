import numpy as np


def intrinsic_corner_angle(first_tangent, second_tangent, tol=1e-12):
    """Return the intrinsic angle between two face-boundary tangents.

    Both tangents must point away from the common vertex along the two
    boundary edges of the same smooth surface patch.
    """
    first = np.asarray(first_tangent, dtype=float)
    second = np.asarray(second_tangent, dtype=float)

    if first.shape != (3,) or second.shape != (3,):
        raise ValueError("Corner tangents must be three-vectors.")

    if not np.all(np.isfinite(first)) or not np.all(np.isfinite(second)):
        raise ValueError("Corner tangents must be finite.")

    first_norm = float(np.linalg.norm(first))
    second_norm = float(np.linalg.norm(second))

    if first_norm < tol or second_norm < tol:
        raise ValueError("Corner tangent is undefined.")

    first /= first_norm
    second /= second_norm

    cosine = float(np.clip(
        np.dot(first, second),
        -1.0,
        1.0,
    ))

    return float(np.arccos(cosine))


def angular_defect(face_angles):
    """Return Gaussian-curvature angular defect at a closed cell vertex."""
    angles = np.asarray(face_angles, dtype=float)

    if angles.ndim != 1 or len(angles) < 3:
        raise ValueError(
            "A closed three-dimensional cell vertex requires "
            "at least three incident face angles."
        )

    if not np.all(np.isfinite(angles)):
        raise ValueError("Face angles must be finite.")

    if np.any(angles <= 0.0) or np.any(angles >= np.pi):
        raise ValueError(
            "Interior face angles must lie strictly between 0 and pi."
        )

    return float(2.0 * np.pi - np.sum(angles))


def outward_edge_tangent_at_vertex(
        edge_geometry,
        vertex_position,
        tol=1e-12):
    """Return the unit tangent pointing away from one endpoint of an edge."""

    if vertex_position == 0:
        tangent = np.asarray(
            edge_geometry.unit_tangent(edge_geometry.t_min, tol=tol),
            dtype=float,
        )
    elif vertex_position == 1:
        tangent = -np.asarray(
            edge_geometry.unit_tangent(edge_geometry.t_max, tol=tol),
            dtype=float,
        )
    else:
        raise ValueError(
            "vertex_position must be 0 for the first endpoint or 1 for the second."
        )

    magnitude = float(np.linalg.norm(tangent))

    if not np.isfinite(magnitude) or magnitude < tol:
        raise ValueError("Edge tangent is undefined at the vertex.")

    return tangent / magnitude


def analytic_edge_tangent_away_from_network_vertex(
        edge_geometry,
        vertex_location,
        tolerance=1e-7):
    """Return the analytic edge tangent pointing away from a network vertex.

    Endpoint identity is determined geometrically rather than from stored
    vertex ordering because analytic edge resolution may canonicalize the
    parameter direction.
    """
    vertex_location = np.asarray(vertex_location, dtype=float)

    if vertex_location.shape != (3,):
        raise ValueError("Vertex location must be a three-vector.")

    start = np.asarray(
        edge_geometry.point(edge_geometry.t_min),
        dtype=float,
    )
    end = np.asarray(
        edge_geometry.point(edge_geometry.t_max),
        dtype=float,
    )

    start_error = float(np.linalg.norm(start - vertex_location))
    end_error = float(np.linalg.norm(end - vertex_location))

    if start_error <= tolerance and end_error <= tolerance:
        raise ValueError(
            "Analytic edge endpoints are indistinguishable at this vertex."
        )

    if start_error <= tolerance:
        return outward_edge_tangent_at_vertex(
            edge_geometry,
            0,
        )

    if end_error <= tolerance:
        return outward_edge_tangent_at_vertex(
            edge_geometry,
            1,
        )

    raise ValueError(
        "Network vertex does not match either endpoint of analytic edge: "
        f"start error={start_error:.6g}, end error={end_error:.6g}"
    )


def aw_vertex_face_angle(
        net,
        vertex_index,
        cell_index,
        surface_index,
        resolved_edges=None,
        tolerance=1e-7):
    """Return one intrinsic face angle for one AW cell at one network vertex.

    The selected pairwise surface must belong to the selected cell and
    exactly two of its boundary edges must meet at the selected vertex.
    """
    vertex_index = int(vertex_index)
    cell_index = int(cell_index)
    surface_index = int(surface_index)

    vertex = net.verts.loc[vertex_index]
    surface = net.surfs.loc[surface_index]

    surface_balls = tuple(int(value) for value in surface["balls"])

    if len(surface_balls) != 2:
        raise ValueError(
            f"Surface {surface_index} does not have exactly two generators."
        )

    if cell_index not in surface_balls:
        raise ValueError(
            f"Surface {surface_index} does not bound cell {cell_index}."
        )

    vertex_edges = {
        int(value)
        for value in vertex["edges"]
    }

    surface_edges = {
        int(value)
        for value in surface["edges"]
    }

    incident_edges = sorted(
        vertex_edges.intersection(surface_edges)
    )

    if len(incident_edges) != 2:
        raise ValueError(
            f"Surface {surface_index} at vertex {vertex_index} has "
            f"{len(incident_edges)} incident boundary edges; expected 2."
        )

    vertex_location = np.asarray(
        vertex["loc"],
        dtype=float,
    )

    tangents = []

    for edge_index in incident_edges:
        if resolved_edges is not None and edge_index in resolved_edges:
            resolved = resolved_edges[edge_index]
        else:
            # Local import avoids coupling vertex_geometry.py to Network
            # during package initialization.
            from vorpy.src.network.edge_geometry_diagnostics import (
                resolve_aw_network_edge,
            )

            resolved = resolve_aw_network_edge(
                net,
                edge_index,
                tolerance=tolerance,
            )

        tangent = analytic_edge_tangent_away_from_network_vertex(
            resolved.geometry,
            vertex_location,
            tolerance=tolerance,
        )

        tangents.append(tangent)

    return intrinsic_corner_angle(
        tangents[0],
        tangents[1],
    )


def aw_cell_vertex_angular_defect(
        net,
        vertex_index,
        cell_index,
        resolved_edges=None,
        tolerance=1e-7):
    """Return the Gaussian-curvature angular defect for one AW cell vertex."""
    vertex_index = int(vertex_index)
    cell_index = int(cell_index)

    vertex = net.verts.loc[vertex_index]

    vertex_balls = {
        int(value)
        for value in vertex["balls"]
    }

    if cell_index not in vertex_balls:
        raise ValueError(
            f"Cell {cell_index} does not meet vertex {vertex_index}."
        )

    incident_surfaces = []

    for surface_index in vertex["surfs"]:
        surface_index = int(surface_index)
        surface_balls = {
            int(value)
            for value in net.surfs.loc[surface_index, "balls"]
        }

        if cell_index in surface_balls:
            incident_surfaces.append(surface_index)

    if len(incident_surfaces) < 3:
        raise ValueError(
            f"Cell {cell_index} has only {len(incident_surfaces)} "
            f"incident surfaces at vertex {vertex_index}; expected at least 3."
        )

    face_angles = [
        aw_vertex_face_angle(
            net=net,
            vertex_index=vertex_index,
            cell_index=cell_index,
            surface_index=surface_index,
            resolved_edges=resolved_edges,
            tolerance=tolerance,
        )
        for surface_index in incident_surfaces
    ]

    return angular_defect(face_angles)


def spherical_triangle_area(first_normal, second_normal, third_normal, tol=1e-12):
    """Return the oriented area of a spherical triangle on the unit sphere.

    The inputs are three nonzero surface normals. After normalization, the
    signed solid-angle formula gives the spherical area enclosed by their
    shorter great-circle arcs.
    """
    normals = []

    for normal in (first_normal, second_normal, third_normal):
        normal = np.asarray(normal, dtype=float)

        if normal.shape != (3,) or not np.all(np.isfinite(normal)):
            raise ValueError("Surface normals must be finite three-vectors.")

        magnitude = float(np.linalg.norm(normal))
        if magnitude < tol:
            raise ValueError("Surface normal is undefined.")

        normals.append(normal / magnitude)

    a, b, c = normals

    numerator = float(np.dot(a, np.cross(b, c)))
    denominator = float(
        1.0
        + np.dot(a, b)
        + np.dot(b, c)
        + np.dot(c, a)
    )

    return float(2.0 * np.arctan2(numerator, denominator))


def aw_cell_surface_normal_at_vertex(
        vertex_location,
        cell_index,
        other_index,
        generator_indices,
        generator_locations,
        tol=1e-12):
    """Return the pairwise AW surface normal pointing outward from one cell.

    For the interface between cells i and j, evaluate the two sphere radial
    unit vectors at the actual vertex. The cell-i outward normal is

        n_i->j = normalize(u_i - u_j)

    where u_i points from generator i to the vertex.

    Reversing the cell perspective reverses the normal exactly.
    """
    vertex_location = np.asarray(vertex_location, dtype=float)
    indices = tuple(int(value) for value in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)

    if vertex_location.shape != (3,):
        raise ValueError("Vertex location must be a three-vector.")

    if locations.shape != (len(indices), 3):
        raise ValueError(
            "Generator locations must contain one three-vector per index."
        )

    try:
        cell_pos = indices.index(int(cell_index))
        other_pos = indices.index(int(other_index))
    except ValueError as exc:
        raise ValueError(
            "Both surface generators must appear in generator_indices."
        ) from exc

    if cell_pos == other_pos:
        raise ValueError("A pairwise surface requires two distinct generators.")

    cell_radial = vertex_location - locations[cell_pos]
    other_radial = vertex_location - locations[other_pos]

    cell_distance = float(np.linalg.norm(cell_radial))
    other_distance = float(np.linalg.norm(other_radial))

    if cell_distance < tol or other_distance < tol:
        raise ValueError("Vertex coincides with a surface generator.")

    cell_radial /= cell_distance
    other_radial /= other_distance

    normal = cell_radial - other_radial
    magnitude = float(np.linalg.norm(normal))

    if not np.isfinite(magnitude) or magnitude < tol:
        raise ValueError("AW pairwise surface normal is undefined at vertex.")

    return normal / magnitude


def aw_vertex_cell_surface_normal(
        net,
        vertex_index,
        cell_index,
        surface_index):
    """Return one incident surface normal outward from a selected AW cell."""
    vertex_index = int(vertex_index)
    cell_index = int(cell_index)
    surface_index = int(surface_index)

    vertex = net.verts.loc[vertex_index]
    surface = net.surfs.loc[surface_index]

    surface_balls = tuple(int(value) for value in surface["balls"])

    if len(surface_balls) != 2:
        raise ValueError(
            f"Surface {surface_index} does not have exactly two generators."
        )

    if cell_index not in surface_balls:
        raise ValueError(
            f"Surface {surface_index} does not bound cell {cell_index}."
        )

    other_index = (
        surface_balls[1]
        if surface_balls[0] == cell_index
        else surface_balls[0]
    )

    vertex_location = np.asarray(vertex["loc"], dtype=float)

    # Ball numbers are not guaranteed to equal DataFrame row positions.
    ball_num_to_pos = {
        int(num): pos
        for pos, num in enumerate(net.balls["num"].to_numpy())
    }

    if cell_index not in ball_num_to_pos or other_index not in ball_num_to_pos:
        raise ValueError(
            f"Could not locate generators for surface {surface_index}."
        )

    cell_location = np.asarray(
        net.balls.iloc[ball_num_to_pos[cell_index]]["loc"],
        dtype=float,
    )
    other_location = np.asarray(
        net.balls.iloc[ball_num_to_pos[other_index]]["loc"],
        dtype=float,
    )

    return aw_cell_surface_normal_at_vertex(
        vertex_location=vertex_location,
        cell_index=cell_index,
        other_index=other_index,
        generator_indices=[cell_index, other_index],
        generator_locations=np.asarray(
            [cell_location, other_location],
            dtype=float,
        ),
    )


def aw_cell_vertex_gauss_map_normals(
        net,
        vertex_index,
        cell_index):
    """Return outward normals of all cell faces incident to one AW vertex."""
    vertex_index = int(vertex_index)
    cell_index = int(cell_index)
    vertex = net.verts.loc[vertex_index]

    incident_surfaces = []

    for surface_index in vertex["surfs"]:
        surface_index = int(surface_index)
        surface_balls = {
            int(value)
            for value in net.surfs.loc[surface_index, "balls"]
        }

        if cell_index in surface_balls:
            incident_surfaces.append(surface_index)

    if len(incident_surfaces) < 3:
        raise ValueError(
            f"Cell {cell_index} has only {len(incident_surfaces)} "
            f"incident surfaces at vertex {vertex_index}."
        )

    normals = [
        aw_vertex_cell_surface_normal(
            net,
            vertex_index,
            cell_index,
            surface_index,
        )
        for surface_index in incident_surfaces
    ]

    return incident_surfaces, np.asarray(normals, dtype=float)


def aw_cell_vertex_surface_cycle(net, vertex_index, cell_index):
    """Return incident cell surfaces in their topological cycle around a vertex.

    Two surfaces are adjacent in the cycle when they share an incident edge
    that also belongs to the selected cell. No reliance is placed on the
    storage order of ``vertex['surfs']``.
    """
    vertex_index = int(vertex_index)
    cell_index = int(cell_index)
    vertex = net.verts.loc[vertex_index]

    incident_surfaces = [
        int(surface_index)
        for surface_index in vertex["surfs"]
        if cell_index in {int(v) for v in net.surfs.loc[int(surface_index), "balls"]}
    ]

    if len(incident_surfaces) < 3:
        raise ValueError(
            f"Cell {cell_index} has only {len(incident_surfaces)} "
            f"incident surfaces at vertex {vertex_index}."
        )

    surface_set = set(incident_surfaces)
    adjacency = {surface_index: set() for surface_index in incident_surfaces}
    shared_edges = {}

    for edge_index in vertex["edges"]:
        edge_index = int(edge_index)
        edge = net.edges.loc[edge_index]

        if cell_index not in {int(v) for v in edge["balls"]}:
            continue

        edge_surfaces = [
            int(surface_index)
            for surface_index in edge["surfs"]
            if int(surface_index) in surface_set
        ]

        if len(edge_surfaces) != 2:
            raise ValueError(
                f"Cell {cell_index}, vertex {vertex_index}, edge {edge_index} "
                f"touches {len(edge_surfaces)} incident cell surfaces; expected 2."
            )

        first, second = edge_surfaces
        adjacency[first].add(second)
        adjacency[second].add(first)
        shared_edges[frozenset((first, second))] = edge_index

    # A closed local manifold must form one cycle: every face has two neighbors.
    bad = {
        surface_index: len(neighbors)
        for surface_index, neighbors in adjacency.items()
        if len(neighbors) != 2
    }

    if bad:
        raise ValueError(
            f"Cell {cell_index} at vertex {vertex_index} does not form a "
            f"closed surface cycle; face degrees={bad}."
        )

    # Deterministic starting face. The traversal direction is deliberately
    # unresolved here; orientation is fixed in the next step.
    start = min(incident_surfaces)
    first_neighbor = min(adjacency[start])

    cycle = [start, first_neighbor]
    previous, current = start, first_neighbor

    while True:
        candidates = adjacency[current] - {previous}

        if len(candidates) != 1:
            raise ValueError(
                f"Ambiguous surface cycle for cell {cell_index} "
                f"at vertex {vertex_index}."
            )

        next_surface = next(iter(candidates))

        if next_surface == start:
            break

        if next_surface in cycle:
            raise ValueError(
                f"Surface cycle closes incorrectly for cell {cell_index} "
                f"at vertex {vertex_index}."
            )

        cycle.append(next_surface)
        previous, current = current, next_surface

    if len(cycle) != len(incident_surfaces):
        raise ValueError(
            f"Surface cycle for cell {cell_index} at vertex {vertex_index} "
            f"contains {len(cycle)} of {len(incident_surfaces)} faces."
        )

    return cycle, shared_edges


def orient_aw_cell_vertex_surface_cycle(
        net,
        vertex_index,
        cell_index,
        cycle,
        shared_edges,
        resolved_edges=None,
        tolerance=1e-7):
    """Orient a cell's face cycle consistently with its physical edge geometry."""
    if len(cycle) < 3:
        raise ValueError("A vertex surface cycle requires at least three faces.")

    first_surface = int(cycle[0])
    second_surface = int(cycle[1])

    first_normal = aw_vertex_cell_surface_normal(
        net, vertex_index, cell_index, first_surface
    )
    second_normal = aw_vertex_cell_surface_normal(
        net, vertex_index, cell_index, second_surface
    )

    edge_key = frozenset((first_surface, second_surface))

    if edge_key not in shared_edges:
        raise ValueError(
            f"No shared edge found for surfaces {first_surface} and "
            f"{second_surface} at vertex {vertex_index}."
        )

    edge_index = shared_edges[edge_key]

    if resolved_edges is not None and edge_index in resolved_edges:
        resolved = resolved_edges[edge_index]
    else:
        from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge

        resolved = resolve_aw_network_edge(
            net,
            edge_index,
            tolerance=tolerance,
        )

    vertex_location = np.asarray(net.verts.loc[vertex_index, "loc"], dtype=float)

    edge_tangent = analytic_edge_tangent_away_from_network_vertex(
        resolved.geometry,
        vertex_location,
        tolerance=tolerance,
    )

    orientation = float(
        np.dot(
            np.cross(first_normal, second_normal),
            edge_tangent,
        )
    )

    if abs(orientation) < tolerance:
        raise ValueError(
            f"Cannot orient surface cycle for cell {cell_index} at "
            f"vertex {vertex_index}: degenerate normal/tangent orientation."
        )

    # For outward cell normals, positive Gauss-map orientation has
    # (n_i x n_j) opposite the edge tangent pointing away from the vertex.
    if orientation > 0.0:
        cycle = [cycle[0], *reversed(cycle[1:])]

    return cycle


def aw_cell_vertex_gauss_map_area(
        net,
        vertex_index,
        cell_index,
        resolved_edges=None,
        tolerance=1e-7):
    """Return the signed singular Gaussian curvature at one AW cell vertex."""
    cycle, shared_edges = aw_cell_vertex_surface_cycle(
        net,
        vertex_index,
        cell_index,
    )

    cycle = orient_aw_cell_vertex_surface_cycle(
        net,
        vertex_index,
        cell_index,
        cycle,
        shared_edges,
        resolved_edges=resolved_edges,
        tolerance=tolerance,
    )

    normals = [
        aw_vertex_cell_surface_normal(
            net,
            vertex_index,
            cell_index,
            surface_index,
        )
        for surface_index in cycle
    ]

    if len(normals) == 3:
        area = spherical_triangle_area(
            normals[0],
            normals[1],
            normals[2],
        )
    else:
        # Higher-order degenerate vertices form a spherical polygon. Fan
        # triangulation is valid once the cyclic orientation is established.
        area = 0.0

        for i in range(1, len(normals) - 1):
            area += spherical_triangle_area(
                normals[0],
                normals[i],
                normals[i + 1],
            )

    return float(area)