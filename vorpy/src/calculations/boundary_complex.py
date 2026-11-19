"""Topology of pairwise boundary surfaces that can contain disconnected disks."""
from collections import Counter, defaultdict

import numpy as np


def boundary_complex(net, surface_ids):
    """Count disk patches and solve their edge-orientation constraints.

    Multiple contours are accepted as separate disks only when the saved face
    triangulation independently confirms that topology. Annular/unsupported
    faces remain unresolved; no hole is filled to obtain a passing result.
    """
    incidences = defaultdict(list)
    face_euler = {}
    unresolved = {}
    patches = set()
    for surface_id in sorted(surface_ids):
        surface = net.surfs.iloc[surface_id]
        edge_ids = set(map(int, surface.get('edges', ())))
        by_vertex = defaultdict(list)
        ends = {}
        for edge_id in edge_ids:
            try:
                if not 0 <= edge_id < len(net.edges):
                    raise IndexError(edge_id)
                vertices = tuple(map(int, net.edges.iloc[edge_id].get('verts', ())))
            except (IndexError, TypeError, ValueError):
                unresolved[surface_id] = 'invalid face edge reference'
                continue
            if len(vertices) != 2 or min(vertices) < 0 or vertices[0] == vertices[1]:
                unresolved[surface_id] = 'edge does not have two distinct endpoint vertices'
                continue
            ends[edge_id] = vertices
            for vertex in vertices:
                by_vertex[vertex].append(edge_id)
        if (not edge_ids or len(ends) != len(edge_ids)
                or any(len(edges) != 2 for edges in by_vertex.values())):
            unresolved[surface_id] = 'face contours are not closed edge cycles'
            continue
        loops = []
        unseen = set(edge_ids)
        while unseen:
            first = min(unseen)
            start, current = ends[first]
            loop = [(first, start, current)]
            unseen.remove(first)
            while current != start:
                next_edge = next((edge for edge in by_vertex[current] if edge in unseen), None)
                if next_edge is None:
                    unresolved[surface_id] = 'face contour cannot be traversed'
                    break
                following = next(v for v in ends[next_edge] if v != current)
                loop.append((next_edge, current, following))
                unseen.remove(next_edge)
                current = following
            loops.append(loop)
        if surface_id in unresolved:
            continue
        if len(loops) > 1:
            triangles = np.asarray(surface.get('tris', ()), dtype=int)
            if triangles.ndim != 2 or triangles.shape[1] != 3 or not len(triangles):
                unresolved[surface_id] = 'multiple contours lack face triangulation'
                continue
            mesh_edges = Counter(tuple(sorted((int(a), int(b)))) for tri in triangles
                                 for a, b in zip(tri, np.roll(tri, -1)))
            mesh_vertices = set(map(int, triangles.ravel()))
            adjacency = {v: set() for v in mesh_vertices}
            for a, b in mesh_edges:
                adjacency[a].add(b)
                adjacency[b].add(a)
            components = []
            remaining = set(mesh_vertices)
            while remaining:
                stack = [remaining.pop()]
                component = set(stack)
                while stack:
                    vertex = stack.pop()
                    for neighbor in adjacency[vertex] & remaining:
                        remaining.remove(neighbor)
                        component.add(neighbor)
                        stack.append(neighbor)
                components.append(component)
            disks = len(components) == len(loops) and all(count <= 2 for count in mesh_edges.values())
            for component in components:
                edges = {edge: count for edge, count in mesh_edges.items() if edge[0] in component}
                faces = [tri for tri in triangles if int(tri[0]) in component]
                boundary_degree = Counter(v for edge, count in edges.items() if count == 1 for v in edge)
                disks = disks and (len(component) - len(edges) + len(faces) == 1)
                disks = disks and bool(boundary_degree) and all(count == 2 for count in boundary_degree.values())
                # Euler totals alone can hide a sphere plus an annulus or a pinched vertex.
                for vertex in component:
                    link = defaultdict(set)
                    for triangle in faces:
                        if vertex in triangle:
                            opposite = [int(v) for v in triangle if v != vertex]
                            if len(opposite) != 2:
                                disks = False
                                continue
                            first, second = opposite
                            link[first].add(second)
                            link[second].add(first)
                    unseen_link = set(link)
                    stack = [unseen_link.pop()] if unseen_link else []
                    while stack:
                        for neighbor in link[stack.pop()] & unseen_link:
                            unseen_link.remove(neighbor)
                            stack.append(neighbor)
                    degrees = Counter(map(len, link.values()))
                    disks = disks and not unseen_link and set(degrees).issubset({1, 2})
                    disks = disks and degrees[1] == (2 if vertex in boundary_degree else 0)
            if not disks:
                unresolved[surface_id] = 'multi-contour face is not verified as disjoint disks'
                continue
        face_euler[surface_id] = len(loops)
        for number, loop in enumerate(loops):
            patch = (surface_id, number)
            patches.add(patch)
            for edge, first, second in loop:
                incidences[edge].append((patch, 1 if first < second else -1))
    adjacency = {patch: [] for patch in patches}
    for values in incidences.values():
        if len(values) != 2:
            continue
        (first, sign_a), (second, sign_b) = values
        relative = -sign_a * sign_b
        adjacency[first].append((second, relative))
        adjacency[second].append((first, relative))
    signs = {}
    components = 0
    orientable = not unresolved and bool(patches) and all(len(v) == 2 for v in incidences.values())
    for patch in sorted(patches):
        if patch in signs:
            continue
        components += 1
        signs[patch] = 1
        stack = [patch]
        while stack:
            first = stack.pop()
            for second, relative in adjacency[first]:
                expected = signs[first] * relative
                if second in signs:
                    orientable = orientable and signs[second] == expected
                else:
                    signs[second] = expected
                    stack.append(second)
    return dict(face_euler=face_euler, components=components, orientable=orientable,
                unresolved_faces=unresolved, patch_count=len(patches))
