"""Full bicolor interfaces derived from solved VorPy Voronoi-side features."""

from __future__ import annotations

import math

from .model import (
    InterfaceEdge,
    InterfaceRepresentation,
    InterfaceSurface,
    InterfaceVertex,
)


def _ids(value):
    if value is None:
        return ()
    try:
        return tuple(int(item) for item in value)
    except (TypeError, ValueError):
        return ()


def _finite_area(value):
    try:
        area = float(value)
    except (TypeError, ValueError):
        return None
    return area if math.isfinite(area) and area >= 0.0 else None


def _bicolor(ids, group_a, group_b):
    return (len(ids) == 2 and
            ((ids[0] in group_a and ids[1] in group_b)
             or (ids[1] in group_a and ids[0] in group_b)))


def build_full_interface(network, dual, group_a, group_b):
    """Return all valid bicolor surfaces for one already-solved network.

    Eligibility is geometric: a surface must have a dual pair mapping, finite
    nonnegative area, bounded/complete primal feature geometry, and supported
    incidence. AW additionally requires complete participating generator
    cells because an incomplete AW clipping neighborhood can truncate a patch.
    No alpha or Cazals M filtering is applied.
    """
    group_a, group_b = frozenset(map(int, group_a)), frozenset(map(int, group_b))
    if not group_a or not group_b or group_a & group_b:
        raise ValueError("Interface groups must be nonempty and disjoint.")
    valid = {int(value) for value in network.balls.index}
    if not (group_a | group_b) <= valid:
        raise ValueError("Interface groups reference an unknown network generator.")
    if getattr(dual, "network", None) is not network:
        raise ValueError("DualComplex must reference the same solved network object.")

    surfaces = []
    candidates = {}
    selected = set()
    for surface_id, row in network.surfs.iterrows():
        surface_id = int(surface_id)
        generator_ids = _ids(row.get("balls", ()))
        if not _bicolor(generator_ids, group_a, group_b):
            continue
        ids = tuple(sorted(generator_ids))
        simplex = dual.simplex(1, ids)
        feature = next((ref for ref in (simplex.primal_features if simplex else ())
                        if ref.kind == "surf" and ref.feature_id == surface_id), None)
        bounded = bool(feature and feature.bounded)
        complete = bool(feature and feature.complete)
        cells_complete = bool(feature and feature.generator_cells_complete)
        supported = bool(simplex and simplex.supported and feature)
        area = _finite_area(row.get("sa"))
        if not feature or not supported:
            status, reason = "unresolved", "surface has no supported pair-simplex incidence mapping"
        elif not bounded or not complete:
            status, reason = "excluded_incomplete", "surface primal geometry is unbounded or incomplete"
        elif area is None:
            status, reason = "excluded_area_unavailable", "surface area is missing or non-finite"
        elif dual.scheme == "aw" and not cells_complete:
            status, reason = "excluded_incomplete_cells", "participating AW generator cells are incomplete"
        else:
            status, reason = "selected", ""
        is_selected = status == "selected"
        record = InterfaceSurface(
            surface_id, ids, simplex.simplex_id if simplex else None, area,
            bounded, complete, cells_complete, supported, is_selected,
            status, reason, _ids(row.get("edges", ())), _ids(row.get("verts", ())),
        )
        surfaces.append(record)
        candidates[surface_id] = record
        if is_selected:
            selected.add(surface_id)

    edges = []
    relevant_edge_ids = set()
    for edge_id, row in network.edges.iterrows():
        edge_id = int(edge_id)
        generator_ids = _ids(row.get("balls", ()))
        if len(generator_ids) != 3 or not (set(generator_ids) & group_a and set(generator_ids) & group_b):
            continue
        relevant_edge_ids.add(edge_id)
        incident_set = set(_ids(row.get("surfs", ()))) & set(candidates)
        if not incident_set:
            incident_set = {surface_id for surface_id, surface in candidates.items()
                            if edge_id in surface.boundary_edge_ids}
        incident = tuple(sorted(incident_set))
        selected_incident = tuple(value for value in incident if value in selected)
        simplex = dual.simplex(2, generator_ids)
        feature = next((ref for ref in (simplex.primal_features if simplex else ())
                        if ref.kind == "edge" and ref.feature_id == edge_id), None)
        bounded = bool(feature and feature.bounded)
        complete = bool(feature and feature.complete)
        supported = bool(simplex and simplex.supported and feature)
        if not feature or not supported:
            status, reason = "unresolved", "edge has no supported triple-simplex incidence mapping"
        elif not bounded or not complete:
            status, reason = "excluded_incomplete", "Voronoi edge is unbounded or incomplete"
        elif len(selected_incident) > 2:
            status, reason = "nonmanifold", "more than two selected interface surfaces meet this edge"
        elif len(selected_incident) == 2:
            status, reason = "shared", "edge is shared by two selected interface surfaces"
        elif len(selected_incident) == 1:
            status, reason = "boundary", "edge bounds one selected interface surface"
        else:
            status, reason = "outside_selected_interface", "no selected interface surface is incident"
        edges.append(InterfaceEdge(
            edge_id, tuple(sorted(generator_ids)), simplex.simplex_id if simplex else None,
            incident, selected_incident, bounded, complete, supported, status, reason,
            _ids(row.get("verts", ())),
        ))

    used_edge_ids = {edge.feature_id for edge in edges if edge.incident_selected_surface_ids}
    surface_vertices = {vertex_id for surface in surfaces if surface.selected
                        for vertex_id in surface.boundary_vertex_ids}
    edge_by_id = {edge.feature_id: edge for edge in edges}
    vertices = []
    for vertex_id, row in network.verts.iterrows():
        vertex_id = int(vertex_id)
        incident_surface_ids = set(_ids(row.get("surfs", ()))) & set(candidates)
        if not incident_surface_ids:
            incident_surface_ids = {surface_id for surface_id, surface in candidates.items()
                                    if vertex_id in surface.boundary_vertex_ids}
        incident_surfaces = tuple(sorted(incident_surface_ids))
        selected_incident_ids = set(incident_surfaces) & selected
        incident_edges = set(_ids(row.get("edges", ())))
        if not incident_edges:
            incident_edges = {edge_id for edge_id, edge in edge_by_id.items()
                              if vertex_id in edge.endpoint_vertex_ids}
        for edge_id in incident_edges & used_edge_ids:
            selected_incident_ids.update(edge_by_id[edge_id].incident_selected_surface_ids)
        selected_incident = tuple(sorted(selected_incident_ids))
        if not selected_incident and vertex_id not in surface_vertices and not (incident_edges & used_edge_ids):
            continue
        generator_ids = _ids(row.get("balls", ()))
        simplex = dual.simplex(3, generator_ids) if len(generator_ids) == 4 else None
        feature = next((ref for ref in (simplex.primal_features if simplex else ())
                        if ref.kind == "vert" and ref.feature_id == vertex_id), None)
        bounded = bool(feature and feature.bounded)
        supported = bool(simplex and simplex.supported and feature)
        status = "supported" if bounded and supported else "unresolved_or_incomplete"
        point = row.get("loc")
        try:
            xyz = tuple(float(value) for value in point)
            if len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
                xyz = None
        except (TypeError, ValueError):
            xyz = None
        vertices.append(InterfaceVertex(
            vertex_id, tuple(sorted(generator_ids)), incident_surfaces, selected_incident,
            bounded, supported, status, xyz,
        ))

    return InterfaceRepresentation(
        dual.scheme, "full", None, group_a, group_b, surfaces, edges, vertices,
        metadata={
            "dual_incidence_source": "same solved VorPy network",
            "physical_geometry_source": "network.surfs/network.edges/network.verts",
            "alpha_or_M_filter_applied": False,
            "area_source": "stored network surface area field sa",
            "euler_characteristic_note": "incidence-complex value only; no closed-surface or genus inference",
            "unrepresented_bicolor_edges": sorted(relevant_edge_ids - {edge.feature_id for edge in edges}),
        },
    )
