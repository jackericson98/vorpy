"""Deterministic analytic records for a completed additively weighted network.

This module is a *producer* adapter.  It does not change a ``Network``, does
not solve missing topology, and does not turn sampled meshes into analytic
geometry.  Its purpose is to give downstream measurement code one explicit
Atoms -> Cells -> Surfaces -> Edges -> Vertices incidence envelope.

Lengths are in the coordinate unit used by :class:`~vorpy.src.network.network.Network`
(Å for ordinary molecular inputs).  Areas and curvature integrals are not
calculated here.  GAUSS owns those measurements; this adapter only preserves
the geometry, orientation convention, provenance, and availability state that
they require.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from vorpy.src.calculations.edge_geometry import (
    AdditivelyWeightedTrisectorBranch,
    AdditivelyWeightedTrisectorConic,
    LineEdgeGeometry,
)
from vorpy.src.calculations.surf import calc_surf_func
from vorpy.src.calculations.surface_geometry import QuadraticSurfaceGeometry
from vorpy.src.network.aw_geometry_reuse import (
    AWGeometryReuseError,
    _digest,
    _stable_generator_id,
    build_solve_universe_key,
)
from vorpy.src.network.edge_geometry_diagnostics import (
    AWNetworkEdgeResolutionError,
    resolve_aw_network_edge,
)
from vorpy.src.results.model import Status


SCHEMA_VERSION = 1
LENGTH_UNITS = "Å"


class AWGeometryEnvelopeError(ValueError):
    """Raised when a network cannot be read as an AW geometry producer."""


@dataclass(frozen=True)
class AWAtomRecord:
    """Input atom/generator and its one-to-one AW cell."""

    atom_id: str
    cell_id: str
    source_index: int
    stable_generator_id: str
    center: tuple[float, float, float] | None
    radius: float | None
    status: Status
    identity_is_input_stable: bool


@dataclass(frozen=True)
class AWCellRecord:
    """AW cell identity and the primitives incident to its generator."""

    cell_id: str
    atom_id: str
    source_index: int
    complete: bool | None
    surface_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    vertex_ids: tuple[str, ...]
    status: Status


@dataclass(frozen=True)
class AWJunctionRecord:
    """A four-generator AW vertex, or an explicit unresolved vertex record."""

    vertex_id: str
    source_index: int
    cell_ids: tuple[str, ...]
    location: tuple[float, float, float] | None
    clearance_radius: float | None
    edge_ids: tuple[str, ...]
    surface_ids: tuple[str, ...]
    component_id: str | None
    status: Status
    reason: str | None = None


@dataclass(frozen=True)
class AWAnalyticEdgeRecord:
    """One exact AW edge, or an explicit unresolved non-analytic record."""

    edge_id: str
    source_index: int
    cell_ids: tuple[str, ...]
    vertex_ids: tuple[str, ...]
    surface_ids: tuple[str, ...]
    analytic: Mapping[str, Any] | None
    component_id: str | None
    status: Status
    reason: str | None = None


@dataclass(frozen=True)
class AWSurfaceOwner:
    """One of the two directed cell views of a pairwise AW surface.

    ``normal_convention`` evaluates the exact clearance-gradient normal from
    ``cell_id`` toward ``neighbor_cell_id``.  A normal is intentionally not
    sampled here because it varies over a non-planar AW surface.
    """

    cell_id: str
    neighbor_cell_id: str
    normal_convention: str


@dataclass(frozen=True)
class AWAnalyticSurfaceRecord:
    """Pairwise AW quadratic surface with explicit shared ownership."""

    surface_id: str
    source_index: int
    cell_ids: tuple[str, ...]
    vertex_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    coefficients: tuple[float, ...] | None
    geometry_kind: str | None
    coefficient_source: str | None
    owners: tuple[AWSurfaceOwner, ...]
    component_id: str | None
    status: Status
    # AW pairwise boundaries are planes or quadrics, not atom-carried
    # spherical molecular patches.  The explicit state prevents callers from
    # silently reinterpreting an AW face as a spherical patch.
    spherical_patch_status: Status = Status.NOT_SUPPORTED
    reason: str | None = None


@dataclass(frozen=True)
class AWGeometryComponent:
    """Stable connected component identity for primitive incidence records."""

    component_id: str
    vertex_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    surface_ids: tuple[str, ...]
    boundary_incidence_status: Status
    boundary_issues: tuple[str, ...]


@dataclass(frozen=True)
class AWGeometryCompleteness:
    """Producer evidence only; it never certifies GAUSS measurements."""

    status: Status
    cell_completeness_known: bool
    globally_complete: bool | None
    solver_completeness_proven: bool
    counts: Mapping[str, int]
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class AWGeometryEnvelope:
    """Immutable analytic AW construction envelope for downstream consumers."""

    schema_version: int
    coordinate_units: str
    solve_key_digest: str
    atoms: tuple[AWAtomRecord, ...]
    cells: tuple[AWCellRecord, ...]
    surfaces: tuple[AWAnalyticSurfaceRecord, ...]
    edges: tuple[AWAnalyticEdgeRecord, ...]
    vertices: tuple[AWJunctionRecord, ...]
    components: tuple[AWGeometryComponent, ...]
    completeness: AWGeometryCompleteness
    provenance: Mapping[str, Any]


def _frozen_mapping(value: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(value))


def _as_vector(value: Any) -> tuple[float, float, float] | None:
    try:
        vector = tuple(float(item) for item in value)
    except (TypeError, ValueError):
        return None
    if len(vector) != 3 or not all(np.isfinite(vector)):
        return None
    return vector


def _as_indices(value: Any) -> tuple[int, ...] | None:
    try:
        result = tuple(int(item) for item in value)
    except (TypeError, ValueError):
        return None
    return result


def _input_stable_identity(row: Mapping[str, Any]) -> bool:
    return any(row.get(field) is not None and str(row.get(field))
               for field in ("stable_id", "atom_stable_id", "generator_id"))


def _record_id(kind: str, solve_key_digest: str, payload: Mapping[str, Any]) -> str:
    return f"aw-envelope:{kind}:{_digest({'solve': solve_key_digest, 'payload': payload})}"


def _require_aw_network(network: Any) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame | None, pd.DataFrame | None]:
    settings = getattr(network, "settings", None) or {}
    if str(settings.get("net_type", "aw")).lower() != "aw":
        raise AWGeometryEnvelopeError("The analytic geometry envelope supports AW networks only.")
    balls = getattr(network, "balls", None)
    if not isinstance(balls, pd.DataFrame) or not {"loc", "rad"}.issubset(balls.columns):
        raise AWGeometryEnvelopeError("AW geometry requires ball locations and radii.")

    def table(name: str) -> pd.DataFrame | None:
        value = getattr(network, name, None)
        if value is None:
            return None
        if not isinstance(value, pd.DataFrame):
            raise AWGeometryEnvelopeError(f"AW {name} table must be a DataFrame.")
        return value

    return balls, table("verts"), table("edges"), table("surfs")


def _edge_analytic_payload(geometry: Any) -> Mapping[str, Any]:
    """Return only reconstructible parameters, never opaque callbacks."""
    if isinstance(geometry, LineEdgeGeometry):
        return _frozen_mapping({
            "kind": "line",
            "parameter_bounds": (float(geometry.t_min), float(geometry.t_max)),
            "start": tuple(float(value) for value in geometry.start),
            "end": tuple(float(value) for value in geometry.end),
        })
    common = {
        "parameter_bounds": (float(geometry.t_min), float(geometry.t_max)),
        "conic_type": str(geometry.conic_type),
        "generator_locations": tuple(
            tuple(float(value) for value in point)
            for point in geometry.generator_locations
        ),
        "generator_radii": tuple(float(value) for value in geometry.generator_radii),
        "origin": tuple(float(value) for value in geometry.origin),
        "basis_u": tuple(float(value) for value in geometry.basis_u),
        "basis_v": tuple(float(value) for value in geometry.basis_v),
        "q_coefficients": tuple(float(value) for value in geometry.q_coefficients),
    }
    if isinstance(geometry, AdditivelyWeightedTrisectorBranch):
        return _frozen_mapping({
            "kind": "aw_trisector_branch",
            **common,
            "branch": int(geometry.branch),
        })
    if isinstance(geometry, AdditivelyWeightedTrisectorConic):
        return _frozen_mapping({
            "kind": "aw_trisector_conic",
            **common,
            "component": geometry.component,
            "rho_center": geometry.rho_center,
            "rho_scale": geometry.rho_scale,
            "normal_scale": geometry.normal_scale,
            "hyperbola_orientation": geometry.hyperbola_orientation,
        })
    raise TypeError(f"Unsupported exact edge geometry: {type(geometry).__name__}")


def _vertex_status(
        row: Mapping[str, Any],
        balls: pd.DataFrame,
        tolerance: float,
) -> tuple[Status, tuple[float, float, float] | None, float | None, str | None]:
    ball_indices = _as_indices(row.get("balls"))
    location = _as_vector(row.get("loc"))
    try:
        radius = float(row.get("rad"))
    except (TypeError, ValueError):
        radius = None
    if ball_indices is None or len(ball_indices) != 4 or len(set(ball_indices)) != 4:
        return Status.UNRESOLVED, location, radius, "regular AW junction requires four distinct cells"
    if location is None or radius is None or not np.isfinite(radius):
        return Status.UNRESOLVED, location, radius, "vertex location or clearance radius is non-finite"
    if any(index < 0 or index >= len(balls) for index in ball_indices):
        return Status.UNRESOLVED, location, radius, "vertex references an unknown cell"
    clearances = []
    for index in ball_indices:
        center = _as_vector(balls.iloc[index]["loc"])
        try:
            atom_radius = float(balls.iloc[index]["rad"])
        except (TypeError, ValueError):
            return Status.UNRESOLVED, location, radius, "cell radius is invalid"
        if center is None or not np.isfinite(atom_radius):
            return Status.UNRESOLVED, location, radius, "cell geometry is invalid"
        clearances.append(float(np.linalg.norm(np.asarray(location) - center) - atom_radius))
    residual = max(abs(value - radius) for value in clearances)
    if residual > tolerance:
        return Status.UNRESOLVED, location, radius, (
            f"vertex clearance residual {residual:.3g} exceeds {tolerance:.3g}"
        )
    return Status.CERTIFIED, location, radius, None


def _component_records(
        vertices: Sequence[AWJunctionRecord],
        edges: Sequence[AWAnalyticEdgeRecord],
        surfaces: Sequence[AWAnalyticSurfaceRecord],
) -> tuple[dict[str, str], tuple[AWGeometryComponent, ...]]:
    """Find primitive components from explicit incidence, independent of row order."""
    graph: dict[str, set[str]] = {}

    def add(node: str) -> None:
        graph.setdefault(node, set())

    def join(left: str, right: str) -> None:
        add(left)
        add(right)
        graph[left].add(right)
        graph[right].add(left)

    for vertex in vertices:
        add(f"v:{vertex.vertex_id}")
        for edge_id in vertex.edge_ids:
            join(f"v:{vertex.vertex_id}", f"e:{edge_id}")
        for surface_id in vertex.surface_ids:
            join(f"v:{vertex.vertex_id}", f"s:{surface_id}")
    for edge in edges:
        add(f"e:{edge.edge_id}")
        for vertex_id in edge.vertex_ids:
            join(f"e:{edge.edge_id}", f"v:{vertex_id}")
        for surface_id in edge.surface_ids:
            join(f"e:{edge.edge_id}", f"s:{surface_id}")
    for surface in surfaces:
        add(f"s:{surface.surface_id}")
        for vertex_id in surface.vertex_ids:
            join(f"s:{surface.surface_id}", f"v:{vertex_id}")
        for edge_id in surface.edge_ids:
            join(f"s:{surface.surface_id}", f"e:{edge_id}")

    component_for_node: dict[str, str] = {}
    components = []
    pending = set(graph)
    surface_by_id = {record.surface_id: record for record in surfaces}
    edge_by_id = {record.edge_id: record for record in edges}
    while pending:
        seed = min(pending)
        stack = [seed]
        nodes = set()
        while stack:
            node = stack.pop()
            if node in nodes:
                continue
            nodes.add(node)
            stack.extend(graph[node] - nodes)
        pending.difference_update(nodes)
        vertex_ids = tuple(sorted(node[2:] for node in nodes if node.startswith("v:")))
        edge_ids = tuple(sorted(node[2:] for node in nodes if node.startswith("e:")))
        surface_ids = tuple(sorted(node[2:] for node in nodes if node.startswith("s:")))
        component_id = _record_id("component", "component", {
            "vertices": vertex_ids, "edges": edge_ids, "surfaces": surface_ids,
        })
        component_for_node.update({node: component_id for node in nodes})

        issues = []
        for surface_id in surface_ids:
            surface = surface_by_id[surface_id]
            boundary_vertices = {vertex_id: 0 for vertex_id in surface.vertex_ids}
            for edge_id in surface.edge_ids:
                edge = edge_by_id.get(edge_id)
                if edge is None or len(edge.vertex_ids) != 2:
                    issues.append(f"surface:{surface_id}:missing_or_invalid_boundary_edge")
                    continue
                if not set(edge.vertex_ids).issubset(boundary_vertices):
                    issues.append(f"surface:{surface_id}:boundary_edge_outside_surface")
                    continue
                for vertex_id in edge.vertex_ids:
                    boundary_vertices[vertex_id] += 1
            if len(surface.edge_ids) != len(surface.vertex_ids):
                issues.append(f"surface:{surface_id}:edge_vertex_count_mismatch")
            if any(degree != 2 for degree in boundary_vertices.values()):
                issues.append(f"surface:{surface_id}:open_boundary")
        components.append(AWGeometryComponent(
            component_id=component_id,
            vertex_ids=vertex_ids,
            edge_ids=edge_ids,
            surface_ids=surface_ids,
            boundary_incidence_status=(Status.CERTIFIED if not issues else Status.PARTIAL),
            boundary_issues=tuple(sorted(set(issues))),
        ))
    return component_for_node, tuple(sorted(components, key=lambda item: item.component_id))


def extract_aw_geometry_envelope(
        network: Any,
        *,
        vertex_tolerance: float = 1.0e-7,
) -> AWGeometryEnvelope:
    """Extract exact AW primitives and explicit unavailable states read-only.

    ``vertex_tolerance`` is only a validation threshold for the already solved
    vertex equations; it does not change solver tolerances or any reference
    value.  An edge resolver failure is represented as ``UNRESOLVED`` rather
    than by a sampled-polyline substitute.
    """
    if vertex_tolerance < 0.0 or not np.isfinite(vertex_tolerance):
        raise ValueError("vertex_tolerance must be finite and nonnegative")
    balls, vertex_table, edge_table, surface_table = _require_aw_network(network)
    try:
        solve_key = build_solve_universe_key(network)
    except AWGeometryReuseError as error:
        raise AWGeometryEnvelopeError(str(error)) from error

    atom_records = []
    atom_ids: list[str] = []
    cell_ids: list[str] = []
    input_stable: list[bool] = []
    for position, (_, row) in enumerate(balls.iterrows()):
        mapping = row.to_dict()
        stable_id = _stable_generator_id(mapping, position)
        center = _as_vector(mapping.get("loc"))
        try:
            radius = float(mapping.get("rad"))
        except (TypeError, ValueError):
            radius = None
        stable = _input_stable_identity(mapping)
        atom_id = _record_id("atom", "input", {"generator": stable_id})
        cell_id = _record_id("cell", solve_key.digest, {
            "generator": stable_id,
            "center": center,
            "radius": radius,
        })
        status = (Status.CERTIFIED if center is not None and radius is not None and np.isfinite(radius)
                  else Status.UNRESOLVED)
        atom_records.append(AWAtomRecord(
            atom_id=atom_id,
            cell_id=cell_id,
            source_index=position,
            stable_generator_id=stable_id,
            center=center,
            radius=radius,
            status=status,
            identity_is_input_stable=stable,
        ))
        atom_ids.append(atom_id)
        cell_ids.append(cell_id)
        input_stable.append(stable)

    vertex_records = []
    vertex_id_by_index: dict[int, str] = {}
    if vertex_table is not None:
        for position, (_, row) in enumerate(vertex_table.iterrows()):
            mapping = row.to_dict()
            raw_cells = _as_indices(mapping.get("balls"))
            referenced_cells = tuple(sorted(
                cell_ids[index] for index in raw_cells
                if 0 <= index < len(cell_ids)
            )) if raw_cells is not None else ()
            status, location, clearance_radius, reason = _vertex_status(
                mapping, balls, vertex_tolerance,
            )
            vertex_id = _record_id("vertex", solve_key.digest, {
                "cells": referenced_cells,
                "location": location,
                "clearance_radius": clearance_radius,
                "doublet": mapping.get("dub", mapping.get("vdub")),
            })
            vertex_id_by_index[position] = vertex_id
            vertex_records.append(AWJunctionRecord(
                vertex_id=vertex_id,
                source_index=position,
                cell_ids=referenced_cells,
                location=location,
                clearance_radius=clearance_radius,
                edge_ids=(),
                surface_ids=(),
                component_id=None,
                status=status,
                reason=reason,
            ))

    edge_records = []
    edge_id_by_index: dict[int, str] = {}
    if edge_table is not None:
        for position, (_, row) in enumerate(edge_table.iterrows()):
            mapping = row.to_dict()
            raw_cells = _as_indices(mapping.get("balls"))
            raw_vertices = _as_indices(mapping.get("verts"))
            record_cells = tuple(sorted(
                cell_ids[index] for index in raw_cells
                if 0 <= index < len(cell_ids)
            )) if raw_cells is not None else ()
            record_vertices = tuple(
                vertex_id_by_index[index] for index in raw_vertices
                if index in vertex_id_by_index
            ) if raw_vertices is not None else ()
            edge_id = _record_id("edge", solve_key.digest, {
                "cells": record_cells,
                "vertices": tuple(sorted(record_vertices)),
            })
            edge_id_by_index[position] = edge_id
            try:
                resolved = resolve_aw_network_edge(network, position)
                analytic = _edge_analytic_payload(resolved.geometry)
                status = Status.CERTIFIED
                reason = None
            except AWNetworkEdgeResolutionError as error:
                analytic = None
                status = Status.UNRESOLVED
                reason = f"{error.status}: {error.reason}"
            except (TypeError, ValueError, np.linalg.LinAlgError) as error:
                analytic = None
                status = Status.UNRESOLVED
                reason = f"{type(error).__name__}: {error}"
            edge_records.append(AWAnalyticEdgeRecord(
                edge_id=edge_id,
                source_index=position,
                cell_ids=record_cells,
                vertex_ids=record_vertices,
                surface_ids=(),
                analytic=analytic,
                component_id=None,
                status=status,
                reason=reason,
            ))

    surface_records = []
    if surface_table is not None:
        for position, (_, row) in enumerate(surface_table.iterrows()):
            mapping = row.to_dict()
            raw_cells = _as_indices(mapping.get("balls"))
            raw_vertices = _as_indices(mapping.get("verts"))
            raw_edges = _as_indices(mapping.get("edges"))
            record_cells = tuple(sorted(
                cell_ids[index] for index in raw_cells
                if 0 <= index < len(cell_ids)
            )) if raw_cells is not None else ()
            record_vertices = tuple(
                vertex_id_by_index[index] for index in raw_vertices
                if index in vertex_id_by_index
            ) if raw_vertices is not None else ()
            record_edges = tuple(
                edge_id_by_index[index] for index in raw_edges
                if index in edge_id_by_index
            ) if raw_edges is not None else ()
            surface_id = _record_id("surface", solve_key.digest, {
                "cells": record_cells,
                "vertices": tuple(sorted(record_vertices)),
                "edges": tuple(sorted(record_edges)),
            })
            coefficients = None
            geometry_kind = None
            coefficient_source = None
            status = Status.CERTIFIED
            reason = None
            if raw_cells is None or len(raw_cells) != 2 or len(set(raw_cells)) != 2:
                status = Status.UNRESOLVED
                reason = "regular AW surface requires two distinct cells"
            elif any(index < 0 or index >= len(balls) for index in raw_cells):
                status = Status.UNRESOLVED
                reason = "surface references an unknown cell"
            else:
                raw_coefficients = mapping.get("func")
                try:
                    values = np.asarray(raw_coefficients, dtype=float)
                    if values.shape != (14,) or not np.all(np.isfinite(values)):
                        raise ValueError("surface function is unavailable")
                    coefficient_source = "network_surface_func"
                except (TypeError, ValueError):
                    try:
                        values = np.asarray(calc_surf_func(
                            balls.iloc[raw_cells[0]]["loc"], balls.iloc[raw_cells[0]]["rad"],
                            balls.iloc[raw_cells[1]]["loc"], balls.iloc[raw_cells[1]]["rad"],
                        ), dtype=float)
                        if values.shape != (14,) or not np.all(np.isfinite(values)):
                            raise ValueError("generator reconstruction was non-finite")
                        coefficient_source = "exact_generator_reconstruction"
                    except (TypeError, ValueError, np.linalg.LinAlgError) as error:
                        values = None
                        status = Status.UNRESOLVED
                        reason = f"analytic surface unavailable: {error}"
                if values is not None:
                    try:
                        QuadraticSurfaceGeometry(values)
                        coefficients = tuple(float(value) for value in values)
                        geometry_kind = "plane" if bool(mapping.get("flat", False)) else "quadratic"
                    except (TypeError, ValueError) as error:
                        status = Status.UNRESOLVED
                        reason = f"invalid quadratic surface: {error}"
            owners = ()
            if len(record_cells) == 2:
                owners = tuple(
                    AWSurfaceOwner(
                        cell_id=owner,
                        neighbor_cell_id=neighbor,
                        normal_convention="aw_clearance_gradient(owner_cell, neighbor_cell)",
                    )
                    for owner, neighbor in (
                        (record_cells[0], record_cells[1]),
                        (record_cells[1], record_cells[0]),
                    )
                )
            surface_records.append(AWAnalyticSurfaceRecord(
                surface_id=surface_id,
                source_index=position,
                cell_ids=record_cells,
                vertex_ids=record_vertices,
                edge_ids=record_edges,
                coefficients=coefficients,
                geometry_kind=geometry_kind,
                coefficient_source=coefficient_source,
                owners=owners,
                component_id=None,
                status=status,
                reason=reason,
            ))

    # Cross-reference the producer's authoritative row-index incidence after
    # all deterministic record IDs have been assigned.
    surface_ids_by_edge: dict[str, list[str]] = {record.edge_id: [] for record in edge_records}
    edge_ids_by_vertex: dict[str, list[str]] = {record.vertex_id: [] for record in vertex_records}
    surface_ids_by_vertex: dict[str, list[str]] = {record.vertex_id: [] for record in vertex_records}
    for surface in surface_records:
        for edge_id in surface.edge_ids:
            surface_ids_by_edge.setdefault(edge_id, []).append(surface.surface_id)
        for vertex_id in surface.vertex_ids:
            surface_ids_by_vertex.setdefault(vertex_id, []).append(surface.surface_id)
    for edge in edge_records:
        for vertex_id in edge.vertex_ids:
            edge_ids_by_vertex.setdefault(vertex_id, []).append(edge.edge_id)
    edge_records = [replace(
        record,
        surface_ids=tuple(sorted(surface_ids_by_edge.get(record.edge_id, ()))),
    ) for record in edge_records]
    vertex_records = [replace(
        record,
        edge_ids=tuple(sorted(edge_ids_by_vertex.get(record.vertex_id, ()))),
        surface_ids=tuple(sorted(surface_ids_by_vertex.get(record.vertex_id, ()))),
    ) for record in vertex_records]

    component_by_node, components = _component_records(vertex_records, edge_records, surface_records)
    vertex_records = [replace(
        record, component_id=component_by_node.get(f"v:{record.vertex_id}"),
    ) for record in vertex_records]
    edge_records = [replace(
        record, component_id=component_by_node.get(f"e:{record.edge_id}"),
    ) for record in edge_records]
    surface_records = [replace(
        record, component_id=component_by_node.get(f"s:{record.surface_id}"),
    ) for record in surface_records]

    cell_surfaces = {cell_id: [] for cell_id in cell_ids}
    cell_edges = {cell_id: [] for cell_id in cell_ids}
    cell_vertices = {cell_id: [] for cell_id in cell_ids}
    for record in surface_records:
        for cell_id in record.cell_ids:
            cell_surfaces.setdefault(cell_id, []).append(record.surface_id)
    for record in edge_records:
        for cell_id in record.cell_ids:
            cell_edges.setdefault(cell_id, []).append(record.edge_id)
    for record in vertex_records:
        for cell_id in record.cell_ids:
            cell_vertices.setdefault(cell_id, []).append(record.vertex_id)

    complete_column_known = "complete" in balls.columns
    complete_values = (
        tuple(bool(value) for value in balls["complete"])
        if complete_column_known else ()
    )
    cell_records = []
    for atom, cell_id in zip(atom_records, cell_ids):
        complete = complete_values[atom.source_index] if complete_column_known else None
        status = atom.status
        if complete is False and status is Status.CERTIFIED:
            status = Status.PARTIAL
        cell_records.append(AWCellRecord(
            cell_id=cell_id,
            atom_id=atom.atom_id,
            source_index=atom.source_index,
            complete=complete,
            surface_ids=tuple(sorted(cell_surfaces[cell_id])),
            edge_ids=tuple(sorted(cell_edges[cell_id])),
            vertex_ids=tuple(sorted(cell_vertices[cell_id])),
            status=status,
        ))

    all_primitives = tuple(surface_records) + tuple(edge_records) + tuple(vertex_records)
    reasons = []
    if not complete_column_known:
        reasons.append("cell_completeness_not_recorded")
    elif not all(complete_values):
        reasons.append("incomplete_cells")
    if not all(input_stable):
        reasons.append("input_stable_generator_ids_not_recorded")
    if vertex_table is None or edge_table is None or surface_table is None:
        reasons.append("topology_tables_not_available")
    if not all_primitives:
        reasons.append("no_geometric_primitives")
    unresolved = [record for record in all_primitives if record.status is not Status.CERTIFIED]
    if unresolved:
        reasons.append("unresolved_analytic_primitives")
    provenance_payload = getattr(network, "aw_solver_provenance", None)
    solver_complete = bool(
        isinstance(provenance_payload, Mapping)
        and isinstance(provenance_payload.get("completeness"), Mapping)
        and provenance_payload["completeness"].get("completeness_proven", False)
    )
    if not solver_complete:
        reasons.append("solver_completeness_not_proven")
    if not all_primitives:
        envelope_status = Status.NOT_CALCULATED
    elif not unresolved and complete_column_known and all(complete_values) and solver_complete and all(input_stable):
        envelope_status = Status.CERTIFIED
    elif any(record.status is Status.CERTIFIED for record in all_primitives):
        envelope_status = Status.PARTIAL
    else:
        envelope_status = Status.UNRESOLVED

    counts = _frozen_mapping({
        "atoms": len(atom_records),
        "cells": len(cell_records),
        "surfaces": len(surface_records),
        "edges": len(edge_records),
        "vertices": len(vertex_records),
        "components": len(components),
        "certified_primitives": sum(record.status is Status.CERTIFIED for record in all_primitives),
        "unresolved_primitives": len(unresolved),
    })
    completeness = AWGeometryCompleteness(
        status=envelope_status,
        cell_completeness_known=complete_column_known,
        globally_complete=(all(complete_values) if complete_column_known else None),
        solver_completeness_proven=solver_complete,
        counts=counts,
        reasons=tuple(sorted(set(reasons))),
    )
    provenance = _frozen_mapping({
        "producer": "vorpy.aw_geometry_envelope",
        "schema_version": SCHEMA_VERSION,
        "coordinate_units": LENGTH_UNITS,
        "solve_key_digest": solve_key.digest,
        "solver_provenance_present": isinstance(provenance_payload, Mapping),
        "analytic_surface_policy": "exact quadratic coefficients; no mesh fitting",
        "analytic_edge_policy": "exact resolver only; unresolved edges have no substitute",
        "spherical_patch_policy": "not supported for pairwise AW boundaries",
    })
    return AWGeometryEnvelope(
        schema_version=SCHEMA_VERSION,
        coordinate_units=LENGTH_UNITS,
        solve_key_digest=solve_key.digest,
        atoms=tuple(atom_records),
        cells=tuple(cell_records),
        surfaces=tuple(surface_records),
        edges=tuple(edge_records),
        vertices=tuple(vertex_records),
        components=components,
        completeness=completeness,
        provenance=provenance,
    )


__all__ = [
    "AWAnalyticEdgeRecord",
    "AWAnalyticSurfaceRecord",
    "AWAtomRecord",
    "AWCellRecord",
    "AWGeometryCompleteness",
    "AWGeometryComponent",
    "AWGeometryEnvelope",
    "AWGeometryEnvelopeError",
    "AWJunctionRecord",
    "AWSurfaceOwner",
    "LENGTH_UNITS",
    "SCHEMA_VERSION",
    "extract_aw_geometry_envelope",
]
