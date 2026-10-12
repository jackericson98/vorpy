"""Immutable selections over an already-solved VorPy network.

This module deliberately contains no network construction path.  A
``NetworkView`` snapshots selected rows and incidence from one source Network;
it never calls vertex search, ``connect``, edge construction, or surface
construction.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np

from vorpy.src.network.aw_geometry_reuse import (
    AWGeometryReuseError,
    SolveContextKey,
    _freeze,
    build_solve_context_key,
)


VIEW_SCHEMA_VERSION = 1


class NetworkViewError(ValueError):
    """Raised when a selection cannot be proven safe as a read-only view."""


@dataclass(frozen=True)
class ViewCompleteness:
    """Closure and source-completeness evidence for one view."""

    proven: bool
    scope: str
    reasons: tuple[str, ...] = ()
    counts: tuple[tuple[str, int], ...] = ()

    @property
    def status(self) -> str:
        return "COMPLETE" if self.proven else "INCOMPLETE"


@dataclass(frozen=True)
class ViewFeature:
    """Immutable source-row snapshot with preserved topology identity."""

    kind: str
    source_id: Any
    topology_index: int
    generator_local_indices: tuple[int, ...]
    generator_ids: tuple[str, ...]
    generator_system_ids: tuple[int | None, ...]
    incidence: Mapping[str, tuple[Any, ...]]
    data: Mapping[str, Any]


@dataclass(frozen=True)
class NetworkView:
    """Read-only group or ordered A/B interface selection.

    The source Network remains the owner of the solve and its DataFrames.  The
    view stores immutable row snapshots so later source-table mutation cannot
    change the selection identity or inspection results.
    """

    solve_context: SolveContextKey
    source_network_id: str
    selection_kind: str
    selection_id: str
    system_ids: tuple[int, ...]
    group_a_ids: tuple[int, ...]
    group_b_ids: tuple[int, ...]
    cell_features: tuple[ViewFeature, ...]
    surface_features: tuple[ViewFeature, ...]
    edge_features: tuple[ViewFeature, ...]
    vertex_features: tuple[ViewFeature, ...]
    support_generator_local_indices: tuple[int, ...]
    boundary_generator_ids: tuple[str, ...]
    completeness: ViewCompleteness

    @classmethod
    def for_group(
        cls,
        source_network: Any,
        system_generator_ids: Sequence[int],
        *,
        source_network_id: str | None = None,
    ) -> "NetworkView":
        """Create a view for explicit system-global physical generator IDs."""
        group = _normalize_system_ids(system_generator_ids, "group")
        return cls._build(
            source_network, "group", group, (), (), source_network_id
        )

    @classmethod
    def for_interface(
        cls,
        source_network: Any,
        group_a_system_ids: Sequence[int],
        group_b_system_ids: Sequence[int],
        *,
        source_network_id: str | None = None,
    ) -> "NetworkView":
        """Create an ordered A -> B interface view."""
        group_a = _normalize_system_ids(group_a_system_ids, "interface group A")
        group_b = _normalize_system_ids(group_b_system_ids, "interface group B")
        if not group_a or not group_b:
            raise NetworkViewError("Interface groups must be non-empty.")
        if set(group_a).intersection(group_b):
            raise NetworkViewError("Interface groups must be disjoint.")
        return cls._build(
            source_network,
            "interface",
            (),
            group_a,
            group_b,
            source_network_id,
        )

    @classmethod
    def from_group(
        cls,
        source_network: Any,
        system_generator_ids: Sequence[int],
        *,
        source_network_id: str | None = None,
    ):
        return cls.for_group(
            source_network,
            system_generator_ids,
            source_network_id=source_network_id,
        )

    @classmethod
    def from_interface(
        cls,
        source_network: Any,
        group_a_system_ids: Sequence[int],
        group_b_system_ids: Sequence[int],
        *,
        source_network_id: str | None = None,
    ):
        return cls.for_interface(
            source_network,
            group_a_system_ids,
            group_b_system_ids,
            source_network_id=source_network_id,
        )

    @classmethod
    def _build(
        cls,
        source_network: Any,
        selection_kind: str,
        system_ids: tuple[int, ...],
        group_a_ids: tuple[int, ...],
        group_b_ids: tuple[int, ...],
        source_network_id: str | None = None,
    ) -> "NetworkView":
        if selection_kind not in {"group", "interface"}:
            raise NetworkViewError(f"Unsupported selection kind: {selection_kind!r}")

        try:
            context = build_solve_context_key(source_network)
        except (AWGeometryReuseError, TypeError, ValueError) as error:
            raise NetworkViewError(f"Unsupported source solve context: {error}") from error

        balls = getattr(source_network, "balls", None)
        verts = getattr(source_network, "verts", None)
        edges = getattr(source_network, "edges", None)
        surfs = getattr(source_network, "surfs", None)
        if any(table is None for table in (balls, verts, edges, surfs)):
            raise NetworkViewError("A view requires solved balls, vertices, edges, and surfaces.")
        if "complete" not in balls.columns:
            raise NetworkViewError("Source cell completeness is unavailable.")

        local_to_system = tuple(context.local_to_system)
        if len(local_to_system) != len(balls):
            raise NetworkViewError("Solve-context mapping does not match source balls.")
        system_to_local = {}
        for local_index, system_index in local_to_system:
            if system_index is None:
                continue
            if system_index in system_to_local:
                raise NetworkViewError("Duplicate system-global generator ID.")
            system_to_local[system_index] = local_index

        boundary_local = _boundary_local_indices(source_network, balls)
        boundary_system_ids = {
            local_to_system[index][1]
            for index in boundary_local
            if 0 <= index < len(local_to_system)
            and local_to_system[index][1] is not None
        }
        physical_system_ids = set(system_to_local).difference(boundary_system_ids)
        requested_ids = (
            set(system_ids)
            if selection_kind == "group"
            else set(group_a_ids).union(group_b_ids)
        )
        missing = requested_ids.difference(physical_system_ids)
        if missing:
            raise NetworkViewError(
                "Selection references unknown or non-physical system IDs: "
                + ", ".join(str(value) for value in sorted(missing))
            )

        selected_local = {
            system_to_local[system_id] for system_id in requested_ids
        }
        group_a_local = {system_to_local[value] for value in group_a_ids}
        group_b_local = {system_to_local[value] for value in group_b_ids}

        surface_candidates = []
        for position in range(len(surfs)):
            row = surfs.iloc[position]
            balls_ref = _row_indices(row, "balls", "surface", position)
            _validate_indices(balls_ref, len(balls), "surface balls", position)
            if selection_kind == "group":
                selected = bool(set(balls_ref).intersection(selected_local))
            else:
                selected = bool(
                    set(balls_ref).intersection(group_a_local)
                    and set(balls_ref).intersection(group_b_local)
                )
            if selected:
                surface_candidates.append(position)

        if not surface_candidates:
            raise NetworkViewError("Selection has no supported source surfaces.")

        selected_edges = set()
        selected_vertices = set()
        for surface_position in surface_candidates:
            surface = surfs.iloc[surface_position]
            edge_refs = _row_indices(surface, "edges", "surface", surface_position)
            vertex_refs = _row_indices(surface, "verts", "surface", surface_position)
            _validate_indices(edge_refs, len(edges), "surface edges", surface_position)
            _validate_indices(vertex_refs, len(verts), "surface vertices", surface_position)
            surface_vertices = set(vertex_refs)
            for edge_position in edge_refs:
                edge_vertices = _row_indices(
                    edges.iloc[edge_position], "verts", "edge", edge_position
                )
                if not set(edge_vertices).issubset(surface_vertices):
                    raise NetworkViewError(
                        f"Surface {surface_position} has an open edge/vertex boundary."
                    )
            selected_edges.update(edge_refs)
            selected_vertices.update(vertex_refs)

        # Downward closure: retain every vertex referenced by a retained edge,
        # even when that vertex is a support feature rather than a side cell.
        for edge_position in tuple(sorted(selected_edges)):
            edge = edges.iloc[edge_position]
            vertex_refs = _row_indices(edge, "verts", "edge", edge_position)
            ball_refs = _row_indices(edge, "balls", "edge", edge_position)
            _validate_indices(vertex_refs, len(verts), "edge vertices", edge_position)
            _validate_indices(ball_refs, len(balls), "edge balls", edge_position)
            selected_vertices.update(vertex_refs)

        relevant_generators = set(selected_local)
        for vertex_position in sorted(selected_vertices):
            vertex = verts.iloc[vertex_position]
            ball_refs = _row_indices(vertex, "balls", "vertex", vertex_position)
            _validate_indices(ball_refs, len(balls), "vertex balls", vertex_position)
            relevant_generators.update(ball_refs)
        for edge_position in sorted(selected_edges):
            relevant_generators.update(
                _row_indices(edges.iloc[edge_position], "balls", "edge", edge_position)
            )
        for surface_position in surface_candidates:
            relevant_generators.update(
                _row_indices(surfs.iloc[surface_position], "balls", "surface", surface_position)
            )

        boundary_completeness_exemptions = (
            boundary_local
            if _has_authoritative_solvent_provenance(source_network)
            else frozenset()
        )
        _validate_complete_cells(
            balls,
            relevant_generators,
            boundary_completeness_exemptions,
        )
        _validate_complete_features(
            verts, selected_vertices, "vertex"
        )
        _validate_complete_features(edges, selected_edges, "edge")
        _validate_complete_features(surfs, surface_candidates, "surface")

        generator_ids, generator_system_ids = _generator_maps(
            source_network, context, balls, boundary_local
        )
        feature_args = {
            "balls": balls,
            "verts": verts,
            "edges": edges,
            "surfs": surfs,
            "generator_ids": generator_ids,
            "generator_system_ids": generator_system_ids,
        }
        cell_positions = tuple(sorted(selected_local))
        surface_positions = tuple(sorted(surface_candidates))
        edge_positions = tuple(sorted(selected_edges))
        vertex_positions = tuple(sorted(selected_vertices))
        cell_features = _snapshot_features("cell", cell_positions, feature_args)
        surface_features = _snapshot_features("surface", surface_positions, feature_args)
        edge_features = _snapshot_features("edge", edge_positions, feature_args)
        vertex_features = _snapshot_features("vertex", vertex_positions, feature_args)

        source_network_id = _resolve_source_network_id(
            source_network,
            context,
            source_network_id,
        )
        selection_payload = {
            "schema_version": VIEW_SCHEMA_VERSION,
            "source_network_id": source_network_id,
            "source_context_id": context.digest,
            "selection_kind": selection_kind,
            "group_ids": system_ids,
            "group_a_ids": group_a_ids,
            "group_b_ids": group_b_ids,
        }
        selection_id = "view:" + _digest(selection_payload)
        boundary_ids = tuple(
            generator_ids[index]
            for index in sorted(relevant_generators.intersection(boundary_local))
        )
        completeness = ViewCompleteness(
            proven=True,
            scope=selection_kind,
            counts=(
                ("cells", len(cell_features)),
                ("surfaces", len(surface_features)),
                ("edges", len(edge_features)),
                ("vertices", len(vertex_features)),
                ("support_generators", len(relevant_generators)),
            ),
        )
        return cls(
            solve_context=context,
            source_network_id=source_network_id,
            selection_kind=selection_kind,
            selection_id=selection_id,
            system_ids=system_ids if selection_kind == "group" else tuple(sorted(requested_ids)),
            group_a_ids=group_a_ids,
            group_b_ids=group_b_ids,
            cell_features=cell_features,
            surface_features=surface_features,
            edge_features=edge_features,
            vertex_features=vertex_features,
            support_generator_local_indices=tuple(sorted(relevant_generators)),
            boundary_generator_ids=boundary_ids,
            completeness=completeness,
        )

    @property
    def solve_universe_id(self) -> str:
        return self.solve_context.solve_universe_id

    @property
    def context_id(self) -> str:
        return self.solve_context.context_id

    @property
    def network_id(self) -> str:
        """Stable source-network namespace used by local topology indices."""
        return self.source_network_id

    @property
    def source_network_id_is_provisional(self) -> bool:
        """Whether the source namespace is the non-certifiable F7 fallback."""
        return self.source_network_id == f"network:{self.solve_context.digest}"

    @property
    def source_network_id_is_authoritative(self) -> bool:
        """Whether the source namespace is producer-attested rather than provisional."""
        return not self.source_network_id_is_provisional

    @property
    def scheme(self) -> str:
        return self.solve_context.scheme

    @property
    def frame_id(self) -> str:
        return self.solve_context.frame_id

    @property
    def local_to_system(self) -> tuple[tuple[int, int | None], ...]:
        return self.solve_context.local_to_system

    @property
    def cells(self) -> tuple[ViewFeature, ...]:
        return self.cell_features

    @property
    def surfaces(self) -> tuple[ViewFeature, ...]:
        return self.surface_features

    @property
    def edges(self) -> tuple[ViewFeature, ...]:
        return self.edge_features

    @property
    def vertices(self) -> tuple[ViewFeature, ...]:
        return self.vertex_features

    @property
    def generator_identity(self) -> tuple[tuple[int, str, int | None], ...]:
        return tuple(
            (
                feature.topology_index,
                feature.generator_ids[0],
                feature.generator_system_ids[0],
            )
            for feature in self.cell_features
        )

    def features(self, kind: str) -> tuple[ViewFeature, ...]:
        tables = {
            "cell": self.cells,
            "surface": self.surfaces,
            "edge": self.edges,
            "vertex": self.vertices,
        }
        try:
            return tables[kind]
        except KeyError as error:
            raise NetworkViewError(f"Unsupported view feature kind: {kind!r}") from error

    def source_ids(self, kind: str) -> tuple[Any, ...]:
        return tuple(feature.source_id for feature in self.features(kind))


def _normalize_system_ids(values: Sequence[int], label: str) -> tuple[int, ...]:
    if isinstance(values, Integral) and not isinstance(values, bool):
        values = (int(values),)
    try:
        values = tuple(values)
    except TypeError as error:
        raise NetworkViewError(f"{label} must contain integer system IDs.") from error
    try:
        normalized = tuple(sorted({int(value) for value in values}))
    except (TypeError, ValueError, OverflowError) as error:
        raise NetworkViewError(f"{label} must contain integer system IDs.") from error
    if any(isinstance(value, bool) for value in values):
        raise NetworkViewError(f"{label} must contain integer system IDs.")
    if any(value < 0 for value in normalized):
        raise NetworkViewError(f"{label} cannot contain negative system IDs.")
    if not normalized:
        raise NetworkViewError(f"{label} must be non-empty.")
    return normalized


def _row_indices(row: Any, field_name: str, kind: str, position: int) -> tuple[int, ...]:
    value = row.get(field_name, ())
    if value is None:
        return ()
    if isinstance(value, (float, np.floating)) and np.isnan(value):
        return ()
    try:
        return tuple(int(item) for item in value)
    except (TypeError, ValueError, OverflowError) as error:
        raise NetworkViewError(
            f"Invalid {kind} {field_name} incidence at source row {position}."
        ) from error


def _validate_indices(values: Sequence[int], size: int, label: str, position: int) -> None:
    if any(value < 0 or value >= size for value in values):
        raise NetworkViewError(f"Invalid {label} incidence at source row {position}.")


def _boundary_local_indices(source_network: Any, balls: Any) -> set[int]:
    indices = set(int(value) for value in (getattr(source_network, "boundary_indices", ()) or ()))
    if "is_boundary_generator" in balls.columns:
        for position in range(len(balls)):
            value = balls.iloc[position].get("is_boundary_generator")
            if value is not None and not (
                isinstance(value, (float, np.floating)) and np.isnan(value)
            ) and bool(value):
                indices.add(position)
    return indices


def _validate_complete_cells(
    balls: Any,
    positions: set[int],
    boundary_positions: set[int] = frozenset(),
) -> None:
    # Numerical boundary generators support finite surfaces but do not own
    # physical cells; only physical generator rows require cell completeness.
    for position in positions:
        if position in boundary_positions:
            continue
        value = balls.iloc[position].get("complete")
        if value is None or (
            isinstance(value, (float, np.floating)) and np.isnan(value)
        ) or not bool(value):
            raise NetworkViewError(
                f"Source cell completeness is not proven for generator {position}."
            )


def _has_authoritative_solvent_provenance(source_network: Any) -> bool:
    """Recognize only explicit F10 producer evidence for boundary support."""
    mode = getattr(source_network, "solvent_mode", None)
    return mode in {"MOLECULE_ONLY", "SOLVENT_CONTEXT", "FULL_SOLVENT"}


def _validate_complete_features(frame: Any, positions: set[int], kind: str) -> None:
    if "complete" not in frame.columns:
        return
    for position in positions:
        value = frame.iloc[position].get("complete")
        if value is None or (
            isinstance(value, (float, np.floating)) and np.isnan(value)
        ) or not bool(value):
            raise NetworkViewError(
                f"Source {kind} completeness is not proven at row {position}."
            )


def _generator_maps(source_network: Any, context: SolveContextKey, balls: Any, boundary_local: set[int]):
    physical_ids = iter(context.physical_generator_ids)
    boundary_ids = iter(context.boundary_generator_ids)
    generator_ids = {}
    generator_system_ids = {}
    for local_index, (_, system_index) in enumerate(context.local_to_system):
        if local_index in boundary_local:
            try:
                stable_id = next(boundary_ids)
            except StopIteration as error:
                raise NetworkViewError("Boundary generator identity is incomplete.") from error
        else:
            try:
                stable_id = next(physical_ids)
            except StopIteration as error:
                raise NetworkViewError("Physical generator identity is incomplete.") from error
        generator_ids[local_index] = stable_id
        generator_system_ids[local_index] = system_index
    try:
        next(physical_ids)
        raise NetworkViewError("Physical generator identity has extra IDs.")
    except StopIteration:
        pass
    try:
        next(boundary_ids)
        raise NetworkViewError("Boundary generator identity has extra IDs.")
    except StopIteration:
        pass
    return generator_ids, generator_system_ids


def _source_id(frame: Any, position: int) -> Any:
    value = frame.index[position]
    return value.item() if hasattr(value, "item") else value


def _resolve_source_network_id(
    source_network: Any,
    context: SolveContextKey,
    explicit_source_network_id: str | None = None,
) -> str:
    if explicit_source_network_id is not None:
        resolved = _validate_source_network_id(explicit_source_network_id)
        producer_id = getattr(source_network, "source_network_id", None)
        if (
            producer_id is not None
            and _validate_source_network_id(str(producer_id)) != resolved
        ):
            raise NetworkViewError(
                "Explicit source_network_id conflicts with the source Network."
            )
        return resolved
    value = getattr(source_network, "source_network_id", None)
    if value is not None:
        return _validate_source_network_id(str(value))
    # archive_id is an artifact identity, never the source-network namespace.
    # Keep the F7 bookkeeping ID for inspection, but expose it as provisional.
    return f"network:{context.digest}"


def _validate_source_network_id(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise NetworkViewError("source_network_id must be a non-empty string.")
    return value


def _snapshot_features(kind: str, positions: Sequence[int], tables: Mapping[str, Any]):
    frame_name = {
        "cell": "balls",
        "vertex": "verts",
        "edge": "edges",
        "surface": "surfs",
    }[kind]
    frame = tables[frame_name]
    features = []
    for position in positions:
        row = frame.iloc[position]
        if kind == "cell":
            generator_local_indices = (position,)
        else:
            generator_local_indices = _row_indices(row, "balls", kind, position)
        generator_ids = tuple(tables["generator_ids"][index] for index in generator_local_indices)
        generator_system_ids = tuple(
            tables["generator_system_ids"][index]
            for index in generator_local_indices
        )
        incidence = {}
        for field_name, target_name in (
            ("verts", "verts"),
            ("edges", "edges"),
            ("surfs", "surfs"),
        ):
            if field_name in row.index:
                refs = _row_indices(row, field_name, kind, position)
                target = tables[target_name]
                incidence[field_name] = tuple(_source_id(target, ref) for ref in refs)
        incidence["balls"] = generator_ids
        features.append(ViewFeature(
            kind=kind,
            source_id=_source_id(frame, position),
            topology_index=int(position),
            generator_local_indices=tuple(generator_local_indices),
            generator_ids=generator_ids,
            generator_system_ids=generator_system_ids,
            incidence=MappingProxyType({
                key: tuple(value) for key, value in incidence.items()
            }),
            data=_freeze(row.to_dict()),
        ))
    return tuple(features)


def _digest(value: Any) -> str:
    encoded = json.dumps(
        _canonical(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def _canonical(value: Any) -> Any:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, Mapping):
        return {
            str(key): _canonical(item)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_canonical(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_canonical(item) for item in value), key=lambda item: json.dumps(item, sort_keys=True))
    if isinstance(value, float):
        if not np.isfinite(value):
            raise NetworkViewError("Non-finite selection identity value.")
        return {"__float_hex__": value.hex()}
    if value is None or isinstance(value, (bool, int, str)):
        return value
    raise NetworkViewError(
        f"Unsupported selection identity value: {type(value).__name__}"
    )


__all__ = [
    "NetworkView",
    "NetworkViewError",
    "ViewCompleteness",
    "ViewFeature",
]
