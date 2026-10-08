"""Read-only AW geometry extraction for the phase-one reuse prototype.

This module deliberately stops at extraction and materialisation.  It does not
alter :class:`Network`, call the vertex solver, or participate in CLI planning.
The returned snapshot is a content-addressed description of one completed AW
network; an interface view only selects records from that snapshot.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
import json
from math import isclose, isfinite
import os
from time import perf_counter
import tracemalloc
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd


SCHEMA_VERSION = 1
DEFAULT_SOLVER_VERSION = "aw-geometry-reuse-v1"
_PRESENTATION_SETTINGS = frozenset(
    {"verbose", "print_actions", "print_vert_metrics", "print_metrics"}
)


class AWGeometryReuseError(ValueError):
    """Raised when a snapshot cannot be proven safe for extraction."""


def _canonical(value: Any) -> Any:
    """Return deterministic JSON-compatible data without using object ids."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, np.ndarray):
        return [_canonical(item) for item in value.tolist()]
    if isinstance(value, Mapping):
        return {
            str(key): _canonical(item)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_canonical(item) for item in value]
    if isinstance(value, (set, frozenset)):
        values = [_canonical(item) for item in value]
        return sorted(values, key=lambda item: json.dumps(item, sort_keys=True))
    if isinstance(value, float):
        if not isfinite(value):
            raise AWGeometryReuseError("Non-finite values cannot identify AW geometry.")
        return {"__float_hex__": value.hex()}
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if hasattr(value, "value") and isinstance(value.value, (str, int, float)):
        return {"__enum__": _canonical(value.value)}
    if hasattr(value, "__dict__"):
        return {"__object__": type(value).__name__, "value": _canonical(value.__dict__)}
    return {"__text__": str(value)}


def _digest(value: Any) -> str:
    encoded = json.dumps(
        _canonical(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def _freeze(value: Any) -> Any:
    """Freeze nested row values so the source DataFrames are never shared."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, np.ndarray):
        return tuple(_freeze(item) for item in value.tolist())
    if isinstance(value, Mapping):
        return MappingProxyType(
            {key: _freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(item) for item in value)
    return value


def _thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_thaw(item) for item in value]
    if isinstance(value, frozenset):
        return sorted(_thaw(item) for item in value)
    return value


def _stable_generator_id(row: Mapping[str, Any], position: int) -> str:
    for field in ("stable_id", "atom_stable_id", "generator_id"):
        value = row.get(field)
        if value is not None and str(value):
            return str(value)
    # Network adjacency is system-global in the current implementation.  This
    # fallback is intentionally explicit; production archives should provide a
    # stable atom/generator id rather than relying on local row position.
    return f"ball:{position}"


def _scope_mode(network: Any) -> str:
    balls = getattr(network, "balls", None)
    count = 0 if balls is None else len(balls)
    if getattr(network, "iface_grps", None) is not None:
        return "interface"
    group = getattr(network, "group", None)
    if group is None or set(group) == set(range(count)):
        return "full"
    return "group"


def _settings_for_identity(settings: Mapping[str, Any] | None) -> dict[str, Any]:
    settings = settings or {}
    return {
        str(key): value
        for key, value in settings.items()
        if key not in _PRESENTATION_SETTINGS
    }


@dataclass(frozen=True)
class SolveUniverseKey:
    """Deterministic identity of one AW geometry construction universe."""

    schema_version: int
    scheme: str
    scope_mode: str
    frame_id: str
    partition_scheme: str
    environment_digest: str
    generator_digest: str
    geometry_digest: str
    resolution: Any
    solver_settings_digest: str
    solver_version: str
    digest: str

    @property
    def competitor_digest(self) -> str:
        """Identity of the full coordinate/radius competitor universe."""
        return _digest(
            {
                "generator_digest": self.generator_digest,
                "geometry_digest": self.geometry_digest,
                "environment_digest": self.environment_digest,
                "partition_scheme": self.partition_scheme,
            }
        )


@dataclass(frozen=True)
class ReuseMetrics:
    operation: str
    seconds: float
    peak_tracemalloc_bytes: int | None = None


@dataclass(frozen=True)
class AWGenerator:
    stable_id: str
    source_index: int
    loc: tuple[float, ...]
    rad: float
    data: Mapping[str, Any]


@dataclass(frozen=True)
class AWVertexPrimitive:
    primitive_id: str
    source_index: int
    balls: tuple[int, ...]
    data: Mapping[str, Any]


@dataclass(frozen=True)
class AWEdgePrimitive:
    primitive_id: str
    source_index: int
    balls: tuple[int, ...]
    vertices: tuple[int, ...]
    data: Mapping[str, Any]


@dataclass(frozen=True)
class AWSurfacePrimitive:
    primitive_id: str
    source_index: int
    balls: tuple[int, ...]
    vertices: tuple[int, ...]
    edges: tuple[int, ...]
    data: Mapping[str, Any]


@dataclass(frozen=True)
class AWGeometryStore:
    """Immutable AW records extracted from one completed Network."""

    key: SolveUniverseKey
    generators: tuple[AWGenerator, ...]
    vertices: tuple[AWVertexPrimitive, ...]
    edges: tuple[AWEdgePrimitive, ...]
    surfaces: tuple[AWSurfacePrimitive, ...]
    table_columns: Mapping[str, tuple[str, ...]]
    global_complete: bool
    metrics: ReuseMetrics


@dataclass(frozen=True)
class AWMaterializedInterface:
    """Detached Network-shaped tables produced without a geometry solve."""

    key: SolveUniverseKey
    group_a: tuple[int, ...]
    group_b: tuple[int, ...]
    balls: pd.DataFrame
    verts: pd.DataFrame
    edges: pd.DataFrame
    surfs: pd.DataFrame
    metrics: ReuseMetrics


@dataclass(frozen=True)
class AWInterfaceView:
    source: AWGeometryStore
    group_a: tuple[int, ...]
    group_b: tuple[int, ...]
    vertex_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    surface_ids: tuple[str, ...]
    discarded: Mapping[str, int]

    def materialize(self, *, measure_memory: bool = False) -> AWMaterializedInterface:
        return materialize_aw_interface(self, measure_memory=measure_memory)


@dataclass(frozen=True)
class AWInterfaceComparison:
    equivalent: bool
    same_competitor_universe: bool
    compatible_geometry_context: bool
    source_scope: str
    direct_scope: str
    counts: Mapping[str, Mapping[str, int]]
    mismatches: Mapping[str, tuple[str, ...]]


class AWReuseDecision(str, Enum):
    EXACT_REUSE = "EXACT_REUSE"
    REJECT_INCOMPLETE = "REJECT_INCOMPLETE"
    REJECT_ENVIRONMENT = "REJECT_ENVIRONMENT"
    REJECT_TOPOLOGY = "REJECT_TOPOLOGY"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True)
class AWLocalCompletenessCertificate:
    """Conservative certificate for one selected interface topology."""

    proven: bool
    globally_complete: bool
    relevant_generators: tuple[int, ...]
    relevant_vertices: tuple[int, ...]
    relevant_edges: tuple[int, ...]
    relevant_surfaces: tuple[int, ...]
    boundary_edges: int
    reasons: tuple[str, ...]
    counts: Mapping[str, int]


@dataclass(frozen=True)
class AWReuseAssessment:
    decision: AWReuseDecision
    certificate: AWLocalCompletenessCertificate
    comparison: AWInterfaceComparison | None
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class AWSolverProvenance:
    """Opt-in evidence emitted by, or observed from, one AW solver run.

    ``completeness_proven`` is intentionally separate from observed counts:
    an exhausted traversal or an empty table is not proof that no geometry was
    omitted.
    """

    schema_version: int
    enabled: bool
    solve_key: SolveUniverseKey
    canonical_generators: tuple[str, ...]
    active_competitor_neighborhood: Mapping[str, tuple[int, ...]]
    vertex_search: Mapping[str, Any]
    edge_construction: Mapping[str, Any]
    surface_construction: Mapping[str, Any]
    discarded_primitives: Mapping[str, Any]
    completeness: Mapping[str, Any]

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-safe diagnostic payload without exposing live tables."""
        return _thaw({
            "schema_version": self.schema_version,
            "enabled": self.enabled,
            "solve_key": self.solve_key.__dict__,
            "canonical_generators": self.canonical_generators,
            "active_competitor_neighborhood": self.active_competitor_neighborhood,
            "vertex_search": self.vertex_search,
            "edge_construction": self.edge_construction,
            "surface_construction": self.surface_construction,
            "discarded_primitives": self.discarded_primitives,
            "completeness": self.completeness,
        })


@dataclass(frozen=True)
class AWSolvedNetworkComparison:
    """Read-only comparison result for two independently solved networks."""

    view: AWInterfaceView
    assessment: AWReuseAssessment
    full_provenance: AWSolverProvenance
    interface_provenance: AWSolverProvenance


def build_solve_universe_key(
    network: Any,
    *,
    frame_id: Any = None,
    environment: Any = None,
    partition_scheme: Any = None,
    solver_version: str = DEFAULT_SOLVER_VERSION,
) -> SolveUniverseKey:
    """Build a deterministic key from network inputs and scientific settings."""
    balls = getattr(network, "balls", None)
    if balls is None or not {"loc", "rad"}.issubset(balls.columns):
        raise AWGeometryReuseError("AW geometry requires ball coordinates and radii.")

    settings = _settings_for_identity(getattr(network, "settings", None))
    scheme = str(settings.get("net_type", "aw")).lower()
    if scheme != "aw":
        raise AWGeometryReuseError(
            f"AW reuse received non-AW network type {scheme!r}."
        )

    generators = []
    seen_ids = set()
    for position, (_, row) in enumerate(balls.iterrows()):
        row_data = row.to_dict()
        stable_id = _stable_generator_id(row_data, position)
        if stable_id in seen_ids:
            raise AWGeometryReuseError(f"Duplicate generator identity: {stable_id}")
        seen_ids.add(stable_id)
        generators.append(
            {
                "stable_id": stable_id,
                "source_index": position,
                "loc": tuple(float(value) for value in row_data["loc"]),
                "rad": float(row_data["rad"]),
            }
        )
    generators.sort(key=lambda item: item["stable_id"])

    frame = frame_id
    if frame is None:
        frame = getattr(network, "frame_id", None)
    if frame is None:
        frame = getattr(network, "frame", "default")

    partition = partition_scheme
    if partition is None:
        partition = settings.get("partition_scheme", settings.get("partition", "default"))

    environment_value = environment
    if environment_value is None:
        environment_value = {
            "network_environment": getattr(network, "environment", None),
            "boundary_indices": getattr(network, "boundary_indices", ()),
            "boundary_config": getattr(network, "boundary_config", None),
            "boundary_mode": getattr(network, "boundary_mode", None),
            "boundary_generators": getattr(network, "boundary_generators", ()),
        }

    resolution = settings.get("surf_res", settings.get("surface_resolution"))
    solver_settings = {
        key: value
        for key, value in settings.items()
        if key not in {"net_type", "surf_res", "surface_resolution"}
    }
    identity_payload = {
        "schema_version": SCHEMA_VERSION,
        "scheme": scheme,
        "scope_mode": _scope_mode(network),
        "frame_id": str(frame),
        "partition_scheme": partition,
        "environment": environment_value,
        "generators": generators,
        "resolution": resolution,
        "solver_settings": solver_settings,
        "solver_version": solver_version,
    }

    geometry_digest = _digest(generators)
    generator_digest = _digest([item["stable_id"] for item in generators])
    environment_digest = _digest(environment_value)
    solver_settings_digest = _digest(solver_settings)
    digest = _digest(identity_payload)
    return SolveUniverseKey(
        schema_version=SCHEMA_VERSION,
        scheme=scheme,
        scope_mode=_scope_mode(network),
        frame_id=str(frame),
        partition_scheme=str(partition),
        environment_digest=environment_digest,
        generator_digest=generator_digest,
        geometry_digest=geometry_digest,
        resolution=resolution,
        solver_settings_digest=solver_settings_digest,
        solver_version=str(solver_version),
        digest=digest,
    )


def _aw_provenance_enabled(network: Any) -> bool:
    settings = getattr(network, "settings", None) or {}
    value = os.environ.get("VORPY_AW_PROVENANCE", "").strip().lower()
    return bool(
        settings.get("aw_provenance", False)
        or value in {"1", "true", "yes", "on"}
    )


def _json_mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(item) for key, item in value.items()})
    return MappingProxyType({})


def _active_vertex_competitors(network: Any, tolerance: float = 1e-9) -> Mapping[str, tuple[int, ...]]:
    """Record generators that can compete at each stored vertex location.

    This is deliberately vertex-local evidence.  It does not certify edge or
    surface neighborhoods, which require solver instrumentation along their
    full geometry.
    """
    balls = getattr(network, "balls", None)
    verts = getattr(network, "verts", None)
    if balls is None or verts is None or not {"loc", "rad"}.issubset(balls.columns):
        return MappingProxyType({})

    locs = np.asarray(balls["loc"].tolist(), dtype=float)
    rads = np.asarray(balls["rad"].tolist(), dtype=float)
    result = {}
    for position, (_, row) in enumerate(verts.iterrows()):
        location = row.get("loc")
        radius = row.get("rad")
        if location is None or radius is None:
            continue
        location = np.asarray(location, dtype=float)
        radius = float(radius)
        if location.shape != (3,) or not np.isfinite(radius):
            continue
        clearance = np.linalg.norm(locs - location, axis=1) - rads
        result[str(position)] = tuple(
            int(index) for index in np.flatnonzero(clearance <= radius + tolerance)
        )
    return MappingProxyType(result)


def collect_aw_solver_provenance(
    network: Any,
    *,
    tolerance: float = 1e-9,
) -> AWSolverProvenance:
    """Collect opt-in, fail-closed provenance from a completed AW Network.

    The function is read-only.  It consumes metadata recorded by the
    production vertex/topology builders when enabled and marks any missing
    evidence as unproven instead of inferring completeness from counts.
    """
    key = build_solve_universe_key(network)
    balls = getattr(network, "balls", None)
    generator_ids = tuple(
        _stable_generator_id(row.to_dict(), position)
        for position, (_, row) in enumerate(balls.iterrows())
    )
    vertex_search = _json_mapping(getattr(network, "aw_vertex_search_provenance", None))
    topology = _json_mapping(getattr(network, "aw_topology_provenance", None))
    edge_construction = _json_mapping(topology.get("edge_construction", {}))
    surface_construction = _json_mapping(topology.get("surface_construction", {}))
    discarded = _json_mapping(topology.get("discarded_primitives", {}))

    if not vertex_search:
        vertex_search = MappingProxyType({
            "status": "unobserved",
            "completeness_proven": False,
            "reason": "vertex builder provenance was not enabled",
        })
    if not edge_construction:
        edge_construction = MappingProxyType({
            "status": "observed_tables_only",
            "completeness_proven": False,
            "reason": "edge builder provenance was not enabled",
        })
    if not surface_construction:
        surface_construction = MappingProxyType({
            "status": "observed_tables_only",
            "completeness_proven": False,
            "reason": "surface builder provenance was not enabled",
        })

    completeness = MappingProxyType({
        "completeness_proven": False,
        "generator_universe_complete": bool(
            balls is not None and {"loc", "rad"}.issubset(balls.columns)
        ),
        "vertex_completeness_proven": bool(
            vertex_search.get("completeness_proven", False)
        ),
        "edge_completeness_proven": bool(
            edge_construction.get("completeness_proven", False)
        ),
        "surface_completeness_proven": bool(
            surface_construction.get("completeness_proven", False)
        ),
        "reason": "counts and exhausted traversal are observational evidence only",
    })
    return AWSolverProvenance(
        schema_version=1,
        enabled=_aw_provenance_enabled(network),
        solve_key=key,
        canonical_generators=generator_ids,
        active_competitor_neighborhood=_active_vertex_competitors(network, tolerance),
        vertex_search=vertex_search,
        edge_construction=edge_construction,
        surface_construction=surface_construction,
        discarded_primitives=discarded,
        completeness=completeness,
    )


def _require_frame(network: Any, field_names: Sequence[str], label: str) -> pd.DataFrame:
    frame = getattr(network, label, None)
    if frame is None or not isinstance(frame, pd.DataFrame):
        raise AWGeometryReuseError(f"AW geometry is missing the {label} table.")
    missing = set(field_names).difference(frame.columns)
    if missing:
        raise AWGeometryReuseError(
            f"AW {label} table is missing required fields: {sorted(missing)}"
        )
    return frame


def _row_mapping(row: pd.Series) -> Mapping[str, Any]:
    return MappingProxyType({str(key): _freeze(value) for key, value in row.to_dict().items()})


def _primitive_id(kind: str, key: SolveUniverseKey, payload: Any) -> str:
    # Scope belongs to the solve identity, not the primitive identity.  This
    # permits a compatible full-system primitive and interface view primitive
    # to share IDs while still keeping their Network build contexts distinct.
    source_identity = {
        "schema_version": key.schema_version,
        "scheme": key.scheme,
        "frame_id": key.frame_id,
        "partition_scheme": key.partition_scheme,
        "environment_digest": key.environment_digest,
        "generator_digest": key.generator_digest,
        "geometry_digest": key.geometry_digest,
        "resolution": key.resolution,
        "solver_settings_digest": key.solver_settings_digest,
        "solver_version": key.solver_version,
    }
    return f"aw:{kind}:{_digest({'source': source_identity, 'payload': payload})}"


def _as_int_tuple(value: Any, label: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item) for item in value)
    except (TypeError, ValueError) as error:
        raise AWGeometryReuseError(f"Invalid AW {label} incidence.") from error
    return values


def extract_aw_geometry(
    network: Any,
    *,
    frame_id: Any = None,
    environment: Any = None,
    partition_scheme: Any = None,
    solver_version: str = DEFAULT_SOLVER_VERSION,
    allow_incomplete: bool = False,
    measure_memory: bool = False,
) -> AWGeometryStore:
    """Extract completed AW tables without mutating or solving ``network``."""
    started = perf_counter()
    if measure_memory:
        tracemalloc.start()
    try:
        key = build_solve_universe_key(
            network,
            frame_id=frame_id,
            environment=environment,
            partition_scheme=partition_scheme,
            solver_version=solver_version,
        )
        balls = _require_frame(network, ("loc", "rad", "complete"), "balls")
        verts = _require_frame(network, ("balls", "loc", "rad"), "verts")
        edges = _require_frame(network, ("balls", "verts"), "edges")
        surfs = _require_frame(network, ("balls", "verts", "edges"), "surfs")

        global_complete = bool(balls["complete"].astype(bool).all())
        if not global_complete and not allow_incomplete:
            raise AWGeometryReuseError(
                "AW geometry is incomplete; exact interface reuse is disabled."
            )

        generator_rows = []
        for position, (_, row) in enumerate(balls.iterrows()):
            data = _row_mapping(row)
            generator_rows.append(
                AWGenerator(
                    stable_id=_stable_generator_id(row.to_dict(), position),
                    source_index=position,
                    loc=tuple(float(value) for value in row["loc"]),
                    rad=float(row["rad"]),
                    data=data,
                )
            )

        vertex_rows = []
        for position, (_, row) in enumerate(verts.iterrows()):
            balls_value = _as_int_tuple(row["balls"], "vertex balls")
            if any(ball < 0 or ball >= len(generator_rows) for ball in balls_value):
                raise AWGeometryReuseError("Vertex references an unknown generator.")
            payload = {
                "balls": balls_value,
                "doublet": row.get("dub", row.get("vdub")),
                "loc": row["loc"],
                "rad": row["rad"],
            }
            vertex_rows.append(
                AWVertexPrimitive(
                    primitive_id=_primitive_id("vertex", key, payload),
                    source_index=position,
                    balls=balls_value,
                    data=_row_mapping(row),
                )
            )

        edge_rows = []
        for position, (_, row) in enumerate(edges.iterrows()):
            balls_value = _as_int_tuple(row["balls"], "edge balls")
            vertices_value = _as_int_tuple(row["verts"], "edge vertices")
            if any(ball < 0 or ball >= len(generator_rows) for ball in balls_value):
                raise AWGeometryReuseError("Edge references an unknown generator.")
            if any(vertex < 0 or vertex >= len(vertex_rows) for vertex in vertices_value):
                raise AWGeometryReuseError("Edge references an unknown vertex.")
            payload = {"balls": balls_value, "vertices": vertices_value}
            edge_rows.append(
                AWEdgePrimitive(
                    primitive_id=_primitive_id("edge", key, payload),
                    source_index=position,
                    balls=balls_value,
                    vertices=vertices_value,
                    data=_row_mapping(row),
                )
            )

        surface_rows = []
        for position, (_, row) in enumerate(surfs.iterrows()):
            balls_value = _as_int_tuple(row["balls"], "surface balls")
            vertices_value = _as_int_tuple(row["verts"], "surface vertices")
            edges_value = _as_int_tuple(row["edges"], "surface edges")
            if any(ball < 0 or ball >= len(generator_rows) for ball in balls_value):
                raise AWGeometryReuseError("Surface references an unknown generator.")
            if any(vertex < 0 or vertex >= len(vertex_rows) for vertex in vertices_value):
                raise AWGeometryReuseError("Surface references an unknown vertex.")
            if any(edge < 0 or edge >= len(edge_rows) for edge in edges_value):
                raise AWGeometryReuseError("Surface references an unknown edge.")
            payload = {
                "balls": balls_value,
                "vertices": vertices_value,
                "edges": edges_value,
            }
            surface_rows.append(
                AWSurfacePrimitive(
                    primitive_id=_primitive_id("surface", key, payload),
                    source_index=position,
                    balls=balls_value,
                    vertices=vertices_value,
                    edges=edges_value,
                    data=_row_mapping(row),
                )
            )

        peak = None
        if measure_memory:
            _, peak = tracemalloc.get_traced_memory()
        return AWGeometryStore(
            key=key,
            generators=tuple(generator_rows),
            vertices=tuple(vertex_rows),
            edges=tuple(edge_rows),
            surfaces=tuple(surface_rows),
            table_columns=MappingProxyType({
                "balls": tuple(str(column) for column in balls.columns),
                "verts": tuple(str(column) for column in verts.columns),
                "edges": tuple(str(column) for column in edges.columns),
                "surfs": tuple(str(column) for column in surfs.columns),
            }),
            global_complete=global_complete,
            metrics=ReuseMetrics(
                operation="extract_aw_geometry",
                seconds=perf_counter() - started,
                peak_tracemalloc_bytes=peak,
            ),
        )
    finally:
        if measure_memory:
            tracemalloc.stop()


def _spans(groups: tuple[set[int], set[int]], balls: Sequence[int]) -> bool:
    return bool(set(balls) & groups[0]) and bool(set(balls) & groups[1])


def build_aw_interface_view(
    source: AWGeometryStore,
    group_a: Sequence[int],
    group_b: Sequence[int],
) -> AWInterfaceView:
    """Select an interface from immutable full-system AW primitives."""
    left = set(int(item) for item in group_a)
    right = set(int(item) for item in group_b)
    if not left or not right or left.intersection(right):
        raise AWGeometryReuseError("Interface groups must be non-empty and disjoint.")
    count = len(source.generators)
    if any(item < 0 or item >= count for item in left | right):
        raise AWGeometryReuseError("Interface group references an unknown generator.")
    groups = (left, right)

    selected_vertices = {
        vertex.source_index
        for vertex in source.vertices
        if _spans(groups, vertex.balls)
    }
    selected_vertex_ids = {
        vertex.source_index: vertex.primitive_id for vertex in source.vertices
        if vertex.source_index in selected_vertices
    }

    selected_edges = {
        edge.source_index
        for edge in source.edges
        if _spans(groups, edge.balls)
        and set(edge.vertices).issubset(selected_vertices)
    }
    selected_edge_ids = {
        edge.source_index: edge.primitive_id for edge in source.edges
        if edge.source_index in selected_edges
    }

    discarded: dict[str, int] = {}
    selected_surfaces = set()
    for surface in source.surfaces:
        if not _spans(groups, surface.balls):
            continue
        if not set(surface.vertices).issubset(selected_vertices):
            discarded["surface_vertex_not_selected"] = discarded.get(
                "surface_vertex_not_selected", 0
            ) + 1
            continue
        if not set(surface.edges).issubset(selected_edges):
            discarded["surface_edge_not_selected"] = discarded.get(
                "surface_edge_not_selected", 0
            ) + 1
            continue
        selected_surfaces.add(surface.source_index)

    selected_surface_ids = {
        surface.source_index: surface.primitive_id for surface in source.surfaces
        if surface.source_index in selected_surfaces
    }
    return AWInterfaceView(
        source=source,
        group_a=tuple(sorted(left)),
        group_b=tuple(sorted(right)),
        vertex_ids=tuple(
            selected_vertex_ids[index]
            for index in sorted(selected_vertex_ids)
        ),
        edge_ids=tuple(
            selected_edge_ids[index]
            for index in sorted(selected_edge_ids)
        ),
        surface_ids=tuple(
            selected_surface_ids[index]
            for index in sorted(selected_surface_ids)
        ),
        discarded=MappingProxyType(dict(sorted(discarded.items()))),
    )


def _source_index_map(records: Sequence[Any]) -> dict[int, Any]:
    return {record.source_index: record for record in records}


def _sample_points(record: Any) -> tuple[np.ndarray, ...]:
    values = record.data.get("points")
    if values is None:
        values = record.data.get("loc", record.data.get("center"))
        if values is None:
            return ()
        values = (values,)
    try:
        return tuple(np.asarray(value, dtype=float) for value in values)
    except (TypeError, ValueError):
        return ()


def _surface_boundary_issues(
    surface: AWSurfacePrimitive,
    edges: Mapping[int, AWEdgePrimitive],
) -> tuple[tuple[str, ...], int]:
    issues = []
    vertices = set(surface.vertices)
    if len(surface.vertices) < 3:
        issues.append("surface_too_small")
    if len(surface.edges) != len(surface.vertices):
        issues.append("edge_vertex_count_mismatch")

    boundary_degree = {vertex: 0 for vertex in vertices}
    adjacency = {vertex: set() for vertex in vertices}
    for edge_index in surface.edges:
        edge = edges.get(edge_index)
        if edge is None or len(edge.vertices) != 2:
            issues.append("missing_or_invalid_boundary_edge")
            continue
        endpoints = set(edge.vertices)
        if not endpoints.issubset(vertices):
            issues.append("boundary_edge_outside_surface")
            continue
        left, right = edge.vertices
        boundary_degree[left] += 1
        boundary_degree[right] += 1
        adjacency[left].add(right)
        adjacency[right].add(left)

    boundary_edges = sum(value != 2 for value in boundary_degree.values())
    if boundary_edges:
        issues.append("open_boundary")

    if adjacency:
        seen = set()
        stack = [next(iter(adjacency))]
        while stack:
            vertex = stack.pop()
            if vertex in seen:
                continue
            seen.add(vertex)
            stack.extend(adjacency[vertex] - seen)
        if len(seen) != len(adjacency):
            issues.append("disconnected_boundary")
    return tuple(sorted(set(issues))), boundary_edges


def certify_aw_interface_completeness(
    view: AWInterfaceView,
    *,
    distance_tolerance: float = 1e-9,
) -> AWLocalCompletenessCertificate:
    """Prove only the local conditions represented by the extracted tables.

    The relevant topology is the full vertex/edge/surface component touching
    the selected interface.  Generator completeness is checked for that
    component and for a conservative geometric competitor neighborhood around
    its stored samples.  Missing samples prevent certification rather than
    being inferred away.
    """
    source = view.source
    vertices = _source_index_map(source.vertices)
    edges = _source_index_map(source.edges)
    surfaces = _source_index_map(source.surfaces)
    selected_vertices = {
        record.source_index
        for record in source.vertices
        if record.primitive_id in set(view.vertex_ids)
    }
    selected_edges = {
        record.source_index
        for record in source.edges
        if record.primitive_id in set(view.edge_ids)
    }
    selected_surfaces = {
        record.source_index
        for record in source.surfaces
        if record.primitive_id in set(view.surface_ids)
    }
    reasons = set()
    if not selected_vertices:
        reasons.add("missing_interface_vertex")
    if not selected_edges:
        reasons.add("missing_interface_edge")
    if not selected_surfaces:
        reasons.add("missing_interface_surface")
    if "surface_vertex_not_selected" in view.discarded:
        reasons.add("missing_interface_vertex")
    if "surface_edge_not_selected" in view.discarded:
        reasons.add("missing_interface_edge")

    # Find the complete topological component touching the interface. This
    # includes non-interface support primitives, which prevents an incomplete
    # local competitor from being mistaken for an unrelated cell.
    vertex_neighbors = {index: set() for index in vertices}
    for edge in source.edges:
        for left in edge.vertices:
            for right in edge.vertices:
                if left != right:
                    vertex_neighbors.setdefault(left, set()).add(right)
    relevant_vertices = set(selected_vertices)
    stack = list(selected_vertices)
    while stack:
        current = stack.pop()
        for neighbor in vertex_neighbors.get(current, ()):
            if neighbor not in relevant_vertices:
                relevant_vertices.add(neighbor)
                stack.append(neighbor)

    relevant_edges = {
        edge.source_index
        for edge in source.edges
        if set(edge.vertices) & relevant_vertices
    }
    relevant_surfaces = {
        surface.source_index
        for surface in source.surfaces
        if set(surface.vertices) & relevant_vertices
        or set(surface.edges) & relevant_edges
    }

    relevant_generators = set(view.group_a) | set(view.group_b)
    for index in relevant_vertices:
        record = vertices.get(index)
        if record is not None:
            relevant_generators.update(record.balls)
    for index in relevant_edges:
        record = edges.get(index)
        if record is not None:
            relevant_generators.update(record.balls)
    for index in relevant_surfaces:
        record = surfaces.get(index)
        if record is not None:
            relevant_generators.update(record.balls)

    # Every local primitive needs stored samples to establish a conservative
    # competitor neighborhood. Without them, topology alone is insufficient.
    local_records = [vertices[index] for index in sorted(relevant_vertices) if index in vertices]
    local_records.extend(edges[index] for index in sorted(relevant_edges) if index in edges)
    local_records.extend(surfaces[index] for index in sorted(relevant_surfaces) if index in surfaces)
    local_samples = []
    max_clearance = 0.0
    for record in local_records:
        samples = _sample_points(record)
        if not samples:
            reasons.add("competitor_relationships_unavailable")
            continue
        local_samples.extend(samples)
        if isinstance(record, AWVertexPrimitive):
            max_clearance = max(max_clearance, abs(float(record.data.get("rad", 0.0))))
        else:
            endpoint_radii = [
                abs(float(vertices[vertex].data.get("rad", 0.0)))
                for vertex in getattr(record, "vertices", ())
                if vertex in vertices
            ]
            if endpoint_radii:
                max_clearance = max(max_clearance, max(endpoint_radii))

    generator_map = _source_index_map(source.generators)
    for generator in source.generators:
        if not local_samples:
            break
        location = np.asarray(generator.loc, dtype=float)
        if any(
            np.linalg.norm(location - sample)
            <= max_clearance + abs(float(generator.rad)) + distance_tolerance
            for sample in local_samples
        ):
            relevant_generators.add(generator.source_index)

    incomplete = [
        index
        for index in sorted(relevant_generators)
        if index not in generator_map
        or not bool(generator_map[index].data.get("complete", False))
    ]
    if incomplete:
        reasons.add("incomplete_relevant_competitor")

    boundary_edges = 0
    for surface_index in sorted(relevant_surfaces):
        surface = surfaces[surface_index]
        issues, count = _surface_boundary_issues(surface, edges)
        boundary_edges += count
        if issues:
            if surface_index in selected_surfaces:
                reasons.add("open_interface_surface_boundary")
            else:
                reasons.add("open_local_surface_boundary")

    counts = MappingProxyType({
        "selected_vertices": len(selected_vertices),
        "selected_edges": len(selected_edges),
        "selected_surfaces": len(selected_surfaces),
        "relevant_vertices": len(relevant_vertices),
        "relevant_edges": len(relevant_edges),
        "relevant_surfaces": len(relevant_surfaces),
        "relevant_generators": len(relevant_generators),
    })
    return AWLocalCompletenessCertificate(
        proven=not reasons,
        globally_complete=source.global_complete,
        relevant_generators=tuple(sorted(relevant_generators)),
        relevant_vertices=tuple(sorted(relevant_vertices)),
        relevant_edges=tuple(sorted(relevant_edges)),
        relevant_surfaces=tuple(sorted(relevant_surfaces)),
        boundary_edges=boundary_edges,
        reasons=tuple(sorted(reasons)),
        counts=counts,
    )


def _lookup(records: Sequence[Any]) -> dict[str, Any]:
    return {record.primitive_id: record for record in records}


def materialize_aw_interface(
    view: AWInterfaceView,
    *,
    measure_memory: bool = False,
) -> AWMaterializedInterface:
    """Materialize detached Network-shaped tables; never mutate the source."""
    started = perf_counter()
    if measure_memory:
        tracemalloc.start()
    try:
        source = view.source
        vertices_by_id = _lookup(source.vertices)
        edges_by_id = _lookup(source.edges)
        surfaces_by_id = _lookup(source.surfaces)
        vertices = [vertices_by_id[item] for item in view.vertex_ids]
        edges = [edges_by_id[item] for item in view.edge_ids]
        surfaces = [surfaces_by_id[item] for item in view.surface_ids]

        vertex_position = {record.source_index: index for index, record in enumerate(vertices)}
        edge_position = {record.source_index: index for index, record in enumerate(edges)}
        surface_position = {record.source_index: index for index, record in enumerate(surfaces)}

        vertex_rows = []
        for record in vertices:
            row = _thaw(record.data)
            row["edges"] = [edge_position[index] for index in row.get("edges", ()) if index in edge_position]
            row["surfs"] = [surface_position[index] for index in row.get("surfs", ()) if index in surface_position]
            vertex_rows.append(row)

        edge_rows = []
        for record in edges:
            row = _thaw(record.data)
            row["verts"] = [vertex_position[index] for index in row.get("verts", ()) if index in vertex_position]
            row["surfs"] = [surface_position[index] for index in row.get("surfs", ()) if index in surface_position]
            edge_rows.append(row)

        surface_rows = []
        for record in surfaces:
            row = _thaw(record.data)
            row["verts"] = [vertex_position[index] for index in row.get("verts", ()) if index in vertex_position]
            row["edges"] = [edge_position[index] for index in row.get("edges", ()) if index in edge_position]
            surface_rows.append(row)

        ball_rows = [_thaw(record.data) for record in source.generators]
        for position, row in enumerate(ball_rows):
            row["verts"] = [
                index for index, record in enumerate(vertices)
                if position in record.balls
            ]
            row["edges"] = [
                index for index, record in enumerate(edges)
                if position in record.balls
            ]
            row["surfs"] = [
                index for index, record in enumerate(surfaces)
                if position in record.balls
            ]

        peak = None
        if measure_memory:
            _, peak = tracemalloc.get_traced_memory()
        return AWMaterializedInterface(
            key=source.key,
            group_a=view.group_a,
            group_b=view.group_b,
            balls=pd.DataFrame(ball_rows, columns=source.table_columns["balls"]),
            verts=pd.DataFrame(vertex_rows, columns=source.table_columns["verts"]),
            edges=pd.DataFrame(edge_rows, columns=source.table_columns["edges"]),
            surfs=pd.DataFrame(surface_rows, columns=source.table_columns["surfs"]),
            metrics=ReuseMetrics(
                operation="materialize_aw_interface",
                seconds=perf_counter() - started,
                peak_tracemalloc_bytes=peak,
            ),
        )
    finally:
        if measure_memory:
            tracemalloc.stop()


def _record_signature(record: Mapping[str, Any], *, exclude: set[str] | None = None) -> str:
    data = {
        key: value
        for key, value in record.items()
        if key not in (exclude or set())
    }
    return _digest(data)


def _table_signatures(network: Any) -> dict[str, set[str]]:
    """Compare topology by content, not by network-local row numbers."""
    def rows(name: str) -> list[dict[str, Any]]:
        table = getattr(network, name, None)
        if table is None or not isinstance(table, pd.DataFrame):
            return []
        return [row.to_dict() for _, row in table.iterrows()]

    vertex_rows = rows("verts")
    vertex_signatures = [
        _record_signature(row, exclude={"edges", "surfs"})
        for row in vertex_rows
    ]
    edge_rows = rows("edges")
    edge_signatures = []
    for row in edge_rows:
        payload = {
            key: value
            for key, value in row.items()
            if key not in {"verts", "surfs"}
        }
        payload["vertex_incidence"] = tuple(
            vertex_signatures[int(index)] for index in row.get("verts", ())
        )
        edge_signatures.append(_digest(payload))

    surface_rows = rows("surfs")
    surface_signatures = []
    for row in surface_rows:
        payload = {
            key: value
            for key, value in row.items()
            if key not in {"verts", "edges"}
        }
        payload["vertex_incidence"] = tuple(
            vertex_signatures[int(index)] for index in row.get("verts", ())
        )
        payload["edge_incidence"] = tuple(
            edge_signatures[int(index)] for index in row.get("edges", ())
        )
        surface_signatures.append(_digest(payload))

    return {
        "vertices": set(vertex_signatures),
        "edges": set(edge_signatures),
        "surfaces": set(surface_signatures),
    }


def _report_value(value: Any) -> Any:
    return _thaw(_freeze(value))


def _canonical_cycle(values: Sequence[Any]) -> tuple[Any, ...]:
    values = tuple(values)
    if not values:
        return ()
    candidates = []
    for sequence in (values, tuple(reversed(values))):
        candidates.extend(
            sequence[offset:] + sequence[:offset]
            for offset in range(len(sequence))
        )
    return min(candidates, key=repr)


def _close_values(left: Any, right: Any, *, tolerance: float = 1.0e-10) -> bool:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return (
            set(left) == set(right)
            and all(_close_values(left[key], right[key], tolerance=tolerance) for key in left)
        )
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        return (
            len(left) == len(right)
            and all(_close_values(a, b, tolerance=tolerance) for a, b in zip(left, right))
        )
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)
    return left == right


def _tolerance_key(value: Any, *, tolerance: float = 1.0e-10) -> Any:
    if isinstance(value, Mapping):
        return tuple(
            (str(key), _tolerance_key(item, tolerance=tolerance))
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        )
    if isinstance(value, (list, tuple)):
        return tuple(_tolerance_key(item, tolerance=tolerance) for item in value)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return round(float(value) / tolerance)
    return value


def _close_multiset(
    left: Sequence[Any],
    right: Sequence[Any],
    *,
    tolerance: float = 1.0e-10,
) -> bool:
    return sorted(
        (_tolerance_key(value, tolerance=tolerance) for value in left),
        key=repr,
    ) == sorted(
        (_tolerance_key(value, tolerance=tolerance) for value in right),
        key=repr,
    )


def _surface_triangle_coordinates(record: Mapping[str, Any]) -> list[Any]:
    points = record.get("points", ())
    triangles = record.get("tris", ())
    result = []
    for triangle in triangles:
        try:
            result.append(tuple(sorted((points[int(index)] for index in triangle), key=repr)))
        except (IndexError, TypeError, ValueError):
            result.append(triangle)
    return result


def _physical_geometry_close(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    *,
    kind: str,
) -> bool:
    if left.get("balls") != right.get("balls"):
        return False
    if kind == "edges":
        left_points = left.get("points", ())
        right_points = right.get("points", ())
        return (
            _close_values(left.get("vals"), right.get("vals"))
            and _close_values(left.get("length"), right.get("length"))
            and (
                _close_values(left_points, right_points)
                or _close_values(left_points, list(reversed(right_points)))
            )
        )
    if kind == "surfaces":
        return (
            _close_values(left.get("sa"), right.get("sa"))
            and _close_values(left.get("orientation"), right.get("orientation"))
            and _close_multiset(left.get("points", ()), right.get("points", ()))
            and _close_multiset(
                _surface_triangle_coordinates(left),
                _surface_triangle_coordinates(right),
            )
        )
    return _close_values(left.get("loc"), right.get("loc")) and _close_values(
        left.get("rad"), right.get("rad")
    )


def _canonical_table_records(
    network: Any,
    *,
    vertex_source_ids: Sequence[str] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """Return deterministic, row-index-independent records for diagnostics."""
    def rows(name: str) -> list[dict[str, Any]]:
        table = getattr(network, name, None)
        if table is None or not isinstance(table, pd.DataFrame):
            return []
        return [row.to_dict() for _, row in table.iterrows()]

    vertex_rows = rows("verts")
    vertices = []
    vertex_ids = []
    for position, row in enumerate(vertex_rows):
        balls = tuple(int(value) for value in row.get("balls", ()))
        payload = {
            "balls": balls,
            "loc": _canonical(row.get("loc")),
            "rad": _canonical(row.get("rad")),
            "dub": _canonical(row.get("dub")),
        }
        identity = _digest(payload)
        vertex_ids.append(identity)
        vertices.append({
            "source_id": (
                vertex_source_ids[position]
                if vertex_source_ids is not None and position < len(vertex_source_ids)
                else None
            ),
            "identity": identity,
            "topology_key": _digest({"generator_set": tuple(sorted(balls))}),
            "geometry_key": _digest({"loc": payload["loc"], "rad": payload["rad"]}),
            "incidence_key": _digest({"balls": tuple(sorted(balls))}),
            "physical_geometry_key": _digest({"loc": payload["loc"], "rad": payload["rad"]}),
            "physical_incidence_key": _digest({"balls": tuple(sorted(balls))}),
            "balls": _report_value(balls),
            "loc": _report_value(row.get("loc")),
            "rad": _report_value(row.get("rad")),
            "dub": _report_value(row.get("dub")),
        })

    edge_rows = rows("edges")
    edges = []
    edge_topology_ids = []
    edge_geometry_ids = []
    for row in edge_rows:
        balls = tuple(int(value) for value in row.get("balls", ()))
        endpoints = tuple(
            vertex_ids[int(index)] if 0 <= int(index) < len(vertex_ids)
            else f"missing-vertex:{index}"
            for index in row.get("verts", ())
        )
        geometry_payload = {
            "balls": tuple(sorted(balls)),
            "points": _canonical(row.get("points", ())),
            "vals": _canonical(row.get("vals", ())),
            "length": _canonical(row.get("length")),
        }
        curve_key = min(
            _digest(geometry_payload["points"]),
            _digest(tuple(reversed(geometry_payload["points"]))),
        )
        incidence_payload = {
            "balls": tuple(sorted(balls)),
            "endpoints": endpoints,
        }
        geometry_key = _digest(geometry_payload)
        incidence_key = _digest(incidence_payload)
        physical_geometry_key = _digest({
            "balls": tuple(sorted(balls)),
            "curve": curve_key,
            "vals": geometry_payload["vals"],
            "length": geometry_payload["length"],
        })
        physical_incidence_key = _digest({
            "balls": tuple(sorted(balls)),
            "endpoints": tuple(sorted(endpoints)),
        })
        topology_key = _digest({"geometry_generators": tuple(sorted(balls)), "incidence": physical_incidence_key})
        identity = _digest({"geometry": geometry_key, "incidence": incidence_key})
        edge_topology_ids.append(topology_key)
        edge_geometry_ids.append(geometry_key)
        edges.append({
            "identity": identity,
            "topology_key": topology_key,
            "geometry_key": geometry_key,
            "physical_geometry_key": physical_geometry_key,
            "incidence_key": incidence_key,
            "physical_incidence_key": physical_incidence_key,
            "balls": _report_value(balls),
            "verts": _report_value(endpoints),
            "points": _report_value(row.get("points", ())),
            "vals": _report_value(row.get("vals", ())),
            "length": _report_value(row.get("length")),
        })

    surface_rows = rows("surfs")
    surfaces = []
    for row in surface_rows:
        balls = tuple(int(value) for value in row.get("balls", ()))
        vertices_incidence = tuple(
            vertex_ids[int(index)] if 0 <= int(index) < len(vertex_ids)
            else f"missing-vertex:{index}"
            for index in row.get("verts", ())
        )
        edges_incidence = tuple(
            edge_topology_ids[int(index)] if 0 <= int(index) < len(edge_topology_ids)
            else f"missing-edge:{index}"
            for index in row.get("edges", ())
        )
        geometry_payload = {
            "balls": tuple(sorted(balls)),
            "points": _canonical(row.get("points", ())),
            "tris": _canonical(row.get("tris", ())),
            "orientation": _canonical(row.get("orientation")),
            "sa": _canonical(row.get("sa")),
        }
        point_values = row.get("points", ()) or ()
        point_keys = tuple(_digest(_canonical(point)) for point in point_values)
        triangle_keys = []
        for triangle in row.get("tris", ()) or ():
            try:
                triangle_keys.append(tuple(sorted(
                    point_keys[int(index)] for index in triangle
                )))
            except (IndexError, TypeError, ValueError):
                triangle_keys.append(_canonical(triangle))
        incidence_payload = {
            "balls": tuple(sorted(balls)),
            "verts": vertices_incidence,
            "edges": edges_incidence,
        }
        geometry_key = _digest(geometry_payload)
        incidence_key = _digest(incidence_payload)
        physical_geometry_key = _digest({
            "balls": tuple(sorted(balls)),
            "points": tuple(sorted(point_keys)),
            "triangles": tuple(sorted(triangle_keys, key=repr)),
            "sa": geometry_payload["sa"],
        })
        physical_incidence_key = _digest({
            "balls": tuple(sorted(balls)),
            "verts": tuple(sorted(vertices_incidence)),
            "edges": tuple(sorted(edges_incidence)),
        })
        topology_key = _digest({"geometry_generators": tuple(sorted(balls)), "incidence": physical_incidence_key})
        identity = _digest({"geometry": geometry_key, "incidence": incidence_key})
        surfaces.append({
            "identity": identity,
            "topology_key": topology_key,
            "geometry_key": geometry_key,
            "physical_geometry_key": physical_geometry_key,
            "incidence_key": incidence_key,
            "physical_incidence_key": physical_incidence_key,
            "balls": _report_value(balls),
            "verts": _report_value(vertices_incidence),
            "edges": _report_value(edges_incidence),
            "points": _report_value(row.get("points", ())),
            "tris": _report_value(row.get("tris", ())),
            "orientation": _report_value(row.get("orientation")),
            "sa": _report_value(row.get("sa")),
        })

    return {"vertices": vertices, "edges": edges, "surfaces": surfaces}


def _compare_canonical_records(
    extracted: Sequence[Mapping[str, Any]],
    direct: Sequence[Mapping[str, Any]],
    *,
    kind: str,
) -> dict[str, Any]:
    def grouped(records: Sequence[Mapping[str, Any]], key: str) -> dict[str, list[Mapping[str, Any]]]:
        result: dict[str, list[Mapping[str, Any]]] = {}
        for record in records:
            result.setdefault(str(record[key]), []).append(record)
        for values in result.values():
            values.sort(key=lambda item: str(item["identity"]))
        return result

    left = grouped(extracted, "topology_key")
    right = grouped(direct, "topology_key")
    only_left = sorted(set(left).difference(right))
    only_right = sorted(set(right).difference(left))
    shared_differences = []
    for key in sorted(set(left).intersection(right)):
        left_records = left[key]
        right_records = right[key]
        if [record["identity"] for record in left_records] == [record["identity"] for record in right_records]:
            continue
        shared_differences.append({
            "topology_key": key,
            "extracted": _report_value(left_records),
            "direct": _report_value(right_records),
            "geometry_equal": [record["geometry_key"] for record in left_records] == [record["geometry_key"] for record in right_records],
            "incidence_equal": [record["incidence_key"] for record in left_records] == [record["incidence_key"] for record in right_records],
        })
    topology_equivalent = not only_left and not only_right
    geometry_equivalent = not shared_differences or all(
        item["geometry_equal"] for item in shared_differences
    )
    incidence_equivalent = not shared_differences or all(
        item["incidence_equal"] for item in shared_differences
    )
    physical_geometry_equivalent = topology_equivalent and all(
        [record["physical_geometry_key"] for record in left[key]]
        == [record["physical_geometry_key"] for record in right[key]]
        for key in sorted(set(left).intersection(right))
    )
    physical_geometry_within_tolerance = topology_equivalent and all(
        len(left[key]) == len(right[key])
        and all(
            _physical_geometry_close(left_record, right_record, kind=kind)
            for left_record, right_record in zip(left[key], right[key])
        )
        for key in sorted(set(left).intersection(right))
    )
    physical_incidence_equivalent = topology_equivalent and all(
        [record["physical_incidence_key"] for record in left[key]]
        == [record["physical_incidence_key"] for record in right[key]]
        for key in sorted(set(left).intersection(right))
    )
    return {
        "extracted_count": len(extracted),
        "direct_count": len(direct),
        "extracted_records": _report_value(extracted),
        "direct_records": _report_value(direct),
        "topology_equivalent": topology_equivalent,
        "geometry_equivalent": geometry_equivalent and topology_equivalent,
        "incidence_equivalent": incidence_equivalent and topology_equivalent,
        "physical_geometry_equivalent": physical_geometry_equivalent,
        "physical_geometry_within_tolerance": physical_geometry_within_tolerance,
        "physical_incidence_equivalent": physical_incidence_equivalent,
        "canonical_equivalent": topology_equivalent and not shared_differences,
        "representation_only": (
            topology_equivalent
            and physical_geometry_within_tolerance
            and physical_incidence_equivalent
            and bool(shared_differences)
        ),
        "only_extracted": [_report_value(left[key]) for key in only_left],
        "only_direct": [_report_value(right[key]) for key in only_right],
        "shared_differences": shared_differences,
    }


def diagnose_aw_interface_difference(
    view: AWInterfaceView,
    direct_network: Any,
) -> dict[str, Any]:
    """Return a deterministic first-difference report without mutating either network."""
    materialized = view.materialize()
    extracted_records = _canonical_table_records(
        materialized,
        vertex_source_ids=view.vertex_ids,
    )
    direct_records = _canonical_table_records(direct_network)
    generator_records = {
        "extracted": [
            {
                "stable_id": generator.stable_id,
                "source_index": generator.source_index,
                "loc": _report_value(generator.loc),
                "rad": generator.rad,
            }
            for generator in view.source.generators
        ],
        "direct": [
            {
                "stable_id": _stable_generator_id(row.to_dict(), position),
                "source_index": position,
                "loc": _report_value(row.get("loc")),
                "rad": _report_value(row.get("rad")),
            }
            for position, (_, row) in enumerate(direct_network.balls.iterrows())
        ],
    }
    stages = {
        kind: _compare_canonical_records(
            extracted_records[kind], direct_records[kind], kind=kind
        )
        for kind in ("vertices", "edges", "surfaces")
    }
    direct_key = build_solve_universe_key(direct_network)
    if view.source.key.competitor_digest != direct_key.competitor_digest:
        first_divergence = "generator_universe"
    elif not stages["vertices"]["canonical_equivalent"]:
        first_divergence = (
            "vertex_orientation_or_canonical_order"
            if stages["vertices"]["representation_only"]
            else "vertex_discovery_or_group_filtering"
        )
    elif not stages["edges"]["canonical_equivalent"]:
        first_divergence = (
            "edge_sample_roundoff"
            if (
                stages["edges"]["representation_only"]
                and stages["edges"]["incidence_equivalent"]
            )
            else "edge_orientation_or_canonical_order"
            if stages["edges"]["representation_only"]
            else "edge_construction"
        )
    elif not stages["surfaces"]["canonical_equivalent"]:
        first_divergence = (
            "surface_orientation_or_canonical_order"
            if stages["surfaces"]["representation_only"]
            else "surface_construction_or_mesh"
        )
    elif view.source.key.digest != direct_key.digest:
        first_divergence = "solver_context"
    else:
        first_divergence = "none"
    first_physical_geometry_difference = next(
        (
            kind
            for kind in ("vertices", "edges", "surfaces")
            if not stages[kind]["physical_geometry_within_tolerance"]
        ),
        "none",
    )

    direct_diagnostics = getattr(direct_network, "interface_topology_diagnostics", None)
    direct_provenance = getattr(direct_network, "aw_solver_provenance", None)
    return {
        "first_divergence": first_divergence,
        "first_physical_geometry_difference": first_physical_geometry_difference,
        "generator_identities": generator_records,
        "stages": stages,
        "discarded_candidates": {
            "extracted_view": _report_value(view.discarded),
            "direct_interface_diagnostics": _report_value(direct_diagnostics or {}),
            "direct_solver_provenance": _report_value(direct_provenance or {}),
        },
        "source_scope": view.source.key.scope_mode,
        "direct_scope": direct_key.scope_mode,
        "same_competitor_universe": view.source.key.competitor_digest == direct_key.competitor_digest,
        "compatible_geometry_context": all((
            view.source.key.scheme == direct_key.scheme,
            view.source.key.frame_id == direct_key.frame_id,
            view.source.key.partition_scheme == direct_key.partition_scheme,
            view.source.key.resolution == direct_key.resolution,
            view.source.key.solver_settings_digest == direct_key.solver_settings_digest,
            view.source.key.solver_version == direct_key.solver_version,
        )),
    }


def compare_aw_interface_view(
    view: AWInterfaceView,
    direct_network: Any,
) -> AWInterfaceComparison:
    """Compare canonical extracted records with a directly built interface."""
    materialized = view.materialize()
    extracted = _table_signatures(materialized)
    direct = _table_signatures(direct_network)
    mismatches = {
        kind: tuple(sorted(extracted[kind] ^ direct[kind]))
        for kind in extracted
        if extracted[kind] != direct[kind]
    }
    direct_key = build_solve_universe_key(direct_network)
    same_competitor = view.source.key.competitor_digest == direct_key.competitor_digest
    compatible_context = same_competitor and all(
        (
            view.source.key.scheme == direct_key.scheme,
            view.source.key.frame_id == direct_key.frame_id,
            view.source.key.partition_scheme == direct_key.partition_scheme,
            view.source.key.resolution == direct_key.resolution,
            view.source.key.solver_settings_digest == direct_key.solver_settings_digest,
            view.source.key.solver_version == direct_key.solver_version,
        )
    )
    if not same_competitor:
        mismatches["competitor_universe"] = (
            view.source.key.competitor_digest,
            direct_key.competitor_digest,
        )
    if not compatible_context:
        mismatches["solve_context"] = (
            view.source.key.digest,
            direct_key.digest,
        )
    counts = {
        kind: {"extracted": len(extracted[kind]), "direct": len(direct[kind])}
        for kind in extracted
    }
    return AWInterfaceComparison(
        equivalent=same_competitor and compatible_context and not mismatches,
        same_competitor_universe=same_competitor,
        compatible_geometry_context=compatible_context,
        source_scope=view.source.key.scope_mode,
        direct_scope=direct_key.scope_mode,
        counts=MappingProxyType(counts),
        mismatches=MappingProxyType(mismatches),
    )


def compare_solved_aw_networks(
    full_network: Any,
    interface_network: Any,
    *,
    group_a: Sequence[int] | None = None,
    group_b: Sequence[int] | None = None,
    tolerance: float = 1e-9,
) -> AWSolvedNetworkComparison:
    """Compare independently solved full and interface AW Networks.

    The full network is extracted into a detached interface view.  Neither
    source network is mutated, and no production build method is invoked.
    Missing provenance remains visible in the returned assessment rather than
    being treated as evidence of equivalence.
    """
    iface_groups = getattr(interface_network, "iface_grps", None)
    if group_a is None or group_b is None:
        if iface_groups is None or len(iface_groups) != 2:
            raise AWGeometryReuseError(
                "Independent interface comparison requires two interface groups."
            )
        group_a, group_b = iface_groups

    source = extract_aw_geometry(full_network, allow_incomplete=True)
    view = build_aw_interface_view(source, group_a, group_b)
    certificate = certify_aw_interface_completeness(view, distance_tolerance=tolerance)
    assessment = assess_aw_reuse(
        view,
        direct_network=interface_network,
        certificate=certificate,
    )
    return AWSolvedNetworkComparison(
        view=view,
        assessment=assessment,
        full_provenance=collect_aw_solver_provenance(full_network),
        interface_provenance=collect_aw_solver_provenance(interface_network),
    )


def assess_aw_reuse(
    view: AWInterfaceView,
    *,
    direct_network: Any = None,
    certificate: AWLocalCompletenessCertificate | None = None,
) -> AWReuseAssessment:
    """Return an explicit, fail-closed decision for prototype reuse."""
    certificate = certificate or certify_aw_interface_completeness(view)
    comparison = compare_aw_interface_view(view, direct_network) if direct_network is not None else None

    if not certificate.proven:
        topology_reasons = {
            "missing_interface_vertex",
            "missing_interface_edge",
        }
        comparison_reasons = set(comparison.mismatches) if comparison is not None else set()
        reasons = tuple(sorted(set(certificate.reasons).union(comparison_reasons)))
        if comparison is not None and (
            not comparison.same_competitor_universe
            or not comparison.compatible_geometry_context
        ):
            decision = AWReuseDecision.REJECT_ENVIRONMENT
        elif topology_reasons.intersection(certificate.reasons) or any(
            key in comparison_reasons for key in ("vertices", "edges", "surfaces")
        ):
            decision = AWReuseDecision.REJECT_TOPOLOGY
        elif "competitor_relationships_unavailable" in certificate.reasons:
            decision = AWReuseDecision.UNRESOLVED
        else:
            decision = AWReuseDecision.REJECT_INCOMPLETE
        return AWReuseAssessment(
            decision=decision,
            certificate=certificate,
            comparison=comparison,
            reasons=reasons,
        )

    if direct_network is None:
        return AWReuseAssessment(
            decision=AWReuseDecision.UNRESOLVED,
            certificate=certificate,
            comparison=None,
            reasons=("direct_interface_equivalence_not_checked",),
        )

    if not comparison.same_competitor_universe or not comparison.compatible_geometry_context:
        decision = AWReuseDecision.REJECT_ENVIRONMENT
    elif any(
        key in comparison.mismatches
        for key in ("vertices", "edges", "surfaces")
    ):
        decision = AWReuseDecision.REJECT_TOPOLOGY
    elif comparison.equivalent:
        decision = AWReuseDecision.EXACT_REUSE
    else:
        decision = AWReuseDecision.UNRESOLVED
    return AWReuseAssessment(
        decision=decision,
        certificate=certificate,
        comparison=comparison,
        reasons=tuple(sorted(comparison.mismatches)),
    )


__all__ = [
    "AWGeometryReuseError",
    "AWGeometryStore",
    "AWGenerator",
    "AWInterfaceComparison",
    "AWInterfaceView",
    "AWMaterializedInterface",
    "AWLocalCompletenessCertificate",
    "AWReuseAssessment",
    "AWReuseDecision",
    "AWSolvedNetworkComparison",
    "AWSolverProvenance",
    "AWEdgePrimitive",
    "AWSurfacePrimitive",
    "AWVertexPrimitive",
    "DEFAULT_SOLVER_VERSION",
    "ReuseMetrics",
    "SolveUniverseKey",
    "assess_aw_reuse",
    "build_aw_interface_view",
    "build_solve_universe_key",
    "certify_aw_interface_completeness",
    "compare_aw_interface_view",
    "compare_solved_aw_networks",
    "collect_aw_solver_provenance",
    "diagnose_aw_interface_difference",
    "extract_aw_geometry",
    "materialize_aw_interface",
]
