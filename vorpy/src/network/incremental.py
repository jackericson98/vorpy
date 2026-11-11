"""Bounded, fail-closed incremental group completion.

This module is deliberately independent of CLI, export, and archive code.  It
supports two exactness levels:

* a complete source network is exposed as an immutable :class:`NetworkView`;
* a small compatible source is enumerated and materialised only after its
  geometry is compared with an independently clean target network.

For larger or unproven inputs the API raises an explicit normal-build error.
It never turns a partial candidate set into a certified group result.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from vorpy.src.calculations import calc_flat_vert, calc_vert, verify_site
from vorpy.src.network.network import Network
from vorpy.src.network.views import NetworkView, NetworkViewError


class IncrementalCompletionError(ValueError):
    """The requested target cannot be certified by the bounded API."""


EXACT_REUSE = "EXACT_REUSE"
VERIFIED_INCREMENTAL = "VERIFIED_INCREMENTAL"
NORMAL_BUILD_REQUIRED = "NORMAL_BUILD_REQUIRED"


@dataclass(frozen=True)
class CandidateVertex:
    """A bounded, independently validated four-generator candidate."""

    local_generators: tuple[int, ...]
    system_generators: tuple[int | None, ...]
    stable_generators: tuple[str, ...]
    location: tuple[float, float, float]
    radius: float
    doublet: int = 0


@dataclass(frozen=True)
class GeometryComparison:
    """Stable-identity comparison between two target networks."""

    equivalent: bool
    mismatches: tuple[str, ...] = ()


@dataclass(frozen=True)
class IncrementalGroupPlan:
    """Read-only plan and evidence for one requested group."""

    decision: str
    group_system_ids: tuple[int, ...]
    target_local_indices: tuple[int, ...]
    source_context_id: str
    source_network_id: str
    scheme: str
    solvent_mode: str | None
    solvent_generator_ids: tuple[str, ...]
    boundary_indices: tuple[int, ...]
    reused_vertex_ids: tuple[Any, ...]
    reused_edge_ids: tuple[Any, ...]
    reused_surface_ids: tuple[Any, ...]
    missing_candidates: tuple[CandidateVertex, ...]
    affected_feature_ids: tuple[Any, ...]
    invalidated_caches: tuple[str, ...]
    candidate_search_complete: bool
    reasons: tuple[str, ...] = ()
    view: NetworkView | None = None


@dataclass(frozen=True)
class IncrementalGroupResult:
    """Result of exact reuse or verified bounded materialisation."""

    plan: IncrementalGroupPlan
    network: Network | None = None
    view: NetworkView | None = None
    comparison: GeometryComparison | None = None


def _values(row: Any, field: str) -> tuple[Any, ...]:
    value = row.get(field, ())
    if value is None:
        return ()
    if isinstance(value, float) and np.isnan(value):
        return ()
    return tuple(value)


def _generator_maps(network: Any) -> tuple[tuple[str, ...], tuple[int | None, ...]]:
    stable = getattr(network, "generator_stable_ids", None)
    system = getattr(network, "generator_system_indices", None)
    balls = getattr(network, "balls", None)
    if balls is None:
        raise IncrementalCompletionError("Source network lacks a balls table.")
    if stable is None and "stable_id" in balls.columns:
        stable = tuple(str(value) for value in balls["stable_id"])
    if system is None and "system_num" in balls.columns:
        system = tuple(
            None if pd.isna(value) else int(value) for value in balls["system_num"]
        )
    if stable is None or system is None:
        raise IncrementalCompletionError(
            "Incremental completion requires stable generator and system identities."
        )
    stable = tuple(str(value) for value in stable)
    system = tuple(None if value is None else int(value) for value in system)
    if len(stable) != len(balls) or len(system) != len(balls):
        raise IncrementalCompletionError("Generator identity tables disagree with balls.")
    return stable, system


def _feature_complete(row: Any) -> bool:
    if "complete" not in row.index:
        return True
    value = row.get("complete")
    return value is not None and not (isinstance(value, float) and np.isnan(value)) and bool(value)


def _generator_ids(network: Any, local_indices: Iterable[int]) -> tuple[str, ...]:
    stable, _ = _generator_maps(network)
    return tuple(sorted(stable[int(index)] for index in local_indices))


def _vertex_id(network: Any, row: Any) -> tuple[Any, ...]:
    balls = _values(row, "balls")
    return (
        _generator_ids(network, balls),
        tuple(float(value) for value in row.get("loc", ())),
    )


def _edge_ids(network: Any, frame: Any) -> tuple[tuple[Any, ...], ...]:
    if frame is None:
        return ()
    vertices = getattr(network, "verts", None)
    vertex_ids = {
        int(position): _vertex_id(network, row)
        for position, row in vertices.iterrows()
    } if vertices is not None else {}
    records = []
    for _, row in frame.iterrows():
        missing = [int(index) for index in _values(row, "verts") if int(index) not in vertex_ids]
        if missing:
            raise IncrementalCompletionError(
                "Edge incidence references unavailable vertices: " + ", ".join(map(str, missing))
            )
        records.append((
            _generator_ids(network, _values(row, "balls")),
            tuple(sorted(vertex_ids[int(index)] for index in _values(row, "verts"))),
        ))
    return tuple(records)


def _surface_ids(network: Any, frame: Any, vertices: Any, edges: Any) -> tuple[tuple[Any, ...], ...]:
    if frame is None:
        return ()
    vertex_ids = {
        int(position): _vertex_id(network, row)
        for position, row in vertices.iterrows()
    } if vertices is not None else {}
    edge_ids = {
        int(position): edge_id
        for position, edge_id in enumerate(_edge_ids(network, edges))
    } if edges is not None else {}
    return tuple(
        (
            _generator_ids(network, _values(row, "balls")),
            tuple(sorted(vertex_ids[int(index)] for index in _values(row, "verts"))),
            tuple(sorted(edge_ids[int(index)] for index in _values(row, "edges"))),
        )
        for _, row in frame.iterrows()
    )


def _context(network: Any):
    try:
        return network.get_solve_context_key()
    except (AttributeError, TypeError, ValueError) as error:
        raise IncrementalCompletionError(
            f"Source solve context is unavailable: {error}"
        ) from error


def _compatible_context(left: Any, right: Any) -> bool:
    """Compare solve inputs while allowing the requested group to differ."""
    left_key = _context(left)
    right_key = _context(right)
    fields = (
        "schema_version", "scheme", "scope_mode", "network_scope", "frame_id",
        "partition_scheme", "boundary_mode", "boundary_config_digest",
        "boundary_generator_digest", "environment_digest", "generator_digest",
        "geometry_digest", "physical_generator_ids", "local_to_system",
        "boundary_generator_ids", "resolution", "vertex_settings_digest",
        "build_options_digest", "solver_settings_digest", "solver_version",
    )
    return all(getattr(left_key, field) == getattr(right_key, field) for field in fields)


def _target_indices(network: Any, group_system_ids: Sequence[int]) -> tuple[int, ...]:
    _, system = _generator_maps(network)
    lookup = {value: index for index, value in enumerate(system) if value is not None}
    requested = tuple(sorted({int(value) for value in group_system_ids}))
    if not requested:
        raise IncrementalCompletionError("Group C must contain at least one generator.")
    missing = [value for value in requested if value not in lookup]
    if missing:
        raise IncrementalCompletionError(
            "Requested generators are outside the solved universe: "
            + ", ".join(str(value) for value in missing)
        )
    return tuple(lookup[value] for value in requested)


def _inside_box(network: Any, location: Sequence[float], tolerance: float = 1e-9) -> bool:
    box = getattr(network, "box", None)
    if not isinstance(box, Mapping) or "verts" not in box:
        return True
    bounds = np.asarray(box["verts"], dtype=float)
    point = np.asarray(location, dtype=float)
    return bool(
        bounds.shape == (2, 3)
        and np.all(point >= bounds[0] - tolerance)
        and np.all(point <= bounds[1] + tolerance)
    )


def _valid_candidate(network: Any, combo: tuple[int, ...], location: Any, radius: Any) -> bool:
    if location is None or radius is None:
        return False
    try:
        location = np.asarray(location, dtype=float)
        radius = float(radius)
    except (TypeError, ValueError, OverflowError):
        return False
    if location.shape != (3,) or not np.all(np.isfinite(location)):
        return False
    balls = network.balls
    settings = getattr(network, "settings", {}) or {}
    max_vert = float(settings.get("max_vert", np.inf))
    scheme = str(settings.get("net_type", "aw")).lower()
    if scheme == "pow":
        # Power vertices store squared power distance; mirror find_site_pow.
        limit = max_vert ** 2 - min(
            float(balls.iloc[index]["rad"]) ** 2 for index in combo
        )
        oversized = radius >= limit
    else:
        limit = max_vert
        oversized = abs(radius) >= limit
    if not np.isfinite(radius) or oversized:
        return False
    if not _inside_box(network, location):
        return False

    other = [index for index in range(len(balls)) if index not in combo]
    if not other:
        return True
    test_locs = np.asarray([balls.iloc[index]["loc"] for index in other], dtype=float)
    test_rads = np.asarray([balls.iloc[index]["rad"] for index in other], dtype=float)
    try:
        return bool(verify_site(location, radius, test_locs, test_rads, net_type=scheme))
    except (TypeError, ValueError, FloatingPointError):
        return False


def discover_missing_vertex_candidates(
    source_network: Any,
    target_local_indices: Sequence[int],
    *,
    max_generators: int = 24,
    max_combinations: int = 100_000,
) -> tuple[tuple[CandidateVertex, ...], bool, tuple[str, ...]]:
    """Enumerate missing four-generator vertices for a bounded source.

    The enumeration is exhaustive only when the combination bound is met.  A
    false ``complete`` result is an explicit reason to use a normal build.
    """
    balls = getattr(source_network, "balls", None)
    verts = getattr(source_network, "verts", None)
    if balls is None or verts is None:
        return (), False, ("source vertex/ball tables are unavailable",)
    count = len(balls)
    if count > max_generators:
        return (), False, (f"generator bound exceeded: {count} > {max_generators}",)
    total = comb(count, 4) if count >= 4 else 0
    if total > max_combinations:
        return (), False, (
            f"candidate combination bound exceeded: {total} > {max_combinations}",
        )

    stable, system = _generator_maps(source_network)
    existing = {
        _vertex_id(source_network, row)
        for _, row in verts.iterrows()
        if _feature_complete(row)
    }
    candidates = []
    candidate_ids = set()
    locations = [balls.iloc[index]["loc"] for index in range(count)]
    radii = [balls.iloc[index]["rad"] for index in range(count)]
    scheme = str((getattr(source_network, "settings", {}) or {}).get("net_type", "aw")).lower()
    if scheme not in {"aw", "pow", "prm"}:
        return (), False, (f"unsupported scheme: {scheme}",)

    for combo in combinations(range(count), 4):
        combo = tuple(int(index) for index in combo)
        if not set(combo).intersection(target_local_indices):
            continue
        combo_locs = np.asarray([locations[index] for index in combo], dtype=float)
        combo_rads = np.asarray([radii[index] for index in combo], dtype=float)
        try:
            if scheme == "aw":
                primary, primary_radius, secondary, secondary_radius = calc_vert(
                    combo_locs,
                    combo_rads,
                )
                roots = ((primary, primary_radius, 2), (secondary, secondary_radius, 1))
            else:
                location, radius = calc_flat_vert(
                    combo_locs,
                    combo_rads,
                    power=scheme == "pow",
                )
                roots = ((location, radius, 0),)
        except (TypeError, ValueError, FloatingPointError, ZeroDivisionError):
            continue

        for location, radius, doublet in roots:
            if not _valid_candidate(source_network, combo, location, radius):
                continue
            row = {
                "balls": list(combo),
                "loc": list(np.asarray(location, dtype=float)),
                "rad": float(radius),
            }
            identity = _vertex_id(source_network, row)
            if identity in existing or identity in candidate_ids:
                continue
            candidate_ids.add(identity)
            candidates.append(CandidateVertex(
                local_generators=combo,
                system_generators=tuple(system[index] for index in combo),
                stable_generators=tuple(stable[index] for index in combo),
                location=tuple(float(value) for value in row["loc"]),
                radius=float(radius),
                doublet=doublet,
            ))

    return tuple(candidates), True, ()


def _reusable_vertex_rows(source_network: Any, target_local: set[int]) -> tuple[list[dict[str, Any]], tuple[Any, ...]]:
    verts = getattr(source_network, "verts", None)
    if verts is None:
        return [], ()
    rows = []
    identities = []
    for _, row in verts.iterrows():
        balls = tuple(int(value) for value in _values(row, "balls"))
        if not set(balls).intersection(target_local):
            continue
        if not _feature_complete(row):
            continue
        if len(balls) != 4 or len(set(balls)) != 4:
            continue
        data = row.to_dict()
        data["balls"] = list(balls)
        data["loc"] = list(np.asarray(data["loc"], dtype=float))
        data["rad"] = float(data["rad"])
        data.setdefault("dub", 0)
        rows.append(data)
        identities.append(_vertex_id(source_network, data))
    return rows, tuple(identities)


def _plan_ids(source_network: Any, target_local: set[int]):
    vertices = getattr(source_network, "verts", None)
    edges = getattr(source_network, "edges", None)
    surfs = getattr(source_network, "surfs", None)
    reused_vertices = tuple(
        _vertex_id(source_network, row)
        for _, row in vertices.iterrows()
        if set(_values(row, "balls")).intersection(target_local) and _feature_complete(row)
    ) if vertices is not None else ()
    reused_edges = ()
    if edges is not None and vertices is not None:
        relevant_edge_indices = {
            int(edge_index)
            for _, row in surfs.iterrows()
            if set(_values(row, "balls")).intersection(target_local)
            for edge_index in _values(row, "edges")
        } if surfs is not None else {
            int(index)
            for index, row in edges.iterrows()
            if set(_values(row, "balls")).intersection(target_local)
        }
        edge_ids = _edge_ids(source_network, edges)
        reused_edges = tuple(edge_ids[index] for index in sorted(relevant_edge_indices) if index < len(edge_ids))
    reused_surfaces = ()
    if surfs is not None:
        surface_ids = _surface_ids(source_network, surfs, vertices, edges)
        reused_surfaces = tuple(
            surface_ids[position]
            for position, (_, row) in enumerate(surfs.iterrows())
            if position < len(surface_ids)
            and set(_values(row, "balls")).intersection(target_local)
        )
    return reused_vertices, reused_edges, reused_surfaces


def assess_incremental_group(
    source_network: Any,
    group_system_ids: Sequence[int],
    *,
    max_generators: int = 24,
    max_combinations: int = 100_000,
) -> IncrementalGroupPlan:
    """Plan exact reuse or bounded completion without mutating the source."""
    context = _context(source_network)
    stable, system = _generator_maps(source_network)
    target_local = _target_indices(source_network, group_system_ids)
    target_local_set = set(target_local)
    requested = tuple(sorted({int(value) for value in group_system_ids}))
    source_id = str(getattr(source_network, "source_network_id", "") or "")
    scheme = str((getattr(source_network, "settings", {}) or {}).get("net_type", "aw")).lower()
    solvent_mode = getattr(source_network, "solvent_mode", None)
    solvent_mode = getattr(solvent_mode, "value", solvent_mode)

    try:
        view = NetworkView.for_group(source_network, requested)
    except NetworkViewError as error:
        view = None
        view_reason = str(error)
    else:
        view_reason = ""

    try:
        reused_vertices, reused_edges, reused_surfaces = _plan_ids(
            source_network, target_local_set
        )
        plan_feature_reason = ()
    except (KeyError, IncrementalCompletionError, TypeError, ValueError) as error:
        reused_vertices, reused_edges, reused_surfaces = (), (), ()
        plan_feature_reason = (f"source incidence is not reusable: {error}",)
    if view is not None:
        return IncrementalGroupPlan(
            decision=EXACT_REUSE,
            group_system_ids=requested,
            target_local_indices=target_local,
            source_context_id=context.context_id,
            source_network_id=source_id,
            scheme=scheme,
            solvent_mode=None if solvent_mode is None else str(solvent_mode),
            solvent_generator_ids=tuple(getattr(source_network, "solvent_generator_ids", ())),
            boundary_indices=tuple(getattr(source_network, "boundary_indices", ()) or ()),
            reused_vertex_ids=reused_vertices,
            reused_edge_ids=reused_edges,
            reused_surface_ids=reused_surfaces,
            missing_candidates=(),
            affected_feature_ids=(),
            invalidated_caches=(),
            candidate_search_complete=True,
            reasons=(),
            view=view,
        )

    candidates, search_complete, reasons = discover_missing_vertex_candidates(
        source_network,
        target_local,
        max_generators=max_generators,
        max_combinations=max_combinations,
    )
    reasons = tuple(dict.fromkeys(
        (view_reason,) + tuple(plan_feature_reason) + tuple(reasons)
    ))
    affected = tuple(
        ("vertex", candidate.local_generators, candidate.location)
        for candidate in candidates
    )
    reusable_rows, _ = _reusable_vertex_rows(source_network, target_local_set)
    if not search_complete:
        decision = NORMAL_BUILD_REQUIRED
        invalidated = ()
    else:
        decision = VERIFIED_INCREMENTAL if (candidates or reusable_rows) else NORMAL_BUILD_REQUIRED
        invalidated = (
            "group_results",
            "curvature",
            "topology",
            "alpha",
            "dual",
            "interface_views",
        ) if (candidates or reusable_rows) else ()
    return IncrementalGroupPlan(
        decision=decision,
        group_system_ids=requested,
        target_local_indices=target_local,
        source_context_id=context.context_id,
        source_network_id=source_id,
        scheme=scheme,
        solvent_mode=None if solvent_mode is None else str(solvent_mode),
        solvent_generator_ids=tuple(getattr(source_network, "solvent_generator_ids", ())),
        boundary_indices=tuple(getattr(source_network, "boundary_indices", ()) or ()),
        reused_vertex_ids=reused_vertices,
        reused_edge_ids=reused_edges,
        reused_surface_ids=reused_surfaces,
        missing_candidates=candidates,
        affected_feature_ids=affected,
        invalidated_caches=invalidated,
        candidate_search_complete=search_complete,
        reasons=reasons,
        view=None,
    )


def _target_network(source_network: Any, plan: IncrementalGroupPlan, rows: Sequence[Mapping[str, Any]]) -> Network:
    source_balls = source_network.balls
    settings = deepcopy(getattr(source_network, "settings", {}) or {})
    settings["build_type"] = "all"
    target = Network(
        locs=[value for value in source_balls["loc"]],
        rads=[float(value) for value in source_balls["rad"]],
        group=list(plan.target_local_indices),
        group_name="incremental:" + ",".join(str(value) for value in plan.group_system_ids),
        settings=settings,
        balls=source_balls.copy(deep=True),
        verts=pd.DataFrame(list(rows)),
        system=getattr(source_network, "sys", None),
        box=deepcopy(getattr(source_network, "box", None)),
        source_network_id=getattr(source_network, "source_network_id", None),
        generator_system_indices=getattr(source_network, "generator_system_indices", None),
        generator_stable_ids=getattr(source_network, "generator_stable_ids", None),
        solvent_local_indices=getattr(source_network, "solvent_local_indices", ()),
        solvent_generator_ids=getattr(source_network, "solvent_generator_ids", ()),
        solvent_mode=getattr(source_network, "solvent_mode", None),
        retained_system_indices=getattr(source_network, "retained_system_indices", ()),
        retained_generator_ids=getattr(source_network, "retained_generator_ids", ()),
    )
    target.boundary_indices = tuple(getattr(source_network, "boundary_indices", ()) or ())
    target.boundary_generators = deepcopy(getattr(source_network, "boundary_generators", ()))
    target.boundary_config = deepcopy(getattr(source_network, "boundary_config", None))
    target.boundary_mode = getattr(source_network, "boundary_mode", None)
    target.incremental_lineage = {
        "source_network_id": getattr(source_network, "source_network_id", None),
        "source_context_id": plan.source_context_id,
        "group_system_ids": plan.group_system_ids,
        "reused_feature_counts": {
            "vertices": len(plan.reused_vertex_ids),
            "edges": len(plan.reused_edge_ids),
            "surfaces": len(plan.reused_surface_ids),
        },
        "invalidated_caches": plan.invalidated_caches,
        "status": "UNVERIFIED_INCREMENTAL",
    }
    return target


def _compare_values(left: Any, right: Any, tolerance: float) -> bool:
    try:
        left = np.asarray(left, dtype=float)
        right = np.asarray(right, dtype=float)
        return left.shape == right.shape and bool(np.allclose(left, right, rtol=tolerance, atol=tolerance))
    except (TypeError, ValueError):
        return left == right


def compare_incremental_geometry(
    incremental_network: Network,
    clean_network: Network,
    *,
    tolerance: float = 1e-9,
) -> GeometryComparison:
    """Compare stable features and selected cell measurements."""
    mismatches = []
    if not _compatible_context(incremental_network, clean_network):
        mismatches.append("solve context")

    def compare_frame(kind: str, left: Any, right: Any, key_builder: Callable[[Any, Any], Any], fields: tuple[str, ...]):
        left_map = {key_builder(incremental_network, row): row for _, row in left.iterrows()}
        right_map = {key_builder(clean_network, row): row for _, row in right.iterrows()}
        if set(left_map) != set(right_map):
            mismatches.append(f"{kind} identities")
            return
        for key in left_map:
            for field in fields:
                if field in left_map[key].index and field in right_map[key].index:
                    if not _compare_values(left_map[key].get(field), right_map[key].get(field), tolerance):
                        mismatches.append(f"{kind} {field} {key!r}")

    if incremental_network.verts is None or clean_network.verts is None:
        mismatches.append("vertices unavailable")
    else:
        compare_frame("vertex", incremental_network.verts, clean_network.verts, _vertex_id, ("loc", "rad"))
    if incremental_network.edges is None or clean_network.edges is None:
        mismatches.append("edges unavailable")
    else:
        compare_frame("edge", incremental_network.edges, clean_network.edges,
                      lambda net, row: (
                          _generator_ids(net, _values(row, "balls")),
                          tuple(sorted(_vertex_id(net, net.verts.iloc[int(index)]) for index in _values(row, "verts"))),
                      ), ("length",))
    if incremental_network.surfs is None or clean_network.surfs is None:
        mismatches.append("surfaces unavailable")
    else:
        compare_frame("surface", incremental_network.surfs, clean_network.surfs,
                      lambda net, row: (
                          _generator_ids(net, _values(row, "balls")),
                          tuple(sorted(_vertex_id(net, net.verts.iloc[int(index)]) for index in _values(row, "verts"))),
                      ), ("sa",))

    selected = set(getattr(incremental_network, "group", ()) or ())
    for field in ("complete", "vol", "int_mean_curv", "int_gauss_curv"):
        if field not in incremental_network.balls.columns or field not in clean_network.balls.columns:
            continue
        for index in sorted(selected):
            if not _compare_values(
                incremental_network.balls.iloc[index].get(field),
                clean_network.balls.iloc[index].get(field),
                tolerance,
            ):
                mismatches.append(f"cell {field} {index}")

    return GeometryComparison(not mismatches, tuple(dict.fromkeys(mismatches)))


def complete_group_incrementally(
    source_network: Any,
    group_system_ids: Sequence[int],
    *,
    clean_reference: Network | None = None,
    max_generators: int = 24,
    max_combinations: int = 100_000,
    calculate_curvature: bool = True,
    normal_build: Callable[[IncrementalGroupPlan], Any] | None = None,
) -> IncrementalGroupResult:
    """Reuse or build a bounded Group C, never silently returning partial data.

    ``clean_reference`` is required for the materialisation path.  A caller
    that cannot provide an independent clean comparison receives an explicit
    normal-build request, or may supply ``normal_build`` as that explicit
    fallback.  The source network is never mutated.
    """
    plan = assess_incremental_group(
        source_network,
        group_system_ids,
        max_generators=max_generators,
        max_combinations=max_combinations,
    )
    if plan.decision == EXACT_REUSE:
        return IncrementalGroupResult(plan=plan, view=plan.view)

    if plan.decision == NORMAL_BUILD_REQUIRED or not plan.candidate_search_complete:
        if normal_build is not None:
            fallback = normal_build(plan)
            return IncrementalGroupResult(plan=plan, network=fallback)
        raise IncrementalCompletionError(
            "Exact incremental completion is unavailable; explicit normal build required. "
            + "; ".join(plan.reasons)
        )

    if clean_reference is None:
        if normal_build is not None:
            fallback = normal_build(plan)
            return IncrementalGroupResult(plan=plan, network=fallback)
        raise IncrementalCompletionError(
            "Incremental candidates are unverified; provide clean_reference or explicit normal build."
        )

    if not _compatible_context(source_network, clean_reference):
        raise IncrementalCompletionError(
            "Clean reference solve context does not match the source context."
        )

    reusable_rows, _ = _reusable_vertex_rows(source_network, set(plan.target_local_indices))
    candidate_rows = [
        {
            "balls": list(candidate.local_generators),
            "loc": list(candidate.location),
            "rad": candidate.radius,
            "dub": candidate.doublet,
        }
        for candidate in plan.missing_candidates
    ]
    target = _target_network(source_network, plan, reusable_rows + candidate_rows)
    target.build(calculate_curvature=calculate_curvature)
    comparison = compare_incremental_geometry(target, clean_reference)
    if not comparison.equivalent:
        raise IncrementalCompletionError(
            "Incremental geometry differs from clean reference: "
            + "; ".join(comparison.mismatches)
        )
    target.incremental_lineage = dict(target.incremental_lineage)
    target.incremental_lineage["status"] = "VERIFIED_INCREMENTAL"
    return IncrementalGroupResult(
        plan=plan,
        network=target,
        comparison=comparison,
    )


__all__ = [
    "CandidateVertex",
    "EXACT_REUSE",
    "GeometryComparison",
    "IncrementalCompletionError",
    "IncrementalGroupPlan",
    "IncrementalGroupResult",
    "NORMAL_BUILD_REQUIRED",
    "VERIFIED_INCREMENTAL",
    "assess_incremental_group",
    "compare_incremental_geometry",
    "complete_group_incrementally",
    "discover_missing_vertex_candidates",
]
