"""Birth candidates on existing VorPy AW primal features."""

from __future__ import annotations

from math import isfinite

import numpy as np

from vorpy.src.analyze.apollonius import _edge_alpha_minimum
from vorpy.src.calculations.edge_geometry import (
    AdditivelyWeightedTrisectorBranch,
    AdditivelyWeightedTrisectorConic,
    AffineEdgeGeometry,
    LineEdgeGeometry,
)
from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge


def _clearances(point, locations, radii):
    return np.linalg.norm(np.asarray(point, dtype=float) - locations, axis=1) - radii


def generator_birth(radius):
    """Formal union-of-balls birth, retained even when negative."""
    return -float(radius)


def _edge_candidate_points(network, edge_id, tolerance):
    """Return the exact minimum and diagnostics on one finite AW edge.

    Geometry is resolved by VorPy's existing AW edge resolver. The stationary
    candidates use its analytic branch/conic parameterization; equal-radius
    line edges use the exact projection onto the finite segment.
    """
    resolved = resolve_aw_network_edge(network, edge_id, tolerance=tolerance)
    geometry = resolved.geometry
    if isinstance(geometry, AffineEdgeGeometry):
        geometry = geometry.source
        lower, upper = sorted(
            (float(resolved.geometry.source_start), float(resolved.geometry.source_end))
        )
    else:
        lower, upper = float(geometry.t_min), float(geometry.t_max)
    if not isfinite(lower) or not isfinite(upper) or upper < lower:
        raise ValueError("Analytic AW edge has invalid finite parameter bounds.")

    candidates = [(lower, "endpoint"), (upper, "endpoint")]
    if isinstance(geometry, AdditivelyWeightedTrisectorBranch):
        # Its parameter is the common clearance rho, so only endpoints can
        # minimize on this monotone interval.
        pass
    elif isinstance(geometry, AdditivelyWeightedTrisectorConic):
        if geometry.conic_type == "ellipse":
            candidates.extend(
                (float(k * np.pi), "stationary_analytic")
                for k in range(
                    int(np.ceil(lower / np.pi)), int(np.floor(upper / np.pi)) + 1
                )
            )
        elif (
            geometry.conic_type == "hyperbola"
            and geometry.hyperbola_orientation == "rho"
        ):
            if lower <= 0.0 <= upper:
                candidates.append((0.0, "stationary_analytic"))
        elif geometry.conic_type == "parabola" and lower <= 0.0 <= upper:
            candidates.append((0.0, "stationary_analytic"))
    elif isinstance(geometry, LineEdgeGeometry):
        direction = geometry.end - geometry.start
        norm2 = float(direction @ direction)
        center = np.asarray(resolved.locations[0], dtype=float)
        projection = float((center - geometry.start) @ direction / norm2)
        parameter = min(1.0, max(0.0, projection))
        candidates = [
            (
                parameter,
                "endpoint" if parameter in (0.0, 1.0) else "stationary_analytic",
            )
        ]
    else:
        # Keep the existing analytic implementation as a supported fallback
        # for future geometry adapters; do not approximate unsupported types.
        value, diagnostic = _edge_alpha_minimum(network, edge_id, tolerance)
        return value, {
            "birth_source": "other",
            "residual_A": None,
            "tolerance_A": tolerance,
            "notes": "Used existing AW edge minimizer; selected point residual unavailable.",
            "diagnostics": diagnostic,
        }

    points = np.asarray(resolved.locations, dtype=float)
    radii = np.asarray(resolved.radii, dtype=float)
    evaluated = []
    for parameter, source in candidates:
        point = np.asarray(geometry.point(parameter), dtype=float)
        clearances = _clearances(point, points, radii)
        residual = float(np.max(clearances) - np.min(clearances))
        if residual > tolerance:
            raise ValueError(
                f"Edge {edge_id} equal-clearance residual {residual:.6g} exceeds {tolerance:g} A."
            )
        if isinstance(geometry, AdditivelyWeightedTrisectorBranch):
            expected = float(parameter)
        elif isinstance(geometry, AdditivelyWeightedTrisectorConic):
            expected = float(geometry.rho(parameter))
        else:
            expected = float(np.mean(clearances))
        consistency = float(np.max(np.abs(clearances - expected)))
        if consistency > tolerance:
            raise ValueError(
                f"Edge {edge_id} analytic clearance residual {consistency:.6g} exceeds {tolerance:g} A."
            )
        evaluated.append(
            {
                "alpha": expected,
                "parameter": float(parameter),
                "point": point,
                "source": source,
                "residual_A": max(residual, consistency),
            }
        )
    best = min(evaluated, key=lambda item: item["alpha"])
    return best["alpha"], {
        "birth_source": best["source"],
        "residual_A": best["residual_A"],
        "tolerance_A": tolerance,
        "parameter": best["parameter"],
        "parameter_bounds": (lower, upper),
        "candidate_count": len(evaluated),
        "candidate_clearances_A": tuple(item["alpha"] for item in evaluated),
        "minimizer_xyz": tuple(map(float, best["point"])),
        "edge_resolution_status": resolved.status,
    }


def edge_birth(network, edge_id, tolerance=1e-6):
    return _edge_candidate_points(network, edge_id, tolerance)


def _pair_contact_candidate(network, pair, tolerance):
    balls = network.balls
    p = np.asarray([balls.loc[index, "loc"] for index in pair], dtype=float)
    r = np.asarray([float(balls.loc[index, "rad"]) for index in pair], dtype=float)
    delta = p[1] - p[0]
    distance = float(np.linalg.norm(delta))
    if distance <= tolerance:
        return None, {"status": "coincident_pair_centers"}
    if distance < abs(float(r[0] - r[1])) - tolerance:
        return None, {"status": "no_pairwise_AW_bisector"}
    a = (distance + r[0] - r[1]) / 2.0
    b = (distance - r[0] + r[1]) / 2.0
    if a < -tolerance or b < -tolerance:
        return None, {"status": "no_nonnegative_offset_contact"}
    point = p[0] + (a / distance) * delta
    values = _clearances(point, p, r)
    residual = float(abs(values[0] - values[1]))
    return {"point": point, "alpha": float(np.mean(values)), "residual_A": residual}, {
        "status": "analytic_pair_contact",
        "center_distance_A": distance,
        "offset_radii_A": (a, b),
        "pair_clearances_A": tuple(map(float, values)),
        "residual_A": residual,
    }


def _vertex_candidate(network, vertex_id, pair, tolerance):
    row = network.verts.loc[vertex_id]
    point = np.asarray(row["loc"], dtype=float)
    radii = np.asarray([network.balls.loc[i, "rad"] for i in pair], dtype=float)
    locations = np.asarray([network.balls.loc[i, "loc"] for i in pair], dtype=float)
    clearances = _clearances(point, locations, radii)
    residual = float(np.ptp(clearances))
    if residual > tolerance:
        raise ValueError(
            f"AW surface boundary vertex {vertex_id} pair residual {residual:g} A exceeds tolerance."
        )
    all_locations = np.asarray(
        [row["loc"] for _, row in network.balls.iterrows()], dtype=float
    )
    all_radii = np.asarray(
        [row["rad"] for _, row in network.balls.iterrows()], dtype=float
    )
    all_clearances = _clearances(point, all_locations, all_radii)
    pair_alpha = float(np.mean(clearances))
    if float(np.min(all_clearances)) < pair_alpha - tolerance:
        raise ValueError(
            f"AW surface boundary vertex {vertex_id} is occluded by another generator."
        )
    return pair_alpha, {
        "birth_source": "boundary_vertex",
        "feature_id": int(vertex_id),
        "residual_A": residual,
        "tolerance_A": tolerance,
        "candidate_xyz": tuple(map(float, point)),
    }


def surface_birth(network, surface_simplex, tolerance=1e-6, *,
                  edge_birth_cache=None, vertex_birth_cache=None):
    """Conservative minimum on the actual network-clipped pair surface.

    The unrestricted pair bisector has one analytic global minimum at the
    external-offset contact point. It is accepted only when that exact point
    satisfies all AW cell inequalities. Otherwise the closed patch minimum
    is sought on its recorded AW edge/vertex boundary. Missing or unresolved
    boundary geometry leaves the pair birth unresolved.
    """
    pair = surface_simplex.generator_ids
    interior, pair_diagnostic = _pair_contact_candidate(network, pair, tolerance)
    if interior is not None and interior["residual_A"] <= tolerance:
        point = interior["point"]
        all_locations = np.asarray(
            [row["loc"] for _, row in network.balls.iterrows()], dtype=float
        )
        all_radii = np.asarray(
            [row["rad"] for _, row in network.balls.iterrows()], dtype=float
        )
        all_clearances = _clearances(point, all_locations, all_radii)
        alpha = interior["alpha"]
        if float(np.min(all_clearances)) >= alpha - tolerance:
            tie_count = int(np.sum(np.abs(all_clearances - alpha) <= tolerance))
            if tie_count == 2:
                return alpha, {
                    "birth_source": "interior_analytic",
                    "residual_A": interior["residual_A"],
                    "tolerance_A": tolerance,
                    "minimizer_xyz": tuple(map(float, point)),
                    "candidate_status": pair_diagnostic["status"],
                    "notes": "Analytic global minimum on pair AW bisector lies in the actual AW surface patch.",
                }

    edge_ids, vertex_ids = set(), set()
    for feature in surface_simplex.primal_features:
        if not (feature.bounded and feature.complete):
            return None, {
                "birth_source": "unresolved",
                "residual_A": None,
                "reason": f"Surface feature {feature.key} is unbounded or incomplete.",
            }
        row = network.surfs.loc[feature.feature_id]
        edge_ids.update(
            int(i) for i in row.get("edges", ()) if i in network.edges.index
        )
        vertex_ids.update(
            int(i) for i in row.get("verts", ()) if i in network.verts.index
        )
    candidates = []
    failures = []
    for edge_id in sorted(edge_ids):
        edge_row = network.edges.loc[edge_id]
        edge_ids_for_feature = tuple(sorted(int(i) for i in edge_row.get("balls", ())))
        if not set(pair).issubset(edge_ids_for_feature):
            failures.append(f"Boundary edge {edge_id} does not contain pair {pair}.")
            continue
        cached = edge_birth_cache.get(int(edge_id)) if edge_birth_cache is not None else None
        if cached is not None:
            alpha, diagnostic, error_text = cached
            if error_text:
                failures.append(f"Boundary edge {edge_id}: {error_text}")
                continue
        else:
            try:
                alpha, diagnostic = edge_birth(network, edge_id, tolerance)
                error_text = None
            except (
                TypeError,
                ValueError,
                KeyError,
                IndexError,
                np.linalg.LinAlgError,
            ) as error:
                alpha, diagnostic, error_text = None, None, str(error)
            if edge_birth_cache is not None:
                edge_birth_cache[int(edge_id)] = (alpha, diagnostic, error_text)
            if error_text:
                failures.append(f"Boundary edge {edge_id}: {error_text}")
                continue
        candidates.append(
            (
                alpha,
                {**diagnostic, "birth_source": "boundary_edge", "feature_id": edge_id},
            )
        )
    for vertex_id in sorted(vertex_ids):
        vertex_row = network.verts.loc[vertex_id]
        if not set(pair).issubset({int(i) for i in vertex_row.get("balls", ())}):
            failures.append(
                f"Boundary vertex {vertex_id} does not contain pair {pair}."
            )
            continue
        cache_key = (int(vertex_id), tuple(sorted(pair)))
        cached = vertex_birth_cache.get(cache_key) if vertex_birth_cache is not None else None
        if cached is not None:
            alpha, diagnostic, error_text = cached
            if error_text:
                failures.append(error_text)
                continue
        else:
            try:
                alpha, diagnostic = _vertex_candidate(network, vertex_id, pair, tolerance)
                error_text = None
            except (TypeError, ValueError, KeyError, IndexError) as error:
                alpha, diagnostic, error_text = None, None, str(error)
            if vertex_birth_cache is not None:
                vertex_birth_cache[cache_key] = (alpha, diagnostic, error_text)
            if error_text:
                failures.append(error_text)
                continue
        candidates.append((alpha, diagnostic))
    if failures:
        return None, {
            "birth_source": "unresolved",
            "residual_A": None,
            "reason": "; ".join(failures[:5]),
            "boundary_candidate_count": len(candidates),
            "pair_candidate": pair_diagnostic,
        }
    if not candidates:
        return None, {
            "birth_source": "unresolved",
            "residual_A": None,
            "reason": "No certified analytic interior minimum or usable actual surface boundary.",
            "pair_candidate": pair_diagnostic,
        }
    alpha, diagnostic = min(candidates, key=lambda item: item[0])
    return float(alpha), {
        **diagnostic,
        "pair_candidate": pair_diagnostic,
        "boundary_candidate_count": len(candidates),
        "notes": "Minimum certified from the recorded closed AW surface boundary.",
    }


def vertex_birth(network, vertex_id, generator_ids, tolerance=1e-6):
    point = np.asarray(network.verts.loc[vertex_id, "loc"], dtype=float)
    locations = np.asarray(
        [network.balls.loc[i, "loc"] for i in generator_ids], dtype=float
    )
    radii = np.asarray(
        [network.balls.loc[i, "rad"] for i in generator_ids], dtype=float
    )
    clearances = _clearances(point, locations, radii)
    residual = float(np.max(clearances) - np.min(clearances))
    if residual > tolerance:
        raise ValueError(
            f"AW vertex {vertex_id} clearance residual {residual:.6g} exceeds {tolerance:g} A."
        )
    if len(generator_ids) != 4:
        raise ValueError("AW alpha vertex birth requires exactly four generators.")
    return float(np.mean(clearances)), {
        "birth_source": "apollonius_vertex",
        "feature_id": int(vertex_id),
        "individual_clearances_A": tuple(map(float, clearances)),
        "residual_A": residual,
        "tolerance_A": tolerance,
    }
