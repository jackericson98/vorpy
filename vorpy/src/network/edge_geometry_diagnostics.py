"""Read-only diagnostics for analytic geometry on constructed network edges."""

from collections import Counter
from dataclasses import dataclass

import numpy as np

from vorpy.src.calculations.edge_geometry import (
    line_geometry_from_endpoints,
    match_aw_trisector_conic_to_samples,
    match_aw_trisector_to_endpoints,
    validate_edge_samples,
)


@dataclass(frozen=True)
class AWNetworkEdgeRecord:
    """Diagnostic result for one stored network edge."""

    edge_index: object
    status: str
    reason: str
    ball_indices: tuple
    vertex_indices: tuple
    sample_count: int
    endpoint_error: float
    clearance_error: float
    maximum_sample_error: float
    rms_sample_error: float
    sampled_length: float
    analytic_length: float
    length_difference: float


@dataclass(frozen=True)
class AWNetworkEdgeDiagnostics:
    """Aggregate and per-edge results from a read-only AW geometry audit."""

    records: tuple

    @property
    def total_edges(self):
        return len(self.records)

    @property
    def status_counts(self):
        return dict(Counter(record.status for record in self.records))

    @property
    def matched_edges(self):
        return sum(record.status.startswith("matched") for record in self.records)

    @property
    def matched_fraction(self):
        if not self.records:
            return 0.0
        return self.matched_edges/len(self.records)

    @property
    def maximum_sample_error(self):
        values = [
            record.maximum_sample_error for record in self.records
            if np.isfinite(record.maximum_sample_error)
        ]
        return max(values, default=float("nan"))

    @property
    def rms_sample_error(self):
        values = [
            record.rms_sample_error for record in self.records
            if np.isfinite(record.rms_sample_error)
        ]
        return (
            float(np.sqrt(np.mean(np.square(values))))
            if values else float("nan")
        )

    @property
    def maximum_absolute_length_difference(self):
        values = [
            abs(record.length_difference) for record in self.records
            if np.isfinite(record.length_difference)
        ]
        return max(values, default=float("nan"))

    def summary(self):
        """Return aggregate values suitable for printing or later logging."""
        return {
            "total_edges": self.total_edges,
            "matched_edges": self.matched_edges,
            "matched_fraction": self.matched_fraction,
            "status_counts": self.status_counts,
            "maximum_sample_error": self.maximum_sample_error,
            "rms_sample_error": self.rms_sample_error,
            "maximum_absolute_length_difference":
                self.maximum_absolute_length_difference,
        }


def _polyline_length(points):
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return 0.0
    return float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=1)))


def _line_sample_errors(geometry, points):
    points = np.asarray(points, dtype=float)
    direction = geometry.end - geometry.start
    denominator = float(direction @ direction)
    parameters = ((points - geometry.start) @ direction)/denominator
    reconstructions = geometry.start + parameters[:, None]*direction
    errors = np.linalg.norm(points - reconstructions, axis=1)
    return float(np.max(errors)), float(np.sqrt(np.mean(errors**2)))


def _failure_record(edge_index, status, reason, balls, vertices,
                    sample_count=0):
    nan = float("nan")
    return AWNetworkEdgeRecord(
        edge_index, status, reason, balls, vertices, sample_count,
        nan, nan, nan, nan, nan, nan, nan,
    )

class AWNetworkEdgeResolutionError(ValueError):
    """Structured failure while resolving one stored AW edge analytically."""

    def __init__(
            self,
            edge_index,
            status,
            reason,
            ball_indices=(),
            vertex_indices=(),
            sample_count=0):
        super().__init__(
            f"Unable to resolve AW edge {edge_index}: "
            f"{status}: {reason}"
        )
        self.edge_index = edge_index
        self.status = status
        self.reason = reason
        self.ball_indices = tuple(ball_indices)
        self.vertex_indices = tuple(vertex_indices)
        self.sample_count = int(sample_count)


@dataclass(frozen=True)
class AWResolvedNetworkEdge:
    """Exact analytic representation of one stored regular AW edge."""

    edge_index: object
    geometry: object
    ball_indices: tuple
    vertex_indices: tuple
    locations: np.ndarray
    radii: np.ndarray
    status: str
    reason: str
    endpoint_error: float
    clearance_error: float


def resolve_aw_network_edge(net, edge_index, tolerance=1e-6):
    """Resolve one stored regular AW edge to exact analytic geometry.

    This is the canonical network-edge resolver. Diagnostics and production
    geometry calculations should both use this function so that they cannot
    silently diverge.

    Parameters
    ----------
    net
        Constructed VorPy Network.
    edge_index
        Index of the edge row in ``net.edges``.
    tolerance : float, default=1e-6
        Numerical tolerance used by the AW analytic matchers.

    Returns
    -------
    AWResolvedNetworkEdge
        Exact line, regular AW branch, or nonsingular AW conic geometry.

    Raises
    ------
    AWNetworkEdgeResolutionError
        If the stored edge cannot be resolved exactly.
    """
    if getattr(net, "edges", None) is None:
        raise ValueError("The network has no edge table.")
    if getattr(net, "verts", None) is None:
        raise ValueError("The network has no vertex table.")
    if getattr(net, "balls", None) is None:
        raise ValueError("The network has no ball table.")

    settings = getattr(net, "settings", None) or {}
    if settings.get("net_type", "aw") != "aw":
        raise ValueError(
            "Analytic AW edge resolution requires an AW network."
        )

    edge = net.edges.loc[edge_index]

    balls = tuple(
        int(value)
        for value in edge.get("balls", ())
    )
    vertices = tuple(
        int(value)
        for value in edge.get("verts", ())
    )

    raw_points = edge.get("points", ())
    points = np.asarray(raw_points, dtype=float)
    sample_count = len(points) if points.ndim == 2 else 0

    if len(balls) != 3 or len(vertices) != 2:
        raise AWNetworkEdgeResolutionError(
            edge_index=edge_index,
            status="unsupported_topology",
            reason=(
                "An analytic regular edge requires "
                "three balls and two vertices."
            ),
            ball_indices=balls,
            vertex_indices=vertices,
            sample_count=sample_count,
        )

    if (
        points.ndim != 2
        or points.shape[1:] != (3,)
        or len(points) == 0
    ):
        raise AWNetworkEdgeResolutionError(
            edge_index=edge_index,
            status="missing_samples",
            reason="The stored edge has no valid sampled points.",
            ball_indices=balls,
            vertex_indices=vertices,
            sample_count=sample_count,
        )

    locations = np.asarray(
        [
            net.balls.loc[index, "loc"]
            for index in balls
        ],
        dtype=float,
    )
    radii = np.asarray(
        [
            net.balls.loc[index, "rad"]
            for index in balls
        ],
        dtype=float,
    )
    vertex_locations = [
        np.asarray(
            net.verts.loc[index, "loc"],
            dtype=float,
        )
        for index in vertices
    ]

    # Equal-radius AW generators give an exact straight edge.
    if np.max(radii) - np.min(radii) <= tolerance:
        geometry = line_geometry_from_endpoints(
            vertices,
            vertex_locations,
        )

        return AWResolvedNetworkEdge(
            edge_index=edge_index,
            geometry=geometry,
            ball_indices=balls,
            vertex_indices=vertices,
            locations=locations,
            radii=radii,
            status="matched_line",
            reason=(
                "Equal-radius AW generators define a straight edge."
            ),
            endpoint_error=0.0,
            clearance_error=0.0,
        )

    # First attempt: regular rho/branch representation.
    try:
        match = match_aw_trisector_to_endpoints(
            locations,
            radii,
            vertices,
            vertex_locations,
            tolerance=tolerance,
        )
    except (TypeError, ValueError, np.linalg.LinAlgError) as error:
        raise AWNetworkEdgeResolutionError(
            edge_index=edge_index,
            status="invalid",
            reason=str(error),
            ball_indices=balls,
            vertex_indices=vertices,
            sample_count=sample_count,
        ) from error

    # Second attempt: nonsingular conic representation for intervals
    # that cross the rho turning point or branch boundary.
    if (
        match.geometry is None
        and match.status in {
            "cross_branch",
            "turning_point",
        }
    ):
        try:
            match = match_aw_trisector_conic_to_samples(
                locations,
                radii,
                vertices,
                vertex_locations,
                points,
                tolerance=tolerance,
            )
        except (TypeError, ValueError, np.linalg.LinAlgError) as error:
            raise AWNetworkEdgeResolutionError(
                edge_index=edge_index,
                status="nonsingular_match_failed",
                reason=str(error),
                ball_indices=balls,
                vertex_indices=vertices,
                sample_count=sample_count,
            ) from error

    if match.geometry is None:
        raise AWNetworkEdgeResolutionError(
            edge_index=edge_index,
            status=match.status,
            reason=match.reason,
            ball_indices=balls,
            vertex_indices=vertices,
            sample_count=sample_count,
        )

    status = (
        "matched_nonsingular"
        if match.status == "matched_nonsingular"
        else "matched_curve"
    )

    return AWResolvedNetworkEdge(
        edge_index=edge_index,
        geometry=match.geometry,
        ball_indices=balls,
        vertex_indices=vertices,
        locations=locations,
        radii=radii,
        status=status,
        reason=match.reason,
        endpoint_error=max(
            match.endpoint_errors,
            default=0.0,
        ),
        clearance_error=match.clearance_error,
    )


def diagnose_aw_edge_geometry(
        net,
        tolerance=1e-6,
        quadrature_order=32):
    """Audit analytic geometry against every edge in a constructed AW network.

    This diagnostic uses the same canonical analytic-edge resolver intended
    for production edge-curvature calculations.

    The function remains read-only. It neither replaces sampled edge points
    nor adds fields to ``net.edges``.
    """
    if getattr(net, "edges", None) is None:
        raise ValueError("The network has no edge table.")
    if getattr(net, "verts", None) is None:
        raise ValueError("The network has no vertex table.")
    if getattr(net, "balls", None) is None:
        raise ValueError("The network has no ball table.")

    settings = getattr(net, "settings", None) or {}
    if settings.get("net_type", "aw") != "aw":
        raise ValueError(
            "Analytic AW diagnostics require an AW network."
        )

    records = []

    for edge_index, edge in net.edges.iterrows():
        balls = tuple(
            int(value)
            for value in edge.get("balls", ())
        )
        vertices = tuple(
            int(value)
            for value in edge.get("verts", ())
        )

        raw_points = edge.get("points", ())
        points = np.asarray(raw_points, dtype=float)
        sample_count = len(points) if points.ndim == 2 else 0

        # Keep sampled-length behavior identical to the old diagnostic.
        sampled_length = (
            _polyline_length(points)
            if (
                points.ndim == 2
                and points.shape[1:] == (3,)
                and len(points) > 0
            )
            else float("nan")
        )

        try:
            resolved = resolve_aw_network_edge(
                net,
                edge_index,
                tolerance=tolerance,
            )
        except AWNetworkEdgeResolutionError as error:
            records.append(
                _failure_record(
                    edge_index=error.edge_index,
                    status=error.status,
                    reason=error.reason,
                    balls=error.ball_indices,
                    vertices=error.vertex_indices,
                    sample_count=error.sample_count,
                )
            )
            continue

        geometry = resolved.geometry

        # Straight line validation retains the historical orthogonal-distance
        # error calculation.
        if resolved.status == "matched_line":
            try:
                maximum_error, rms_error = _line_sample_errors(
                    geometry,
                    points,
                )
                analytic_length = geometry.arc_length(
                    order=quadrature_order
                )
            except (
                TypeError,
                ValueError,
                np.linalg.LinAlgError,
            ) as error:
                records.append(
                    _failure_record(
                        edge_index,
                        "sample_validation_failed",
                        str(error),
                        balls,
                        vertices,
                        sample_count,
                    )
                )
                continue

            records.append(
                AWNetworkEdgeRecord(
                    edge_index=edge_index,
                    status=resolved.status,
                    reason=resolved.reason,
                    ball_indices=resolved.ball_indices,
                    vertex_indices=resolved.vertex_indices,
                    sample_count=sample_count,
                    endpoint_error=0.0,
                    clearance_error=0.0,
                    maximum_sample_error=maximum_error,
                    rms_sample_error=rms_error,
                    sampled_length=sampled_length,
                    analytic_length=analytic_length,
                    length_difference=(
                        analytic_length - sampled_length
                    ),
                )
            )
            continue

        # Curved AW edges use the existing sample-validation machinery.
        try:
            validation = validate_edge_samples(
                geometry,
                points,
            )
            analytic_length = geometry.arc_length(
                order=quadrature_order
            )
        except (
            TypeError,
            ValueError,
            np.linalg.LinAlgError,
        ) as error:
            records.append(
                _failure_record(
                    edge_index,
                    "sample_validation_failed",
                    str(error),
                    balls,
                    vertices,
                    sample_count,
                )
            )
            continue

        records.append(
            AWNetworkEdgeRecord(
                edge_index=edge_index,
                status=resolved.status,
                reason=resolved.reason,
                ball_indices=resolved.ball_indices,
                vertex_indices=resolved.vertex_indices,
                sample_count=sample_count,
                endpoint_error=resolved.endpoint_error,
                clearance_error=resolved.clearance_error,
                maximum_sample_error=validation.maximum_error,
                rms_sample_error=validation.rms_error,
                sampled_length=sampled_length,
                analytic_length=analytic_length,
                length_difference=(
                    analytic_length - sampled_length
                ),
            )
        )

    return AWNetworkEdgeDiagnostics(tuple(records))
