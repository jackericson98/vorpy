"""Read-only diagnostics for analytic geometry on constructed network edges."""

from collections import Counter
from dataclasses import dataclass

import numpy as np

from vorpy.src.calculations.edge_geometry import (
    line_geometry_from_endpoints,
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


def diagnose_aw_edge_geometry(net, tolerance=1e-6, quadrature_order=32):
    """Audit analytic geometry against every edge in a constructed AW network.

    The function is deliberately read-only. It neither replaces sampled edge
    points nor adds fields to ``net.edges``.
    """
    if getattr(net, "edges", None) is None:
        raise ValueError("The network has no edge table.")
    if getattr(net, "verts", None) is None:
        raise ValueError("The network has no vertex table.")
    if getattr(net, "balls", None) is None:
        raise ValueError("The network has no ball table.")
    settings = getattr(net, "settings", None) or {}
    if settings.get("net_type", "aw") != "aw":
        raise ValueError("Analytic AW diagnostics require an AW network.")

    records = []
    for edge_index, edge in net.edges.iterrows():
        balls = tuple(int(value) for value in edge.get("balls", ()))
        vertices = tuple(int(value) for value in edge.get("verts", ()))
        raw_points = edge.get("points", ())
        points = np.asarray(raw_points, dtype=float)
        sample_count = len(points) if points.ndim == 2 else 0

        if len(balls) != 3 or len(vertices) != 2:
            records.append(_failure_record(
                edge_index, "unsupported_topology",
                "An analytic regular edge requires three balls and two vertices.",
                balls, vertices, sample_count,
            ))
            continue
        if points.ndim != 2 or points.shape[1:] != (3,) or len(points) == 0:
            records.append(_failure_record(
                edge_index, "missing_samples",
                "The stored edge has no valid sampled points.",
                balls, vertices, sample_count,
            ))
            continue

        locations = np.asarray([net.balls.loc[index, "loc"] for index in balls])
        radii = np.asarray([net.balls.loc[index, "rad"] for index in balls])
        vertex_locations = [
            np.asarray(net.verts.loc[index, "loc"]) for index in vertices
        ]
        sampled_length = _polyline_length(points)

        if np.max(radii) - np.min(radii) <= tolerance:
            geometry = line_geometry_from_endpoints(vertices, vertex_locations)
            maximum_error, rms_error = _line_sample_errors(geometry, points)
            analytic_length = geometry.arc_length()
            records.append(AWNetworkEdgeRecord(
                edge_index, "matched_line",
                "Equal-radius AW generators define a straight edge.",
                balls, vertices, sample_count, 0.0, 0.0,
                maximum_error, rms_error, sampled_length, analytic_length,
                analytic_length - sampled_length,
            ))
            continue

        try:
            match = match_aw_trisector_to_endpoints(
                locations, radii, vertices, vertex_locations,
                tolerance=tolerance,
            )
        except (TypeError, ValueError, np.linalg.LinAlgError) as error:
            records.append(_failure_record(
                edge_index, "invalid", str(error), balls, vertices,
                sample_count,
            ))
            continue
        if match.geometry is None:
            records.append(_failure_record(
                edge_index, match.status, match.reason, balls, vertices,
                sample_count,
            ))
            continue

        try:
            validation = validate_edge_samples(match.geometry, points)
            analytic_length = match.geometry.arc_length(order=quadrature_order)
        except (TypeError, ValueError, np.linalg.LinAlgError) as error:
            records.append(_failure_record(
                edge_index, "sample_validation_failed", str(error),
                balls, vertices, sample_count,
            ))
            continue
        records.append(AWNetworkEdgeRecord(
            edge_index, "matched_curve", match.reason, balls, vertices,
            sample_count, max(match.endpoint_errors, default=0.0),
            match.clearance_error, validation.maximum_error,
            validation.rms_error, sampled_length, analytic_length,
            analytic_length - sampled_length,
        ))

    return AWNetworkEdgeDiagnostics(tuple(records))
