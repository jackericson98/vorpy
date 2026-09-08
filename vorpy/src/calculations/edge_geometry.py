"""Geometry interfaces and normalized construction data for VorPy edges.

The interface is independent of network topology and storage.  It lets later
boundary code integrate intrinsic curve descriptors without depending on the
legacy sampled ``edge.points`` representation.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


def _optional_vector(value, label):
    if value is None:
        return None
    vector = np.asarray(value, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{label} must be a finite three-vector or None.")
    return vector


def _optional_float(value, label):
    if value is None:
        return None
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"{label} must be finite or None.")
    return number


@dataclass(frozen=True)
class EdgeConstructionMetadata:
    """Stable view of the legacy dictionary stored in ``net.edges['vals']``."""

    case: str | None
    center: np.ndarray | None
    radius: float | None
    secondary_center: np.ndarray | None
    secondary_radius: float | None
    midpoint: np.ndarray | None
    plane_normal: np.ndarray | None
    projection_direction: np.ndarray | None
    projection_point: np.ndarray | None
    outside: bool


@dataclass(frozen=True)
class AWEdgeGeometryMatch:
    """Result of matching stored edge endpoints to one analytic AW branch."""

    status: str
    reason: str
    geometry: object | None
    endpoint_errors: tuple
    clearance_error: float


@dataclass(frozen=True)
class EdgeSampleValidation:
    """Geometric residuals between stored samples and an analytic edge."""

    count: int
    maximum_error: float
    rms_error: float


def normalize_edge_construction(values):
    """Normalize old, current, straight, and missing edge metadata dictionaries.

    The legacy dictionary remains the serialized source of truth.  This helper
    provides explicit names and validation without changing that storage format.
    """
    values = {} if values is None else values
    if not hasattr(values, "get"):
        raise TypeError("Edge construction metadata must be a mapping or None.")
    case = values.get("case")
    return EdgeConstructionMetadata(
        case=None if case is None else str(case),
        center=_optional_vector(values.get("loc"), "loc"),
        radius=_optional_float(values.get("rad"), "rad"),
        secondary_center=_optional_vector(values.get("loc2"), "loc2"),
        secondary_radius=_optional_float(values.get("rad2"), "rad2"),
        midpoint=_optional_vector(values.get("vmid"), "vmid"),
        plane_normal=_optional_vector(values.get("pnorm"), "pnorm"),
        projection_direction=_optional_vector(values.get("dnorm"), "dnorm"),
        projection_point=_optional_vector(values.get("pa"), "pa"),
        outside=bool(values.get("outside", False)),
    )


def canonical_edge_endpoints(vertex_indices, vertex_locations):
    """Return endpoints ordered by stable vertex index.

    Geometry construction may discover the same edge in either direction.  A
    stable index order makes tangent-dependent signed quantities repeatable.
    """
    indices = tuple(int(index) for index in vertex_indices)
    locations = tuple(np.asarray(location, dtype=float) for location in vertex_locations)
    if len(indices) != 2 or len(locations) != 2:
        raise ValueError("An edge must have exactly two endpoint vertices.")
    if indices[0] == indices[1]:
        raise ValueError("Edge endpoint vertex indices must be distinct.")
    if any(location.shape != (3,) or not np.all(np.isfinite(location))
           for location in locations):
        raise ValueError("Edge endpoint locations must be finite three-vectors.")
    order = np.argsort(indices)
    return (
        (indices[int(order[0])], locations[int(order[0])]),
        (indices[int(order[1])], locations[int(order[1])]),
    )


def line_geometry_from_endpoints(vertex_indices, vertex_locations):
    """Construct exact line geometry using canonical endpoint orientation."""
    endpoints = canonical_edge_endpoints(vertex_indices, vertex_locations)
    return LineEdgeGeometry(endpoints[0][1], endpoints[1][1])


class EdgeGeometry(ABC):
    """Parametric edge curve ``r(t)`` on a closed parameter interval."""

    def __init__(self, t_min, t_max):
        self.t_min = float(t_min)
        self.t_max = float(t_max)
        if not np.isfinite(self.t_min) or not np.isfinite(self.t_max):
            raise ValueError("Edge parameter bounds must be finite.")
        if self.t_max <= self.t_min:
            raise ValueError("Edge parameter bounds must satisfy t_max > t_min.")

    @abstractmethod
    def point(self, t):
        """Return ``r(t)``."""

    @abstractmethod
    def tangent(self, t):
        """Return ``r'(t)`` (not necessarily unit length)."""

    @abstractmethod
    def second_derivative(self, t):
        """Return ``r''(t)``."""

    def speed(self, t):
        return float(np.linalg.norm(self.tangent(t)))

    def unit_tangent(self, t, tol=1e-12):
        tangent = np.asarray(self.tangent(t), dtype=float)
        speed = float(np.linalg.norm(tangent))
        if not np.isfinite(speed) or speed < tol:
            raise ValueError("Edge tangent is undefined at this parameter.")
        return tangent/speed

    def curvature(self, t, tol=1e-12):
        first = np.asarray(self.tangent(t), dtype=float)
        second = np.asarray(self.second_derivative(t), dtype=float)
        speed = float(np.linalg.norm(first))
        if not np.isfinite(speed) or speed < tol:
            raise ValueError("Edge curvature is undefined at this parameter.")
        value = float(np.linalg.norm(np.cross(first, second)))/(speed**3)
        if not np.isfinite(value):
            raise ValueError("Edge curvature is non-finite at this parameter.")
        return value

    def integrate(self, integrand, t0=None, t1=None, order=32):
        """Integrate a scalar parameter-space function by Gauss--Legendre."""
        lower = self.t_min if t0 is None else float(t0)
        upper = self.t_max if t1 is None else float(t1)
        if lower < self.t_min or upper > self.t_max or upper < lower:
            raise ValueError("Integration bounds fall outside the edge interval.")
        if upper == lower:
            return 0.0
        nodes, weights = np.polynomial.legendre.leggauss(int(order))
        half = 0.5*(upper - lower)
        center = 0.5*(upper + lower)
        return float(half*sum(
            weight*float(integrand(center + half*node))
            for node, weight in zip(nodes, weights)
        ))

    def arc_length(self, t0=None, t1=None, order=32):
        return self.integrate(self.speed, t0=t0, t1=t1, order=order)

    def integrated_curvature(self, order=32):
        return self.integrate(
            lambda t: self.curvature(t)*self.speed(t), order=order
        )

    def integrated_curvature_squared(self, order=32):
        return self.integrate(
            lambda t: self.curvature(t)**2*self.speed(t), order=order
        )


class LineEdgeGeometry(EdgeGeometry):
    """Exact straight edge between two endpoints."""

    conic_type = "line"

    def __init__(self, start, end):
        self.start = np.asarray(start, dtype=float)
        self.end = np.asarray(end, dtype=float)
        if self.start.shape != (3,) or self.end.shape != (3,):
            raise ValueError("Line endpoints must be three-vectors.")
        self._direction = self.end - self.start
        if not np.all(np.isfinite(self._direction)):
            raise ValueError("Line endpoints must be finite.")
        if np.linalg.norm(self._direction) < 1e-12:
            raise ValueError("Line endpoints must be distinct.")
        super().__init__(0.0, 1.0)

    def point(self, t):
        return self.start + float(t)*self._direction

    def tangent(self, t):
        return self._direction.copy()

    def second_derivative(self, t):
        return np.zeros(3, dtype=float)

    def arc_length(self, t0=None, t1=None, order=32):
        lower = self.t_min if t0 is None else float(t0)
        upper = self.t_max if t1 is None else float(t1)
        if lower < self.t_min or upper > self.t_max or upper < lower:
            raise ValueError("Integration bounds fall outside the edge interval.")
        return float(np.linalg.norm(self._direction))*(upper - lower)


class ParametricEdgeGeometry(EdgeGeometry):
    """Validated adapter for an analytic parameterization and derivatives."""

    def __init__(self, point, tangent, second_derivative, t_min, t_max,
                 conic_type="conic"):
        self._point = point
        self._tangent = tangent
        self._second_derivative = second_derivative
        self.conic_type = str(conic_type)
        super().__init__(t_min, t_max)

    @staticmethod
    def _vector(value, label):
        vector = np.asarray(value, dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError(f"{label} must return a finite three-vector.")
        return vector

    def point(self, t):
        return self._vector(self._point(float(t)), "point")

    def tangent(self, t):
        return self._vector(self._tangent(float(t)), "tangent")

    def second_derivative(self, t):
        return self._vector(
            self._second_derivative(float(t)), "second_derivative"
        )


class AdditivelyWeightedTrisectorBranch(EdgeGeometry):
    """One analytic branch of a three-sphere AW trisector.

    The parameter ``t`` is the common signed clearance (equivalently, tangent
    sphere radius) ``rho`` satisfying ``|x-p_i| - r_i = rho`` for all three
    generators.  Eliminating pairs of these equations gives

    ``x(rho) = x0 + u*rho + branch*n*sqrt(q(rho))``.

    A branch interval may not include a root of ``q`` because this particular
    parameterization becomes singular at a conic turning point. Such an edge
    must be represented by multiple branches or a later angle/hyperbolic
    parameterization.
    """

    def __init__(self, generator_locations, generator_radii, rho_min, rho_max,
                 branch=1, tol=1e-12):
        locations = np.asarray(generator_locations, dtype=float)
        radii = np.asarray(generator_radii, dtype=float)
        if locations.shape != (3, 3) or radii.shape != (3,):
            raise ValueError("AW trisectors require three 3D locations and radii.")
        if not np.all(np.isfinite(locations)) or not np.all(np.isfinite(radii)):
            raise ValueError("AW generator locations and radii must be finite.")
        if branch not in (-1, 1):
            raise ValueError("AW trisector branch must be -1 or +1.")

        rows = locations[1:] - locations[0]
        gram = rows @ rows.T
        if abs(float(np.linalg.det(gram))) < tol:
            raise ValueError("AW trisector generators must not be collinear.")

        radius_differences = radii[1:] - radii[0]
        constants = 0.5*(
            np.sum(locations[1:]**2, axis=1)
            - np.sum(locations[0]**2)
            - (radii[1:]**2 - radii[0]**2)
        )
        inverse_actions = np.linalg.solve(
            gram, np.column_stack((constants, radius_differences))
        )
        origin = rows.T @ inverse_actions[:, 0]
        linear = -(rows.T @ inverse_actions[:, 1])
        plane_normal = np.cross(rows[0], rows[1])
        plane_normal /= np.linalg.norm(plane_normal)

        # Choose the affine origin in the generator plane. Adding a nullspace
        # component does not alter the two eliminated linear equations.
        origin += plane_normal*float(plane_normal @ locations[0])
        relative_origin = origin - locations[0]
        q2 = 1.0 - float(linear @ linear)
        q1 = 2.0*radii[0] - 2.0*float(relative_origin @ linear)
        q0 = radii[0]**2 - float(relative_origin @ relative_origin)

        self.generator_locations = locations
        self.generator_radii = radii
        self.origin = origin
        self.basis_u = linear
        self.basis_v = plane_normal
        self.branch = int(branch)
        self.q_coefficients = np.array([q2, q1, q0], dtype=float)
        self.tol = float(tol)
        if q2 < -tol:
            self.conic_type = "ellipse"
        elif q2 > tol:
            self.conic_type = "hyperbola"
        else:
            self.conic_type = "parabola"

        super().__init__(rho_min, rho_max)
        candidates = [self.t_min, self.t_max]
        if abs(q2) > tol:
            critical = -q1/(2.0*q2)
            if self.t_min <= critical <= self.t_max:
                candidates.append(critical)
        if min(self._q(value) for value in candidates) <= tol:
            raise ValueError(
                "AW rho interval reaches or crosses a conic turning point."
            )

    def _q(self, rho):
        q2, q1, q0 = self.q_coefficients
        return float((q2*rho + q1)*rho + q0)

    def _q_prime(self, rho):
        q2, q1, _ = self.q_coefficients
        return float(2.0*q2*rho + q1)

    def point(self, t):
        rho = float(t)
        q = self._q(rho)
        if q < -self.tol:
            raise ValueError("AW trisector parameter lies outside the real conic.")
        return (
            self.origin + self.basis_u*rho
            + self.branch*self.basis_v*np.sqrt(max(q, 0.0))
        )

    def tangent(self, t):
        rho = float(t)
        q = self._q(rho)
        if q <= self.tol:
            raise ValueError("AW rho tangent is singular at a conic turning point.")
        return (
            self.basis_u
            + self.branch*self.basis_v*self._q_prime(rho)/(2.0*np.sqrt(q))
        )

    def second_derivative(self, t):
        rho = float(t)
        q = self._q(rho)
        if q <= self.tol:
            raise ValueError("AW rho derivative is singular at a conic turning point.")
        q_prime = self._q_prime(rho)
        q_second = 2.0*self.q_coefficients[0]
        scalar = q_second/(2.0*np.sqrt(q)) - q_prime**2/(4.0*q**1.5)
        return self.branch*self.basis_v*scalar

    def clearance_residuals(self, t):
        """Return each generator clearance minus the parameter value."""
        point = self.point(t)
        clearances = np.linalg.norm(
            point - self.generator_locations, axis=1
        ) - self.generator_radii
        return clearances - float(t)

    def parameter_for_point(self, point):
        """Recover rho from a point expected to lie on this trisector."""
        point = np.asarray(point, dtype=float)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("AW trisector point must be a finite three-vector.")
        return float(np.linalg.norm(point - self.generator_locations[0])
                     - self.generator_radii[0])


class AffineEdgeGeometry(EdgeGeometry):
    """Reparameterize another edge onto [0, 1] with affine orientation."""

    def __init__(self, source, source_start, source_end):
        if not isinstance(source, EdgeGeometry):
            raise TypeError("source must implement EdgeGeometry.")
        self.source = source
        self.source_start = float(source_start)
        self.source_end = float(source_end)
        if self.source_start == self.source_end:
            raise ValueError("Affine edge parameter endpoints must be distinct.")
        lower = min(self.source_start, self.source_end)
        upper = max(self.source_start, self.source_end)
        if lower < source.t_min or upper > source.t_max:
            raise ValueError("Affine parameters fall outside the source interval.")
        self._scale = self.source_end - self.source_start
        self.conic_type = getattr(source, "conic_type", "curve")
        super().__init__(0.0, 1.0)

    def _source_parameter(self, t):
        return self.source_start + float(t)*self._scale

    def point(self, t):
        return self.source.point(self._source_parameter(t))

    def tangent(self, t):
        return self.source.tangent(self._source_parameter(t))*self._scale

    def second_derivative(self, t):
        return self.source.second_derivative(
            self._source_parameter(t)
        )*self._scale**2


def match_aw_trisector_to_endpoints(generator_locations, generator_radii,
                                    vertex_indices, vertex_locations,
                                    tolerance=1e-8):
    """Match two stored vertices to a single analytic AW trisector branch.

    A successful geometry runs from the lower canonical vertex index at t=0
    to the higher index at t=1. Failure statuses are diagnostic: crossing a
    turning point or switching branches needs a different parameterization.
    """
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive.")
    endpoints = canonical_edge_endpoints(vertex_indices, vertex_locations)
    start, end = endpoints[0][1], endpoints[1][1]
    locations = np.asarray(generator_locations, dtype=float)
    radii = np.asarray(generator_radii, dtype=float)
    if locations.shape != (3, 3) or radii.shape != (3,):
        raise ValueError("AW trisectors require three 3D locations and radii.")

    rhos = np.array([
        np.linalg.norm(start - locations[0]) - radii[0],
        np.linalg.norm(end - locations[0]) - radii[0],
    ])
    if abs(float(rhos[1] - rhos[0])) <= tolerance:
        return AWEdgeGeometryMatch(
            "turning_point",
            "Endpoint clearances coincide; rho cannot parameterize this edge.",
            None, (), float("inf"),
        )

    try:
        probe = AdditivelyWeightedTrisectorBranch(
            locations, radii, float(np.min(rhos)), float(np.max(rhos)),
            branch=1, tol=tolerance,
        )
    except ValueError as error:
        status = "turning_point" if "turning point" in str(error) else "invalid"
        return AWEdgeGeometryMatch(status, str(error), None, (), float("inf"))

    branch_coordinates = []
    clearance_error = 0.0
    for point, rho in zip((start, end), rhos):
        affine_point = probe.origin + probe.basis_u*rho
        branch_coordinates.append(float(probe.basis_v @ (point - affine_point)))
        clearances = np.linalg.norm(point - locations, axis=1) - radii
        clearance_error = max(
            clearance_error, float(np.max(np.abs(clearances - rho)))
        )

    if min(abs(value) for value in branch_coordinates) <= tolerance:
        return AWEdgeGeometryMatch(
            "turning_point",
            "An endpoint lies at or too near the rho turning point.",
            None, (), clearance_error,
        )
    branches = tuple(1 if value > 0 else -1 for value in branch_coordinates)
    if branches[0] != branches[1]:
        return AWEdgeGeometryMatch(
            "cross_branch",
            "Endpoints lie on opposite analytic conic branches.",
            None, (), clearance_error,
        )

    branch = AdditivelyWeightedTrisectorBranch(
        locations, radii, float(np.min(rhos)), float(np.max(rhos)),
        branch=branches[0], tol=tolerance,
    )
    geometry = AffineEdgeGeometry(branch, rhos[0], rhos[1])
    endpoint_errors = (
        float(np.linalg.norm(geometry.point(0.0) - start)),
        float(np.linalg.norm(geometry.point(1.0) - end)),
    )
    maximum_error = max(*endpoint_errors, clearance_error)
    if maximum_error > tolerance:
        return AWEdgeGeometryMatch(
            "residual_too_large",
            f"Analytic endpoint residual {maximum_error:.6g} exceeds tolerance.",
            None, endpoint_errors, clearance_error,
        )
    return AWEdgeGeometryMatch(
        "matched", "Endpoints match one analytic AW branch.",
        geometry, endpoint_errors, clearance_error,
    )


def validate_edge_samples(geometry, points):
    """Compare stored AW samples with their rho-based reconstruction."""
    if not isinstance(geometry, AffineEdgeGeometry) or not isinstance(
            geometry.source, AdditivelyWeightedTrisectorBranch):
        raise TypeError("Sample validation requires a matched analytic AW edge.")
    samples = np.asarray(points, dtype=float)
    if samples.ndim != 2 or samples.shape[1:] != (3,) or len(samples) == 0:
        raise ValueError("points must be a non-empty array of three-vectors.")
    if not np.all(np.isfinite(samples)):
        raise ValueError("points must contain only finite values.")
    branch = geometry.source
    residuals = []
    for point in samples:
        rho = branch.parameter_for_point(point)
        residuals.append(float(np.linalg.norm(branch.point(rho) - point)))
    residuals = np.asarray(residuals)
    return EdgeSampleValidation(
        len(samples), float(np.max(residuals)),
        float(np.sqrt(np.mean(residuals**2))),
    )
