"""Geometry interfaces and normalized construction data for VorPy edges.

The interface is independent of network topology and storage.  It lets later
boundary code integrate intrinsic curve descriptors without depending on the
legacy sampled ``edge.points`` representation.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass

from time import perf_counter

import numpy as np

from vorpy.src.calculations import aw_outward_normal


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


class AdditivelyWeightedTrisectorConic(EdgeGeometry):
    """Nonsingular parameterization of a complete three-sphere AW conic.

    Unlike the clearance-parameterized branch, this representation remains
    regular where the normal coordinate crosses zero. Ellipses use an angle,
    hyperbolas use a hyperbolic angle, and parabolas use the signed normal
    coordinate.
    """

    def __init__(self, generator_locations, generator_radii, t_min, t_max,
                 component=None, tol=1e-12):
        locations = np.asarray(generator_locations, dtype=float)
        radii = np.asarray(generator_radii, dtype=float)
        if locations.shape != (3, 3) or radii.shape != (3,):
            raise ValueError("AW trisectors require three 3D locations and radii.")
        if not np.all(np.isfinite(locations)) or not np.all(np.isfinite(radii)):
            raise ValueError("AW generator locations and radii must be finite.")

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
        self.q_coefficients = np.array([q2, q1, q0], dtype=float)
        self.tol = float(tol)
        self.component = component
        self.rho_center = None
        self.rho_scale = None
        self.normal_scale = None
        self.hyperbola_orientation = None

        if q2 < -tol:
            self.conic_type = "ellipse"
            self.rho_center = -q1/(2.0*q2)
            q_center = self._q(self.rho_center)
            if q_center <= tol:
                raise ValueError("Degenerate or non-real AW ellipse.")
            self.rho_scale = np.sqrt(q_center/(-q2))
            self.normal_scale = np.sqrt(q_center)
        elif q2 > tol:
            self.conic_type = "hyperbola"
            if component not in (-1, 1):
                raise ValueError("AW hyperbola component must be -1 or +1.")
            self.component = int(component)
            self.rho_center = -q1/(2.0*q2)
            q_center = self._q(self.rho_center)
            if q_center < -tol:
                self.hyperbola_orientation = "rho"
                self.rho_scale = np.sqrt((-q_center)/q2)
                self.normal_scale = np.sqrt(-q_center)
            elif q_center > tol:
                self.hyperbola_orientation = "normal"
                self.rho_scale = np.sqrt(q_center/q2)
                self.normal_scale = np.sqrt(q_center)
            else:
                raise ValueError("Degenerate AW hyperbola is not supported.")
        else:
            self.conic_type = "parabola"
            if abs(q1) <= tol:
                raise ValueError("Degenerate AW parabola is not supported.")

        super().__init__(t_min, t_max)

    def _q(self, rho):
        q2, q1, q0 = self.q_coefficients
        return float((q2*rho + q1)*rho + q0)

    def rho(self, t):
        value = float(t)
        if self.conic_type == "ellipse":
            return self.rho_center + self.rho_scale*np.cos(value)
        if self.conic_type == "hyperbola":
            if self.hyperbola_orientation == "rho":
                return (
                    self.rho_center
                    + self.component*self.rho_scale*np.cosh(value)
                )
            return self.rho_center + self.rho_scale*np.sinh(value)
        _, q1, q0 = self.q_coefficients
        return (value**2 - q0)/q1

    def normal_coordinate(self, t):
        value = float(t)
        if self.conic_type == "ellipse":
            return self.normal_scale*np.sin(value)
        if self.conic_type == "hyperbola":
            if self.hyperbola_orientation == "rho":
                return self.normal_scale*np.sinh(value)
            return self.component*self.normal_scale*np.cosh(value)
        return value

    def point(self, t):
        return (
            self.origin + self.basis_u*self.rho(t)
            + self.basis_v*self.normal_coordinate(t)
        )

    def tangent(self, t):
        value = float(t)
        if self.conic_type == "ellipse":
            rho_prime = -self.rho_scale*np.sin(value)
            normal_prime = self.normal_scale*np.cos(value)
        elif self.conic_type == "hyperbola":
            if self.hyperbola_orientation == "rho":
                rho_prime = (
                    self.component*self.rho_scale*np.sinh(value)
                )
                normal_prime = self.normal_scale*np.cosh(value)
            else:
                rho_prime = self.rho_scale*np.cosh(value)
                normal_prime = (
                    self.component*self.normal_scale*np.sinh(value)
                )
        else:
            q1 = self.q_coefficients[1]
            rho_prime = 2.0*value/q1
            normal_prime = 1.0
        return self.basis_u*rho_prime + self.basis_v*normal_prime

    def second_derivative(self, t):
        value = float(t)
        if self.conic_type == "ellipse":
            rho_second = -self.rho_scale*np.cos(value)
            normal_second = -self.normal_scale*np.sin(value)
        elif self.conic_type == "hyperbola":
            if self.hyperbola_orientation == "rho":
                rho_second = (
                    self.component*self.rho_scale*np.cosh(value)
                )
                normal_second = self.normal_scale*np.sinh(value)
            else:
                rho_second = self.rho_scale*np.sinh(value)
                normal_second = (
                    self.component*self.normal_scale*np.cosh(value)
                )
        else:
            rho_second = 2.0/self.q_coefficients[1]
            normal_second = 0.0
        return self.basis_u*rho_second + self.basis_v*normal_second

    def clearance_residuals(self, t):
        point = self.point(t)
        clearances = np.linalg.norm(
            point - self.generator_locations, axis=1
        ) - self.generator_radii
        return clearances - self.rho(t)

    def parameter_for_point(self, point):
        point = np.asarray(point, dtype=float)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("AW trisector point must be a finite three-vector.")
        rho = float(
            np.linalg.norm(point - self.generator_locations[0])
            - self.generator_radii[0]
        )
        affine = self.origin + self.basis_u*rho
        normal = float(self.basis_v @ (point - affine))

        if self.conic_type == "ellipse":
            return float(np.arctan2(
                normal/self.normal_scale,
                (rho - self.rho_center)/self.rho_scale,
            ))
        if self.conic_type == "hyperbola":
            if self.hyperbola_orientation == "rho":
                if self.component*(rho - self.rho_center) < -self.tol:
                    raise ValueError("Point lies on the other hyperbola component.")
                return float(np.arcsinh(normal/self.normal_scale))
            if self.component*normal < -self.tol:
                raise ValueError("Point lies on the other hyperbola component.")
            return float(np.arcsinh(
                (rho - self.rho_center)/self.rho_scale
            ))
        return normal


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


def match_aw_trisector_conic_to_samples(
        generator_locations, generator_radii, vertex_indices, vertex_locations,
        points, tolerance=1e-8):
    """Match a sampled edge to a nonsingular complete AW conic.

    Samples select the intended elliptical arc and the connected hyperbola
    component. The result retains canonical vertex-index orientation.
    """
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive.")
    endpoints = canonical_edge_endpoints(vertex_indices, vertex_locations)
    start, end = endpoints[0][1], endpoints[1][1]
    samples = np.asarray(points, dtype=float)
    if samples.ndim != 2 or samples.shape[1:] != (3,) or len(samples) < 2:
        raise ValueError("points must contain at least two finite three-vectors.")
    if not np.all(np.isfinite(samples)):
        raise ValueError("points must contain only finite values.")

    forward_cost = (
        np.linalg.norm(samples[0] - start)
        + np.linalg.norm(samples[-1] - end)
    )
    reverse_cost = (
        np.linalg.norm(samples[-1] - start)
        + np.linalg.norm(samples[0] - end)
    )
    if reverse_cost < forward_cost:
        samples = samples[::-1]

    locations = np.asarray(generator_locations, dtype=float)
    radii = np.asarray(generator_radii, dtype=float)
    candidates = []
    construction_errors = []
    for component in (None, 1, -1):
        try:
            probe = AdditivelyWeightedTrisectorConic(
                locations, radii, 0.0, 1.0, component=component,
                tol=tolerance,
            )
        except ValueError as error:
            construction_errors.append(str(error))
            continue

        # Ellipses and parabolas have no disconnected component, so avoid
        # evaluating duplicate candidates after the component-free probe.
        if component is not None and probe.conic_type != "hyperbola":
            continue
        try:
            sample_parameters = np.array([
                probe.parameter_for_point(point) for point in samples
            ])
            start_parameter = probe.parameter_for_point(start)
            end_parameter = probe.parameter_for_point(end)
        except ValueError:
            continue

        if probe.conic_type == "ellipse":
            sample_parameters = np.unwrap(sample_parameters)
            two_pi = 2.0*np.pi
            start_parameter += two_pi*round(
                (sample_parameters[0] - start_parameter)/two_pi
            )
            end_parameter += two_pi*round(
                (sample_parameters[-1] - end_parameter)/two_pi
            )

        if abs(float(end_parameter - start_parameter)) <= tolerance:
            continue
        source = AdditivelyWeightedTrisectorConic(
            locations, radii,
            min(start_parameter, end_parameter),
            max(start_parameter, end_parameter),
            component=probe.component,
            tol=tolerance,
        )
        geometry = AffineEdgeGeometry(
            source, start_parameter, end_parameter
        )
        endpoint_errors = (
            float(np.linalg.norm(geometry.point(0.0) - start)),
            float(np.linalg.norm(geometry.point(1.0) - end)),
        )
        clearances = []
        for point in (start, end):
            values = np.linalg.norm(point - locations, axis=1) - radii
            clearances.append(float(np.max(np.abs(values - values[0]))))
        clearance_error = max(clearances)

        sample_errors = []
        for point, parameter in zip(samples, sample_parameters):
            sample_errors.append(float(
                np.linalg.norm(source.point(parameter) - point)
            ))
        candidates.append((
            max(*endpoint_errors, clearance_error, max(sample_errors)),
            geometry,
            endpoint_errors,
            clearance_error,
        ))

        if probe.conic_type != "hyperbola":
            break

    if not candidates:
        reason = (
            "Stored samples do not lie on one connected nonsingular AW conic."
        )
        if construction_errors:
            reason += f" Construction detail: {construction_errors[-1]}"
        return AWEdgeGeometryMatch(
            "disconnected", reason, None, (), float("inf")
        )

    maximum_error, geometry, endpoint_errors, clearance_error = min(
        candidates, key=lambda candidate: candidate[0]
    )
    if maximum_error > tolerance:
        return AWEdgeGeometryMatch(
            "residual_too_large",
            f"Nonsingular conic residual {maximum_error:.6g} exceeds tolerance.",
            None, endpoint_errors, clearance_error,
        )
    return AWEdgeGeometryMatch(
        "matched_nonsingular",
        "Samples match one connected nonsingular AW conic.",
        geometry, endpoint_errors, clearance_error,
    )


def validate_edge_samples(geometry, points):
    """Compare stored AW samples with their analytic reconstruction."""
    if not isinstance(geometry, AffineEdgeGeometry) or not isinstance(
            geometry.source, (
                AdditivelyWeightedTrisectorBranch,
                AdditivelyWeightedTrisectorConic,
            )):
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


@dataclass(frozen=True)
class AWEdgePointGeometry:
    """Local differential geometry of one AW cell along a regular edge."""

    point: np.ndarray
    tangent: np.ndarray
    first_normal: np.ndarray
    second_normal: np.ndarray
    cell_index: int
    neighbor_indices: tuple


def aw_edge_point_geometry(
        edge_geometry,
        t,
        generator_indices,
        generator_locations,
        cell_index,
        tol=1e-12):
    """Return exact local boundary geometry for one AW cell along an edge.

    A regular AW edge is generated by exactly three balls. For the selected
    cell, the two incident boundary faces are the pairwise AW surfaces between
    that cell and each of the other two generators.

    Parameters
    ----------
    edge_geometry : EdgeGeometry
        Exact analytic geometry of the edge.
    t : float
        Edge parameter.
    generator_indices : iterable of int
        The three generating ball indices.
    generator_locations : array-like, shape (3, 3)
        Locations corresponding to generator_indices.
    cell_index : int
        Generator whose cell boundary is being evaluated.
    tol : float
        Numerical tolerance.

    Returns
    -------
    AWEdgePointGeometry
        Exact point, unit tangent, and the two outward incident normals.
    """
    if not isinstance(edge_geometry, EdgeGeometry):
        raise TypeError("edge_geometry must implement EdgeGeometry.")

    indices = tuple(int(index) for index in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)

    if len(indices) != 3 or locations.shape != (3, 3):
        raise ValueError(
            "A regular AW edge requires exactly three generator indices "
            "and three generator locations."
        )

    if len(set(indices)) != 3:
        raise ValueError("AW edge generator indices must be distinct.")

    if not np.all(np.isfinite(locations)):
        raise ValueError("AW edge generator locations must be finite.")

    cell_index = int(cell_index)

    if cell_index not in indices:
        raise ValueError(
            "cell_index must be one of the three AW edge generators."
        )

    point = np.asarray(edge_geometry.point(t), dtype=float)
    tangent = np.asarray(
        edge_geometry.unit_tangent(t, tol=tol),
        dtype=float,
    )

    cell_position = indices.index(cell_index)

    neighbor_positions = [
        position
        for position in range(3)
        if position != cell_position
    ]

    neighbor_indices = tuple(
        indices[position]
        for position in neighbor_positions
    )

    cell_location = locations[cell_position]

    first_normal = aw_outward_normal(
        point=point,
        cell_location=cell_location,
        neighbor_location=locations[neighbor_positions[0]],
        tol=tol,
    )

    second_normal = aw_outward_normal(
        point=point,
        cell_location=cell_location,
        neighbor_location=locations[neighbor_positions[1]],
        tol=tol,
    )

    return AWEdgePointGeometry(
        point=point,
        tangent=tangent,
        first_normal=np.asarray(first_normal, dtype=float),
        second_normal=np.asarray(second_normal, dtype=float),
        cell_index=cell_index,
        neighbor_indices=neighbor_indices,
    )


def boundary_turning_angle(first_normal, second_normal, tol=1e-12):
    """Return the unsigned turning angle between two outward surface normals.

    The result lies in [0, pi]. It represents the local normal jump across
    a regular edge of a single piecewise-smooth cell boundary.

    This quantity is independent of edge orientation and of the ordering
    of the two incident surfaces.
    """
    first = np.asarray(first_normal, dtype=float)
    second = np.asarray(second_normal, dtype=float)

    if first.shape != (3,) or second.shape != (3,):
        raise ValueError("Boundary normals must be three-vectors.")

    if not np.all(np.isfinite(first)) or not np.all(np.isfinite(second)):
        raise ValueError("Boundary normals must be finite.")

    first_mag = float(np.linalg.norm(first))
    second_mag = float(np.linalg.norm(second))

    if first_mag < tol or second_mag < tol:
        raise ValueError("Boundary normal is undefined.")

    first = first / first_mag
    second = second / second_mag

    cosine = float(np.dot(first, second))

    # Protect arccos from roundoff such as 1.0000000000000002.
    cosine = float(np.clip(cosine, -1.0, 1.0))

    return float(np.arccos(cosine))


def aw_cell_turning_angle(
        edge_geometry,
        t,
        generator_indices,
        generator_locations,
        cell_index,
        tol=1e-12):
    """Return the local boundary turning angle for one AW cell at an edge point."""

    local = aw_edge_point_geometry(
        edge_geometry=edge_geometry,
        t=t,
        generator_indices=generator_indices,
        generator_locations=generator_locations,
        cell_index=cell_index,
        tol=tol,
    )

    return boundary_turning_angle(
        local.first_normal,
        local.second_normal,
        tol=tol,
    )


def aw_cell_edge_mean_curvature(
        edge_geometry,
        generator_indices,
        generator_locations,
        cell_index,
        order=32,
        tol=1e-12):
    """Return one edge's contribution to integrated mean curvature for one AW cell.

    Computes

        M_edge = 1/2 * integral(theta(s) ds)

    where theta is the local turning angle between the two outward incident
    surface normals of the selected AW cell.

    This is boundary mean curvature concentrated on a piecewise-smooth edge.
    It is not the intrinsic curvature integral of the edge curve itself.
    """
    if not isinstance(edge_geometry, EdgeGeometry):
        raise TypeError("edge_geometry must implement EdgeGeometry.")

    order = int(order)
    if order < 1:
        raise ValueError("Quadrature order must be positive.")

    def integrand(t):
        angle = aw_cell_turning_angle(
            edge_geometry=edge_geometry,
            t=t,
            generator_indices=generator_indices,
            generator_locations=generator_locations,
            cell_index=cell_index,
            tol=tol,
        )

        speed = edge_geometry.speed(t)

        if not np.isfinite(angle) or not np.isfinite(speed):
            raise ValueError(
                "AW edge mean-curvature integrand is non-finite."
            )

        return angle * speed

    return 0.5 * edge_geometry.integrate(
        integrand,
        order=order,
    )


def aw_edge_mean_curvatures(
        edge_geometry,
        generator_indices,
        generator_locations,
        order=32,
        quadrature=None,
        tol=1e-12):
    """Integrate edge mean curvature for all three AW cells in one pass.

    For an AW edge generated by cells i, j, k, each quadrature point is shared
    by all three cell perspectives. The edge point, speed, radial unit vectors,
    and pairwise surface normals are therefore evaluated once.

    Returns
    -------
    dict
        {cell_index: 0.5 * integral(theta_cell ds)}
    """
    indices = tuple(int(value) for value in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)

    if len(indices) != 3 or locations.shape != (3, 3):
        raise ValueError("An AW edge requires exactly three generators.")

    if quadrature is None:
        nodes, weights = np.polynomial.legendre.leggauss(int(order))
    else:
        nodes, weights = quadrature

    lower, upper = edge_geometry.t_min, edge_geometry.t_max
    half = 0.5 * (upper - lower)
    center = 0.5 * (upper + lower)

    totals = np.zeros(3, dtype=float)

    for node, weight in zip(nodes, weights):
        t = center + half * node

        # Geometry of the edge itself is common to all three cells.
        point = np.asarray(edge_geometry.point(t), dtype=float)
        tangent = np.asarray(edge_geometry.tangent(t), dtype=float)
        speed = float(np.linalg.norm(tangent))

        if not np.isfinite(speed) or speed < tol:
            raise ValueError("Edge tangent is undefined during AW curvature integration.")

        # Radial unit vectors from the three generators to this edge point.
        radial = point - locations
        distances = np.linalg.norm(radial, axis=1)

        if np.any(~np.isfinite(distances)) or np.any(distances < tol):
            raise ValueError("AW edge point coincides with a generator.")

        radial /= distances[:, None]

        # For cell i, outward normal across its boundary with j is
        #
        #     n_i->j = normalize(u_i - u_j)
        #
        # where u_i is the unit radial vector from generator i.
        n01 = radial[0] - radial[1]
        n02 = radial[0] - radial[2]
        n12 = radial[1] - radial[2]

        norm01 = float(np.linalg.norm(n01))
        norm02 = float(np.linalg.norm(n02))
        norm12 = float(np.linalg.norm(n12))

        if min(norm01, norm02, norm12) < tol:
            raise ValueError("AW pairwise surface normal is undefined on edge.")

        n01 /= norm01
        n02 /= norm02
        n12 /= norm12

        # Cell 0: n_0->1 and n_0->2
        # Cell 1: n_1->0 = -n01 and n_1->2 = n12
        # Cell 2: n_2->0 = -n02 and n_2->1 = -n12
        cos0 = np.clip(np.dot(n01, n02), -1.0, 1.0)
        cos1 = np.clip(np.dot(-n01, n12), -1.0, 1.0)
        cos2 = np.clip(np.dot(-n02, -n12), -1.0, 1.0)

        angles = np.arccos([cos0, cos1, cos2])
        totals += weight * angles * speed

    totals *= 0.5 * half

    if not np.all(np.isfinite(totals)):
        raise ValueError("Non-finite AW edge mean-curvature integral.")

    return {
        indices[0]: float(totals[0]),
        indices[1]: float(totals[1]),
        indices[2]: float(totals[2]),
    }


def signed_geodesic_curvature(edge_geometry, t, surface_normal, tol=1e-12):
    """Return signed geodesic curvature of an edge on an oriented surface.

    The sign convention is

        k_g = (dT/ds) · (n × T),

    where T is the unit edge tangent and n is the oriented surface normal.
    For a positively oriented CCW circle in the xy-plane with n = +z,
    k_g is positive.
    """
    first = np.asarray(edge_geometry.tangent(t), dtype=float)
    second = np.asarray(edge_geometry.second_derivative(t), dtype=float)
    normal = np.asarray(surface_normal, dtype=float)

    if first.shape != (3,) or second.shape != (3,) or normal.shape != (3,):
        raise ValueError("Edge derivatives and surface normal must be three-vectors.")

    if not np.all(np.isfinite(first)) or not np.all(np.isfinite(second)):
        raise ValueError("Edge derivatives must be finite.")
    if not np.all(np.isfinite(normal)):
        raise ValueError("Surface normal must be finite.")

    speed = float(np.linalg.norm(first))
    normal_norm = float(np.linalg.norm(normal))

    if speed < tol:
        raise ValueError("Edge tangent is undefined.")
    if normal_norm < tol:
        raise ValueError("Surface normal is undefined.")

    normal /= normal_norm

    value = float(np.dot(second, np.cross(normal, first)) / speed**3)

    if not np.isfinite(value):
        raise ValueError("Geodesic curvature is non-finite.")

    return value


def integrated_geodesic_curvature(
        edge_geometry,
        surface_normal,
        order=32):
    """Integrate signed geodesic curvature along an analytic edge.

    ``surface_normal`` is a callable accepting the edge parameter t and
    returning the correctly oriented local surface normal.
    """
    return edge_geometry.integrate(
        lambda t: signed_geodesic_curvature(
            edge_geometry=edge_geometry,
            t=t,
            surface_normal=surface_normal(t),
        ) * edge_geometry.speed(t),
        order=order,
    )


def aw_face_edge_geodesic_curvature(
        edge_geometry,
        t,
        generator_indices,
        generator_locations,
        cell_index,
        other_index,
        tol=1e-12):
    """Return signed geodesic curvature for one AW face-edge incidence.

    The selected face is the interface between ``cell_index`` and
    ``other_index``. Its normal is oriented outward from ``cell_index``.

    This function uses the edge's current parameter direction. Boundary-cycle
    orientation is handled separately when integrating around a complete face.
    """
    indices = tuple(int(value) for value in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)

    if len(indices) != 3 or locations.shape != (3, 3):
        raise ValueError("An AW edge requires exactly three generators.")

    cell_index = int(cell_index)
    other_index = int(other_index)

    if cell_index not in indices or other_index not in indices:
        raise ValueError("Selected face generators must belong to the AW edge.")

    if cell_index == other_index:
        raise ValueError("An AW face requires two distinct generators.")

    point = np.asarray(edge_geometry.point(t), dtype=float)

    cell_pos = indices.index(cell_index)
    other_pos = indices.index(other_index)

    cell_radial = point - locations[cell_pos]
    other_radial = point - locations[other_pos]

    cell_norm = float(np.linalg.norm(cell_radial))
    other_norm = float(np.linalg.norm(other_radial))

    if cell_norm < tol or other_norm < tol:
        raise ValueError("AW edge point coincides with a generator.")

    cell_radial /= cell_norm
    other_radial /= other_norm

    surface_normal = cell_radial - other_radial
    normal_norm = float(np.linalg.norm(surface_normal))

    if normal_norm < tol:
        raise ValueError("AW face normal is undefined on the edge.")

    surface_normal /= normal_norm

    return signed_geodesic_curvature(
        edge_geometry=edge_geometry,
        t=t,
        surface_normal=surface_normal,
        tol=tol,
    )


def integrated_aw_face_edge_geodesic_curvature(
        edge_geometry,
        generator_indices,
        generator_locations,
        cell_index,
        other_index,
        order=32):
    """Integrate signed geodesic curvature along one AW face boundary edge."""
    return edge_geometry.integrate(
        lambda t: aw_face_edge_geodesic_curvature(
            edge_geometry=edge_geometry,
            t=t,
            generator_indices=generator_indices,
            generator_locations=generator_locations,
            cell_index=cell_index,
            other_index=other_index,
        ) * edge_geometry.speed(t),
        order=order,
    )


def network_surface_boundary_cycle(net, surface_index):
    """Return a surface's edges as one exact topological boundary cycle.

    Returns
    -------
    list[tuple[int, int, int]]
        ``(edge_index, start_vertex, end_vertex)`` for each boundary edge.

    Notes
    -----
    This uses edge/vertex incidence only. It does not depend on sampled edge
    points, mesh triangulation, or the storage order of ``surface['edges']``.
    The returned cycle has a deterministic start/direction but is not yet
    oriented relative to either incident cell.
    """
    surface_index = int(surface_index)
    edge_indices = [int(i) for i in net.surfs.loc[surface_index, "edges"]]

    if len(edge_indices) < 2:
        raise ValueError(
            f"Surface {surface_index} has only {len(edge_indices)} edges."
        )

    edge_vertices = {}

    for edge_index in edge_indices:
        vertices = [int(v) for v in net.edges.loc[edge_index, "verts"]]

        if len(vertices) != 2:
            raise ValueError(
                f"Surface {surface_index}, edge {edge_index} has "
                f"{len(vertices)} vertices; expected 2."
            )

        if vertices[0] == vertices[1]:
            raise ValueError(
                f"Surface {surface_index}, edge {edge_index} has identical "
                f"endpoint vertices."
            )

        edge_vertices[edge_index] = tuple(vertices)

    # Each boundary vertex of a closed simple face must touch exactly two
    # boundary edges.
    vertex_edges = {}

    for edge_index, (v0, v1) in edge_vertices.items():
        vertex_edges.setdefault(v0, []).append(edge_index)
        vertex_edges.setdefault(v1, []).append(edge_index)

    bad_vertices = {
        vertex: len(edges)
        for vertex, edges in vertex_edges.items()
        if len(edges) != 2
    }

    if bad_vertices:
        raise ValueError(
            f"Surface {surface_index} is not a simple closed edge cycle; "
            f"boundary vertex degrees={bad_vertices}."
        )

    # Deterministic starting edge and direction. Physical cell orientation
    # will be imposed separately.
    start_edge = min(edge_indices)
    start_v0, start_v1 = edge_vertices[start_edge]
    start_vertex, next_vertex = sorted((start_v0, start_v1))

    cycle = [(start_edge, start_vertex, next_vertex)]
    used_edges = {start_edge}
    current_vertex = next_vertex

    while len(used_edges) < len(edge_indices):
        candidates = [
            edge_index
            for edge_index in vertex_edges[current_vertex]
            if edge_index not in used_edges
        ]

        if len(candidates) != 1:
            raise ValueError(
                f"Surface {surface_index} boundary is ambiguous at "
                f"vertex {current_vertex}; unused edges={candidates}."
            )

        edge_index = candidates[0]
        v0, v1 = edge_vertices[edge_index]

        if v0 == current_vertex:
            next_vertex = v1
        elif v1 == current_vertex:
            next_vertex = v0
        else:
            raise ValueError(
                f"Surface {surface_index}, edge {edge_index} is not incident "
                f"to expected vertex {current_vertex}."
            )

        cycle.append((edge_index, current_vertex, next_vertex))
        used_edges.add(edge_index)
        current_vertex = next_vertex

    if current_vertex != start_vertex:
        raise ValueError(
            f"Surface {surface_index} boundary does not close: ended at "
            f"vertex {current_vertex}, expected {start_vertex}."
        )

    return cycle


def aw_face_edge_boundary_orientation(
        edge_geometry,
        generator_indices,
        generator_locations,
        cell_index,
        face_other_index,
        tol=1e-12):
    """Return +1/-1 for the analytic edge direction on an oriented AW face.

    The selected face is the interface between ``cell_index`` and
    ``face_other_index``. The third generator determines which side of the
    edge is the interior of that face patch.

    Returns
    -------
    int
        +1 if increasing edge parameter follows the positive boundary
        orientation of the selected face; -1 otherwise.
    """
    indices = tuple(int(v) for v in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)
    cell_index = int(cell_index)
    face_other_index = int(face_other_index)

    if len(indices) != 3 or locations.shape != (3, 3):
        raise ValueError("An AW edge requires exactly three generators.")
    if cell_index not in indices or face_other_index not in indices:
        raise ValueError("Selected face generators must belong to the edge.")
    if cell_index == face_other_index:
        raise ValueError("An AW face requires two distinct generators.")

    third_candidates = [
        index for index in indices
        if index not in {cell_index, face_other_index}
    ]

    if len(third_candidates) != 1:
        raise ValueError("Could not identify the third AW edge generator.")

    third_index = third_candidates[0]
    t = 0.5 * (edge_geometry.t_min + edge_geometry.t_max)
    point = np.asarray(edge_geometry.point(t), dtype=float)

    def outward_normal(other_index):
        cell_pos = indices.index(cell_index)
        other_pos = indices.index(other_index)

        u_cell = point - locations[cell_pos]
        u_other = point - locations[other_pos]

        d_cell = float(np.linalg.norm(u_cell))
        d_other = float(np.linalg.norm(u_other))

        if d_cell < tol or d_other < tol:
            raise ValueError("AW edge point coincides with a generator.")

        u_cell /= d_cell
        u_other /= d_other

        normal = u_cell - u_other
        magnitude = float(np.linalg.norm(normal))

        if not np.isfinite(magnitude) or magnitude < tol:
            raise ValueError("AW face normal is undefined on the edge.")

        return normal / magnitude

    face_normal = outward_normal(face_other_index)
    third_normal = outward_normal(third_index)

    tangent = np.asarray(edge_geometry.unit_tangent(t, tol=tol), dtype=float)

    # Cell A lies opposite its outward A->C normal. Project that direction
    # into the tangent plane of face AB to obtain the local face interior.
    inward = -third_normal + np.dot(third_normal, face_normal) * face_normal
    inward_norm = float(np.linalg.norm(inward))

    if not np.isfinite(inward_norm) or inward_norm < tol:
        raise ValueError("AW face interior direction is undefined at the edge.")

    inward /= inward_norm

    # Positive oriented-boundary convention:
    #
    #     n_face × T_boundary -> into the face.
    #
    orientation = float(np.dot(np.cross(face_normal, tangent), inward))

    if not np.isfinite(orientation) or abs(orientation) < tol:
        raise ValueError("AW face-edge boundary orientation is degenerate.")

    return 1 if orientation > 0.0 else -1


def oriented_aw_face_edge_geodesic_curvature(
        edge_geometry,
        generator_indices,
        generator_locations,
        cell_index,
        face_other_index,
        order=32):
    """Return ∫kg ds using the positive boundary orientation of an AW face."""
    sign = aw_face_edge_boundary_orientation(
        edge_geometry=edge_geometry,
        generator_indices=generator_indices,
        generator_locations=generator_locations,
        cell_index=cell_index,
        face_other_index=face_other_index,
    )

    value = integrated_aw_face_edge_geodesic_curvature(
        edge_geometry=edge_geometry,
        generator_indices=generator_indices,
        generator_locations=generator_locations,
        cell_index=cell_index,
        other_index=face_other_index,
        order=order,
    )

    return float(sign * value)

def aw_edge_geodesic_curvatures(
        edge_geometry,
        generator_indices,
        generator_locations,
        face_pairs,
        order=32,
        quadrature=None,
        tol=1e-12):
    """Integrate oriented geodesic curvature for many AW face incidences at once.

    Parameters
    ----------
    edge_geometry : EdgeGeometry
        Exact analytic geometry of one physical AW edge.
    generator_indices : iterable of int
        The three generators defining the edge.
    generator_locations : array-like, shape (3, 3)
        Generator positions corresponding to ``generator_indices``.
    face_pairs : iterable of (cell_index, other_index)
        Directed cell/face incidences to evaluate. ``(i, j)`` means the
        boundary face between cell ``i`` and neighbor ``j``, oriented as the
        positive boundary of cell ``i``.
    order : int
        Gauss-Legendre quadrature order.
    quadrature : tuple(array, array), optional
        Precomputed Gauss-Legendre nodes and weights. Supplying this avoids
        rebuilding the same rule for every network edge.
    tol : float
        Numerical tolerance.

    Returns
    -------
    dict
        ``{(cell_index, other_index): integral(k_g ds), ...}``

    Notes
    -----
    Edge position, first/second derivatives, radial unit vectors, and the
    three pairwise AW normals are evaluated once per quadrature node and
    reused for every requested face incidence.
    """
    if not isinstance(edge_geometry, EdgeGeometry):
        raise TypeError("edge_geometry must implement EdgeGeometry.")

    indices = tuple(int(value) for value in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)
    pairs = tuple((int(cell), int(other)) for cell, other in face_pairs)

    if len(indices) != 3 or locations.shape != (3, 3):
        raise ValueError("An AW edge requires exactly three generators.")
    if len(set(indices)) != 3:
        raise ValueError("AW edge generator indices must be distinct.")
    if not np.all(np.isfinite(locations)):
        raise ValueError("AW edge generator locations must be finite.")

    if len(set(pairs)) != len(pairs):
        raise ValueError("AW face incidences must be unique.")

    positions = {index: pos for pos, index in enumerate(indices)}
    for cell_index, other_index in pairs:
        if cell_index == other_index:
            raise ValueError("An AW face requires two distinct generators.")
        if cell_index not in positions or other_index not in positions:
            raise ValueError("Selected face generators must belong to the AW edge.")

    if not pairs:
        return {}

    # Straight physical edges have identically zero geodesic curvature.
    if isinstance(edge_geometry, LineEdgeGeometry):
        return {pair: 0.0 for pair in pairs}

    if quadrature is None:
        nodes, weights = np.polynomial.legendre.leggauss(int(order))
    else:
        nodes, weights = quadrature
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)

    if len(nodes) != len(weights):
        raise ValueError("Quadrature nodes and weights must have equal length.")

    lower, upper = edge_geometry.t_min, edge_geometry.t_max
    half = 0.5 * (upper - lower)
    center = 0.5 * (upper + lower)

    # Boundary orientation is constant on a regular edge. Determine it once
    # per directed face incidence using the already validated scalar helper.
    orientation = {
        pair: aw_face_edge_boundary_orientation(
            edge_geometry=edge_geometry,
            generator_indices=indices,
            generator_locations=locations,
            cell_index=pair[0],
            face_other_index=pair[1],
            tol=tol,
        )
        for pair in pairs
    }

    totals = {pair: 0.0 for pair in pairs}

    for node, weight in zip(nodes, weights):
        t = center + half * float(node)

        # Shared edge geometry: one evaluation per quadrature node.
        point = np.asarray(edge_geometry.point(t), dtype=float)
        first = np.asarray(edge_geometry.tangent(t), dtype=float)
        second = np.asarray(edge_geometry.second_derivative(t), dtype=float)

        speed_sq = float(np.dot(first, first))
        if not np.isfinite(speed_sq) or speed_sq < tol**2:
            raise ValueError("Edge tangent is undefined during geodesic integration.")
        if not np.all(np.isfinite(second)):
            raise ValueError("Edge second derivative is non-finite.")

        # Shared radial unit vectors from all three generators.
        radial = point - locations
        distances = np.linalg.norm(radial, axis=1)

        if np.any(~np.isfinite(distances)) or np.any(distances < tol):
            raise ValueError("AW edge point coincides with a generator.")

        radial /= distances[:, None]

        # Canonical pairwise outward normals.
        n01 = radial[0] - radial[1]
        n02 = radial[0] - radial[2]
        n12 = radial[1] - radial[2]

        norm01 = float(np.linalg.norm(n01))
        norm02 = float(np.linalg.norm(n02))
        norm12 = float(np.linalg.norm(n12))

        if min(norm01, norm02, norm12) < tol:
            raise ValueError("AW pairwise surface normal is undefined on edge.")

        n01 /= norm01
        n02 /= norm02
        n12 /= norm12

        normals = {
            (0, 1): n01, (1, 0): -n01,
            (0, 2): n02, (2, 0): -n02,
            (1, 2): n12, (2, 1): -n12,
        }

        # Because ds = |r'|dt,
        #
        #   k_g ds = r'' · (n × r') / |r'|^2 dt.
        #
        # This avoids separately recomputing speed and signed k_g.
        for pair in pairs:
            cell_pos = positions[pair[0]]
            other_pos = positions[pair[1]]
            normal = normals[(cell_pos, other_pos)]
            integrand = float(np.dot(second, np.cross(normal, first)) / speed_sq)

            if not np.isfinite(integrand):
                raise ValueError("AW geodesic-curvature integrand is non-finite.")

            totals[pair] += float(weight) * integrand

    for pair in pairs:
        totals[pair] = float(orientation[pair] * half * totals[pair])

    if not all(np.isfinite(value) for value in totals.values()):
        raise ValueError("Non-finite AW geodesic-curvature integral.")

    return totals

def aw_edge_curvature_measures(
        edge_geometry,
        generator_indices,
        generator_locations,
        face_pairs=(),
        order=32,
        quadrature=None,
        tol=1e-12,
        timing=None):
    """Integrate AW edge mean- and Gaussian-curvature terms in one traversal.

    ``mean`` contains the three cell-relative edge-M contributions. ``gaussian``
    contains oriented geodesic-curvature integrals for the requested directed
    face incidences ``(cell_index, other_index)``.

    Edge position, tangent, radial unit vectors, pairwise AW normals, and the
    second derivative are shared between M and G. If ``timing`` is supplied,
    cumulative kernel timings are added to that dictionary without changing
    the numerical path.
    """
    if not isinstance(edge_geometry, EdgeGeometry):
        raise TypeError("edge_geometry must implement EdgeGeometry.")

    profile = timing is not None
    if profile:
        for key in (
            "orientation",
            "point_tangent",
            "second_derivative",
            "radial",
            "normals",
            "mean",
            "gaussian",
            "other",
        ):
            timing.setdefault(key, 0.0)

    function_start = perf_counter() if profile else None

    indices = tuple(int(value) for value in generator_indices)
    locations = np.asarray(generator_locations, dtype=float)
    pairs = tuple((int(cell), int(other)) for cell, other in face_pairs)

    if len(indices) != 3 or locations.shape != (3, 3):
        raise ValueError("An AW edge requires exactly three generators.")
    if len(set(indices)) != 3:
        raise ValueError("AW edge generator indices must be distinct.")
    if not np.all(np.isfinite(locations)):
        raise ValueError("AW edge generator locations must be finite.")
    if len(set(pairs)) != len(pairs):
        raise ValueError("AW face incidences must be unique.")

    positions = {index: pos for pos, index in enumerate(indices)}
    for cell_index, other_index in pairs:
        if cell_index == other_index:
            raise ValueError("An AW face requires two distinct generators.")
        if cell_index not in positions or other_index not in positions:
            raise ValueError("Selected face generators must belong to the AW edge.")

    if quadrature is None:
        nodes, weights = np.polynomial.legendre.leggauss(int(order))
    else:
        nodes, weights = quadrature
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)

    if len(nodes) != len(weights):
        raise ValueError("Quadrature nodes and weights must have equal length.")

    lower, upper = edge_geometry.t_min, edge_geometry.t_max
    half = 0.5 * (upper - lower)
    center = 0.5 * (upper + lower)

    t0 = perf_counter() if profile else None
    orientations = {
        pair: aw_face_edge_boundary_orientation(
            edge_geometry=edge_geometry,
            generator_indices=indices,
            generator_locations=locations,
            cell_index=pair[0],
            face_other_index=pair[1],
            tol=tol,
        )
        for pair in pairs
    }
    if profile:
        timing["orientation"] += perf_counter() - t0

    mean_totals = np.zeros(3, dtype=float)
    gauss_totals = {pair: 0.0 for pair in pairs}
    curved_for_g = bool(pairs) and not isinstance(edge_geometry, LineEdgeGeometry)

    for node, weight in zip(nodes, weights):
        t = center + half * float(node)

        # Analytic edge point and first derivative.
        t0 = perf_counter() if profile else None
        point = np.asarray(edge_geometry.point(t), dtype=float)
        first = np.asarray(edge_geometry.tangent(t), dtype=float)
        speed_sq = float(np.dot(first, first))
        if profile:
            timing["point_tangent"] += perf_counter() - t0

        if not np.isfinite(speed_sq) or speed_sq < tol**2:
            raise ValueError("Edge tangent is undefined during AW curvature integration.")

        speed = float(np.sqrt(speed_sq))

        # Radial vectors and normalization from all three generators.
        t0 = perf_counter() if profile else None
        radial = point - locations
        distances = np.linalg.norm(radial, axis=1)
        if np.any(~np.isfinite(distances)) or np.any(distances < tol):
            raise ValueError("AW edge point coincides with a generator.")
        radial /= distances[:, None]
        if profile:
            timing["radial"] += perf_counter() - t0

        # Pairwise AW face normals.
        t0 = perf_counter() if profile else None
        n01 = radial[0] - radial[1]
        n02 = radial[0] - radial[2]
        n12 = radial[1] - radial[2]
        norm01 = float(np.linalg.norm(n01))
        norm02 = float(np.linalg.norm(n02))
        norm12 = float(np.linalg.norm(n12))

        if min(norm01, norm02, norm12) < tol:
            raise ValueError("AW pairwise surface normal is undefined on edge.")

        n01 /= norm01
        n02 /= norm02
        n12 /= norm12
        if profile:
            timing["normals"] += perf_counter() - t0

        # Edge M for all three cells.
        t0 = perf_counter() if profile else None
        cos0 = np.clip(np.dot(n01, n02), -1.0, 1.0)
        cos1 = np.clip(np.dot(-n01, n12), -1.0, 1.0)
        cos2 = np.clip(np.dot(-n02, -n12), -1.0, 1.0)
        mean_totals += float(weight) * np.arccos([cos0, cos1, cos2]) * speed
        if profile:
            timing["mean"] += perf_counter() - t0

        # Straight edges have zero geodesic curvature; avoid r'' entirely.
        if curved_for_g:
            t0 = perf_counter() if profile else None
            second = np.asarray(edge_geometry.second_derivative(t), dtype=float)
            if profile:
                timing["second_derivative"] += perf_counter() - t0

            if not np.all(np.isfinite(second)):
                raise ValueError("Edge second derivative is non-finite.")

            normals = {
                (0, 1): n01, (1, 0): -n01,
                (0, 2): n02, (2, 0): -n02,
                (1, 2): n12, (2, 1): -n12,
            }

            t0 = perf_counter() if profile else None
            for pair in pairs:
                normal = normals[(positions[pair[0]], positions[pair[1]])]
                value = float(np.dot(second, np.cross(normal, first)) / speed_sq)
                if not np.isfinite(value):
                    raise ValueError("AW geodesic-curvature integrand is non-finite.")
                gauss_totals[pair] += float(weight) * value
            if profile:
                timing["gaussian"] += perf_counter() - t0

    mean_totals *= 0.5 * half
    for pair in pairs:
        gauss_totals[pair] = float(orientations[pair] * half * gauss_totals[pair])

    if not np.all(np.isfinite(mean_totals)):
        raise ValueError("Non-finite AW edge mean-curvature integral.")
    if not all(np.isfinite(value) for value in gauss_totals.values()):
        raise ValueError("Non-finite AW edge geodesic-curvature integral.")

    if profile:
        measured = sum(
            timing[key]
            for key in (
                "orientation",
                "point_tangent",
                "second_derivative",
                "radial",
                "normals",
                "mean",
                "gaussian",
            )
        )
        timing["other"] += max(perf_counter() - function_start - measured, 0.0)

    return {
        "mean": {
            indices[0]: float(mean_totals[0]),
            indices[1]: float(mean_totals[1]),
            indices[2]: float(mean_totals[2]),
        },
        "gaussian": gauss_totals,
    }

